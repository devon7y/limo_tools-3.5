function results = limo_likelihood(test_file, H0_file, varargin)
% LIMO_LIKELIHOOD - Bayes Factor analysis using Kang et al. (2015) approach
%
% Computes likelihood ratios (Bayes Factors) for mass univariate EEG analysis
% using the method of Kang et al. (2015) adapted to use LIMO's empirical
% bootstrap null distribution instead of parametric assumptions.
%
% FORMAT: results = limo_likelihood(test_file, H0_file)
%         results = limo_likelihood(test_file, H0_file, 'tau', 1.0, 'save', true)
%
% INPUTS:
%   test_file - Path to LIMO test file (e.g., 'paired_samples_ttest_parameter_35.mat')
%   H0_file   - Path to H0 bootstrap file (e.g., 'H0_paired_samples_ttest_parameter_35.mat')
%
% OPTIONAL PARAMETERS (Name-Value pairs):
%   'tau'     - Prior SD on standardized effect size (default: 1.0)
%               Kang et al. recommend: 0.5 (conservative), 1.0 (medium), sqrt(2) (unit info)
%   'save'    - Save augmented file with BF values (default: true)
%   'plot'    - Generate diagnostic plots (default: false)
%
% OUTPUTS:
%   results structure with fields:
%     .BF              - Bayes Factors [channels × time] or [channels × freq × time]
%     .logBF           - log(BF) for numerical stability
%     .BF_H0_dist      - BF distribution from bootstrap [ch × time × n_boots]
%     .evidence_labels - Evidence categories ('Strong H1', 'Inconclusive', etc.)
%     .params          - Analysis parameters (tau, n, formula, etc.)
%     .test_file_path  - Path to augmented test file (if saved)
%
% REFERENCE:
%   Kang, I., Cohen, A. L., & Luck, S. J. (2015).
%   A Bayesian approach to measuring the strength of evidence.
%   Adapted for mass univariate EEG analysis with empirical H0.
%
% METHOD:
%   For each channel×time point, computes:
%   BF₁₀ = p(data|H1) / p(data|H0)
%
%   Where H1: effect size δ ~ Normal(0, τ²) [composite alternative]
%         H0: effect size δ = 0 [point null]
%
%   Using closed-form Bayes Factor (normal-normal conjugacy):
%   BF = √(1/(1+n×τ²)) × exp[t²/2 × n×τ²/(1+n×τ²)]
%
%   The null distribution p(data|H0) is estimated empirically using
%   LIMO's bootstrap samples rather than assuming a parametric form.
%
% INTERPRETATION:
%   BF > 10:   Strong evidence for alternative hypothesis (effect exists)
%   BF > 3:    Moderate evidence for H1
%   1/3 < BF < 3: Inconclusive (data don't favor either hypothesis strongly)
%   BF < 1/3:  Moderate evidence for null (no effect)
%   BF < 1/10: Strong evidence for H0
%
% EXAMPLE:
%   % Compute Bayes Factors for paired t-test
%   results = limo_likelihood('paired_samples_ttest_parameter_1.mat', ...
%                             'H0/H0_paired_samples_ttest_parameter_1.mat', ...
%                             'tau', 1.0, 'plot', true);
%
%   % Access results
%   fprintf('Number of strong H1 evidence: %d\n', results.n_strong_H1);
%   imagesc(results.logBF);  % Plot log Bayes Factors
%
% See also: LIMO_RANDOM_ROBUST, LIMO_STAT_VALUES, LIMO_MAX_CORRECTION
%
% Devon Yanitski & Claude (Anthropic) 2025
% Adapted from: Kang et al. (2015)
% ------------------------------
%  Copyright (C) LIMO Team 2025

%% Parse inputs
p = inputParser;
addRequired(p, 'test_file', @ischar);
addRequired(p, 'H0_file', @ischar);
addParameter(p, 'tau', 1.0, @(x) isnumeric(x) && x > 0);
addParameter(p, 'save', true, @islogical);
addParameter(p, 'plot', false, @islogical);
parse(p, test_file, H0_file, varargin{:});

tau = p.Results.tau;
save_output = p.Results.save;
do_plot = p.Results.plot;

%% Validate files exist
if ~exist(test_file, 'file')
    error('Test file not found: %s', test_file);
end

if ~exist(H0_file, 'file')
    error('H0 file not found: %s', H0_file);
end

fprintf('\n=== LIMO LIKELIHOOD RATIO ANALYSIS ===\n');
fprintf('Method: Kang et al. (2015) with empirical H0\n');
fprintf('Test file: %s\n', test_file);
fprintf('H0 file: %s\n', H0_file);
fprintf('Prior τ (standardized effect SD): %.3f\n', tau);
fprintf('======================================\n\n');

%% Load test file and identify type
fprintf('Loading test file...\n');
test_data_struct = load(test_file);
field_names = fieldnames(test_data_struct);
test_data = test_data_struct.(field_names{1});

% Determine test type from filename
[~, filename, ~] = fileparts(test_file);
if contains(filename, 'one_sample')
    test_type = 'one_sample';
    var_name = 'one_sample';
elseif contains(filename, 'two_samples')
    test_type = 'two_samples';
    var_name = 'two_samples';
elseif contains(filename, 'paired_samples')
    test_type = 'paired_samples';
    var_name = 'paired_samples';
else
    error('Cannot determine test type from filename. Expected: one_sample, two_samples, or paired_samples');
end

fprintf('Test type detected: %s\n', test_type);

%% Extract t-statistics and degrees of freedom
% Data structure: [channels × time × 5] or [channels × freq × time × 5]
% Dimension 3 (or 4): [mean, SE, df, t, p]

data_dims = ndims(test_data);
if data_dims == 3
    % Time or frequency domain: [channels × time × 5]
    t_obs = test_data(:,:,4);  % t-statistic
    df = test_data(:,:,3);     % degrees of freedom
    n = df + 1;                % effective sample size
    is_tf = false;
elseif data_dims == 4
    % Time-frequency domain: [channels × freq × time × 5]
    t_obs = test_data(:,:,:,4);  % t-statistic
    df = test_data(:,:,:,3);     % degrees of freedom
    n = df + 1;                  % effective sample size
    is_tf = true;
else
    error('Unexpected data dimensions: %d. Expected 3 or 4.', data_dims);
end

% Check for all NaN
if all(isnan(t_obs(:)))
    error('All t-statistics are NaN. Cannot compute Bayes Factors.');
end

fprintf('Data dimensions: ');
if is_tf
    fprintf('[%d channels × %d frequencies × %d times]\n', size(t_obs,1), size(t_obs,2), size(t_obs,3));
else
    fprintf('[%d channels × %d times]\n', size(t_obs,1), size(t_obs,2));
end

%% Load H0 bootstrap file
fprintf('\nLoading H0 bootstrap file...\n');
H0_data_struct = load(H0_file);
H0_field_names = fieldnames(H0_data_struct);
H0_data = H0_data_struct.(H0_field_names{1});

% Extract t-statistics from H0
% H0 structure: [channels × time × 2 × n_bootstraps] or [ch × freq × time × 2 × n_boots]
% Dimension 1: t-statistics under H0
% Dimension 2: p-values under H0 (not used)

if is_tf
    t_H0 = squeeze(H0_data(:,:,:,1,:));  % [ch × freq × time × n_boots]
    n_boots = size(t_H0, 4);
else
    t_H0 = squeeze(H0_data(:,:,1,:));  % [ch × time × n_boots]
    n_boots = size(t_H0, 3);
end

fprintf('Number of bootstrap samples: %d\n', n_boots);

% Check dimension compatibility
if is_tf
    if ~isequal([size(t_obs,1), size(t_obs,2), size(t_obs,3)], [size(t_H0,1), size(t_H0,2), size(t_H0,3)])
        error('Dimension mismatch between test file and H0 file');
    end
else
    if ~isequal(size(t_obs), [size(t_H0,1), size(t_H0,2)])
        error('Dimension mismatch between test file and H0 file');
    end
end

%% Compute observed Bayes Factors using Kang et al. (2015) formula
fprintf('\nComputing observed Bayes Factors...\n');

% Kang formula: BF = √(1/(1+n×τ²)) × exp[t²/2 × n×τ²/(1+n×τ²)]
% log(BF) = -0.5×log(1+n×τ²) + (t²/2) × [n×τ²/(1+n×τ²)]

n_tau_sq = n .* (tau^2);  % n×τ²
scaling_factor = n_tau_sq ./ (1 + n_tau_sq);  % n×τ²/(1+n×τ²)

% Compute log(BF) for numerical stability
logBF_obs = -0.5 * log(1 + n_tau_sq) + (t_obs.^2 / 2) .* scaling_factor;

% Convert to BF (for those who prefer non-log scale)
BF_obs = exp(logBF_obs);

% Handle NaN propagation
BF_obs(isnan(t_obs)) = NaN;
logBF_obs(isnan(t_obs)) = NaN;

fprintf('Observed BF range: [%.2e, %.2e]\n', min(BF_obs(:)), max(BF_obs(:)));
fprintf('Observed log(BF) range: [%.2f, %.2f]\n', min(logBF_obs(:)), max(logBF_obs(:)));

%% Compute BF for bootstrap H0 samples
fprintf('\nComputing BF for %d bootstrap samples...\n', n_boots);

% Expand n and scaling_factor to match bootstrap dimension
if is_tf
    n_expanded = repmat(n, [1 1 1 n_boots]);
    n_tau_sq_H0 = n_expanded .* (tau^2);
else
    n_expanded = repmat(n, [1 1 n_boots]);
    n_tau_sq_H0 = n_expanded .* (tau^2);
end

scaling_factor_H0 = n_tau_sq_H0 ./ (1 + n_tau_sq_H0);

% Compute log(BF) for H0 samples
logBF_H0 = -0.5 * log(1 + n_tau_sq_H0) + (t_H0.^2 / 2) .* scaling_factor_H0;
BF_H0 = exp(logBF_H0);

fprintf('H0 BF range: [%.2e, %.2e]\n', min(BF_H0(:)), max(BF_H0(:)));
fprintf('H0 log(BF) range: [%.2f, %.2f]\n', min(logBF_H0(:)), max(logBF_H0(:)));

%% Categorize evidence
fprintf('\nCategorizing evidence...\n');

% Following Kang et al. interpretation thresholds
evidence_labels = cell(size(BF_obs));
evidence_numeric = zeros(size(BF_obs));

% Strong evidence for H1 (BF > 10)
strong_H1 = BF_obs > 10;
evidence_labels(strong_H1) = {'Strong H1'};
evidence_numeric(strong_H1) = 2;

% Moderate evidence for H1 (3 < BF <= 10)
moderate_H1 = BF_obs > 3 & BF_obs <= 10;
evidence_labels(moderate_H1) = {'Moderate H1'};
evidence_numeric(moderate_H1) = 1;

% Inconclusive (1/3 <= BF <= 3)
inconclusive = BF_obs >= 1/3 & BF_obs <= 3;
evidence_labels(inconclusive) = {'Inconclusive'};
evidence_numeric(inconclusive) = 0;

% Moderate evidence for H0 (1/10 <= BF < 1/3)
moderate_H0 = BF_obs < 1/3 & BF_obs >= 1/10;
evidence_labels(moderate_H0) = {'Moderate H0'};
evidence_numeric(moderate_H0) = -1;

% Strong evidence for H0 (BF < 1/10)
strong_H0 = BF_obs < 1/10;
evidence_labels(strong_H0) = {'Strong H0'};
evidence_numeric(strong_H0) = -2;

% NaN handling
evidence_labels(isnan(BF_obs)) = {'NaN'};
evidence_numeric(isnan(BF_obs)) = NaN;

% Summary statistics
n_strong_H1 = sum(strong_H1(:));
n_moderate_H1 = sum(moderate_H1(:));
n_inconclusive = sum(inconclusive(:));
n_moderate_H0 = sum(moderate_H0(:));
n_strong_H0 = sum(strong_H0(:));
n_total = numel(BF_obs) - sum(isnan(BF_obs(:)));

fprintf('\nEvidence Summary:\n');
fprintf('  Strong H1 (BF>10):        %6d (%.1f%%)\n', n_strong_H1, 100*n_strong_H1/n_total);
fprintf('  Moderate H1 (3<BF<=10):   %6d (%.1f%%)\n', n_moderate_H1, 100*n_moderate_H1/n_total);
fprintf('  Inconclusive (1/3<BF<3):  %6d (%.1f%%)\n', n_inconclusive, 100*n_inconclusive/n_total);
fprintf('  Moderate H0 (1/10<BF<1/3):%6d (%.1f%%)\n', n_moderate_H0, 100*n_moderate_H0/n_total);
fprintf('  Strong H0 (BF<1/10):      %6d (%.1f%%)\n', n_strong_H0, 100*n_strong_H0/n_total);

%% Create results structure
results = struct();
results.BF = BF_obs;
results.logBF = logBF_obs;
results.BF_H0_dist = BF_H0;
results.logBF_H0_dist = logBF_H0;
results.evidence_labels = evidence_labels;
results.evidence_numeric = evidence_numeric;

% Summary statistics
results.n_strong_H1 = n_strong_H1;
results.n_moderate_H1 = n_moderate_H1;
results.n_inconclusive = n_inconclusive;
results.n_moderate_H0 = n_moderate_H0;
results.n_strong_H0 = n_strong_H0;
results.proportion_strong_H1 = n_strong_H1 / n_total;
results.proportion_inconclusive = n_inconclusive / n_total;

% Parameters
results.params.tau = tau;
results.params.n = n;
results.params.formula = 'Kang et al. (2015): BF = sqrt(1/(1+n*tau^2)) * exp(t^2/2 * n*tau^2/(1+n*tau^2))';
results.params.test_type = test_type;
results.params.n_bootstraps = n_boots;
results.params.is_time_frequency = is_tf;

%% Save augmented file
if save_output
    fprintf('\nSaving augmented test file...\n');

    % Create augmented data with BF in 6th dimension, log(BF) in 7th
    if is_tf
        % Time-frequency: [ch × freq × time × 7]
        augmented_data = cat(4, test_data, BF_obs, logBF_obs);
    else
        % Time or frequency: [ch × time × 7]
        augmented_data = cat(3, test_data, BF_obs, logBF_obs);
    end

    % Create new filename
    [filepath, filename_only, ext] = fileparts(test_file);
    new_filename = [filename_only '_likelihood' ext];
    if isempty(filepath)
        new_filepath = new_filename;
    else
        new_filepath = fullfile(filepath, new_filename);
    end

    % Save with same variable name as original
    eval([var_name ' = augmented_data;']);
    save(new_filepath, var_name, '-v7.3');

    results.test_file_path = new_filepath;
    fprintf('Saved: %s\n', new_filepath);
    fprintf('  Dimension 6: Bayes Factors (BF)\n');
    fprintf('  Dimension 7: log(BF)\n');
end

%% Diagnostic plots
if do_plot
    fprintf('\nGenerating diagnostic plots...\n');

    figure('Name', 'LIMO Likelihood Ratio Analysis', 'Position', [100 100 1400 900]);

    % For plotting, reduce to 2D if time-frequency
    if is_tf
        % Average over frequencies for visualization
        logBF_plot = squeeze(mean(logBF_obs, 2, 'omitnan'));
        evidence_plot = squeeze(mean(evidence_numeric, 2, 'omitnan'));
        t_obs_plot = squeeze(mean(t_obs, 2, 'omitnan'));
    else
        logBF_plot = logBF_obs;
        evidence_plot = evidence_numeric;
        t_obs_plot = t_obs;
    end

    % Subplot 1: Observed log(BF) heatmap
    subplot(2,3,1);
    imagesc(logBF_plot');
    colorbar;
    title(sprintf('Observed log(BF) - τ=%.2f', tau));
    xlabel('Channel'); ylabel('Time');
    caxis([-5 5]);  % Symmetric color scale
    colormap(subplot(2,3,1), redblue);

    % Subplot 2: Evidence categories
    subplot(2,3,2);
    imagesc(evidence_plot');
    colorbar;
    title('Evidence Categories');
    xlabel('Channel'); ylabel('Time');
    caxis([-2 2]);
    colormap(subplot(2,3,2), redblue);

    % Subplot 3: Distribution of log(BF) - observed vs H0
    subplot(2,3,3);
    edges = linspace(-10, 10, 50);
    histogram(logBF_obs(:), edges, 'FaceAlpha', 0.5, 'Normalization', 'probability');
    hold on;
    histogram(logBF_H0(:), edges, 'FaceAlpha', 0.5, 'Normalization', 'probability');
    legend({'Observed', 'H0 Bootstrap'}, 'Location', 'best');
    xlabel('log(BF)'); ylabel('Probability');
    title('Distribution of log(BF)');
    xline(0, 'k--', 'LineWidth', 1.5);  % BF = 1
    xline(log(10), 'r--', 'Strong H1');
    xline(log(1/10), 'r--', 'Strong H0');
    grid on;

    % Subplot 4: BF vs t-statistic (scatterplot)
    subplot(2,3,4);
    % Sample for efficiency
    n_sample = min(5000, numel(t_obs_plot));
    idx_sample = randsample(numel(t_obs_plot), n_sample);
    scatter(t_obs_plot(idx_sample), logBF_plot(idx_sample), 10, 'filled', 'MarkerFaceAlpha', 0.3);
    xlabel('t-statistic'); ylabel('log(BF)');
    title('Relationship: log(BF) vs t-statistic');
    grid on;
    hold on;
    % Theoretical curve
    t_range = linspace(-8, 8, 100);
    n_mean = nanmean(n(:));
    logBF_theory = -0.5*log(1 + n_mean*tau^2) + (t_range.^2/2) * (n_mean*tau^2/(1 + n_mean*tau^2));
    plot(t_range, logBF_theory, 'r-', 'LineWidth', 2);
    legend('Data', sprintf('Theory (n=%.0f)', n_mean), 'Location', 'best');

    % Subplot 5: Empirical null at max |logBF| location
    subplot(2,3,5);
    [~, max_idx] = max(abs(logBF_plot(:)));
    [ch_max, t_max] = ind2sub(size(logBF_plot), max_idx);

    if is_tf
        % Average over frequency for this channel-time
        logBF_H0_slice = squeeze(mean(logBF_H0(ch_max, :, t_max, :), 2));
    else
        logBF_H0_slice = squeeze(logBF_H0(ch_max, t_max, :));
    end

    histogram(logBF_H0_slice, 30, 'Normalization', 'probability');
    hold on;
    xline(logBF_plot(ch_max, t_max), 'r-', 'LineWidth', 2.5, 'Label', 'Observed');
    xline(log(10), 'k--', 'Strong H1');
    xline(log(1/10), 'k--', 'Strong H0');
    title(sprintf('H0 dist at max |logBF| (ch %d, t %d)', ch_max, t_max));
    xlabel('log(BF)'); ylabel('Probability');
    grid on;

    % Subplot 6: Evidence summary bar chart
    subplot(2,3,6);
    evidence_counts = [n_strong_H0, n_moderate_H0, n_inconclusive, n_moderate_H1, n_strong_H1];
    bar(evidence_counts);
    set(gca, 'XTickLabel', {'Strong H0', 'Mod H0', 'Inconcl.', 'Mod H1', 'Strong H1'});
    ylabel('Count');
    title('Evidence Summary');
    xtickangle(45);
    grid on;

    % Add percentage labels
    for i = 1:length(evidence_counts)
        text(i, evidence_counts(i), sprintf('%.1f%%', 100*evidence_counts(i)/n_total), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
end

fprintf('\n=== ANALYSIS COMPLETE ===\n');

end

%% Helper function: Red-Blue colormap
function cmap = redblue(m)
    % REDBLUE Creates a red-white-blue colormap
    if nargin < 1
        m = 256;
    end

    % Create diverging colormap: blue (negative) -> white (zero) -> red (positive)
    r = [0*ones(m/2,1); linspace(0,1,m/2)'];
    g = [linspace(0,1,m/2)'; linspace(1,0,m/2)'];
    b = [linspace(1,0,m/2)'; 0*ones(m/2,1)];

    cmap = [r g b];
end
