function limo_mua_yuen_ttest_y2r_y1r(y1r_file, y2r_file, electrode, time_window, varargin)

% Perform robust paired t-tests on ERP data from Y1r.mat and Y2r.mat files at each time point
%
% FORMAT
% limo_mua_yuen_ttest_y2r_y1r(y1r_file, y2r_file, electrode, time_window)
%
% INPUTS
% y1r_file     - path to Y1r.mat file (condition 1 parameter estimates)
% y2r_file     - path to Y2r.mat file (condition 2 parameter estimates)
% electrode    - electrode number to analyze (e.g., 21 for Cz)
% time_window  - time window in milliseconds for analysis (e.g., [400 700])
%
% OUTPUTS
% Prints t-statistics and p-values for robust paired t-tests at each time point
% Uses limo_yuend_ttest.m for statistical computation (Yuen's robust method)
%
% Copyright (C) Devon Yanitski and Claude Code 

% Parse optional inputs
p = inputParser;
parse(p, varargin{:});

% Handle file inputs
if ~ischar(y1r_file) || ~exist(y1r_file, 'file')
    error('Y1r file %s does not exist', y1r_file);
end

if ~ischar(y2r_file) || ~exist(y2r_file, 'file')
    error('Y2r file %s does not exist', y2r_file);
end

% Load Y1r and Y2r data
fprintf('Loading Y1r data from %s...\n', y1r_file);
Y1r_data = load(y1r_file);
Y1r = Y1r_data.(cell2mat(fieldnames(Y1r_data)));

fprintf('Loading Y2r data from %s...\n', y2r_file);
Y2r_data = load(y2r_file);
Y2r = Y2r_data.(cell2mat(fieldnames(Y2r_data)));

% Validate data dimensions
if size(Y1r, 1) ~= size(Y2r, 1) || size(Y1r, 2) ~= size(Y2r, 2) || size(Y1r, 3) ~= size(Y2r, 3)
    error('Y1r and Y2r files have incompatible dimensions');
end

n_subjects = size(Y1r, 3);
fprintf('Processing %d subjects for robust paired t-test between Y1r and Y2r...\n', n_subjects);

% Get time vector - we need to load a LIMO.mat file to get timing info
% Try to find LIMO.mat in the same directory as Y1r
[y1r_path, ~, ~] = fileparts(y1r_file);
limo_file = fullfile(y1r_path, 'LIMO.mat');

if exist(limo_file, 'file')
    LIMO_data = load(limo_file);
    LIMO = LIMO_data.LIMO;
    
    if isfield(LIMO.data, 'timevect')
        time_vector = LIMO.data.timevect;
    else
        time_vector = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
    end
else
    error('Could not find LIMO.mat file in %s to get time vector information', y1r_path);
end

% Convert time window from milliseconds to sample indices
if length(time_window) == 2
    % Time window given as [start_ms end_ms]
    start_ms = time_window(1);
    end_ms = time_window(2);
    
    % Find closest sample indices
    [~, start_idx] = min(abs(time_vector - start_ms));
    [~, end_idx] = min(abs(time_vector - end_ms));
    
    time_window_samples = start_idx:end_idx;
else
    error('Time window must be a 2-element vector [start_ms end_ms]');
end

% Validate time window
if max(time_window_samples) > length(time_vector) || min(time_window_samples) < 1
    error('Time window exceeds data length. Data ranges from %.1f to %.1f ms.', time_vector(1), time_vector(end));
end

fprintf('Time window: %.1f to %.1f ms corresponds to samples %d to %d\n', ...
    start_ms, end_ms, time_window_samples(1), time_window_samples(end));

% Extract data for the specified electrode and time window
% For limo_yuend_ttest, we need to prepare 3D data: 1 electrode x time_points x subjects
data1_3d = Y1r(electrode, time_window_samples, :);  % 1 x time_points x subjects
data2_3d = Y2r(electrode, time_window_samples, :);  % 1 x time_points x subjects

% Reshape to ensure correct dimensions for limo_yuend_ttest
data1_3d = reshape(data1_3d, [1, length(time_window_samples), n_subjects]);
data2_3d = reshape(data2_3d, [1, length(time_window_samples), n_subjects]);

time_window_ms = time_vector(time_window_samples);
n_timepoints = length(time_window_samples);

fprintf('\nPerforming robust paired t-tests using Yuen''s method...\n');
fprintf('Time window: %.1f to %.1f ms\n', time_window_ms(1), time_window_ms(end));
fprintf('Electrode: %d\n', electrode);
fprintf('Using limo_yuend_ttest for robust statistical computation (20%% trimming)\n\n');

% Perform robust paired t-test using limo_yuend_ttest
% This function performs the test across all time points simultaneously
try
    [Ty, diff, se, CI, p_vals, tcrit, df] = limo_yuend_ttest(data1_3d, data2_3d, 20, 0.05);
    
    % Extract results for the single electrode
    t_stats = squeeze(Ty(1, :));
    p_values = squeeze(p_vals(1, :));
    means_diff = squeeze(diff(1, :));
    std_errors = squeeze(se(1, :));
    
catch ME
    error('limo_yuend_ttest failed: %s', ME.message);
end

% Print results
fprintf('Results Summary:\n');
fprintf('================\n');
fprintf('Index\t\tTime (ms)\tt-statistic\tp-value\t\tMean Diff\tSE\n');
fprintf('-----\t\t---------\t-----------\t-------\t\t---------\t--\n');

for t = 1:n_timepoints
    fprintf('%d\t\t%.1f\t\t%.4f\t\t%.6f\t%.4f\t\t%.4f\n', ...
        t, time_window_ms(t), t_stats(t), p_values(t), means_diff(t), std_errors(t));
end

% Print summary statistics
valid_tests = ~isnan(t_stats);
n_valid = sum(valid_tests);

if n_valid > 0
    fprintf('\nSummary Statistics (Yuen''s Robust Method):\n');
    fprintf('==========================================\n');
    fprintf('Valid time points: %d / %d\n', n_valid, n_timepoints);
    fprintf('Mean t-statistic: %.4f\n', mean(t_stats(valid_tests)));
    fprintf('Max |t-statistic|: %.4f\n', max(abs(t_stats(valid_tests))));
    fprintf('Min p-value: %.6f\n', min(p_values(valid_tests)));
    fprintf('Mean difference (Y1r - Y2r): %.4f\n', mean(means_diff(valid_tests)));
    fprintf('Trimming percentage: 20%%\n');
    fprintf('Degrees of freedom: %.0f\n', df);
    fprintf('Significant time points (p < 0.05): %d / %d\n', ...
        sum(p_values(valid_tests) < 0.05), n_valid);
    fprintf('Significant time points (p < 0.01): %d / %d\n', ...
        sum(p_values(valid_tests) < 0.01), n_valid);
end

end