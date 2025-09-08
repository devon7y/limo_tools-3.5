function limo_mua_ttest_y2r_y1r(y1r_file, y2r_file, electrode, time_window, varargin)

% Perform paired t-tests on ERP data from Y1r.mat and Y2r.mat files at each time point
%
% FORMAT
% limo_mua_ttest_y2r_y1r(y1r_file, y2r_file, electrode, time_window)
%
% INPUTS
% y1r_file     - path to Y1r.mat file (condition 1 parameter estimates)
% y2r_file     - path to Y2r.mat file (condition 2 parameter estimates)
% electrode    - electrode number to analyze (e.g., 21 for Cz)
% time_window  - time window in milliseconds for analysis (e.g., [400 700])
%
% OUTPUTS
% Prints t-statistics and p-values for paired t-tests at each time point
% Uses limo_ttest.m for statistical computation
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
fprintf('Processing %d subjects for paired t-test between Y1r and Y2r...\n', n_subjects);

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
data1_window = squeeze(Y1r(electrode, time_window_samples, :));  % Y1r data: time x subjects
data2_window = squeeze(Y2r(electrode, time_window_samples, :));  % Y2r data: time x subjects
time_window_ms = time_vector(time_window_samples);

% Ensure data is oriented correctly (time points x subjects)
if size(data1_window, 2) ~= n_subjects
    data1_window = data1_window';
    data2_window = data2_window';
end

% Perform paired t-tests at each time point using limo_ttest
n_timepoints = length(time_window_samples);
t_stats = zeros(n_timepoints, 1);
p_values = zeros(n_timepoints, 1);
means_diff = zeros(n_timepoints, 1);
std_errors = zeros(n_timepoints, 1);

fprintf('\nPerforming paired t-tests at %d time points...\n', n_timepoints);
fprintf('Time window: %.1f to %.1f ms\n', time_window_ms(1), time_window_ms(end));
fprintf('Electrode: %d\n', electrode);
fprintf('Using limo_ttest for statistical computation\n\n');

for t = 1:n_timepoints
    % Extract data for this time point across all subjects
    x1 = data1_window(t, :);  % Y1r data for this time point
    x2 = data2_window(t, :);  % Y2r data for this time point
    
    % Remove subjects with NaN values
    valid_subjects = ~isnan(x1) & ~isnan(x2);
    x1_clean = x1(valid_subjects);
    x2_clean = x2(valid_subjects);
    
    if length(x1_clean) < 2
        warning('Insufficient valid subjects at time point %d (%.1f ms)', t, time_window_ms(t));
        t_stats(t) = NaN;
        p_values(t) = NaN;
        means_diff(t) = NaN;
        std_errors(t) = NaN;
        continue;
    end
    
    % Perform paired t-test using limo_ttest (type = 1 for paired)
    % limo_ttest expects data as row vectors for the last dimension
    try
        [m, dfe, ci, sd, n, t_val, p_val] = limo_ttest(1, x1_clean, x2_clean, 0.05);
        
        t_stats(t) = t_val;
        p_values(t) = p_val;
        means_diff(t) = m;
        std_errors(t) = sd / sqrt(n);
        
    catch ME
        warning('limo_ttest failed at time point %d (%.1f ms): %s', t, time_window_ms(t), ME.message);
        t_stats(t) = NaN;
        p_values(t) = NaN;
        means_diff(t) = NaN;
        std_errors(t) = NaN;
    end
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
    fprintf('\nSummary Statistics:\n');
    fprintf('==================\n');
    fprintf('Valid time points: %d / %d\n', n_valid, n_timepoints);
    fprintf('Mean t-statistic: %.4f\n', mean(t_stats(valid_tests)));
    fprintf('Max |t-statistic|: %.4f\n', max(abs(t_stats(valid_tests))));
    fprintf('Min p-value: %.6f\n', min(p_values(valid_tests)));
    fprintf('Mean difference (Y1r - Y2r): %.4f\n', mean(means_diff(valid_tests)));
    fprintf('Significant time points (p < 0.05): %d / %d\n', ...
        sum(p_values(valid_tests) < 0.05), n_valid);
    fprintf('Significant time points (p < 0.01): %d / %d\n', ...
        sum(p_values(valid_tests) < 0.01), n_valid);
end

end