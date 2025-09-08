function limo_mua_ttest(file_list, parameter, electrode, time_window, varargin)

% Perform paired t-tests on ERP data from LIMO.mat files at each time point
%
% FORMAT
% limo_mua_ttest(file_list, parameter, electrode, time_window)
%
% INPUTS
% file_list    - path to .txt file listing LIMO.mat files (one per line)
% parameter    - two condition parameters to compare (e.g., [1 2] for first vs second condition)
% electrode    - electrode number to analyze (e.g., 21 for Cz)
% time_window  - time window in milliseconds for analysis (e.g., [400 700])
%
% OUTPUTS
% Prints t-statistics and p-values for paired t-tests at each time point
% Uses same data extraction logic as limo_erp_plot.m
%
% Copyright (C) Devon Yanitski and Claude Code 

% Parse optional inputs
p = inputParser;
parse(p, varargin{:});

% Validate parameter input - must be exactly 2 parameters for paired t-test
if length(parameter) ~= 2
    error('Parameter must contain exactly 2 conditions for paired t-test (e.g., [1 2])');
end

% Handle file list input
if ischar(file_list)
    if ~exist(file_list, 'file')
        error('File list %s does not exist', file_list);
    end
    [Names, Paths, Files] = limo_get_files([], [], [], file_list);
else
    error('file_list must be a path to a .txt file listing LIMO.mat files');
end

if isempty(Names)
    error('No files found in file list');
elseif size(Names, 2) < 1
    error('At least one LIMO.mat file required');
end

% Validate that all files are LIMO.mat files
is_limo = zeros(1, size(Names, 2));
for i = size(Names, 2):-1:1
    if contains(Names{i}, 'LIMO')
        is_limo(i) = 1;
    end
end

if ~all(is_limo)
    error('All files must be LIMO.mat files');
end

% Ensure parameter is a row vector
if iscolumn(parameter)
    parameter = parameter';
end

% Initialize data storage
n_subjects = size(Names, 2);
n_parameters = 2;  % Always 2 for paired t-test
data = cell(n_parameters, 1);
time_vector = [];

fprintf('Processing %d subjects for paired t-test between parameters %d and %d...\n', n_subjects, parameter(1), parameter(2));

% Process all subjects - extract data for both parameters
subject_data = cell(n_subjects, n_parameters);

parfor i = 1:n_subjects
    fprintf('Processing subject %g\n', i);
    
    % Load LIMO structure once per subject
    LIMO = load(fullfile(Paths{i}, 'LIMO.mat'));
    LIMO = LIMO.LIMO;
    
    % Load Yr.mat once per subject
    Yr = load(fullfile(Paths{i}, 'Yr.mat'));
    Yr = Yr.(cell2mat(fieldnames(Yr)));
    
    % Process both parameters for this subject
    subject_means = cell(n_parameters, 1);
    
    for p_idx = 1:n_parameters
        current_param = parameter(p_idx);
        
        % Check if parameter is valid categorical variable
        if current_param <= sum(LIMO.design.nb_conditions + LIMO.design.nb_interactions) || ...
                current_param == size(LIMO.design.X, 2)
            
            % Find trials for this condition
            trial_indices = logical(LIMO.design.X(:, current_param) == 1);
            
            if sum(trial_indices) == 0
                warning('No trials found for parameter %d in subject %d', current_param, i);
                continue;
            end
            
            % Extract electrode data for this parameter
            electrode_data = squeeze(Yr(electrode, :, trial_indices));
            
            % Calculate mean across trials for this subject and parameter
            if ndims(electrode_data) == 1
                % Only one trial - electrode_data is already the time series
                subject_means{p_idx} = electrode_data;
            else
                % Multiple trials - average across trials (dimension 2)
                subject_means{p_idx} = nanmean(electrode_data, 2);
            end
            
        else
            error('Parameter %d appears to be continuous - only categorical parameters supported', current_param);
        end
    end
    
    % Store results for this subject
    for p_idx = 1:n_parameters
        subject_data{i, p_idx} = subject_means{p_idx};
    end
end

% Combine results for each parameter
for p_idx = 1:n_parameters
    param_data = [];
    for i = 1:n_subjects
        if ~isempty(subject_data{i, p_idx})
            if isempty(param_data)
                param_data = zeros(length(subject_data{i, p_idx}), n_subjects);
            end
            param_data(:, i) = subject_data{i, p_idx};
        end
    end
    data{p_idx} = param_data;
end

% Get time vector from first subject  
if isempty(time_vector)
    LIMO_first = load(fullfile(Paths{1}, 'LIMO.mat'));
    LIMO_first = LIMO_first.LIMO;
    if isfield(LIMO_first.data, 'timevect')
        time_vector = LIMO_first.data.timevect;
    else
        time_vector = LIMO_first.data.start:(1000/LIMO_first.data.sampling_rate):LIMO_first.data.end;
    end
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

% Extract data for the specified time window
data1_window = data{1}(time_window_samples, :);  % Parameter 1 data
data2_window = data{2}(time_window_samples, :);  % Parameter 2 data
time_window_ms = time_vector(time_window_samples);

% Perform paired t-tests at each time point
n_timepoints = length(time_window_samples);
t_stats = zeros(n_timepoints, 1);
p_values = zeros(n_timepoints, 1);

fprintf('\nPerforming paired t-tests at %d time points...\n', n_timepoints);
fprintf('Time window: %.1f to %.1f ms\n', time_window_ms(1), time_window_ms(end));
fprintf('Electrode: %d\n', electrode);
fprintf('Parameters compared: %d vs %d\n\n', parameter(1), parameter(2));

for t = 1:n_timepoints
    % Extract data for this time point across all subjects
    x1 = data1_window(t, :)';  % Parameter 1, all subjects
    x2 = data2_window(t, :)';  % Parameter 2, all subjects
    
    % Remove subjects with NaN values
    valid_subjects = ~isnan(x1) & ~isnan(x2);
    x1_clean = x1(valid_subjects);
    x2_clean = x2(valid_subjects);
    
    if length(x1_clean) < 2
        warning('Insufficient valid subjects at time point %d (%.1f ms)', t, time_window_ms(t));
        t_stats(t) = NaN;
        p_values(t) = NaN;
        continue;
    end
    
    % Perform paired t-test
    [~, p, ~, stats] = ttest(x1_clean, x2_clean);
    
    t_stats(t) = stats.tstat;
    p_values(t) = p;
end

% Print results
fprintf('Results Summary:\n');
fprintf('================\n');
fprintf('Index\t\tTime (ms)\tt-statistic\tp-value\n');
fprintf('-----\t\t---------\t-----------\t-------\n');

for t = 1:n_timepoints
    fprintf('%d\t\t%.1f\t\t%.4f\t\t%.6f\n', ...
        t, time_window_ms(t), t_stats(t), p_values(t));
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
    fprintf('Significant time points (p < 0.05): %d / %d\n', ...
        sum(p_values(valid_tests) < 0.05), n_valid);
    fprintf('Significant time points (p < 0.01): %d / %d\n', ...
        sum(p_values(valid_tests) < 0.01), n_valid);
end

end