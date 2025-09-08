function limo_erp_plot(file_list, parameter, electrode, varargin)

% Plot ERPs from LIMO.mat files for specified condition(s) and electrode
%
% FORMAT
% limo_erp_plot(file_list, parameter, electrode)
% limo_erp_plot(file_list, parameter, electrode, 'estimator', 'Mean')
%
% INPUTS
% file_list    - path to .txt file listing LIMO.mat files (one per line)
% parameter    - condition parameter(s) to plot (e.g., 1 for first condition, or [1 2 3] for multiple)
% electrode    - electrode number to plot (e.g., 21 for Cz)
%
% OPTIONAL INPUTS
% 'title'      - custom plot title (default: auto-generated)
% 'show_individual' - plot individual subject traces (default: false)
%
% OUTPUTS
% Creates ERP plot showing voltage over time for specified condition(s) and electrode
% Multiple parameters are plotted as separate lines on the same plot
% Uses mean estimator at both within-subject and across-subject levels
% Uses same logic as limo_central_tendency_and_ci.m and limo_add_plots.m
%
% Copyright (C) Devon Yanitski and Claude Code 

% Parse optional inputs
p = inputParser;
addParameter(p, 'title', '', @ischar);
addParameter(p, 'show_individual', false, @islogical);
parse(p, varargin{:});
plot_title = p.Results.title;
show_individual = p.Results.show_individual;

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
n_parameters = length(parameter);
data = cell(n_parameters, 1);
time_vector = [];

fprintf('Processing %d subjects for %d parameters...\n', n_subjects, n_parameters);

% Process all subjects once - extract data for all parameters simultaneously
subject_data = cell(n_subjects, n_parameters);

parfor i = 1:n_subjects
    fprintf('Processing subject %g\n', i);
    
    % Load LIMO structure once per subject
    LIMO = load(fullfile(Paths{i}, 'LIMO.mat'));
    LIMO = LIMO.LIMO;
    
    % Load Yr.mat once per subject
    Yr = load(fullfile(Paths{i}, 'Yr.mat'));
    Yr = Yr.(cell2mat(fieldnames(Yr)));
    
    % Process all requested parameters for this subject
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

% Calculate grand averages for each parameter
grand_averages = cell(n_parameters, 1);
for p_idx = 1:n_parameters
    if ~isempty(data{p_idx})
        grand_averages{p_idx} = mean(data{p_idx}, 2, 'omitnan');
    end
end


% Create plot following limo_add_plots style
figure('Name', 'ERP Plot', 'color', 'w');
hold on;

% Define colors for different parameters
colors = lines(n_parameters);
legend_labels = cell(n_parameters, 1);

% Plot individual subjects in light gray (if requested)
if show_individual
    for p_idx = 1:n_parameters
        if ~isempty(data{p_idx})
            for i = 1:n_subjects
                if ~all(isnan(data{p_idx}(:, i)))
                    plot(time_vector, data{p_idx}(:, i), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.25);
                end
            end
        end
    end
end

% Plot grand averages for each parameter
for p_idx = 1:n_parameters
    if ~isempty(grand_averages{p_idx})
        plot(time_vector, grand_averages{p_idx}, 'Color', colors(p_idx, :), 'LineWidth', 3);
        legend_labels{p_idx} = sprintf('Parameter %d', parameter(p_idx));
    end
end

% Format plot following limo_add_plots style
grid on; 
axis tight; 
box on;
xlabel('Time (ms)', 'FontSize', 14);
ylabel('Amplitude', 'FontSize', 14);

% Set title
if isempty(plot_title)
    if n_parameters == 1
        plot_title = sprintf('Electrode %d - Parameter %d (Mean estimator)', electrode, parameter);
    else
        plot_title = sprintf('Electrode %d - Parameters %s (Mean estimator)', electrode, mat2str(parameter));
    end
end
title(plot_title, 'FontSize', 16, 'Interpreter', 'none');

% Add legend
if show_individual
    legend(['Individual Subjects', legend_labels'], 'Location', 'best');
else
    legend(legend_labels, 'Location', 'best');
end

hold off;


end