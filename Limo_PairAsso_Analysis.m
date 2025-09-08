% TODO
% See if weighted means have different values
% Are weights justified
% 

%% EEGLAB LIMO light mode
settings.matlab.appearance.figure.GraphicsTheme.PersonalValue = "light";

%% Create a plotting file for each condition

% 1. trial_type = Study_hits 
% 2. trial_type = Study_misses 
% 3. trial_type = Test_hits 
% 4. trial_type = Test_misses
% 5. trial_type = Correct_rejections 
% 6. trial_type = False_alarms

chan_loc = '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/limo_gp_level_chanlocs.mat';
files = '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/LIMO_list.txt';
parameters = 2;
estimator = 'Mean'; % 'Mean' 'Weighted mean' 'Trimmed mean' 'HD' 'Median'
analysis_type = 'Mean'; % 'Mean' 'Trimmed mean' 'HD' 'Median'
savename = [ '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/parameter_' num2str(parameters) '.mat'];
limo_central_tendency_and_ci(files, parameters, chan_loc, estimator, analysis_type, [],savename);

%% Plot averages for each condition

average_files = {'/Volumes/T7/ERP Files/Epoched Files 50/derivatives/parameter_5_Mean_of_Mean.mat'};
% average_files = {'/Volumes/T7/ERP Files/Original Epoched Files/derivatives/parameter_5_Mean_of_Mean.mat', '/Volumes/T7/ERP Files/Original Epoched Files/derivatives/parameter_6_Mean_of_Mean.mat'};
channel = [21];

% Check if average_files is empty and open file picker if needed
if isempty(average_files)
    [filename, pathname] = uigetfile('*.mat', 'Select average files', 'MultiSelect', 'on');
    if isequal(filename, 0) || isequal(pathname, 0)
        disp('User cancelled file selection');
        return;
    end
    if iscell(filename)
        average_files = cellfun(@(x) fullfile(pathname, x), filename, 'UniformOutput', false);
    else
        average_files = {fullfile(pathname, filename)};
    end
end

% Check if channel is empty and open channel picker if needed
if isempty(channel)
    % Load channel locations to get list of available channels
    chanlocs = load(chan_loc);
    channel_names = {chanlocs.expected_chanlocs.labels};
    [selection, ok] = listdlg('PromptString', 'Select a channel:', ...
                              'SelectionMode', 'single', ...
                              'ListString', channel_names);
    if ok
        channel = selection;
    else
        disp('User cancelled channel selection');
        return;
    end
end

limo_add_plots(average_files, [], 'channel', channel, 'restrict', 'Time', 'figure', 'new');

% DEBUG: Let's examine what's actually in the parameter files
fprintf('\n=== DEBUGGING PARAMETER FILES VS RAW DATA ===\n');

% 1. Load parameter file data
fprintf('\n1. PARAMETER FILE DATA:\n');
for i = 1:length(average_files)
    fprintf('\nFile %d: %s\n', i, average_files{i});
    
    % Load and examine the parameter file
    param_data = load(average_files{i});
    field_names = fieldnames(param_data);
    fprintf('  Fields: %s\n', strjoin(field_names, ', '));
    
    % Get the main data field (usually the largest array)
    main_field = '';
    for j = 1:length(field_names)
        field_value = param_data.(field_names{j});
        if isnumeric(field_value)
            fprintf('  %s size: %s\n', field_names{j}, mat2str(size(field_value)));
            if numel(field_value) > 100  % Likely the main data
                main_field = field_names{j};
            end
        elseif isstruct(field_value)
            fprintf('  %s: [struct with fields: %s]\n', field_names{j}, strjoin(fieldnames(field_value), ', '));
            % Check if it's the Data structure
            if strcmp(field_names{j}, 'Data')
                data_fields = fieldnames(field_value);
                for k = 1:length(data_fields)
                    data_field_value = field_value.(data_fields{k});
                    if isnumeric(data_field_value)
                        fprintf('    Data.%s size: %s\n', data_fields{k}, mat2str(size(data_field_value)));
                        if contains(data_fields{k}, 'mean') || contains(data_fields{k}, 'data')
                            main_field = sprintf('Data.%s', data_fields{k});
                        end
                    end
                end
            end
        end
    end
    
    % Examine the main data for this channel
    if ~isempty(main_field)
        if contains(main_field, 'Data.')
            parts = split(main_field, '.');
            main_data = param_data.(parts{1}).(parts{2});
        else
            main_data = param_data.(main_field);
        end
        
        if ndims(main_data) >= 2
            if ndims(main_data) == 4 && size(main_data, 4) == 3
                % Format: [channels x time x conditions x stats]
                % Take middle statistic (mean) from the 4th dimension
                channel_data = squeeze(main_data(channel, :, 1, 2));
                fprintf('  Channel %d data (mean from 4D): range [%.6f, %.6f]\n', channel, min(channel_data), max(channel_data));
            elseif ndims(main_data) == 3 && size(main_data, 3) == 3
                % Format: [channels x time x stats] 
                % Take middle statistic (mean) from the 3rd dimension
                channel_data = squeeze(main_data(channel, :, 2));
                fprintf('  Channel %d data (mean from 3D): range [%.6f, %.6f]\n', channel, min(channel_data), max(channel_data));
            else
                % Regular 2D format
                channel_data = squeeze(main_data(channel, :));
                fprintf('  Channel %d data (2D): range [%.6f, %.6f]\n', channel, min(channel_data), max(channel_data));
            end
            fprintf('  Channel %d first 10 points: %s\n', channel, mat2str(channel_data(1:min(10,end))', 6));
            
            % Store for comparison
            if i == 1
                param_data_for_comparison = channel_data;
            end
        end
    end
end

% 2. Now compute the same condition using our method for comparison
fprintf('\n\n2. RAW ERP DATA (our method):\n');
% Read the subject list with proper handling of spaces in paths  
% Use the CORRECT file list that matches the parameter file dataset
fid = fopen('/Volumes/T7/ERP Files/Epoched Files 50/derivatives/LIMO_list_full.txt', 'r');
if fid == -1
    fprintf('  Could not open LIMO file list - trying alternative path\n');
    % If the file list doesn't exist, construct paths manually for a few subjects
    file_list = {
        '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/sub-1009/PairAsso_PairAsso_GLM_Channels_Time_WLS/LIMO.mat'
        '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/sub-1010/PairAsso_PairAsso_GLM_Channels_Time_WLS/LIMO.mat'
        '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/sub-1012/PairAsso_PairAsso_GLM_Channels_Time_WLS/LIMO.mat'
        '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/sub-1013/PairAsso_PairAsso_GLM_Channels_Time_WLS/LIMO.mat'
        '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/sub-1016/PairAsso_PairAsso_GLM_Channels_Time_WLS/LIMO.mat'
    };
    fprintf('  Using manually constructed file list with %d subjects\n', length(file_list));
else
    file_list = {};
    line_count = 0;
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(line)
            line_count = line_count + 1;
            file_list{line_count} = strtrim(line);  % Remove leading/trailing whitespace
        end
    end
    fclose(fid);
    fprintf('  Found %d subjects in file list\n', length(file_list));
end

if length(file_list) > 0
    fprintf('  First file path: "%s"\n', file_list{1});
end
    
    % Process first 5 subjects to get raw ERP data for parameter 2
    raw_erp_data = [];
    valid_subjects = 0;
    
    for subj = 1:min(5, length(file_list))
        fprintf('  Trying subject %d: %s\n', subj, file_list{subj});
        try
            % Check if file exists
            if ~exist(file_list{subj}, 'file')
                fprintf('    File does not exist\n');
                continue;
            end
            
            % Load LIMO file
            [path, ~, ~] = fileparts(file_list{subj});
            LIMO = load(file_list{subj});
            LIMO = LIMO.LIMO;
            
            % Debug time info for first subject
            if subj == 1
                fprintf('    LIMO time info: start=%.1f, end=%.1f, SR=%.1f Hz\n', ...
                    LIMO.data.start, LIMO.data.end, LIMO.data.sampling_rate);
                if isfield(LIMO.data, 'timevect')
                    fprintf('    Time vector length: %d, range: %.1f to %.1f ms\n', ...
                        length(LIMO.data.timevect), LIMO.data.timevect(1), LIMO.data.timevect(end));
                else
                    expected_length = round((LIMO.data.end - LIMO.data.start) * LIMO.data.sampling_rate / 1000) + 1;
                    fprintf('    Expected time vector length: %d\n', expected_length);
                end
            end
            
            % Load Yr file
            yr_path = fullfile(path, 'Yr.mat');
            if ~exist(yr_path, 'file')
                fprintf('    Yr.mat does not exist at: %s\n', yr_path);
                continue;
            end
            
            Yr = load(yr_path);
            Yr = Yr.(cell2mat(fieldnames(Yr)));
            
            % Find trials for parameter 2
            trial_indices = logical(LIMO.design.X(:, 2) == 1);
            
            if sum(trial_indices) > 0
                % Extract electrode data and average across trials
                electrode_data = squeeze(Yr(channel, :, trial_indices));
                if ndims(electrode_data) == 1
                    subject_mean = electrode_data;
                else
                    subject_mean = nanmean(electrode_data, 2);
                end
                
                if valid_subjects == 0
                    raw_erp_data = zeros(length(subject_mean), min(5, length(file_list)));
                end
                
                valid_subjects = valid_subjects + 1;
                raw_erp_data(:, valid_subjects) = subject_mean;
                
                fprintf('    SUCCESS: %d trials, range [%.6f, %.6f]\n', sum(trial_indices), min(subject_mean), max(subject_mean));
            else
                fprintf('    No trials found for parameter 2\n');
            end
            
        catch ME
            fprintf('    Error: %s\n', ME.message);
        end
    end
    
    if valid_subjects > 0
        % Compute grand average
        raw_grand_avg = mean(raw_erp_data(:, 1:valid_subjects), 2);
        fprintf('  \nRaw ERP Grand Average (first %d subjects):\n', valid_subjects);
        fprintf('    Range: [%.6f, %.6f]\n', min(raw_grand_avg), max(raw_grand_avg));
        fprintf('    First 10 points: %s\n', mat2str(raw_grand_avg(1:min(10,end))', 6));
        
        % Compare with parameter file data
        if exist('param_data_for_comparison', 'var')
            fprintf('  \n3. COMPARISON:\n');
            fprintf('    Parameter file length: %d, Raw ERP length: %d\n', length(param_data_for_comparison), length(raw_grand_avg));
            if length(param_data_for_comparison) == length(raw_grand_avg)
                %correlation = corr(param_data_for_comparison, raw_grand_avg);
                difference = mean(abs(param_data_for_comparison - raw_grand_avg));
                ratio = mean(param_data_for_comparison) / mean(raw_grand_avg);
                %fprintf('    Correlation: %.6f\n', correlation);
                fprintf('    Mean absolute difference: %.6f\n', difference);
                fprintf('    Amplitude ratio (param/raw): %.6f\n', ratio);
                if ratio < 1
                    fprintf('    Parameter file amplitude: %.2fx smaller than raw ERP\n', 1/ratio);
                else
                    fprintf('    Parameter file amplitude: %.2fx larger than raw ERP\n', ratio);
                end
            else
                fprintf('    CANNOT COMPARE: Different time lengths\n');
                fprintf('    This suggests parameter file was created from different data:\n');
                fprintf('      - Different time window (epoch length)\n');
                fprintf('      - Different sampling rate\n');
                fprintf('      - Different file list/dataset\n');
            end
        end
    else
        fprintf('  No valid subjects processed\n');
    end



% Extract parameter numbers from filenames for legend and title
legend_labels = {};
param_numbers = [];
for i = 1:length(average_files)
    [~, filename, ~] = fileparts(average_files{i});
    % Extract parameter number from filename (assumes format: parameter_X_...)
    param_match = regexp(filename, 'parameter_(\d+)', 'tokens');
    if ~isempty(param_match)
        param_num = str2double(param_match{1}{1});
        param_numbers(end+1) = param_num;
        legend_labels{end+1} = sprintf('Condition %d', param_num);
    else
        legend_labels{end+1} = sprintf('File %d', i);
    end
end

% Create dynamic title based on conditions
if length(param_numbers) > 1
    condition_str = sprintf('Conditions %s', strjoin(arrayfun(@num2str, param_numbers, 'UniformOutput', false), ', '));
else
    condition_str = sprintf('Condition %d', param_numbers(1));
end

% Update plot title and subtitle
title_str = sprintf('Channel %d ERP, %s', channel, condition_str);
subtitle_str = sprintf('(Within-subjects %s, Between-subjects %s, 95%% high density interval)', lower(estimator), lower(analysis_type));

% Apply the custom title and subtitle
title(title_str);
subtitle(subtitle_str);

% IMPROVED: More robust line handling approach
% Instead of trying to detect which lines are "main" lines, 
% we'll work with what limo_add_plots actually creates
if length(average_files) > 1
    % Get all line objects
    h_lines = findobj(gca, 'Type', 'line');
    
    if length(h_lines) >= length(average_files)
        % Get default color order
        colorOrder = get(gca, 'ColorOrder');
        
        % Debug output
        fprintf('Found %d line objects for %d conditions\n', length(h_lines), length(average_files));
        
        % Strategy: Modify the most recently plotted lines (typically the main ERP lines)
        % limo_add_plots usually plots main lines last, so they appear first in findobj
        legend_handles = [];
        
        for i = 1:min(length(average_files), length(h_lines))
            line_handle = h_lines(i);
            color_idx = mod(i-1, size(colorOrder, 1)) + 1;
            
            % Set distinctive colors and line width for each condition
            set(line_handle, 'Color', colorOrder(color_idx,:), 'LineWidth', 3);
            legend_handles(end+1) = line_handle;
        end
        
        % Create legend if we have handles
        if ~isempty(legend_handles)
            % Reverse the order to match the typical plotting order
            legend(fliplr(legend_handles), fliplr(legend_labels), 'Location', 'best');
        end
    else
        warning('Fewer line objects (%d) than expected conditions (%d)', length(h_lines), length(average_files));
    end
else
    % Single condition - just update the title, no need for complex legend
    fprintf('Single condition plot - no legend modification needed\n');
end

%% Plot single-subject data for each condition
single_subject_files = ['/Volumes/T7/ERP Files/Epoched Files 50/derivatives/parameter_2_single_subjects_Mean.mat'];
channel = [21];

% Check if single_subject_files is empty and open file picker if needed
if isempty(single_subject_files)
    [filename, pathname] = uigetfile('*.mat', 'Select single-subject files', 'MultiSelect', 'on');
    if isequal(filename, 0) || isequal(pathname, 0)
        disp('User cancelled file selection');
        return;
    end
    if iscell(filename)
        single_subject_files = cellfun(@(x) fullfile(pathname, x), filename, 'UniformOutput', false);
    else
        single_subject_files = {fullfile(pathname, filename)};
    end
end

% Check if channel is empty and open channel picker if needed
if isempty(channel)
    % Load channel locations to get list of available channels
    chanlocs = load(chan_loc);
    channel_names = {chanlocs.expected_chanlocs.labels};
    [selection, ok] = listdlg('PromptString', 'Select a channel:', ...
                              'SelectionMode', 'single', ...
                              'ListString', channel_names);
    if ok
        channel = selection;
    else
        disp('User cancelled channel selection');
        return;
    end
end

limo_add_plots({single_subject_files}, [], 'channel', channel, 'restrict', 'Time', 'figure', 'new');

% Extract parameter number from single_subject_files filename for title
[~, filename, ~] = fileparts(single_subject_files);
param_match = regexp(filename, 'parameter_(\d+)', 'tokens');
if ~isempty(param_match)
    param_num = str2double(param_match{1}{1});
    condition_str = sprintf('Condition %d', param_num);
else
    condition_str = 'Unknown Condition';
end

% Update plot title and subtitle for single-subject data
title_str = sprintf('Channel %d Single-subject ERPs, %s', channel, condition_str);
subtitle_str = sprintf('(Within-subjects %s)', lower(estimator));

% Apply the custom title and subtitle
title(title_str);
subtitle(subtitle_str);

%% Plot the difference between two conditions

% Define the parameters
%data1 = '/Volumes/T7/ERP Files/Epoched Files/derivatives/parameter_1_single_subjects_Weighted mean.mat';
%data2 = '/Volumes/T7/ERP Files/Epoched Files/derivatives/parameter_4_single_subjects_Weighted mean.mat';

data1 = '/Volumes/T7/ERP Files/Epoched Files/derivatives/parameter_1_Mean_of_Mean.mat';
data2 = '/Volumes/T7/ERP Files/Epoched Files/derivatives/parameter_4_Mean_of_Mean.mat';
type = 'independent';
percent = 'mean';
alpha = 5;
channel = 21;
output_filename = fullfile(pwd, 'difference_results.mat');

% Run the limo_plot_difference function
Data = limo_plot_difference(data1, data2, 'type', type, 'percent', percent, 'alpha', alpha/100, 'channel', channel, 'name', output_filename);

% Display the results
disp('Function executed successfully!');
disp('Output structure fields:');
disp(fieldnames(Data));

%% Specify data for plots

load('/Volumes/T7/ERP Files/Epoched Files 50/paired_ttest2/LIMO.mat');
data = load('/Volumes/T7/ERP Files/Epoched Files 50/paired_ttest2/paired_samples_ttest_parameter_34.mat');
one_sample = data.paired_samples;

plot_data = squeeze(one_sample(:,:,4)); % t-values (4th dimension)
%p_values = squeeze(one_sample(:,:,5)); % p-values (5th dimension)
p_values = squeeze(LIMO.cache.fig.pval); % p-values (5th dimension)
df_values = squeeze(one_sample(:,:,3)); % df-values (3rd dimension)

%% Plot the electodes versus time 2D and 3D graphs

% The function has 3 parameters pointing to the test directory, so it's
% easier to go to that directory first and use pwd
cd '/Volumes/T7/ERP Files/Epoched Files 50/paired_ttest2'
limo_display_results(1, ... # electodes x time graph
    'paired_samples_ttest_parameter_34.mat', ... # test statistics mat file
    pwd, ... # path to the test directory
    0.05, ... # significance level
    3, ... # MCC: 1 = none, 2 = cluster, 3 = TFCE, 4 = max
    fullfile(pwd,'LIMO.mat'), ... # LIMO mat file
    0 ... # interactive figure: 0 = false, 1 = true
    );

%% Plot topoplot graphs (uses pop_topoplot)

EEG = [];
EEG.data = plot_data;  % t-values, could be replaced with voltages
EEG.setname = 'Paired T-test Results';
EEG.nbchan = size(EEG.data,1);
EEG.pnts = size(EEG.data,2);
EEG.trials = 1;
EEG.chanlocs = LIMO.data.chanlocs;
EEG.xmin = LIMO.data.start/1000;  % Convert to seconds
EEG.xmax = LIMO.data.end/1000;    % Convert to seconds
EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);

latencies = 0:100:1500; % time plots to display

pop_topoplot(EEG, 1, latencies, 'Paired T-test Topoplots', 0);

%% Custom topoplot graphs with modified titles, size, and colorbar (uses topoplot)

% LIMO-style masking
plot_data_masked = plot_data;
plot_data_masked(p_values >= 0.05) = NaN;

% LIMO-style colormap
cc = limo_color_images(plot_data(~isnan(plot_data)));

% Time setup
times = linspace(LIMO.data.start, LIMO.data.end, size(plot_data,2));
latencies = 0:100:1500;

% Find time indices for each latency
time_indices = arrayfun(@(x) find(abs(times - x) == min(abs(times - x)), 1), latencies);

% Create figure with subplots for 2D topoplots (much taller to give title space)
fig = figure('Position', [100 50 1200 1200]); 
ax_handles = [];  % Store axis handles

% Calculate symmetric color limits like the 3D plots
abs_max = max(abs(plot_data(:)));
colorbar_limits = [-abs_max, abs_max];

for i = 1:length(latencies)
    ax = subplot(4, 4, i);
    ax_handles(end+1) = ax;
    
    % Get data for this time point
    data_timepoint = plot_data(:, time_indices(i));
    
    % Create 2D topoplot using direct topoplot() call
    topoplot(data_timepoint, LIMO.data.chanlocs, ... % DO NOT APPLY THE MASK TO TOPOPLOTS
        'maplimits', colorbar_limits, ...
        'electrodes', 'off', ...
        'style', 'both', ...
        'shading', 'interp', ...
        'numcontour', 6);
    
    % Add small latency title above each plot
    title(sprintf('%d ms', latencies(i)), 'FontSize', 8, 'Color', 'k');
    
    % Simple 10% scaling without repositioning to avoid breaking handles
    pos = get(ax, 'Position');
    center_x = pos(1) + pos(3)/2;
    center_y = pos(2) + pos(4)/2;
    new_width = pos(3) * 1.1;
    new_height = pos(4) * 1.1;
    new_x = center_x - new_width/2;
    new_y = center_y - new_height/2;
    set(ax, 'Position', [new_x, new_y, new_width, new_height]);
end

% Set colormap to match 3D plots
colormap(fig, cc);

% Set color axis limits for all subplots to match 3D plotting style
for i = 1:length(ax_handles)
    set(ax_handles(i), 'CLim', colorbar_limits);
end

% Create invisible reference image for colorbar with proper gradient (same as 3D)
ax_ref = axes('Position', [0.001 0.001 0.001 0.001], 'Visible', 'off');
imagesc(ax_ref, linspace(colorbar_limits(1), colorbar_limits(2), 100)');
set(ax_ref, 'CLim', colorbar_limits);
colormap(ax_ref, cc);

% Create colorbar linked to the reference image with 3D-style formatting
h = colorbar(ax_ref, 'Position', [0.92 0.1 0.03 0.8]);  % [x y width height] - larger like 3D
% Set 5 tick marks like the 3D plots: min, mid-low, 0, mid-high, max
tick_values = [colorbar_limits(1), colorbar_limits(1)/2, 0, colorbar_limits(2)/2, colorbar_limits(2)];
tick_labels = arrayfun(@(x) sprintf('%.1f', x), tick_values, 'UniformOutput', false);
set(h, 'Ticks', tick_values, 'TickLabels', tick_labels, 'Color', 'k', 'Box', 'on');
title(h, 't-values', 'FontSize', 10, 'Color', 'k');

% Make sure the reference axes is completely hidden
set(ax_ref, 'XTick', [], 'YTick', [], 'Box', 'off');

% Add main title matching 3D plotting style with more space
sgtitle([LIMO.cache.fig.title, ' Topoplots'], 'FontSize', 14, 'Color', 'k');

% Adjust subplot positions to move everything down and create title space
for i = 1:length(ax_handles)
    pos = get(ax_handles(i), 'Position');
    set(ax_handles(i), 'Position', [pos(1), pos(2)*0.925, pos(3), pos(4)]);
end

%% Plot headplot (3D topoplot) graphs (uses headplot)

% headplot('setup', LIMO.data.chanlocs, 'paired_ttest_headplot.spl');

% LIMO-style masking
plot_data_masked = plot_data;
plot_data_masked(p_values >= 0.05) = NaN;

% LIMO-style colormap
cc = limo_color_images(plot_data_masked(~isnan(plot_data_masked)));

% Time setup
times = linspace(LIMO.data.start, LIMO.data.end, size(plot_data,2));
latencies = 0:100:1500;

% Find time indices for each latency
time_indices = arrayfun(@(x) find(abs(times - x) == min(abs(times - x)), 1), latencies);

% Create figure with subplots for 3D headplots (much taller to give title space)
fig = figure('Position', [100 50 1200 1200], 'Color', 'w', 'InvertHardcopy', 'off'); 
ax_handles = [];  % Store axis handles for synchronized rotation

% Use same logic as topoplot: symmetric limits based on absolute maximum
abs_max = max(abs(plot_data(:)));  % Get absolute maximum like 'absmax' does
colorbar_limits = [-abs_max, abs_max];  % Symmetric limits like topoplot

for i = 1:length(latencies)
    ax = subplot(4, 4, i);
    ax_handles(end+1) = ax;  % Store axis handle

    % Get data for this time point
    data_timepoint = plot_data(:, time_indices(i));
    data_timepoint(isnan(data_timepoint)) = 0; % Replace NaN with 0 for headplot

    % Create 3D headplot
    headplot(data_timepoint, 'paired_ttest_headplot.spl', ...
        'electrodes', 'off', ...
        'maplimits', colorbar_limits, ...
        'view', [143 18], ...
        'verbose', 'off');
    
    % Add small latency title above each head
    title(sprintf('%d ms', latencies(i)), 'FontSize', 8, 'Color', 'k');
        
    % Simple 10% scaling without repositioning to avoid breaking handles
    pos = get(ax, 'Position');
    center_x = pos(1) + pos(3)/2;
    center_y = pos(2) + pos(4)/2;
    new_width = pos(3) * 1.1;
    new_height = pos(4) * 1.1;
    new_x = center_x - new_width/2;
    new_y = center_y - new_height/2;
    set(ax, 'Position', [new_x, new_y, new_width, new_height]);
end

% Set colormap
colormap(fig, cc);

% Set color axis limits for all subplots to match topoplot style
for i = 1:length(ax_handles)
    set(ax_handles(i), 'CLim', colorbar_limits);
end

% Create invisible reference image for colorbar with proper gradient
ax_ref = axes('Position', [0.001 0.001 0.001 0.001], 'Visible', 'off');  % Tiny invisible axes off-screen
imagesc(ax_ref, linspace(colorbar_limits(1), colorbar_limits(2), 100)');  % Gradient data
set(ax_ref, 'CLim', colorbar_limits);
colormap(ax_ref, cc);  % Use same colormap as topoplots

% Create colorbar linked to the reference image with topoplot-style ticks
h = colorbar(ax_ref, 'Position', [0.92 0.1 0.03 0.8]);  % [x y width height]
% Set 5 tick marks like topoplot: min, mid-low, 0, mid-high, max
tick_values = [colorbar_limits(1), colorbar_limits(1)/2, 0, colorbar_limits(2)/2, colorbar_limits(2)];
tick_labels = arrayfun(@(x) sprintf('%.1f', x), tick_values, 'UniformOutput', false);
set(h, 'Ticks', tick_values, 'TickLabels', tick_labels, 'Color', 'k', 'Box', 'on');
title(h, 't-values', 'FontSize', 10, 'Color', 'k');  % Label the colorbar

% Make sure the reference axes is completely hidden
set(ax_ref, 'XTick', [], 'YTick', [], 'Box', 'off');
sgtitle([LIMO.cache.fig.title, ' Headplots'], 'FontSize', 14, 'Color', 'k');

% Adjust subplot positions to move everything down and create title space
for i = 1:length(ax_handles)
    pos = get(ax_handles(i), 'Position');
    set(ax_handles(i), 'Position', [pos(1), pos(2)*0.925, pos(3), pos(4)]);
end

% Link camera properties for synchronized rotation
hlink = linkprop(ax_handles, {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'});
setappdata(fig, 'rotation_link', hlink); % Store the link to prevent it from being deleted

% Enable figure-level rotation mode
rotate3d(fig, 'on');

%% Time-Channel Plot with Topographic Plots

% Define rectangles for highlighting. Each cell contains a 2x2 matrix for
% a rectangle, defined by its top-left and bottom-right corners in the
% format: [time_ms, channel; time_ms, channel].
highlight_rects = {[400, 21; 800, 21], [400, 87; 800, 87], [400, 101; 800, 101], [400, 153; 800, 153]}; 
% Example: [400, 21; 700, 21], [400, 36; 700, 36], [400, 101; 700, 101], [-100, 224; 1500, 224]

% Define labels for the highlighted rectangles. The number of labels should
% match the number of rectangles.
highlight_labels = {'E21', 'E87', 'E101', 'E153'};
% Example: 'Fz (E21)', 'F3 (E36)', 'F4 (E224)', 'Pz (E101)'

% LIMO-style masking
plot_data_masked = plot_data;
plot_data_masked(p_values >= 0.05) = NaN;

% LIMO-style colormap
cc = limo_color_images(plot_data_masked(~isnan(plot_data_masked)));

% Create time-channel plot data
% Placeholder values from tutorial: [400 8; 350 14; 500 24; 1050 11]
% Format: [time_ms, channel]
time_channel_points_original = [];

% Check if topoplots should be plotted
plot_topoplots = ~isempty(time_channel_points_original);

if plot_topoplots
    % Sort points by time (x-axis) for proper topoplot ordering
    [~, sort_idx] = sort(time_channel_points_original(:, 1));
    time_channel_points = time_channel_points_original(sort_idx, :);
    
    % Define the number of points to plot
    num_tf_points = size(time_channel_points, 1);
else
    time_channel_points = [];
    num_tf_points = 0;
end

% Time setup
times = linspace(LIMO.data.start, LIMO.data.end, size(plot_data,2));

% Create figure for time-channel plot with topoplots above
fig_exp = figure('Position', [100 100 1000 700], 'Color', 'w', 'InvertHardcopy', 'off');

% Apply LIMO colormap to this figure
colormap(fig_exp, cc);

% Create main time-channel plot - adjust position based on whether topoplots are present
if plot_topoplots
    imgax_exp = axes('Position', [0.1 0.08 0.75 0.55]); % Lower position with space for topoplots above
else
    imgax_exp = axes('Position', [0.1 0.15 0.75 0.65]); % Higher position when no topoplots
end
set(imgax_exp, 'Color', [0.9 0.9 0.9], 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1); % Set axes background and and outline to black with thicker border

% Plot the masked data
h_img_exp = imagesc(times, 1:size(plot_data,1), plot_data_masked);

% Set NaN values to be transparent by creating an alpha mask
alpha_mask_exp = ~isnan(plot_data_masked);
set(h_img_exp, 'AlphaData', alpha_mask_exp);

% Set background color to show through transparent areas
set(imgax_exp, 'Color', [0.9 0.9 0.9]);

% Set color limits
abs_max = max(abs(plot_data(:)));
caxis([-abs_max, abs_max]);

% Add labels
xlabel('Time (ms)', 'Color', 'k');
ylabel('Channel', 'Color', 'k');

% Flip y-axis so channel 1 is at the top
set(imgax_exp, 'YDir', 'reverse');
ylim([1 size(plot_data,1)]);
grid(imgax_exp, 'off'); % Turn off all grid lines for the specific axes

% Set tick label colors to black
set(imgax_exp, 'XColor', 'k', 'YColor', 'k');
hold on;

% Draw highlighted rectangles and labels if any are defined
if ~isempty(highlight_rects)
    % Get properties for label positioning and sizing
    axis_font_size = get(imgax_exp, 'FontSize');
    time_range = times(end) - times(1);
    label_offset = 0.004 * time_range; % 2% offset to the left

    % Get time step for pixel adjustment, assuming uniform spacing
    if length(times) > 1
        time_step = times(2) - times(1);
    else
        time_step = 0;
    end

    for i = 1:length(highlight_rects)
        coords = highlight_rects{i};
        time1 = coords(1,1);
        chan1 = coords(1,2);
        time2 = coords(2,1);
        chan2 = coords(2,2);
        
        % Adjust coordinates to cover the full pixels.
        % For imagesc, pixels are centered on the coordinate points.
        % We need to extend the rectangle by half a pixel in each direction.
        
        % Find the time indices corresponding to the time values
        [~, time_idx1] = min(abs(times - time1));
        [~, time_idx2] = min(abs(times - time2));

        % Get the exact time values from the times vector
        exact_time1 = times(time_idx1);
        exact_time2 = times(time_idx2);

        x = min(exact_time1, exact_time2) - time_step/2;
        y = min(chan1, chan2) - 0.5;
        width = abs(exact_time2 - exact_time1) + time_step;
        height = abs(chan2 - chan1) + 1;
        
        rectangle('Position', [x, y, width, height], ...
                  'FaceColor', [0 1 0], ... % Green
                  'FaceAlpha', 0.5, ...      % Set alpha value for transparency
                  'EdgeColor', 'none');

        % Add a label to the left of the rectangle
        if i <= length(highlight_labels) && ~isempty(highlight_labels{i})
            label_text = highlight_labels{i};
            % Position the label vertically in the middle of the rectangle, with a slight upward adjustment
            label_y = y + height/2 - 0.25;
            
            % Position the label horizontally to the left of the plot
            label_x = times(1) - label_offset;
            
            text(label_x, label_y, label_text, ...
                 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'middle', ...
                 'Color', 'k', ...
                 'FontSize', axis_font_size, ...
                 'Clipping', 'off'); % Prevent label from being clipped
        end
    end
end

% Add vertical line at time=0
ylim_vals = get(imgax_exp, 'YLim');
plot(imgax_exp, [0 0], ylim_vals, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1);


if plot_topoplots
    % Pre-calculate topoplot positions for line drawing
    topo_positions = zeros(num_tf_points, 4); % [x, y, width, height] for each topoplot
    for n = 1:num_tf_points
        % Position topoplots evenly across the top of the figure
        % Space them evenly from left to right, regardless of their time points
        topo_x = 0.1 + (n-1) * (0.75 / num_tf_points) + (0.75 / num_tf_points - 0.15) / 2; % Adjusted for new topo_size
        topo_y = 0.69; % Position above main plot with more space from title
        topo_size = 0.15; % Size of each topoplot
        
        topo_positions(n, :) = [topo_x, topo_y, topo_size, topo_size];
    end

    % Store all data for drawing lines and markers in proper order
    marker_data_exp = zeros(num_tf_points, 2);
    line_coords_exp = zeros(num_tf_points, 4);

    % Calculate line coordinates using pre-calculated topoplot positions
    for n = 1:num_tf_points
        time_point = time_channel_points(n, 1);
        channel_point = time_channel_points(n, 2);
        
        % Store marker data for later
        marker_data_exp(n,:) = [time_point, channel_point];
        
        % Get coordinates for the line in figure units
        main_pos = get(imgax_exp, 'Position');
        ylims = get(imgax_exp, 'YLim');
        time_norm = (time_point - min(times)) / (max(times) - min(times));
        channel_norm = (channel_point - ylims(1)) / (ylims(2) - ylims(1));
        from_x = main_pos(1) + time_norm * main_pos(3);
        from_y = main_pos(2) + (1 - channel_norm) * main_pos(4); % Corrected for reversed y-axis
        
        topo_pos = topo_positions(n, :); % Use pre-calculated position
        to_x = topo_pos(1) + topo_pos(3) * 0.5;
        to_y = topo_pos(2) + 0.009;
        
        line_coords_exp(n,:) = [from_x, to_x, from_y, to_y];
    end
end

% Set focus to the main plot
axes(imgax_exp);

if plot_topoplots
    % Draw connecting lines on the main plot by converting figure to data coordinates
    main_pos = get(imgax_exp, 'Position');
    xlims = get(imgax_exp, 'XLim');
    ylims = get(imgax_exp, 'YLim');

    for n = 1:size(line_coords_exp, 1)
        % Line endpoints in figure units
        from_x_fig = line_coords_exp(n, 1);
        to_x_fig = line_coords_exp(n, 2);
        from_y_fig = line_coords_exp(n, 3);
        to_y_fig = line_coords_exp(n, 4);

        % Vector from start to end in figure units
        dx_fig = to_x_fig - from_x_fig;
        dy_fig = to_y_fig - from_y_fig;
        len_fig = sqrt(dx_fig^2 + dy_fig^2);

        % Normalize the vector
        udx_fig = dx_fig / len_fig;
        udy_fig = dy_fig / len_fig;

        % Shorten by a fixed amount in figure units. This value may need tuning.
        shorten_amount = 0.0057; 
        from_x_fig_new = from_x_fig + shorten_amount * udx_fig;
        from_y_fig_new = from_y_fig + shorten_amount * udy_fig;

        % Convert new start point from figure units back to data units
        time_norm_start_new = (from_x_fig_new - main_pos(1)) / main_pos(3);
        time_start_new = time_norm_start_new * (xlims(2) - xlims(1)) + xlims(1);
        channel_norm_start_new = 1 - ((from_y_fig_new - main_pos(2)) / main_pos(4));
        channel_start_new = channel_norm_start_new * (ylims(2) - ylims(1)) + ylims(1);

        % Convert end point from figure units back to data units
        time_norm_end = (to_x_fig - main_pos(1)) / main_pos(3);
        time_end = time_norm_end * (xlims(2) - xlims(1)) + xlims(1);
        channel_norm_end = 1 - ((to_y_fig - main_pos(2)) / main_pos(4));
        channel_end = channel_norm_end * (ylims(2) - ylims(1)) + ylims(1);

        % Draw connecting line in black
        plot([time_start_new, time_end], [channel_start_new, channel_end], 'k', 'LineWidth', 2, 'Clipping', 'off');
    end

    % Then draw all markers in the plot, so they appear on top of the lines
    for n = 1:num_tf_points
        time_point = marker_data_exp(n, 1);
        channel_point = marker_data_exp(n, 2);
        
        % Add marker at the time-channel point with black outline
        plot(time_point, channel_point, 'ko', 'MarkerSize', 11, ...
            'MarkerFaceColor', 'none', 'LineWidth', 2);
    end

    % Store topoplot axes handles
    topoaxes_exp = zeros(1, num_tf_points);

    % Now create topoplots evenly spaced across the top (like tftopo) - they will appear above lines
    for n = 1:num_tf_points
        time_point = time_channel_points(n, 1);
        channel_point = time_channel_points(n, 2);
        
        % Use pre-calculated topoplot position
        topo_pos = topo_positions(n, :);
        
        % Create topoplot axes
        topoaxes_exp(n) = axes('Position', topo_pos);
        
        % Get the time point and find closest time index
        [~, time_idx] = min(abs(times - time_point));
        
        % Get data for this time point (all channels)
        topo_data = plot_data(:, time_idx);
        
        % Get t, p, and df for title
        t_value = plot_data(channel_point, time_idx);
        p_value = p_values(channel_point, time_idx);
        df = df_values(channel_point, time_idx);

        % Create topoplot and then modify white background patches to be transparent
        topoplot(topo_data, LIMO.data.chanlocs, ... % DO NOT APPLY THE MASK TO TOPOPLOTS
            'maplimits', [-abs_max, abs_max], ...
            'colormap', cc, ...
            'electrodes', 'off', ...
            'style', 'both', ...
            'shading', 'interp', ...
            'numcontour', 6, ...
            'hcolor', 'k'); % Keep head outline black
        
        % Mark the specific channel electrode with a solid black dot
        hold on;
        % Create a simple overlay topoplot with just electrode markers
        topoplot([], LIMO.data.chanlocs, ...
            'plotchans', channel_point, ...
            'colormap', cc, ...
            'electrodes', 'on', ...
            'emarker', {'.', 'k', 16, 1}, ...
            'style', 'blank', ...
            'hcolor', 'none');
        hold off;
        
        % Find and make the white background patches transparent
        patch_handles = findobj(gca, 'Type', 'patch');
        for i = 1:length(patch_handles)
            % Check if this patch has the BACKCOLOR (light blue/white)
            face_color = get(patch_handles(i), 'FaceColor');
            if isnumeric(face_color) && length(face_color) == 3
                % Check if it's close to the BACKCOLOR [.93 .96 1]
                if all(abs(face_color - [0.93 0.96 1]) < 0.1)
                    % Make this patch transparent
                    set(patch_handles(i), 'FaceAlpha', 0);
                end
            end
        end
        
        % Add title showing time and channel, and stats
        title_line1 = sprintf('E%d, %d ms', channel_point, time_point);
        if p_value < 0.001
            p_string = 'p < 0.001';
        else
            p_string = sprintf('p = %.3f', p_value);
        end
        title_line2 = sprintf('t(%d) = %.3f, %s', df, t_value, p_string);
        title({title_line1, title_line2}, 'FontSize', 10, 'Color', 'k'); % 8 * 1.2 = 9.6, rounded to 10
        
        % Make topoplots circular
        axis tight;
        axis equal;
        
        % Add circular outline around topoplot with transparent white background
        set(gca, 'XTick', [], 'YTick', [], 'Box', 'off');
        hold on;
        % Draw circular outline (smaller to match color border)
        theta = 0:0.01:2*pi;
        xlims = xlim;
        ylims = ylim;
        radius = min(diff(xlims), diff(ylims))/2 * 0.93; % Make 8% smaller to match color border
        center_x = mean(xlims);
        center_y = mean(ylims);
        circle_x = center_x + radius * cos(theta);
        circle_y = center_y + radius * sin(theta);
        % Draw white background circle first (for typical outline effect)
        h_white = plot(circle_x, circle_y, 'Color', [1 1 1 0], 'LineWidth', 4.5); % Transparent white
        % Draw black outline on top (thinner to match other elements)
        plot(circle_x, circle_y, 'k', 'LineWidth', 1);
        hold off;
    end
end

% Ensure LIMO colormap is applied to the main time-channel axes
colormap(imgax_exp, cc);

% Add colorbar positioned to the right of the main plot - adjust based on plot position
if plot_topoplots
    h = colorbar(imgax_exp, 'Position', [0.87 0.08 0.03 0.55]); % Match lower plot position
else
    h = colorbar(imgax_exp, 'Position', [0.87 0.15 0.03 0.65]); % Match higher plot position
end
title(h, 't-values', 'FontSize', 10, 'Color', 'k');
set(h, 'Color', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', 0); % Set colorbar text and outline to black, remove internal gridlines

% Final formatting: ensure all axes elements are black
set(imgax_exp, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', [0 0]);
% Set custom y-axis ticks. If custom labels are used, only show the
% first and last channel numbers. Otherwise, show 10 equally spaced values.
num_channels = size(plot_data, 1);
has_custom_labels = ~isempty(highlight_labels) && any(cellfun(@(x) ~isempty(x), highlight_labels));

if has_custom_labels
    ytick_values = [1, num_channels];
else
    ytick_values = round(linspace(1, num_channels, 10));
end
set(imgax_exp, 'YTick', ytick_values);

% Set custom x-axis ticks - always include last time point
tick_interval = 250;
min_tick = ceil(times(1) / tick_interval) * tick_interval;
max_tick = floor(times(end) / tick_interval) * tick_interval;

% Create regular ticks
regular_ticks = min_tick:tick_interval:max_tick;

% Always include the actual end time
all_ticks = [regular_ticks, times(end)];

% Remove duplicates and sort
all_ticks = unique(all_ticks);

set(imgax_exp, 'XTick', all_ticks);

% Add main title using annotation - adjust position based on whether topoplots are present
if plot_topoplots
    annotation('textbox', [0, 0.85, 1, 0.15], ...
               'String', [LIMO.cache.fig.title, ' Time-Channel Plot'], ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FontSize', 21, ...
               'Color', 'k', ...
               'FitBoxToText', 'on');
else
    annotation('textbox', [0, 0.85, 1, 0.08], ...
               'String', [LIMO.cache.fig.title, ' Time-Channel Plot'], ...
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FontSize', 21, ...
               'Color', 'k', ...
               'FitBoxToText', 'on');
end

%% New TODO
% Add a feature in my channel-time graph so that you can define rectangles
% that highlight a given range of values, defined by the two corner values.
% These rectangles will have an associated label on the left side of the graph
% Consider plotting TFCE statistics instead of t-values on the channel-time graph
% Look at how topoplots merge values at each electrode point and consider
% changing them to more direct pixel translations like the channel-time graph
% Consider using TFCE stats for topoplots too
% Provide an easy way to update LIMO.cache.fig.title values
% Make a channel-time plot that it is only the plot, and it is shaped such
% that it always has square pixels, so the shape is determined by the
% channels x time lengths
% Reduce filesize of eegmovies

% Verify the six condition numerical codings

% Find a way to apply individual subject correlation analysis to MUA. Other
% papers likely tried this before


%% Experimental
% This section demonstrates an alternative highlighting approach where
% non-selected channels are masked with a semi-transparent black overlay.

% Define rectangles for highlighting. Each cell contains a 2x2 matrix for
% a rectangle, defined by its top-left and bottom-right corners in the
% format: [time_ms, channel; time_ms, channel].
highlight_rects = {[400, 21; 700, 21], [400, 36; 700, 36], [400, 101; 700, 101], [-100, 224; 1500, 224]}; % Example: [200, 20; 400, 40], [600, 80; 800, 100]

% Define labels for the highlighted rectangles. The number of labels should
% match the number of rectangles.
highlight_labels = {'Fz (E21)', 'F3 (E36)', 'F4 (E224)', 'Pz (E101)'};
% 

% LIMO-style masking
plot_data_masked = plot_data;
plot_data_masked(p_values >= 0.05) = NaN;

% LIMO-style colormap
cc = limo_color_images(plot_data_masked(~isnan(plot_data_masked)));

% Create time-channel plot data
% Placeholder values from tutorial: [400 8; 350 14; 500 24; 1050 11]
% Format: [time_ms, channel]
time_channel_points_original = [];

% Check if topoplots should be plotted
plot_topoplots_experimental = ~isempty(time_channel_points_original);

if plot_topoplots_experimental
    % Sort points by time (x-axis) for proper topoplot ordering
    [~, sort_idx] = sort(time_channel_points_original(:, 1));
    time_channel_points = time_channel_points_original(sort_idx, :);
    
    % Define the number of points to plot
    num_tf_points = size(time_channel_points, 1);
else
    time_channel_points = [];
    num_tf_points = 0;
end

% Time setup
times = linspace(LIMO.data.start, LIMO.data.end, size(plot_data,2));

% Create figure for time-channel plot with topoplots above
fig_exp_experimental = figure('Position', [100 100 1000 700], 'Color', 'w', 'InvertHardcopy', 'off');

% Apply LIMO colormap to this figure
colormap(fig_exp_experimental, cc);

% Create main time-channel plot taking up most of the lower area
imgax_exp_experimental = axes('Position', [0.1 0.08 0.75 0.55]); % Adjusted height to give title space
set(imgax_exp_experimental, 'Color', [0.9 0.9 0.9], 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1); % Set axes background and and outline to black with thicker border

% Plot the masked data
h_img_exp_experimental = imagesc(times, 1:size(plot_data,1), plot_data_masked);

% Set NaN values to be transparent by creating an alpha mask
alpha_mask_exp = ~isnan(plot_data_masked);
set(h_img_exp_experimental, 'AlphaData', alpha_mask_exp);

% Set background color to show through transparent areas
set(imgax_exp_experimental, 'Color', [0.9 0.9 0.9]);

% Set color limits
abs_max = max(abs(plot_data(:)));
caxis(imgax_exp_experimental, [-abs_max, abs_max]);

% Add labels
xlabel(imgax_exp_experimental, 'Time (ms)', 'Color', 'k');
ylabel(imgax_exp_experimental, 'Channel', 'Color', 'k');

% Flip y-axis so channel 1 is at the top
set(imgax_exp_experimental, 'YDir', 'reverse');
ylim(imgax_exp_experimental, [1 size(plot_data,1)]);
grid(imgax_exp_experimental, 'off'); % Turn off all grid lines for the specific axes

% Set tick label colors to black
set(imgax_exp_experimental, 'XColor', 'k', 'YColor', 'k');
hold(imgax_exp_experimental, 'on');

% Draw inverted mask and labels if any are defined
if ~isempty(highlight_rects)
    % --- Inverted Masking Logic ---
    num_channels = size(plot_data, 1);
    
    % 1. Get channel ranges to highlight
    highlight_channels = cellfun(@(c) sort([c(1,2), c(2,2)]), highlight_rects, 'UniformOutput', false);
    highlight_channels = vertcat(highlight_channels{:});
    highlight_channels = sortrows(highlight_channels);

    % 2. Merge overlapping highlight intervals
    if ~isempty(highlight_channels)
        merged = highlight_channels(1,:);
        for k = 2:size(highlight_channels, 1)
            if highlight_channels(k,1) <= merged(end,2) + 1
                merged(end,2) = max(merged(end,2), highlight_channels(k,2));
            else
                merged = [merged; highlight_channels(k,:)];
            end
        end
        highlight_channels = merged;
    end
    
    % 3. Generate mask ranges
    mask_ranges = [];
    last_chan = 0;
    for k = 1:size(highlight_channels, 1)
        if highlight_channels(k,1) > last_chan + 1
            mask_ranges = [mask_ranges; last_chan + 1, highlight_channels(k,1) - 1];
        end
        last_chan = highlight_channels(k,2);
    end
    if last_chan < num_channels
        mask_ranges = [mask_ranges; last_chan + 1, num_channels];
    end
    
    % 4. Draw masking rectangles
    if length(times) > 1; time_step = times(2) - times(1); else; time_step = 0; end
    for k = 1:size(mask_ranges, 1)
        y = mask_ranges(k,1) - 0.5;
        height = mask_ranges(k,2) - mask_ranges(k,1) + 1;
        x = times(1) - time_step/2;
        width = times(end) - times(1) + time_step;
        
        rectangle('Position', [x, y, width, height], ...
                  'FaceColor', [0 0 0], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
    end

    % --- Original Label Logic ---
    axis_font_size = get(imgax_exp_experimental, 'FontSize');
    time_range = times(end) - times(1);
    label_offset = 0.02 * time_range;

    for i = 1:length(highlight_rects)
        coords = highlight_rects{i};
        y = min(coords(:,2)) - 0.5;
        height = abs(coords(1,2) - coords(2,2)) + 1;
        if i <= length(highlight_labels) && ~isempty(highlight_labels{i})
            text(times(1) - label_offset, y + height/2 - 0.25, ...
                 highlight_labels{i}, 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'middle', 'Color', 'k', ...
                 'FontSize', axis_font_size, 'Clipping', 'off');
        end
    end
end

% Add vertical line at time=0
ylim_vals = get(imgax_exp_experimental, 'YLim');
plot(imgax_exp_experimental, [0 0], ylim_vals, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1);


if plot_topoplots_experimental
    % Pre-calculate topoplot positions for line drawing
    topo_positions = zeros(num_tf_points, 4); % [x, y, width, height] for each topoplot
    for n = 1:num_tf_points
        % Position topoplots evenly across the top of the figure
        % Space them evenly from left to right, regardless of their time points
        topo_x = 0.1 + (n-1) * (0.75 / num_tf_points) + (0.75 / num_tf_points - 0.15) / 2; % Adjusted for new topo_size
        topo_y = 0.69; % Position above main plot with more space from title
        topo_size = 0.15; % Size of each topoplot
        
        topo_positions(n, :) = [topo_x, topo_y, topo_size, topo_size];
    end

    % Store all data for drawing lines and markers in proper order
    marker_data_exp_experimental = zeros(num_tf_points, 2);
    line_coords_exp_experimental = zeros(num_tf_points, 4);

    % Calculate line coordinates using pre-calculated topoplot positions
    for n = 1:num_tf_points
        time_point = time_channel_points(n, 1);
        channel_point = time_channel_points(n, 2);
        
        % Store marker data for later
        marker_data_exp_experimental(n,:) = [time_point, channel_point];
        
        % Get coordinates for the line in figure units
        main_pos = get(imgax_exp_experimental, 'Position');
        ylims = get(imgax_exp_experimental, 'YLim');
        time_norm = (time_point - min(times)) / (max(times) - min(times));
        channel_norm = (channel_point - ylims(1)) / (ylims(2) - ylims(1));
        from_x = main_pos(1) + time_norm * main_pos(3);
        from_y = main_pos(2) + (1 - channel_norm) * main_pos(4); % Corrected for reversed y-axis
        
        topo_pos = topo_positions(n, :); % Use pre-calculated position
        to_x = topo_pos(1) + topo_pos(3) * 0.5;
        to_y = topo_pos(2) + 0.009;
        
        line_coords_exp_experimental(n,:) = [from_x, to_x, from_y, to_y];
    end
end

% Set focus to the main plot
axes(imgax_exp_experimental);

if plot_topoplots_experimental
    % Draw connecting lines on the main plot by converting figure to data coordinates
    main_pos = get(imgax_exp_experimental, 'Position');
    xlims = get(imgax_exp_experimental, 'XLim');
    ylims = get(imgax_exp_experimental, 'YLim');

    for n = 1:size(line_coords_exp_experimental, 1)
        % Line endpoints in figure units
        from_x_fig = line_coords_exp_experimental(n, 1);
        to_x_fig = line_coords_exp_experimental(n, 2);
        from_y_fig = line_coords_exp_experimental(n, 3);
        to_y_fig = line_coords_exp_experimental(n, 4);

        % Vector from start to end in figure units
        dx_fig = to_x_fig - from_x_fig;
        dy_fig = to_y_fig - from_y_fig;
        len_fig = sqrt(dx_fig^2 + dy_fig^2);

        % Normalize the vector
        udx_fig = dx_fig / len_fig;
        udy_fig = dy_fig / len_fig;

        % Shorten by a fixed amount in figure units. This value may need tuning.
        shorten_amount = 0.0057; 
        from_x_fig_new = from_x_fig + shorten_amount * udx_fig;
        from_y_fig_new = from_y_fig + shorten_amount * udy_fig;

        % Convert new start point from figure units back to data units
        time_norm_start_new = (from_x_fig_new - main_pos(1)) / main_pos(3);
        time_start_new = time_norm_start_new * (xlims(2) - xlims(1)) + xlims(1);
        channel_norm_start_new = 1 - ((from_y_fig_new - main_pos(2)) / main_pos(4));
        channel_start_new = channel_norm_start_new * (ylims(2) - ylims(1)) + ylims(1);

        % Convert end point from figure units back to data units
        time_norm_end = (to_x_fig - main_pos(1)) / main_pos(3);
        time_end = time_norm_end * (xlims(2) - xlims(1)) + xlims(1);
        channel_norm_end = 1 - ((to_y_fig - main_pos(2)) / main_pos(4));
        channel_end = channel_norm_end * (ylims(2) - ylims(1)) + ylims(1);

        % Draw connecting line in black
        plot([time_start_new, time_end], [channel_start_new, channel_end], 'k', 'LineWidth', 2, 'Clipping', 'off');
    end

    % Then draw all markers in the plot, so they appear on top of the lines
    for n = 1:num_tf_points
        time_point = marker_data_exp_experimental(n, 1);
        channel_point = marker_data_exp_experimental(n, 2);
        
        % Add marker at the time-channel point with black outline
        plot(time_point, channel_point, 'ko', 'MarkerSize', 11, ...
            'MarkerFaceColor', 'none', 'LineWidth', 2);
    end

    % Store topoplot axes handles
    topoaxes_exp_experimental = zeros(1, num_tf_points);

    % Now create topoplots evenly spaced across the top (like tftopo) - they will appear above lines
    for n = 1:num_tf_points
    time_point = time_channel_points(n, 1);
    channel_point = time_channel_points(n, 2);
    
    % Use pre-calculated topoplot position
    topo_pos = topo_positions(n, :);
    
    % Create topoplot axes
    topoaxes_exp_experimental(n) = axes('Position', topo_pos);
    
    % Get the time point and find closest time index
    [~, time_idx] = min(abs(times - time_point));
    
    % Get data for this time point (all channels)
    topo_data = plot_data(:, time_idx);
    
    % Get t, p, and df for title
    t_value = plot_data(channel_point, time_idx);
    p_value = p_values(channel_point, time_idx);
    df = df_values(channel_point, time_idx);

    % Create topoplot and then modify white background patches to be transparent
    topoplot(topo_data, LIMO.data.chanlocs, ... % DO NOT APPLY THE MASK TO TOPOPLOTS
        'maplimits', [-abs_max, abs_max], ...
        'colormap', cc, ...
        'electrodes', 'off', ...
        'style', 'both', ...
        'shading', 'interp', ...
        'numcontour', 6, ...
        'hcolor', 'k'); % Keep head outline black
    
    % Mark the specific channel electrode with a solid black dot
    hold on;
    % Create a simple overlay topoplot with just electrode markers
    topoplot([], LIMO.data.chanlocs, ...
        'plotchans', channel_point, ...
        'colormap', cc, ...
        'electrodes', 'on', ...
        'emarker', {'.', 'k', 16, 1}, ...
        'style', 'blank', ...
        'hcolor', 'none');
    hold off;
    
    % Find and make the white background patches transparent
    patch_handles = findobj(gca, 'Type', 'patch');
    for i = 1:length(patch_handles)
        % Check if this patch has the BACKCOLOR (light blue/white)
        face_color = get(patch_handles(i), 'FaceColor');
        if isnumeric(face_color) && length(face_color) == 3
            % Check if it's close to the BACKCOLOR [.93 .96 1]
            if all(abs(face_color - [0.93 0.96 1]) < 0.1)
                % Make this patch transparent
                set(patch_handles(i), 'FaceAlpha', 0);
            end
        end
    end
    
    % Add title showing time and channel, and stats
    title_line1 = sprintf('E%d, %d ms', channel_point, time_point);
    if p_value < 0.001
        p_string = 'p < 0.001';
    else
        p_string = sprintf('p = %.3f', p_value);
    end
    title_line2 = sprintf('t(%d) = %.3f, %s', df, t_value, p_string);
    title({title_line1, title_line2}, 'FontSize', 10, 'Color', 'k'); % 8 * 1.2 = 9.6, rounded to 10
    
    % Make topoplots circular
    axis tight;
    axis equal;
    
    % Add circular outline around topoplot with transparent white background
    set(gca, 'XTick', [], 'YTick', [], 'Box', 'off');
    hold on;
    % Draw circular outline (smaller to match color border)
    theta = 0:0.01:2*pi;
    xlims = xlim;
    ylims = ylim;
    radius = min(diff(xlims), diff(ylims))/2 * 0.93; % Make 8% smaller to match color border
    center_x = mean(xlims);
    center_y = mean(ylims);
    circle_x = center_x + radius * cos(theta);
    circle_y = center_y + radius * sin(theta);
    % Draw white background circle first (for typical outline effect)
    h_white = plot(circle_x, circle_y, 'Color', [1 1 1 0], 'LineWidth', 4.5); % Transparent white
    % Draw black outline on top (thinner to match other elements)
    plot(circle_x, circle_y, 'k', 'LineWidth', 1);
    hold off;
    end
end

% Ensure LIMO colormap is applied to the main time-channel axes
colormap(imgax_exp_experimental, cc);

% Add colorbar positioned to the right of the main plot
h = colorbar(imgax_exp_experimental, 'Position', [0.87 0.08 0.03 0.55]); % Adjusted to match new main plot position
title(h, 't-values', 'FontSize', 10, 'Color', 'k');
set(h, 'Color', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', 0); % Set colorbar text and outline to black, remove internal gridlines

% Final formatting: ensure all axes elements are black
set(imgax_exp_experimental, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', [0 0]);
% Set custom y-axis ticks. If custom labels are used, only show the
% first and last channel numbers. Otherwise, show 10 equally spaced values.
num_channels = size(plot_data, 1);
has_custom_labels = ~isempty(highlight_labels) && any(cellfun(@(x) ~isempty(x), highlight_labels));

if has_custom_labels
    ytick_values = [1, num_channels];
else
    ytick_values = round(linspace(1, num_channels, 10));
end
set(imgax_exp_experimental, 'YTick', ytick_values);

% Set custom x-axis ticks - always include last time point
tick_interval = 250;
min_tick = ceil(times(1) / tick_interval) * tick_interval;
max_tick = floor(times(end) / tick_interval) * tick_interval;

% Create regular ticks
regular_ticks = min_tick:tick_interval:max_tick;

% Always include the actual end time
all_ticks_exp = [regular_ticks, times(end)];

% Remove duplicates and sort
all_ticks_exp = unique(all_ticks_exp);

set(imgax_exp_experimental, 'XTick', all_ticks_exp);

% Add main title using annotation
annotation('textbox', [0, 0.85, 1, 0.15], ...
           'String', [LIMO.cache.fig.title, ' Time-Channel Plot (Experimental)'], ...
           'EdgeColor', 'none', ...
           'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'middle', ...
           'FontSize', 21, ...
           'Color', 'k', ...
           'FitBoxToText', 'on');