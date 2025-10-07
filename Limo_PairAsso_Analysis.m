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
parameters = 5;
estimator = 'Mean'; % 'Mean' 'Weighted mean' 'Trimmed mean' 'HD' 'Median'
analysis_type = 'Mean'; % 'Mean' 'Trimmed mean' 'HD' 'Median'
savename = [ '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/parameter_' num2str(parameters) '.mat'];
limo_central_tendency_and_ci(files, parameters, chan_loc, estimator, analysis_type, [],savename);

%% Plot averages for each condition

average_files = {'/Volumes/T7/ERP Files/Epoched Files 50/derivatives/parameter_3_Mean_of_Mean.mat', '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/parameter_4_Mean_of_Mean.mat'};
% average_files = {'/Volumes/T7/ERP Files/Original Epoched Files/derivatives/parameter_5_Mean_of_Mean.mat', '/Volumes/T7/ERP Files/Original Epoched Files/derivatives/parameter_6_Mean_of_Mean.mat'};
channel = [155];

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

% IMPROVED: Proper color coordination between ERP lines and confidence intervals
% Strategy: Identify main ERP lines by their LineWidth and coordinate colors
if length(average_files) > 1
    % Get all line objects and patch objects (confidence intervals)
    h_lines = findobj(gca, 'Type', 'line');
    h_patches = findobj(gca, 'Type', 'patch');

    fprintf('Found %d line objects and %d patch objects for %d conditions\n', ...
            length(h_lines), length(h_patches), length(average_files));

    % Get default color order
    colorOrder = get(gca, 'ColorOrder');

    % Identify main ERP lines (typically have thicker LineWidth)
    main_lines = [];
    for i = 1:length(h_lines)
        linewidth = get(h_lines(i), 'LineWidth');
        if linewidth > 1 % Main ERP lines are usually thicker
            main_lines(end+1) = h_lines(i);
        end
    end

    % If we don't find thick lines, use the first N lines (most recent)
    if length(main_lines) < length(average_files)
        main_lines = h_lines(1:min(length(average_files), length(h_lines)));
    end

    % Apply consistent colors to patches first, then lines (for proper layering)
    legend_handles = [];

    % First, set colors for all confidence interval patches (background layer)
    for i = 1:min(length(average_files), length(h_patches))
        color_idx = mod(i-1, size(colorOrder, 1)) + 1;
        condition_color = colorOrder(color_idx,:);

        % Set matching color for confidence interval patch
        set(h_patches(i), 'FaceColor', condition_color, 'FaceAlpha', 0.3, ...
            'EdgeColor', condition_color, 'EdgeAlpha', 0.5);
    end

    % Then, set colors for main ERP lines and ensure they're on top
    for i = 1:min(length(average_files), length(main_lines))
        color_idx = mod(i-1, size(colorOrder, 1)) + 1;
        condition_color = colorOrder(color_idx,:);

        % Set color for main ERP line
        set(main_lines(i), 'Color', condition_color, 'LineWidth', 3);
        legend_handles(end+1) = main_lines(i);

        % Bring line to front by moving it to the end of the children list
        uistack(main_lines(i), 'top');
    end

    % Create legend if we have handles
    if ~isempty(legend_handles)
        % Reverse the order to match the typical plotting order
        legend(fliplr(legend_handles), fliplr(legend_labels), 'Location', 'best');
    end
else
    % Single condition - just update the title, no need for complex legend
    fprintf('Single condition plot - no legend modification needed\n');
end

%% Plot single-subject data for each condition
single_subject_files = ['/Volumes/T7/ERP Files/Epoched Files 50/derivatives/parameter_1_single_subjects_Mean.mat'];
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

data1 = '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/parameter_3_single_subjects_Mean.mat';
data2 = '/Volumes/T7/ERP Files/Epoched Files 50/derivatives/parameter_5_single_subjects_Mean.mat';
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

load('/Volumes/T7/ERP Files/Epoched Files 50/paired_ttest4/LIMO.mat');
%data = load('/Volumes/T7/ERP Files/Epoched Files 50/paired_ttest3/paired_samples_ttest_parameter_35.mat');
data = load('/Volumes/T7/ERP Files/Epoched Files 50/paired_ttest4/paired_samples_ttest_parameter_34.mat');
one_sample = data.paired_samples;

%[M, mask, mytitle] = limo_stat_values('paired_samples_ttest_parameter_12.mat', 0.05, 3, LIMO);

%LIMO.cache.fig.pval = M;
%LIMO.cache.fig.mask = mask;
%LIMO.cache.fig.MCC = 3;
%LIMO.cache.fig.threshold = 0.05;

df_values = squeeze(one_sample(:,:,3)); % df-values (3rd dimension)
plot_data = squeeze(one_sample(:,:,4)); % t-values (4th dimension)
%p_values = squeeze(one_sample(:,:,5)); % p-values (5th dimension)
p_values = squeeze(LIMO.cache.fig.pval); % cluster/tfce p-values (LIMO cache)
%likelihood_values = squeeze(one_sample(:,:,6)); % likelihood values (6th dimension)

%% Plot the electodes versus time 2D and 3D graphs

% The function has 3 parameters pointing to the test directory, so it's
% easier to go to that directory first and use pwd
cd '/Volumes/T7/ERP Files/Epoched Files 50/paired_ttest2'
limo_display_results(1, ... # electodes x time graph
    'paired_samples_ttest_parameter_34.mat', ... # test statistics mat file
    pwd, ... # path to the test directory
    0.05, ... # significance level
    2, ... # MCC: 1 = none, 2 = cluster, 3 = TFCE, 4 = max
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

latencies = 0:100:2500; % time plots to display

pop_topoplot(EEG, 1, latencies, 'Paired T-test Topoplots', 0);

%% Custom topoplot graphs with modified titles, size, and colorbar (uses topoplot)

% LIMO-style masking
plot_data_masked = plot_data;
plot_data_masked(p_values >= 0.05) = NaN;

% LIMO-style colormap
cc = limo_color_images(plot_data(~isnan(plot_data)));

% Time setup
times = linspace(LIMO.data.start, LIMO.data.end, size(plot_data,2));
latencies = 200:200:2500;

% Find time indices for each latency
time_indices = arrayfun(@(x) find(abs(times - x) == min(abs(times - x)), 1), latencies);

% Create figure with subplots for 2D topoplots (much taller to give title space)
fig = figure('Position', [100 50 1200 1200]); 
ax_handles = [];  % Store axis handles

% Calculate symmetric color limits like the 3D plots
abs_max = max(abs(plot_data(:)));
colorbar_limits = [-abs_max, abs_max];

for i = 1:length(latencies)
    ax = subplot(3, 4, i);
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
sgtitle('Paired T-test Between Test Hits & Correct Rejections Topoplots', 'FontSize', 14, 'Color', 'k');

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
latencies = 0:100:2500;

% Find time indices for each latency
time_indices = arrayfun(@(x) find(abs(times - x) == min(abs(times - x)), 1), latencies);

% Create figure with subplots for 3D headplots (much taller to give title space)
fig = figure('Position', [100 50 1200 1200], 'Color', 'w', 'InvertHardcopy', 'off'); 
ax_handles = [];  % Store axis handles for synchronized rotation

% Use same logic as topoplot: symmetric limits based on absolute maximum
abs_max = max(abs(plot_data(:)));  % Get absolute maximum like 'absmax' does
colorbar_limits = [-abs_max, abs_max];  % Symmetric limits like topoplot

for i = 1:length(latencies)
    ax = subplot(5, 7, i);
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
highlight_rects = {[400, 21; 800, 21], [400, 36; 800, 36], [400, 101; 800, 101], [400, 153; 800, 153]};
% Example: [400, 21; 700, 21], [400, 36; 700, 36], [400, 101; 700, 101], [-100, 224; 1500, 224]

% Define labels for the highlighted rectangles. The number of labels should
% match the number of rectangles.
highlight_labels = {'E21', 'E36', 'E101', 'E153'};

% Significance threshold (p-value)
alpha_threshold = 0.05;

% Alpha value for non-significant pixels (0 = invisible, 0.3 = 30% opacity, etc.)
alpha_nonsig = 0.4;
% Example: 'Fz (E21)', 'F3 (E36)', 'F4 (E224)', 'Pz (E101)'

% Create signed -log10(p-values) data
% Compute -log10(p-values) and inherit sign from t-values
signed_logp_data = -log10(p_values) .* sign(plot_data);

% Create alpha mask for significance-based transparency
% Significant pixels get full opacity (1), non-significant get reduced opacity
alpha_mask_significance = ones(size(p_values));
alpha_mask_significance(p_values >= alpha_threshold) = alpha_nonsig;

% LIMO-style colormap for signed -log10(p-values), ensuring it's symmetric around 0
abs_max_val = max(abs(signed_logp_data(:)));
if isempty(abs_max_val) || abs_max_val == 0
    abs_max_val = 1; % handle case with no significant data
end
cc = limo_color_images([-abs_max_val, abs_max_val]);

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
    imgax_exp = axes('Position', [0.08 0.10 0.78 0.55]); % Lower position with space for topoplots above
else
    imgax_exp = axes('Position', [0.08 0.12 0.78 0.68]); % Lower the plot to create more space above for title
end
set(imgax_exp, 'Color', [0.9 0.9 0.9], 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1); % Set axes background and and outline to black with thicker border

% Plot the signed -log10(p-value) data with alpha masking
h_img_exp = imagesc(times, 1:size(signed_logp_data,1), signed_logp_data);
xlim([times(1), times(end)]);

% Apply significance-based alpha masking
set(h_img_exp, 'AlphaData', alpha_mask_significance);

% Set background color to show through transparent areas
set(imgax_exp, 'Color', [0.9 0.9 0.9]);

% Set color limits for signed -log10(p-values)
caxis([-abs_max_val, abs_max_val]);

% Add labels
xlabel('Time (ms)', 'Color', 'k');
ylabel('Channel', 'Color', 'k');

% Flip y-axis so channel 1 is at the top
set(imgax_exp, 'YDir', 'reverse');
ylim([1 size(signed_logp_data,1)]);
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

% Update variable name reference for consistency
signed_logp_data_for_limits = signed_logp_data;

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
        
        % Get data for this time point (all channels) - use full data, not masked
        topo_data = signed_logp_data(:, time_idx);
        
        % Get t, p, and df for title
        t_value = plot_data(channel_point, time_idx);
        p_value = p_values(channel_point, time_idx);
        df = df_values(channel_point, time_idx);
        signed_logp_value = signed_logp_data(channel_point, time_idx);

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
        title_line3 = sprintf('-log10(p) = %.3f', signed_logp_value);
        title({title_line1, title_line2, title_line3}, 'FontSize', 9, 'Color', 'k');
        
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
    h = colorbar(imgax_exp, 'Position', [0.88 0.10 0.04 0.55]); % Match lower plot position
else
    h = colorbar(imgax_exp, 'Position', [0.88 0.12 0.04 0.68]); % Match adjusted plot position
end

% Set up colorbar with signed -log10(p-value) labels
title(h, 'signed -log10(p)', 'FontSize', 10, 'Color', 'k');
set(h, 'Color', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', 0); % Set colorbar text and outline to black, remove internal gridlines

% Add significance threshold lines to colorbar
sig_threshold = -log10(alpha_threshold); % Threshold on -log10 scale
% Get colorbar position and limits
cbar_pos = get(h, 'Position'); % [x, y, width, height] in figure normalized coordinates
cbar_limits = get(h, 'Limits'); % [min, max] data values

% Calculate normalized positions of threshold lines within colorbar
y_pos = NaN; % Initialize
y_neg = NaN; % Initialize

% For positive threshold
if sig_threshold >= cbar_limits(1) && sig_threshold <= cbar_limits(2)
    norm_pos_pos = (sig_threshold - cbar_limits(1)) / (cbar_limits(2) - cbar_limits(1));
    y_pos = cbar_pos(2) + norm_pos_pos * cbar_pos(4);
end

% For negative threshold
if -sig_threshold >= cbar_limits(1) && -sig_threshold <= cbar_limits(2)
    norm_pos_neg = (-sig_threshold - cbar_limits(1)) / (cbar_limits(2) - cbar_limits(1));
    y_neg = cbar_pos(2) + norm_pos_neg * cbar_pos(4);
end

% Add semi-transparent overlay for non-significant region (between thresholds)
if ~isnan(y_neg) && ~isnan(y_pos)
    % Calculate the overlay rectangle position with slight inset to avoid covering borders
    border_inset = 0.0006; % Small inset to avoid covering colorbar borders
    overlay_x = cbar_pos(1) + border_inset;
    overlay_y = y_neg;
    overlay_width = cbar_pos(3) - 2*border_inset;
    overlay_height = y_pos - y_neg;

    % Create semi-transparent rectangle to "fade out" non-significant region
    % Color #E9E6E6 in RGB: [233, 230, 230] / 255 = [0.9137, 0.9020, 0.9020]
    annotation('rectangle', [overlay_x, overlay_y, overlay_width, overlay_height], ...
        'FaceColor', [0.9137, 0.9020, 0.9020], ...  % Light gray overlay (#E9E6E6)
        'FaceAlpha', 1 - alpha_nonsig, ...  % Transparency inverse of plot alpha
        'EdgeColor', 'none');
end

% Draw threshold lines on top of the overlay
if ~isnan(y_pos)
    annotation('line', [cbar_pos(1), cbar_pos(1) + cbar_pos(3)], [y_pos, y_pos], ...
        'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
end

if ~isnan(y_neg)
    annotation('line', [cbar_pos(1), cbar_pos(1) + cbar_pos(3)], [y_neg, y_neg], ...
        'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
end

% Final formatting: ensure all axes elements are black
set(imgax_exp, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', [0 0]);
% Set custom y-axis ticks. If custom labels are used, only show the
% first and last channel numbers. Otherwise, show 10 equally spaced values.
num_channels = size(signed_logp_data_for_limits, 1);
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

% Adjust text sizes for high-resolution export
set(imgax_exp, 'FontSize', 7); % Axis tick labels
xlabel(imgax_exp, 'Time (ms)', 'Color', 'k', 'FontSize', 8); % X-axis label
ylabel(imgax_exp, 'Channel', 'Color', 'k', 'FontSize', 8); % Y-axis label

% Adjust axis label positions
h_xlabel = get(imgax_exp, 'XLabel');
h_ylabel = get(imgax_exp, 'YLabel');
xlabel_pos = get(h_xlabel, 'Position');
ylabel_pos = get(h_ylabel, 'Position');
set(h_xlabel, 'Position', [xlabel_pos(1), xlabel_pos(2) + 3, xlabel_pos(3)]); % Move x-label down by 5 units
set(h_ylabel, 'Position', [ylabel_pos(1) - 25, ylabel_pos(2), ylabel_pos(3)]); % Move y-label left by 15 units

% Colorbar text
title(h, 'signed -log10(p)', 'FontSize', 8, 'Color', 'k'); % Colorbar title
set(h, 'FontSize', 7); % Colorbar tick labels

% Adjust highlight labels text size (if any exist)
if ~isempty(highlight_rects)
    text_objs = findobj(imgax_exp, 'Type', 'text');
    for i = 1:length(text_objs)
        if contains(get(text_objs(i), 'String'), {'E', 'Fz', 'F3', 'F4', 'Pz'}) % Common electrode patterns
            set(text_objs(i), 'FontSize', 7); % Adjust highlight label font size
        end
    end
end

% Add title using annotation positioned above the plot
% Determine correction method suffix based on LIMO.cache.fig.MCC
if LIMO.cache.fig.MCC == 1
    correction_suffix = ' Uncorrected';
elseif LIMO.cache.fig.MCC == 2
    correction_suffix = ' With Cluster Correction';
elseif LIMO.cache.fig.MCC == 3
    correction_suffix = ' With TFCE Correction';
else
    correction_suffix = ''; % Default case
end

% Create dynamic title as single string
dynamic_title = strcat(LIMO.cache.fig.title, correction_suffix);

if plot_topoplots
    annotation('textbox', [0, 0.85, 1, 0.15], ...
               'String', 'Subsequent Memory Effect Paired T-test With TFCE Correction', ... % dynamic_title
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FontSize', 12, ...
               'Color', 'k');
else
    annotation('textbox', [0, 0.78, 1, 0.10], ... % Higher position above the plot
               'String', 'Subsequent Memory Effect Paired T-test With TFCE Correction', ... % dynamic_title
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FontSize', 14, ...
               'Color', 'k');
end

% Set figure properties for high-resolution export with scaled text
set(fig_exp, 'PaperPositionMode', 'auto');
set(fig_exp, 'PaperUnits', 'inches');
set(fig_exp, 'InvertHardCopy', 'off');

% Export at high resolution (600 DPI) - text will scale automatically
print(fig_exp, '/Users/devon7y/Downloads/time_channel_plot_high_res.png', '-dpng', '-r600');


%% Experimental Darken for Rectangles
% This section demonstrates an alternative highlighting approach where
% non-selected channels are masked with a semi-transparent black overlay.

% Define rectangles for highlighting. Each cell contains a 2x2 matrix for
% a rectangle, defined by its top-left and bottom-right corners in the
% format: [time_ms, channel; time_ms, channel].
highlight_rects = {[400, 21; 700, 21], [400, 36; 700, 36], [400, 101; 700, 101], [-100, 224; 1500, 224]}; % Example: [200, 20; 400, 40], [600, 80; 800, 100]

% Define labels for the highlighted rectangles. The number of labels should
% match the number of rectangles.
highlight_labels = {'Fz (E21)', 'F3 (E36)', 'F4 (E224)', 'Pz (E101)'};

% Alpha value for non-significant pixels (0 = invisible, 0.3 = 30% opacity, etc.)
% Note: Using same parameter as main section for consistency
% alpha_nonsig = 0.3; % Uncomment if running experimental section independently

% Create alpha mask for experimental section (same approach as main section)
alpha_mask_significance_exp = ones(size(p_values));
alpha_mask_significance_exp(p_values >= 0.05) = alpha_nonsig;

% LIMO-style colormap (use full data range, not masked)
cc = limo_color_images(plot_data(~isnan(plot_data)));

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

% Plot the data with significance-based alpha masking
h_img_exp_experimental = imagesc(times, 1:size(plot_data,1), plot_data);

% Apply significance-based alpha masking (same as main section)
set(h_img_exp_experimental, 'AlphaData', alpha_mask_significance_exp);

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
    
    % Get data for this time point (all channels) - use full data, not masked
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

%% Time-Channel Plot with Likelihood Ratios

% Define rectangles for highlighting. Each cell contains a 2x2 matrix for
% a rectangle, defined by its top-left and bottom-right corners in the
% format: [time_ms, channel; time_ms, channel].
highlight_rects_bf = {[400, 21; 800, 21], [400, 36; 800, 36], [400, 101; 800, 101], [400, 153; 800, 153]};
% Example: [400, 21; 700, 21], [400, 36; 700, 36], [400, 101; 700, 101], [-100, 224; 1500, 224]

% Define labels for the highlighted rectangles. The number of labels should
% match the number of rectangles.
highlight_labels_bf = {'E21', 'E36', 'E101', 'E153'};

% Likelihood Ratio threshold for strong evidence categorization
% User can modify this value - values above this show strong evidence for H1,
% values below 1/threshold show strong evidence for H0
LR_strong_threshold = 10;       % Strong evidence threshold (Kang et al. 2015)

% Alpha value for inconclusive evidence pixels (0 = invisible, 0.4 = 40% opacity, etc.)
alpha_nonsig_lr = 0.4;

% Use log10-transformed Likelihood Ratios for symmetric visualization
% log10(LR) = 0 when LR = 1 (no evidence either way)
% log10(LR) > 0 when LR > 1 (evidence for H1, plotted in red)
% log10(LR) < 0 when LR < 1 (evidence for H0, plotted in blue)
LR_data = log10(likelihood_values);

% Create alpha mask for evidence-based transparency
% Strong evidence pixels get full opacity (1), inconclusive evidence gets reduced opacity
alpha_mask_evidence = ones(size(likelihood_values));
alpha_mask_evidence(likelihood_values < LR_strong_threshold & likelihood_values > 1/LR_strong_threshold) = alpha_nonsig_lr;

% LIMO-style colormap for log10(LR), ensuring it's symmetric around 0
abs_max_lr = max(abs(LR_data(:)));
if isempty(abs_max_lr) || abs_max_lr == 0
    abs_max_lr = 1; % handle case with no strong evidence
end
cc_lr = limo_color_images([-abs_max_lr, abs_max_lr]);

% Create time-channel plot data
% Placeholder values from tutorial: [400 8; 350 14; 500 24; 1050 11]
% Format: [time_ms, channel]
time_channel_points_original_bf = [];

% Check if topoplots should be plotted
plot_topoplots_bf = ~isempty(time_channel_points_original_bf);

if plot_topoplots_bf
    % Sort points by time (x-axis) for proper topoplot ordering
    [~, sort_idx] = sort(time_channel_points_original_bf(:, 1));
    time_channel_points_bf = time_channel_points_original_bf(sort_idx, :);

    % Define the number of points to plot
    num_tf_points_bf = size(time_channel_points_bf, 1);
else
    time_channel_points_bf = [];
    num_tf_points_bf = 0;
end

% Time setup
times_bf = linspace(LIMO.data.start, LIMO.data.end, size(plot_data,2));

% Create figure for time-channel plot with topoplots above
fig_bf = figure('Position', [100 100 1000 700], 'Color', 'w', 'InvertHardcopy', 'off');

% Apply LIMO colormap to this figure
colormap(fig_bf, cc_lr);

% Create main time-channel plot - adjust position based on whether topoplots are present
if plot_topoplots_bf
    imgax_bf = axes('Position', [0.08 0.10 0.78 0.55]); % Lower position with space for topoplots above
else
    imgax_bf = axes('Position', [0.08 0.12 0.78 0.68]); % Lower the plot to create more space above for title
end
set(imgax_bf, 'Color', [0.9 0.9 0.9], 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1); % Set axes background and and outline to black with thicker border

% Plot the log10(LR) data (no alpha masking)
h_img_bf = imagesc(times_bf, 1:size(LR_data,1), LR_data);
xlim([times_bf(1), times_bf(end)]);

% Apply alpha masking (all ones = no masking, full opacity)
set(h_img_bf, 'AlphaData', alpha_mask_evidence);

% Set background color to show through transparent areas
set(imgax_bf, 'Color', [0.9 0.9 0.9]);

% Set color limits for log10(LR) - symmetric around 0
caxis([-abs_max_lr, abs_max_lr]);

% Add labels
xlabel('Time (ms)', 'Color', 'k');
ylabel('Channel', 'Color', 'k');

% Flip y-axis so channel 1 is at the top
set(imgax_bf, 'YDir', 'reverse');
ylim([1 size(LR_data,1)]);
grid(imgax_bf, 'off'); % Turn off all grid lines for the specific axes

% Set tick label colors to black
set(imgax_bf, 'XColor', 'k', 'YColor', 'k');
hold on;

% Draw highlighted rectangles and labels if any are defined
if ~isempty(highlight_rects_bf)
    % Get properties for label positioning and sizing
    axis_font_size = get(imgax_bf, 'FontSize');
    time_range = times_bf(end) - times_bf(1);
    label_offset = 0.004 * time_range; % 2% offset to the left

    % Get time step for pixel adjustment, assuming uniform spacing
    if length(times_bf) > 1
        time_step = times_bf(2) - times_bf(1);
    else
        time_step = 0;
    end

    for i = 1:length(highlight_rects_bf)
        coords = highlight_rects_bf{i};
        time1 = coords(1,1);
        chan1 = coords(1,2);
        time2 = coords(2,1);
        chan2 = coords(2,2);

        % Adjust coordinates to cover the full pixels.
        % For imagesc, pixels are centered on the coordinate points.
        % We need to extend the rectangle by half a pixel in each direction.

        % Find the time indices corresponding to the time values
        [~, time_idx1] = min(abs(times_bf - time1));
        [~, time_idx2] = min(abs(times_bf - time2));

        % Get the exact time values from the times vector
        exact_time1 = times_bf(time_idx1);
        exact_time2 = times_bf(time_idx2);

        x = min(exact_time1, exact_time2) - time_step/2;
        y = min(chan1, chan2) - 0.5;
        width = abs(exact_time2 - exact_time1) + time_step;
        height = abs(chan2 - chan1) + 1;

        rectangle('Position', [x, y, width, height], ...
                  'FaceColor', [0 1 0], ... % Green
                  'FaceAlpha', 0.5, ...      % Set alpha value for transparency
                  'EdgeColor', 'none');

        % Add a label to the left of the rectangle
        if i <= length(highlight_labels_bf) && ~isempty(highlight_labels_bf{i})
            label_text = highlight_labels_bf{i};
            % Position the label vertically in the middle of the rectangle, with a slight upward adjustment
            label_y = y + height/2 - 0.25;

            % Position the label horizontally to the left of the plot
            label_x = times_bf(1) - label_offset;

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
ylim_vals = get(imgax_bf, 'YLim');
plot(imgax_bf, [0 0], ylim_vals, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1);


if plot_topoplots_bf
    % Pre-calculate topoplot positions for line drawing
    topo_positions_bf = zeros(num_tf_points_bf, 4); % [x, y, width, height] for each topoplot
    for n = 1:num_tf_points_bf
        % Position topoplots evenly across the top of the figure
        % Space them evenly from left to right, regardless of their time points
        topo_x = 0.1 + (n-1) * (0.75 / num_tf_points_bf) + (0.75 / num_tf_points_bf - 0.15) / 2; % Adjusted for new topo_size
        topo_y = 0.69; % Position above main plot with more space from title
        topo_size = 0.15; % Size of each topoplot

        topo_positions_bf(n, :) = [topo_x, topo_y, topo_size, topo_size];
    end

    % Store all data for drawing lines and markers in proper order
    marker_data_bf = zeros(num_tf_points_bf, 2);
    line_coords_bf = zeros(num_tf_points_bf, 4);

    % Calculate line coordinates using pre-calculated topoplot positions
    for n = 1:num_tf_points_bf
        time_point = time_channel_points_bf(n, 1);
        channel_point = time_channel_points_bf(n, 2);

        % Store marker data for later
        marker_data_bf(n,:) = [time_point, channel_point];

        % Get coordinates for the line in figure units
        main_pos = get(imgax_bf, 'Position');
        ylims = get(imgax_bf, 'YLim');
        time_norm = (time_point - min(times_bf)) / (max(times_bf) - min(times_bf));
        channel_norm = (channel_point - ylims(1)) / (ylims(2) - ylims(1));
        from_x = main_pos(1) + time_norm * main_pos(3);
        from_y = main_pos(2) + (1 - channel_norm) * main_pos(4); % Corrected for reversed y-axis

        topo_pos = topo_positions_bf(n, :); % Use pre-calculated position
        to_x = topo_pos(1) + topo_pos(3) * 0.5;
        to_y = topo_pos(2) + 0.009;

        line_coords_bf(n,:) = [from_x, to_x, from_y, to_y];
    end
end

% Update variable name reference for consistency
LR_data_for_limits = LR_data;

% Set focus to the main plot
axes(imgax_bf);

if plot_topoplots_bf
    % Draw connecting lines on the main plot by converting figure to data coordinates
    main_pos = get(imgax_bf, 'Position');
    xlims = get(imgax_bf, 'XLim');
    ylims = get(imgax_bf, 'YLim');

    for n = 1:size(line_coords_bf, 1)
        % Line endpoints in figure units
        from_x_fig = line_coords_bf(n, 1);
        to_x_fig = line_coords_bf(n, 2);
        from_y_fig = line_coords_bf(n, 3);
        to_y_fig = line_coords_bf(n, 4);

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
    for n = 1:num_tf_points_bf
        time_point = marker_data_bf(n, 1);
        channel_point = marker_data_bf(n, 2);

        % Add marker at the time-channel point with black outline
        plot(time_point, channel_point, 'ko', 'MarkerSize', 11, ...
            'MarkerFaceColor', 'none', 'LineWidth', 2);
    end

    % Store topoplot axes handles
    topoaxes_bf = zeros(1, num_tf_points_bf);

    % Now create topoplots evenly spaced across the top (like tftopo) - they will appear above lines
    for n = 1:num_tf_points_bf
        time_point = time_channel_points_bf(n, 1);
        channel_point = time_channel_points_bf(n, 2);

        % Use pre-calculated topoplot position
        topo_pos = topo_positions_bf(n, :);

        % Create topoplot axes
        topoaxes_bf(n) = axes('Position', topo_pos);

        % Get the time point and find closest time index
        [~, time_idx] = min(abs(times_bf - time_point));

        % Get data for this time point (all channels) - use full data, not masked
        topo_data_bf = LR_data(:, time_idx);

        % Get t, LR, and df for title
        t_value = plot_data(channel_point, time_idx);
        LR_value = likelihood_values(channel_point, time_idx);
        df = df_values(channel_point, time_idx);
        log10_LR_value = LR_data(channel_point, time_idx);

        % Create topoplot and then modify white background patches to be transparent
        topoplot(topo_data_bf, LIMO.data.chanlocs, ... % DO NOT APPLY THE MASK TO TOPOPLOTS
            'maplimits', [-abs_max_lr, abs_max_lr], ...
            'colormap', cc_lr, ...
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
            'colormap', cc_lr, ...
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
        if LR_value >= 1000
            LR_string = sprintf('LR = %.2e', LR_value);
        elseif LR_value < 0.001
            LR_string = sprintf('LR = %.2e', LR_value);
        else
            LR_string = sprintf('LR = %.3f', LR_value);
        end
        title_line2 = sprintf('t(%d) = %.3f, %s', df, t_value, LR_string);
        title_line3 = sprintf('log10(LR) = %.3f', log10_LR_value);
        title({title_line1, title_line2, title_line3}, 'FontSize', 9, 'Color', 'k');

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
colormap(imgax_bf, cc_lr);

% Add colorbar positioned to the right of the main plot - adjust based on plot position
if plot_topoplots_bf
    h_bf = colorbar(imgax_bf, 'Position', [0.88 0.10 0.04 0.55]); % Match lower plot position
else
    h_bf = colorbar(imgax_bf, 'Position', [0.88 0.12 0.04 0.68]); % Match adjusted plot position
end

% Set up colorbar with log10(LR) labels
title(h_bf, 'log10(LR)', 'FontSize', 10, 'Color', 'k');
set(h_bf, 'Color', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', 0); % Set colorbar text and outline to black, remove internal gridlines

% Add Likelihood Ratio threshold lines to colorbar (in log10 space)
% LR thresholds: only strong evidence (10 for H1, 1/10 for H0)
log10_LR_threshold_strong_H1 = log10(LR_strong_threshold);      % e.g., log10(10) = 1.0
log10_LR_threshold_strong_H0 = log10(1/LR_strong_threshold);    % e.g., log10(1/10) = -1.0

% Get colorbar position and limits
cbar_pos_bf = get(h_bf, 'Position'); % [x, y, width, height] in figure normalized coordinates
cbar_limits_bf = get(h_bf, 'Limits'); % [min, max] data values

% Calculate normalized positions of threshold lines within colorbar
y_strong_H1 = NaN;
y_strong_H0 = NaN;

% For strong H1 threshold
if log10_LR_threshold_strong_H1 >= cbar_limits_bf(1) && log10_LR_threshold_strong_H1 <= cbar_limits_bf(2)
    norm_pos = (log10_LR_threshold_strong_H1 - cbar_limits_bf(1)) / (cbar_limits_bf(2) - cbar_limits_bf(1));
    y_strong_H1 = cbar_pos_bf(2) + norm_pos * cbar_pos_bf(4);
end

% For strong H0 threshold
if log10_LR_threshold_strong_H0 >= cbar_limits_bf(1) && log10_LR_threshold_strong_H0 <= cbar_limits_bf(2)
    norm_pos = (log10_LR_threshold_strong_H0 - cbar_limits_bf(1)) / (cbar_limits_bf(2) - cbar_limits_bf(1));
    y_strong_H0 = cbar_pos_bf(2) + norm_pos * cbar_pos_bf(4);
end

% Add semi-transparent overlay for inconclusive region (between thresholds)
if ~isnan(y_strong_H0) && ~isnan(y_strong_H1)
    % Calculate the overlay rectangle position with slight inset to avoid covering borders
    border_inset = 0.0006; % Small inset to avoid covering colorbar borders
    overlay_x = cbar_pos_bf(1) + border_inset;
    overlay_y = y_strong_H0;
    overlay_width = cbar_pos_bf(3) - 2*border_inset;
    overlay_height = y_strong_H1 - y_strong_H0;

    % Create semi-transparent rectangle to "fade out" inconclusive region
    % Color #E9E6E6 in RGB: [233, 230, 230] / 255 = [0.9137, 0.9020, 0.9020]
    annotation('rectangle', [overlay_x, overlay_y, overlay_width, overlay_height], ...
        'FaceColor', [0.9137, 0.9020, 0.9020], ...  % Light gray overlay (#E9E6E6)
        'FaceAlpha', 1 - alpha_nonsig_lr, ...  % Transparency inverse of plot alpha
        'EdgeColor', 'none');
end

% Draw threshold lines on top of the overlay
if ~isnan(y_strong_H1)
    annotation('line', [cbar_pos_bf(1), cbar_pos_bf(1) + cbar_pos_bf(3)], [y_strong_H1, y_strong_H1], ...
        'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
end

if ~isnan(y_strong_H0)
    annotation('line', [cbar_pos_bf(1), cbar_pos_bf(1) + cbar_pos_bf(3)], [y_strong_H0, y_strong_H0], ...
        'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
end

% Final formatting: ensure all axes elements are black
set(imgax_bf, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', [0 0]);
% Set custom y-axis ticks. If custom labels are used, only show the
% first and last channel numbers. Otherwise, show 10 equally spaced values.
num_channels_bf = size(LR_data, 1);
has_custom_labels_bf = ~isempty(highlight_labels_bf) && any(cellfun(@(x) ~isempty(x), highlight_labels_bf));

if has_custom_labels_bf
    ytick_values = [1, num_channels_bf];
else
    ytick_values = round(linspace(1, num_channels_bf, 10));
end
set(imgax_bf, 'YTick', ytick_values);

% Set custom x-axis ticks - always include last time point
tick_interval = 250;
min_tick = ceil(times_bf(1) / tick_interval) * tick_interval;
max_tick = floor(times_bf(end) / tick_interval) * tick_interval;

% Create regular ticks
regular_ticks = min_tick:tick_interval:max_tick;

% Always include the actual end time
all_ticks_bf = [regular_ticks, times_bf(end)];

% Remove duplicates and sort
all_ticks_bf = unique(all_ticks_bf);

set(imgax_bf, 'XTick', all_ticks_bf);

% Adjust text sizes for high-resolution export
set(imgax_bf, 'FontSize', 7); % Axis tick labels
xlabel(imgax_bf, 'Time (ms)', 'Color', 'k', 'FontSize', 8); % X-axis label
ylabel(imgax_bf, 'Channel', 'Color', 'k', 'FontSize', 8); % Y-axis label

% Adjust axis label positions
h_xlabel_bf = get(imgax_bf, 'XLabel');
h_ylabel_bf = get(imgax_bf, 'YLabel');
xlabel_pos_bf = get(h_xlabel_bf, 'Position');
ylabel_pos_bf = get(h_ylabel_bf, 'Position');
set(h_xlabel_bf, 'Position', [xlabel_pos_bf(1), xlabel_pos_bf(2) + 3, xlabel_pos_bf(3)]); % Move x-label down by 5 units
set(h_ylabel_bf, 'Position', [ylabel_pos_bf(1) - 25, ylabel_pos_bf(2), ylabel_pos_bf(3)]); % Move y-label left by 15 units

% Colorbar text
title(h_bf, 'log10(LR)', 'FontSize', 8, 'Color', 'k'); % Colorbar title
set(h_bf, 'FontSize', 7); % Colorbar tick labels

% Adjust highlight labels text size (if any exist)
if ~isempty(highlight_rects_bf)
    text_objs = findobj(imgax_bf, 'Type', 'text');
    for i = 1:length(text_objs)
        if contains(get(text_objs(i), 'String'), {'E', 'Fz', 'F3', 'F4', 'Pz'}) % Common electrode patterns
            set(text_objs(i), 'FontSize', 7); % Adjust highlight label font size
        end
    end
end

% Add title using annotation positioned above the plot
% Determine correction method suffix based on LIMO.cache.fig.MCC
if LIMO.cache.fig.MCC == 1
    correction_suffix_bf = ' Uncorrected';
elseif LIMO.cache.fig.MCC == 2
    correction_suffix_bf = ' With Cluster Correction';
elseif LIMO.cache.fig.MCC == 3
    correction_suffix_bf = ' With TFCE Correction';
else
    correction_suffix_bf = ''; % Default case
end

% Create dynamic title as single string
dynamic_title_bf = strcat(LIMO.cache.fig.title, ' (Likelihood Ratios)', correction_suffix_bf);

if plot_topoplots_bf
    annotation('textbox', [0, 0.85, 1, 0.15], ...
               'String', 'Subsequent Memory Effect Paired T-test Likelihood Ratios', ... % dynamic_title_bf
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FontSize', 12, ...
               'Color', 'k');
else
    annotation('textbox', [0, 0.78, 1, 0.10], ... % Higher position above the plot
               'String', 'Subsequent Memory Effect Paired T-test Likelihood Ratios', ... % dynamic_title_bf
               'EdgeColor', 'none', ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'middle', ...
               'FontSize', 14, ...
               'Color', 'k');
end

% Set figure properties for high-resolution export with scaled text
set(fig_bf, 'PaperPositionMode', 'auto');
set(fig_bf, 'PaperUnits', 'inches');
set(fig_bf, 'InvertHardCopy', 'off');

% Export at high resolution (600 DPI) - text will scale automatically
print(fig_bf, '/Users/devon7y/Downloads/time_channel_plot_LR_high_res.png', '-dpng', '-r600');
