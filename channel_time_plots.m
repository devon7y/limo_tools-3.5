%% Specify data for plots

% Tamari SME
tamari_sme_limo = load('/home/devon7y/scratch/devon7y/tamari/paired_ttest_shits_smisses_tamari/LIMO.mat');
tamari_sme_data = load('/home/devon7y/scratch/devon7y/tamari/paired_ttest_shits_smisses_tamari/paired_samples_ttest_parameter_12.mat');
[tfce_p_values, mask, title] = limo_stat_values(FileName, 0.05, 3, LIMO);

% Tamari RSE
tamari_rse_limo = load('/home/devon7y/scratch/devon7y/tamari/paired_ttest_thits_tmisses_tamari/LIMO.mat');
tamari_rse_data = load('/home/devon7y/scratch/devon7y/tamari/paired_ttest_thits_tmisses_tamari/paired_samples_ttest_parameter_34.mat');

% Tamari ONE
tamari_one_limo = load('/home/devon7y/scratch/devon7y/tamari/paired_ttest_thits_crs_tamari/LIMO.mat');
tamari_one_data = load('/home/devon7y/scratch/devon7y/tamari/paired_ttest_thits_crs_tamari/paired_samples_ttest_parameter_35.mat');

% Yvonne SME
yvonne_sme_limo = load('/home/devon7y/scratch/devon7y/yvonne/paired_ttest_shits_smisses/LIMO.mat');
yvonne_sme_data = load('/home/devon7y/scratch/devon7y/yvonne/paired_ttest_shits_smisses/paired_samples_ttest_parameter_34.mat');

% Yvonne RSE
yvonne_rse_limo = load('/home/devon7y/scratch/devon7y/yvonne/paired_ttest_thits_tmisses/LIMO.mat');
yvonne_rse_data = load('/home/devon7y/scratch/devon7y/yvonne/paired_ttest_thits_tmisses/paired_samples_ttest_parameter_56.mat');

% Yvonne ONE
yvonne_one_limo = load('/home/devon7y/scratch/devon7y/yvonne/paired_ttest_thits_crs/LIMO.mat');
yvonne_one_data = load('/home/devon7y/scratch/devon7y/yvonne/paired_ttest_thits_crs/paired_samples_ttest_parameter_51.mat');


load('/Volumes/T7/Yvonnes ERP Files/Interpol/Reconstructed/No Reference Electrode/Epoched Files 50/paired_ttest_shits_smisses_mac/LIMO.mat');
%data = load('/Volumes/T7/ERP Files/Epoched Files 50/paired_ttest3/paired_samples_ttest_parameter_35.mat');
data = load('/Volumes/T7/Yvonnes ERP Files/Interpol/Reconstructed/No Reference Electrode/Epoched Files 50/paired_ttest_shits_smisses_mac/paired_samples_ttest_parameter_34_likelihood.mat');
one_sample = data.paired_samples;

%[M, mask, mytitle] = limo_stat_values('paired_samples_ttest_parameter_12.mat', 0.05, 3, LIMO);

%LIMO.cache.fig.pval = M;
%LIMO.cache.fig.mask = mask;
%LIMO.cache.fig.MCC = 3;
%LIMO.cache.fig.threshold = 0.05;

df_values = squeeze(one_sample(:,:,3)); % df-values (3rd dimension)
t_values = squeeze(one_sample(:,:,4)); % t-values (4th dimension)
%p_values = squeeze(one_sample(:,:,5)); % p-values (5th dimension)
p_values = squeeze(LIMO.cache.fig.pval); % cluster/tfce p-values (LIMO cache)
likelihood_values = squeeze(one_sample(:,:,6)); % likelihood values (6th dimension)

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
signed_logp_data = -log10(p_values) .* sign(t_values);

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
times = linspace(LIMO.data.start, LIMO.data.end, size(t_values,2));

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
        t_value = t_values(channel_point, time_idx);
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