%% Parallelized Channel-Time Plotting for Multiple Tests
% This script processes 6 paired t-tests in parallel and generates
% channel-time plots with TFCE-corrected p-values

%% Setup: Define all tests with titles and paths
tests = struct(...
    'name', {}, ...
    'title', {}, ...
    'limo_path', {}, ...
    'data_path', {});

% Tamari SME
tests(1).name = 'tamari_sme';
tests(1).title = 'Tamari Subsequent Memory Effect (Hits vs Misses)';
tests(1).limo_path = '/home/devon7y/scratch/devon7y/tamari/paired_ttest_shits_smisses_tamari/LIMO.mat';
tests(1).data_path = '/home/devon7y/scratch/devon7y/tamari/paired_ttest_shits_smisses_tamari/paired_samples_ttest_parameter_12.mat';

% Tamari RSE
tests(2).name = 'tamari_rse';
tests(2).title = 'Tamari Retrieval Success Effect (Hits vs Misses)';
tests(2).limo_path = '/home/devon7y/scratch/devon7y/tamari/paired_ttest_thits_tmisses_tamari/LIMO.mat';
tests(2).data_path = '/home/devon7y/scratch/devon7y/tamari/paired_ttest_thits_tmisses_tamari/paired_samples_ttest_parameter_34.mat';

% Tamari ONE
tests(3).name = 'tamari_one';
tests(3).title = 'Tamari Old/New Effect (Hits vs Correct Rejections)';
tests(3).limo_path = '/home/devon7y/scratch/devon7y/tamari/paired_ttest_thits_crs_tamari/LIMO.mat';
tests(3).data_path = '/home/devon7y/scratch/devon7y/tamari/paired_ttest_thits_crs_tamari/paired_samples_ttest_parameter_35.mat';

% Yvonne SME
tests(4).name = 'yvonne_sme';
tests(4).title = 'Yvonne Subsequent Memory Effect (Hits vs Misses)';
tests(4).limo_path = '/home/devon7y/scratch/devon7y/yvonne/paired_ttest_shits_smisses/LIMO.mat';
tests(4).data_path = '/home/devon7y/scratch/devon7y/yvonne/paired_ttest_shits_smisses/paired_samples_ttest_parameter_34.mat';

% Yvonne RSE
tests(5).name = 'yvonne_rse';
tests(5).title = 'Yvonne Retrieval Success Effect (Hits vs Misses)';
tests(5).limo_path = '/home/devon7y/scratch/devon7y/yvonne/paired_ttest_thits_tmisses/LIMO.mat';
tests(5).data_path = '/home/devon7y/scratch/devon7y/yvonne/paired_ttest_thits_tmisses/paired_samples_ttest_parameter_56.mat';

% Yvonne ONE
tests(6).name = 'yvonne_one';
tests(6).title = 'Yvonne Old/New Effect (Hits vs Correct Rejections)';
tests(6).limo_path = '/home/devon7y/scratch/devon7y/yvonne/paired_ttest_thits_crs/LIMO.mat';
tests(6).data_path = '/home/devon7y/scratch/devon7y/yvonne/paired_ttest_thits_crs/paired_samples_ttest_parameter_51.mat';

% Output directory
output_dir = '/home/devon7y/scratch/devon7y/plots/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Parallel Processing Setup
% Start parallel pool with maximum available workers
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');  % Uses all available cores
end

fprintf('Processing %d tests in parallel using %d workers...\n', length(tests), pool.NumWorkers);

%% Parallel Loop - Process all tests simultaneously
parfor test_idx = 1:length(tests)
    try
        fprintf('Worker %d: Starting %s\n', test_idx, tests(test_idx).name);

        % Extract file information
        [~, FileName, ext] = fileparts(tests(test_idx).data_path);
        FileName = [FileName ext];

        % Load LIMO and data
        test_limo = load(tests(test_idx).limo_path);
        test_data = load(tests(test_idx).data_path);
        LIMO = test_limo.LIMO;

        % Get TFCE-corrected p-values
        [tfce_p_values, mask, ~] = limo_stat_values(FileName, 0.05, 3, LIMO);

        % Load TFCE scores to get thresholds
        tfce_file = fullfile(LIMO.dir, 'tfce', ['tfce_' FileName]);
        H0_tfce_file = fullfile(LIMO.dir, 'H0', ['tfce_H0_' FileName]);

        if exist(tfce_file, 'file') && exist(H0_tfce_file, 'file')
            tfce_data = load(tfce_file);
            tfce_score = tfce_data.tfce_score;
            H0_tfce_data = load(H0_tfce_file);
            H0_tfce_score = H0_tfce_data.tfce_H0_score;

            % Get thresholds using limo_max_correction
            [~, ~, bootstrap_threshold] = limo_max_correction(tfce_score, H0_tfce_score, 0.05);
            max_observed_tfce = max(tfce_score(:));

            % Always print thresholds for comparison
            if max_observed_tfce < bootstrap_threshold
                fprintf('Worker %d (%s): Max observed TFCE = %.3f, Bootstrap threshold = %.3f [NOT SIGNIFICANT]\n', ...
                    test_idx, tests(test_idx).name, max_observed_tfce, bootstrap_threshold);
            else
                fprintf('Worker %d (%s): Max observed TFCE = %.3f, Bootstrap threshold = %.3f [SIGNIFICANT]\n', ...
                    test_idx, tests(test_idx).name, max_observed_tfce, bootstrap_threshold);
            end
        end

        % Extract data components
        one_sample = test_data.paired_samples;
        df_values = squeeze(one_sample(:,:,3));
        t_values = squeeze(one_sample(:,:,4));
        p_values = tfce_p_values;  % Use TFCE-corrected p-values

        % Define output file path
        output_file = fullfile(output_dir, sprintf('%s_channel_time_plot.png', tests(test_idx).name));

        % Check if file exists (for logging only)
        file_exists = exist(output_file, 'file');

        %% Time-Channel Plot with Topographic Plots

        % Define rectangles for highlighting
        highlight_rects = {[400, 21; 800, 21], [400, 36; 800, 36], [400, 101; 800, 101], [400, 153; 800, 153]};
        highlight_labels = {'E21', 'E36', 'E101', 'E153'};

        % Significance threshold (p-value)
        alpha_threshold = 0.05;

        % Alpha value for non-significant pixels
        alpha_nonsig = 0.4;

        % Create signed -log10(p-values) data
        signed_logp_data = -log10(p_values) .* sign(t_values);

        % Create alpha mask for significance-based transparency
        alpha_mask_significance = ones(size(p_values));
        alpha_mask_significance(p_values >= alpha_threshold) = alpha_nonsig;

        % LIMO-style colormap for signed -log10(p-values)
        abs_max_val = max(abs(signed_logp_data(:)));
        if isempty(abs_max_val) || abs_max_val == 0
            abs_max_val = 1;
        end
        cc = limo_color_images([-abs_max_val, abs_max_val]);

        % Topoplot configuration (empty = no topoplots)
        time_channel_points_original = [];
        plot_topoplots = ~isempty(time_channel_points_original);

        if plot_topoplots
            [~, sort_idx] = sort(time_channel_points_original(:, 1));
            time_channel_points = time_channel_points_original(sort_idx, :);
            num_tf_points = size(time_channel_points, 1);
        else
            time_channel_points = [];
            num_tf_points = 0;
        end

        % Time setup
        times = linspace(LIMO.data.start, LIMO.data.end, size(t_values,2));

        % Create figure (INVISIBLE for parallel processing)
        fig_exp = figure('Visible', 'off', 'Position', [100 100 1000 700], 'Color', 'w', 'InvertHardcopy', 'off');

        % Apply LIMO colormap
        colormap(fig_exp, cc);

        % Create main time-channel plot
        if plot_topoplots
            imgax_exp = axes('Position', [0.08 0.10 0.78 0.55]);
        else
            imgax_exp = axes('Position', [0.08 0.12 0.78 0.68]);
        end
        set(imgax_exp, 'Color', [0.9 0.9 0.9], 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1);

        % Plot the signed -log10(p-value) data
        h_img_exp = imagesc(times, 1:size(signed_logp_data,1), signed_logp_data);
        xlim([times(1), times(end)]);

        % Apply significance-based alpha masking
        set(h_img_exp, 'AlphaData', alpha_mask_significance);
        set(imgax_exp, 'Color', [0.9 0.9 0.9]);

        % Set color limits
        caxis([-abs_max_val, abs_max_val]);

        % Add labels
        xlabel('Time (ms)', 'Color', 'k');
        ylabel('Channel', 'Color', 'k');

        % Flip y-axis
        set(imgax_exp, 'YDir', 'reverse');
        ylim([1 size(signed_logp_data,1)]);
        grid(imgax_exp, 'off');
        set(imgax_exp, 'XColor', 'k', 'YColor', 'k');
        hold on;

        % Draw highlighted rectangles and labels
        if ~isempty(highlight_rects)
            axis_font_size = get(imgax_exp, 'FontSize');
            time_range = times(end) - times(1);
            label_offset = 0.004 * time_range;

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

                [~, time_idx1] = min(abs(times - time1));
                [~, time_idx2] = min(abs(times - time2));

                exact_time1 = times(time_idx1);
                exact_time2 = times(time_idx2);

                x = min(exact_time1, exact_time2) - time_step/2;
                y = min(chan1, chan2) - 0.5;
                width = abs(exact_time2 - exact_time1) + time_step;
                height = abs(chan2 - chan1) + 1;

                % Use patch instead of rectangle for better headless compatibility
                patch([x x+width x+width x], [y y y+height y+height], [0 1 0], ...
                      'FaceAlpha', 0.5, 'EdgeColor', 'none');

                if i <= length(highlight_labels) && ~isempty(highlight_labels{i})
                    label_text = highlight_labels{i};
                    label_y = y + height/2 - 0.25;
                    label_x = times(1) - label_offset;

                    text(label_x, label_y, label_text, ...
                         'HorizontalAlignment', 'right', ...
                         'VerticalAlignment', 'middle', ...
                         'Color', 'k', ...
                         'FontSize', axis_font_size, ...
                         'Clipping', 'off');
                end
            end
        end

        % Add vertical line at time=0
        ylim_vals = get(imgax_exp, 'YLim');
        plot(imgax_exp, [0 0], ylim_vals, 'Color', [0.3 0.3 0.3], 'LineStyle', '--', 'LineWidth', 1);

        % Topoplot section (if enabled)
        if plot_topoplots
            topo_positions = zeros(num_tf_points, 4);
            for n = 1:num_tf_points
                topo_x = 0.1 + (n-1) * (0.75 / num_tf_points) + (0.75 / num_tf_points - 0.15) / 2;
                topo_y = 0.69;
                topo_size = 0.15;
                topo_positions(n, :) = [topo_x, topo_y, topo_size, topo_size];
            end

            marker_data_exp = zeros(num_tf_points, 2);
            line_coords_exp = zeros(num_tf_points, 4);

            for n = 1:num_tf_points
                time_point = time_channel_points(n, 1);
                channel_point = time_channel_points(n, 2);
                marker_data_exp(n,:) = [time_point, channel_point];

                main_pos = get(imgax_exp, 'Position');
                ylims = get(imgax_exp, 'YLim');
                time_norm = (time_point - min(times)) / (max(times) - min(times));
                channel_norm = (channel_point - ylims(1)) / (ylims(2) - ylims(1));
                from_x = main_pos(1) + time_norm * main_pos(3);
                from_y = main_pos(2) + (1 - channel_norm) * main_pos(4);

                topo_pos = topo_positions(n, :);
                to_x = topo_pos(1) + topo_pos(3) * 0.5;
                to_y = topo_pos(2) + 0.009;

                line_coords_exp(n,:) = [from_x, to_x, from_y, to_y];
            end

            axes(imgax_exp);
            main_pos = get(imgax_exp, 'Position');
            xlims = get(imgax_exp, 'XLim');
            ylims = get(imgax_exp, 'YLim');

            for n = 1:size(line_coords_exp, 1)
                from_x_fig = line_coords_exp(n, 1);
                to_x_fig = line_coords_exp(n, 2);
                from_y_fig = line_coords_exp(n, 3);
                to_y_fig = line_coords_exp(n, 4);

                dx_fig = to_x_fig - from_x_fig;
                dy_fig = to_y_fig - from_y_fig;
                len_fig = sqrt(dx_fig^2 + dy_fig^2);

                udx_fig = dx_fig / len_fig;
                udy_fig = dy_fig / len_fig;

                shorten_amount = 0.0057;
                from_x_fig_new = from_x_fig + shorten_amount * udx_fig;
                from_y_fig_new = from_y_fig + shorten_amount * udy_fig;

                time_norm_start_new = (from_x_fig_new - main_pos(1)) / main_pos(3);
                time_start_new = time_norm_start_new * (xlims(2) - xlims(1)) + xlims(1);
                channel_norm_start_new = 1 - ((from_y_fig_new - main_pos(2)) / main_pos(4));
                channel_start_new = channel_norm_start_new * (ylims(2) - ylims(1)) + ylims(1);

                time_norm_end = (to_x_fig - main_pos(1)) / main_pos(3);
                time_end = time_norm_end * (xlims(2) - xlims(1)) + xlims(1);
                channel_norm_end = 1 - ((to_y_fig - main_pos(2)) / main_pos(4));
                channel_end = channel_norm_end * (ylims(2) - ylims(1)) + ylims(1);

                plot([time_start_new, time_end], [channel_start_new, channel_end], 'k', 'LineWidth', 2, 'Clipping', 'off');
            end

            for n = 1:num_tf_points
                time_point = marker_data_exp(n, 1);
                channel_point = marker_data_exp(n, 2);
                plot(time_point, channel_point, 'ko', 'MarkerSize', 11, 'MarkerFaceColor', 'none', 'LineWidth', 2);
            end

            for n = 1:num_tf_points
                time_point = time_channel_points(n, 1);
                channel_point = time_channel_points(n, 2);
                topo_pos = topo_positions(n, :);
                topoaxes_exp_n = axes('Position', topo_pos);

                [~, time_idx] = min(abs(times - time_point));
                topo_data = signed_logp_data(:, time_idx);

                t_value = t_values(channel_point, time_idx);
                p_value = p_values(channel_point, time_idx);
                df = df_values(channel_point, time_idx);
                signed_logp_value = signed_logp_data(channel_point, time_idx);

                topoplot(topo_data, LIMO.data.chanlocs, ...
                    'maplimits', [-abs_max_val, abs_max_val], ...
                    'colormap', cc, ...
                    'electrodes', 'off', ...
                    'style', 'both', ...
                    'shading', 'interp', ...
                    'numcontour', 6, ...
                    'hcolor', 'k');

                hold on;
                topoplot([], LIMO.data.chanlocs, ...
                    'plotchans', channel_point, ...
                    'colormap', cc, ...
                    'electrodes', 'on', ...
                    'emarker', {'.', 'k', 16, 1}, ...
                    'style', 'blank', ...
                    'hcolor', 'none');
                hold off;

                patch_handles = findobj(gca, 'Type', 'patch');
                for patch_i = 1:length(patch_handles)
                    face_color = get(patch_handles(patch_i), 'FaceColor');
                    if isnumeric(face_color) && length(face_color) == 3
                        if all(abs(face_color - [0.93 0.96 1]) < 0.1)
                            set(patch_handles(patch_i), 'FaceAlpha', 0);
                        end
                    end
                end

                title_line1 = sprintf('E%d, %d ms', channel_point, time_point);
                if p_value < 0.001
                    p_string = 'p < 0.001';
                else
                    p_string = sprintf('p = %.3f', p_value);
                end
                title_line2 = sprintf('t(%d) = %.3f, %s', df, t_value, p_string);
                title_line3 = sprintf('-log10(p) = %.3f', signed_logp_value);
                title({title_line1, title_line2, title_line3}, 'FontSize', 9, 'Color', 'k');

                axis tight;
                axis equal;
                set(gca, 'XTick', [], 'YTick', [], 'Box', 'off');
                hold on;
                theta = 0:0.01:2*pi;
                xlims_topo = xlim;
                ylims_topo = ylim;
                radius = min(diff(xlims_topo), diff(ylims_topo))/2 * 0.93;
                center_x = mean(xlims_topo);
                center_y = mean(ylims_topo);
                circle_x = center_x + radius * cos(theta);
                circle_y = center_y + radius * sin(theta);
                plot(circle_x, circle_y, 'Color', [1 1 1 0], 'LineWidth', 4.5);
                plot(circle_x, circle_y, 'k', 'LineWidth', 1);
                hold off;
            end
        end

        % Ensure colormap is applied
        colormap(imgax_exp, cc);

        % Add colorbar
        if plot_topoplots
            h = colorbar(imgax_exp, 'Position', [0.88 0.10 0.04 0.55]);
        else
            h = colorbar(imgax_exp, 'Position', [0.88 0.12 0.04 0.68]);
        end

        title(h, 'signed -log10(p)', 'FontSize', 10, 'Color', 'k');
        set(h, 'Color', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', 0);

        % Add significance threshold lines to colorbar
        sig_threshold = -log10(alpha_threshold);
        cbar_pos = get(h, 'Position');
        cbar_limits = get(h, 'Limits');

        y_pos = NaN;
        y_neg = NaN;

        if sig_threshold >= cbar_limits(1) && sig_threshold <= cbar_limits(2)
            norm_pos_pos = (sig_threshold - cbar_limits(1)) / (cbar_limits(2) - cbar_limits(1));
            y_pos = cbar_pos(2) + norm_pos_pos * cbar_pos(4);
        end

        if -sig_threshold >= cbar_limits(1) && -sig_threshold <= cbar_limits(2)
            norm_pos_neg = (-sig_threshold - cbar_limits(1)) / (cbar_limits(2) - cbar_limits(1));
            y_neg = cbar_pos(2) + norm_pos_neg * cbar_pos(4);
        end

        % If no significant pixels exist, cover entire colorbar
        if isnan(y_neg) && isnan(y_pos)
            y_neg = cbar_pos(2);  % Bottom of colorbar
            y_pos = cbar_pos(2) + cbar_pos(4);  % Top of colorbar
        end

        if ~isnan(y_neg) && ~isnan(y_pos)
            border_inset = 0.0006;
            overlay_x = cbar_pos(1) + border_inset;
            overlay_y = y_neg;
            overlay_width = cbar_pos(3) - 2*border_inset;
            overlay_height = y_pos - y_neg;

            annotation('rectangle', [overlay_x, overlay_y, overlay_width, overlay_height], ...
                'FaceColor', [0.9137, 0.9020, 0.9020], ...
                'FaceAlpha', 1 - alpha_nonsig, ...
                'EdgeColor', 'none');
        end

        if ~isnan(y_pos)
            annotation('line', [cbar_pos(1), cbar_pos(1) + cbar_pos(3)], [y_pos, y_pos], ...
                'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
        end

        if ~isnan(y_neg)
            annotation('line', [cbar_pos(1), cbar_pos(1) + cbar_pos(3)], [y_neg, y_neg], ...
                'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);
        end

        % Final formatting
        set(imgax_exp, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 1, 'TickLength', [0 0]);

        num_channels = size(signed_logp_data, 1);
        has_custom_labels = ~isempty(highlight_labels) && any(cellfun(@(x) ~isempty(x), highlight_labels));

        if has_custom_labels
            ytick_values = [1, num_channels];
        else
            ytick_values = round(linspace(1, num_channels, 10));
        end
        set(imgax_exp, 'YTick', ytick_values);

        tick_interval = 250;
        min_tick = ceil(times(1) / tick_interval) * tick_interval;
        max_tick = floor(times(end) / tick_interval) * tick_interval;
        regular_ticks = min_tick:tick_interval:max_tick;
        all_ticks = unique([regular_ticks, times(end)]);
        set(imgax_exp, 'XTick', all_ticks);

        % Adjust text sizes
        set(imgax_exp, 'FontSize', 7);
        xlabel(imgax_exp, 'Time (ms)', 'Color', 'k', 'FontSize', 8);
        ylabel(imgax_exp, 'Channel', 'Color', 'k', 'FontSize', 8);

        h_xlabel = get(imgax_exp, 'XLabel');
        h_ylabel = get(imgax_exp, 'YLabel');
        xlabel_pos = get(h_xlabel, 'Position');
        ylabel_pos = get(h_ylabel, 'Position');
        set(h_xlabel, 'Position', [xlabel_pos(1), xlabel_pos(2) + 3, xlabel_pos(3)]);
        set(h_ylabel, 'Position', [ylabel_pos(1) - 25, ylabel_pos(2), ylabel_pos(3)]);

        title(h, 'signed -log10(p)', 'FontSize', 8, 'Color', 'k');
        set(h, 'FontSize', 7);

        if ~isempty(highlight_rects)
            text_objs = findobj(imgax_exp, 'Type', 'text');
            for text_i = 1:length(text_objs)
                if contains(get(text_objs(text_i), 'String'), {'E', 'Fz', 'F3', 'F4', 'Pz'})
                    set(text_objs(text_i), 'FontSize', 7);
                end
            end
        end

        % Add title using dynamic test title
        plot_title = sprintf('%s With TFCE Correction', tests(test_idx).title);

        % Add "No Significance" suffix if no pixels are significant
        if ~any(p_values(:) < alpha_threshold)
            plot_title = sprintf('%s - No Significance', plot_title);
        end

        if plot_topoplots
            annotation('textbox', [0, 0.85, 1, 0.15], ...
                       'String', plot_title, ...
                       'EdgeColor', 'none', ...
                       'HorizontalAlignment', 'center', ...
                       'VerticalAlignment', 'middle', ...
                       'FontSize', 12, ...
                       'Color', 'k');
        else
            annotation('textbox', [0, 0.78, 1, 0.10], ...
                       'String', plot_title, ...
                       'EdgeColor', 'none', ...
                       'HorizontalAlignment', 'center', ...
                       'VerticalAlignment', 'middle', ...
                       'FontSize', 14, ...
                       'Color', 'k');
        end

        % Set figure properties for export
        set(fig_exp, 'PaperPositionMode', 'auto');
        set(fig_exp, 'PaperUnits', 'inches');
        set(fig_exp, 'InvertHardCopy', 'off');

        % Export at high resolution (automatically overwrites)
        print(fig_exp, output_file, '-dpng', '-r600');

        % Clean up
        close(fig_exp);

        % Log completion with overwrite status
        if file_exists
            fprintf('Worker %d: Completed %s -> %s (OVERWRITTEN)\n', test_idx, tests(test_idx).name, output_file);
        else
            fprintf('Worker %d: Completed %s -> %s (NEW)\n', test_idx, tests(test_idx).name, output_file);
        end

    catch ME
        fprintf('Worker %d: ERROR in %s: %s\n', test_idx, tests(test_idx).name, ME.message);
        fprintf('  Stack trace:\n');
        for stack_i = 1:length(ME.stack)
            fprintf('    %s (line %d)\n', ME.stack(stack_i).name, ME.stack(stack_i).line);
        end
    end
end

fprintf('\nAll tests completed! Plots saved to: %s\n', output_dir);
