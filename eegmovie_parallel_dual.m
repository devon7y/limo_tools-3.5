function [Movie, Colormap] = eegmovie_parallel_dual(data,srate,eloc_locs,varargin)
% eegmovie_parallel() - Optimized parallel version of eegmovie with dual topoplot support
%                       Generates movie frames in parallel for faster processing
%                       Now supports side-by-side topoplots for comparing two datasets
%
% Usage:
%   [Movie, Colormap] = eegmovie_parallel(data, srate, eloc_locs, 'key', value, ...)
%
% Inputs:
%   data       - EEG data matrix (channels x frames)
%   srate      - Sampling rate in Hz
%   eloc_locs  - Electrode locations structure
%
% Optional parameters:
%   'data2'      - Second dataset for dual topoplot mode (channels x frames)
%   'layout'     - 'single' (default) or 'dual' for side-by-side topoplots
%   'subtitle1'  - Title for left/first topoplot (e.g., 'Voltage (μV)')
%   'subtitle2'  - Title for right/second topoplot (e.g., 't-value')
%   'minmax2'    - Color scale limits for second dataset [min max]
%   'colormap2'  - Colormap for second dataset (n x 3 RGB matrix)
%                  Default: cool(64) for diverging data, hot(64) otherwise
%   'showcolorbars' - 'on' or 'off' to show colorbars for each topoplot (default: 'off')
%   'timecourse_channels' - Array of channel indices to display in timecourse (e.g., [1 5 10 15])
%   'timecourse_maxchans' - Maximum number of channels to display (auto-selected, evenly spaced)
%   'timecourse_scaling' - 'linked' (default) or 'independent' y-axis scaling for multi-channel plots
%   'resolution' - [width height] in pixels. Default: [1200 1000] for single, [1600 800] for dual

set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultAxesColor', 'white');

if nargin<1
    help eegmovie_parallel
    return
end

% Check if Parallel Computing Toolbox is available
if ~license('test', 'Distrib_Computing_Toolbox')
    warning('Parallel Computing Toolbox not available. Falling back to regular eegmovie.');
    % Note: Regular eegmovie doesn't support dual layout, so we'd need to handle this
    if any(strcmpi(varargin, 'layout')) && any(strcmpi(varargin, 'dual'))
        error('Dual layout requires Parallel Computing Toolbox for eegmovie_parallel');
    end
    [Movie, Colormap] = eegmovie(data,srate,eloc_locs,varargin{:});
    return
end

% Suppress the interpolated shading warning
warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');

[chans,frames] = size(data);

% Define DEFAULT_SRATE directly instead of calling icadefs
DEFAULT_SRATE = 256;  % Standard EEG sampling rate default

if nargin <2
    srate = 0;
end
if nargin <3
    eloc_locs = 0;
end

% Verify data and electrode locations compatibility
if ~isequal(eloc_locs, 0) && size(data, 1) ~= length(eloc_locs)
    error('Number of channels in data (%d) does not match number of electrode locations (%d)', size(data, 1), length(eloc_locs));
end

% Parse arguments (same as original)
if nargin > 5 && ~ischar(varargin{3}) || nargin == 4 && ~ischar(varargin{2})
    % legacy mode
    options = {};
    if nargin>=8, options = { options{:} 'topoplotopt' varargin(5:end) }; end
    if nargin>=7, options = { options{:} 'startsec'    varargin{4} }; end
    if nargin>=6, options = { options{:} 'minmax'      varargin{3} }; end
    if nargin>=5, options = { options{:} 'movieframes' varargin{2} }; end
    if nargin>=4, options = { options{:} 'title'       varargin{1} }; end
else
    options = varargin;
end

% Check for layout option early to set appropriate default resolution
layout_idx = find(strcmpi(options, 'layout'));
is_dual_layout = false;
if ~isempty(layout_idx) && layout_idx < length(options)
    is_dual_layout = strcmpi(options{layout_idx+1}, 'dual');
end

% Set default resolution based on layout
if is_dual_layout
    default_resolution = [1600 800];  % Wider for dual layout
else
    default_resolution = [1200 1000];  % Original default
end

opt = finputcheck(options, { 'startsec'    'real'    {}    0;
                             'minmax'      'real'    {}    0;
                             'movieframes' 'integer' {}    0;
                             'title'       'string'  {}    '';
                             'vert'        'real'    {}    [];
                             'mode'        'string'  { '2D' '3D'  }    '2D';
                             'timecourse'  'string'  { 'on' 'off' }    'off';
                             'framenum'    'string'  { 'on' 'off' }    'on';
                             'camerapath'  'real'    []                0;
                             'time'        'string'  { 'on' 'off' }    'off';
                             'resolution'  'real'    []                default_resolution;
                             'layout'      'string'  { 'single' 'dual' } 'single';  % New parameter
                             'data2'       'real'    []                [];          % New parameter
                             'subtitle1'   'string'  {}                '';          % New parameter
                             'subtitle2'   'string'  {}                '';          % New parameter
                             'minmax2'     'real'    {}                0;           % New parameter
                             'colormap1'   'real'    []                [];          % New parameter
                             'colormap2'   'real'    []                [];          % New parameter
                             'showcolorbars' 'string' { 'on' 'off' }   'off';       % New parameter for colorbars
                             'timecourse_channels' 'real' []    [];        % New parameter for channel selection
                             'timecourse_maxchans' 'integer' []    [];     % New parameter for max channels
                             'timecourse_scaling' 'string' {'linked' 'independent'} 'linked'; % New parameter for scaling
                             'headlinewidth' 'real'    {}    [];          % New parameter
                             'topo_linewidth' 'real'    {}    [];          % New parameter
                             'topoplotopt' 'cell'    {}    {};
                             'headplotopt' 'cell'    {}    {} }, 'eegmovie_parallel');
if ischar(opt), error(opt); end

% Validate dual layout requirements
if strcmpi(opt.layout, 'dual')
    if isempty(opt.data2)
        error('Dual layout requires ''data2'' parameter with second dataset');
    end
    
    % Check dimensions match
    if size(opt.data2, 1) ~= chans || size(opt.data2, 2) ~= frames
        error('data2 dimensions (%dx%d) must match data dimensions (%dx%d)', ...
              size(opt.data2, 1), size(opt.data2, 2), chans, frames);
    end
    
    % 3D mode not supported for dual layout (at least initially)
    if strcmpi(opt.mode, '3d')
        error('3D mode is not currently supported with dual layout');
    end
    
    % If timecourse is on with dual layout, adjust resolution for vertical layout
    if strcmpi(opt.timecourse, 'on') && isequal(opt.resolution, default_resolution)
        opt.resolution = [1400 1400];  % More square for vertical stacking
        fprintf('Adjusting resolution to %dx%d for dual layout with timecourse\n', ...
                opt.resolution(1), opt.resolution(2));
    end
end

% Calculate data limits for first dataset
if opt.minmax ==0,
    datamin = min(min(data));
    datamax = max(max(data));
    absmax  = max([abs(datamin), abs(datamax)]);
    fudge   = 0.05*(datamax-datamin);
    datamin = -absmax-fudge;
    datamax =  absmax+fudge;
    opt.minmax = [datamin datamax];
end

% Calculate data limits for second dataset if in dual mode
if strcmpi(opt.layout, 'dual') && isscalar(opt.minmax2) && opt.minmax2 == 0
    datamin2 = min(min(opt.data2));
    datamax2 = max(max(opt.data2));
    absmax2  = max([abs(datamin2), abs(datamax2)]);
    fudge2   = 0.05*(datamax2-datamin2);
    datamin2 = -absmax2-fudge2;
    datamax2 =  absmax2+fudge2;
    opt.minmax2 = [datamin2 datamax2];
end

if opt.movieframes == 0
    opt.movieframes = 1:frames;
end

% Validate movieframes against the actual data size
if any(opt.movieframes < 1) || any(opt.movieframes > frames)
    error('Some movieframes indices are out of range (1 to %d)', frames);
end

if srate ==0,
    srate = DEFAULT_SRATE;
end
if strcmpi(opt.time, 'on'), opt.framenum = 'off'; end

mframes = length(opt.movieframes);

% Debug output to verify dimensions
fprintf('Data dimensions: %d channels x %d frames\n', chans, frames);
if strcmpi(opt.layout, 'dual')
    fprintf('Data2 dimensions: %d channels x %d frames\n', size(opt.data2,1), size(opt.data2,2));
    fprintf('Layout: DUAL (side-by-side topoplots)\n');
end
fprintf('Movieframes: %d frames requested\n', mframes);
fprintf('Frame resolution: %d x %d pixels\n', opt.resolution(1), opt.resolution(2));

%%%%%%%%%%%%%%%%%%%%%%%%% OPTIMIZED PARALLEL SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimize parallel pool - use all available cores
pool = gcp('nocreate');
if isempty(pool)
    numCores = feature('numcores');
    fprintf('Starting parallel pool with %d workers (all available cores)...\n', numCores);
    parpool('local', numCores);
    pool = gcp;
else
    numCores = feature('numcores');
    if pool.NumWorkers < numCores
        fprintf('Expanding pool from %d to %d workers...\n', pool.NumWorkers, numCores);
        delete(pool);
        parpool('local', numCores);
        pool = gcp;
    end
end
fprintf('Using %d workers for parallel processing\n', pool.NumWorkers);
fprintf('Generating %d frames...\n\n', mframes);

% Pre-calculate camera positions for 3D mode (single layout only)
camera_positions = [];
if strcmpi(opt.mode, '3d') && strcmpi(opt.layout, 'single')
    headplot('setup', eloc_locs, 'tmp.spl', opt.headplotopt{:});
    
    if isequal(opt.camerapath, 0)
        opt.camerapath = [-127 0 30 0];
        fprintf('Using default view [-127 0 30 0].\n');
    end
    if size(opt.camerapath,2)~=4
        error('Camerapath parameter must have exact 4 columns');
    end
    camera_positions = calculate_camera_positions(opt.camerapath, opt.movieframes);
end

% Initialize output
MovieFrames = cell(1, mframes);
Colormap = [jet(64); [1 1 1]];

% Setup colormaps for dual layout
if strcmpi(opt.layout, 'dual')
    if ~isempty(opt.colormap1)
        colormap1 = opt.colormap1;
    else
        % Default colormap for first topoplot (typically voltage data)
        colormap1 = jet(64);
    end  
    
    % Colormap for second topoplot
    if ~isempty(opt.colormap2)
        colormap2 = opt.colormap2;
    else
        % Choose default colormap based on data range
        if opt.minmax2(1) < 0 && opt.minmax2(2) > 0
            % Diverging data (e.g., t-values) - use cool colormap
            colormap2 = cool(64);
        else
            % Non-diverging data - use hot colormap for contrast
            colormap2 = hot(64);
        end
        fprintf('Using default colormap for second topoplot (use ''colormap2'' parameter to customize)\n');
    end
else
    if ~isempty(opt.colormap1)
        colormap1 = opt.colormap1;
    else
        colormap1 = jet(64);
    end
    colormap2 = [];
end

% Create shared variables for workers
selected_data = data(:, opt.movieframes);
shared_data = selected_data;
if strcmpi(opt.layout, 'dual')
    selected_data2 = opt.data2(:, opt.movieframes);
    shared_data2 = selected_data2;
else
    shared_data2 = [];
end
shared_eloc_locs = eloc_locs;
shared_opt = opt;
shared_srate = srate;
shared_colormap1 = colormap1;
shared_colormap2 = colormap2;

%%%%%%%%%%%%%%%%%%%%%%%%% PARALLEL FRAME GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup parallel progress tracking using DataQueue
dq = parallel.pool.DataQueue;
frameCount = 0;
startTime = tic;

% Listener for progress updates
afterEach(dq, @(idx) updateProgress(idx));

    function updateProgress(idx)
        frameCount = frameCount + 1;
        fprintf('Frame %d/%d completed (%.1f%%)\n', ...
            frameCount, mframes, (frameCount/mframes)*100);
    end

% OPTIMIZED PARFOR LOOP WITH HIGH RESOLUTION
fprintf('Starting parallel frame generation at %dx%d resolution...\n', opt.resolution(1), opt.resolution(2));

% Pass colormaps to parallel workers
local_colormap1 = colormap1;
local_colormap2 = colormap2;

% Set base font size for subtitles - THIS IS THE VALUE TO MODIFY
subtitle_base_fontsize = 30;
shared_subtitle_base_fontsize = subtitle_base_fontsize;

parfor f = 1:mframes
    try
        % Initialize variables to prevent parfor warnings
        frame_data2 = [];
        ax1 = [];
        ax2 = [];

        % Create figure at exactly the target resolution
        fig = figure('Visible', 'off', ...
                    'Units', 'pixels', ...
                    'Position', [0 0 shared_opt.resolution(1) shared_opt.resolution(2)], ...
                    'Color', 'white', ...
                    'PaperPositionMode', 'auto', ...
                    'Renderer', 'zbuffer', ...
                    'Resize', 'off');
        
        % --- DYNAMIC SCALING ---
        baseline_height = 1000;
        scaling_factor = shared_opt.resolution(2) / baseline_height;
        
        % Scale font size and line widths
        dynamic_fontsize = max(8, round(28 * scaling_factor));
        dynamic_subtitle_fontsize = max(8, round(shared_subtitle_base_fontsize * scaling_factor));  % Smaller for subtitles
        if ~isempty(shared_opt.headlinewidth)
            head_linewidth = shared_opt.headlinewidth;
        else
            head_linewidth = 50; % Fixed large value for testing
        end
        if ~isempty(shared_opt.topo_linewidth)
            topo_linewidth = shared_opt.topo_linewidth;
        else
            topo_linewidth = 10; % Fixed large value for testing
        end
        
        % Get frame data
        frame_data = shared_data(:, f);
        if strcmpi(shared_opt.layout, 'dual')
            frame_data2 = shared_data2(:, f);
        end
        original_frame_idx = shared_opt.movieframes(f);
        
        % Handle layout based on mode
        if strcmpi(shared_opt.layout, 'dual')
            % DUAL LAYOUT - Create two topoplots side by side
            
            % Left topoplot - centered vertically
            ax1 = axes('Parent', fig, ...
                      'Units', 'normalized', ...
                      'Position', [0.025, 0.025, 0.45, 0.8], ...  % Vertically centered
                      'Color', 'white');
            
            % Suppress warning locally
            oldWarning = warning('off', 'all');
            
            % First topoplot with colormap1
            limo_topoplot_dual(frame_data, shared_eloc_locs, ...
                'maplimits', shared_opt.minmax, ...
                'electrodes', 'off', ...
                'shading', 'interp', ...
                'style', 'both', ...
                'numcontour', 6, ...
                'whitebk', 'on', ...
                'contourlinewidth', topo_linewidth, ...
                'colormap', shared_colormap1, ...
                shared_opt.topoplotopt{:});
            
            axis(ax1, 'equal'); % Enforce circular aspect ratio
            set(ax1, 'Color', 'none');
            
            % Set head outline width for first topoplot
            head_lines = findobj(ax1, 'Type', 'Line', 'Color', 'k');
            set(head_lines, 'LineWidth', head_linewidth);

            % Add circular black outline for first topoplot
            hold(ax1, 'on');
            theta = linspace(0, 2*pi, 360);
            xlims = get(ax1, 'XLim');
            ylims = get(ax1, 'YLim');
            radius = min(diff(xlims), diff(ylims))/2 * 0.93;
            center_x = mean(xlims);
            center_y = mean(ylims);
            circle_x = center_x + radius * cos(theta);
            circle_y = center_y + radius * sin(theta);
            patch(ax1, 'XData', circle_x, 'YData', circle_y, 'ZData', ones(size(circle_x))*10, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', head_linewidth / 3, 'Clipping', 'off');
            hold(ax1, 'off');
            
            % Add subtitle for left topoplot - lower position
            if ~isempty(shared_opt.subtitle1)
                text(0.5, 1.02, shared_opt.subtitle1, ...
                    'Units', 'normalized', 'FontSize', dynamic_subtitle_fontsize, ...
                    'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
            end
            
            % Add colorbar for left topoplot if requested
            if strcmpi(shared_opt.showcolorbars, 'on')
                % Create invisible reference image for colorbar with proper gradient
                ax_ref1 = axes('Parent', fig, 'Position', [0.001 0.001 0.001 0.001], 'Visible', 'off');
                imagesc(ax_ref1, linspace(shared_opt.minmax(1), shared_opt.minmax(2), 100)');
                set(ax_ref1, 'CLim', shared_opt.minmax);
                colormap(ax_ref1, shared_colormap1);
                
                % Create colorbar linked to the reference image
                ax1_pos = get(ax1, 'Position');
                cb1_x = ax1_pos(1) + ax1_pos(3) - 0.02;  % Position to the right of topoplot
                cb1_y = ax1_pos(2) + 0.15;  % Vertically centered
                cb1_width = 0.015;
                cb1_height = ax1_pos(4) - 0.3;
                h_cb1 = colorbar(ax_ref1, 'Position', [cb1_x cb1_y cb1_width cb1_height]);
                
                % Set tick marks (5 ticks: min, mid-low, 0, mid-high, max)
                tick_values1 = [shared_opt.minmax(1), shared_opt.minmax(1)/2, 0, ...
                               shared_opt.minmax(2)/2, shared_opt.minmax(2)];
                tick_labels1 = arrayfun(@(x) sprintf('%.1f', x), tick_values1, 'UniformOutput', false);
                set(h_cb1, 'Ticks', tick_values1, 'TickLabels', tick_labels1, ...
                          'Color', 'k', 'Box', 'on', 'FontSize', dynamic_fontsize * 0.6);
                
                % Hide the reference axes completely
                set(ax_ref1, 'XTick', [], 'YTick', [], 'Box', 'off');
            end
            
            % Right topoplot - centered vertically
            ax2 = axes('Parent', fig, ...
                      'Units', 'normalized', ...
                      'Position', [0.525, 0.025, 0.45, 0.8], ...  % Vertically centered
                      'Color', 'white');
            
            % Second topoplot with colormap2
            limo_topoplot_dual(frame_data2, shared_eloc_locs, ...
                'maplimits', shared_opt.minmax2, ...
                'electrodes', 'off', ...
                'shading', 'interp', ...
                'style', 'both', ...
                'numcontour', 6, ...
                'whitebk', 'on', ...
                'contourlinewidth', topo_linewidth, ...
                'colormap', shared_colormap2, ...
                shared_opt.topoplotopt{:});
            
            axis(ax2, 'equal'); % Enforce circular aspect ratio
            set(ax2, 'Color', 'none');
            
            % Set head outline width for second topoplot
            head_lines = findobj(ax2, 'Type', 'Line', 'Color', 'k');
            set(head_lines, 'LineWidth', head_linewidth);

            % Add circular black outline for second topoplot
            hold(ax2, 'on');
            theta = linspace(0, 2*pi, 360);
            xlims = get(ax2, 'XLim');
            ylims = get(ax2, 'YLim');
            radius = min(diff(xlims), diff(ylims))/2 * 0.93;
            center_x = mean(xlims);
            center_y = mean(ylims);
            circle_x = center_x + radius * cos(theta);
            circle_y = center_y + radius * sin(theta);
            patch(ax2, 'XData', circle_x, 'YData', circle_y, 'ZData', ones(size(circle_x))*10, ...
                'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', head_linewidth / 3, 'Clipping', 'off');
            hold(ax2, 'off');
            
            % Add subtitle for right topoplot - lower position
            if ~isempty(shared_opt.subtitle2)
                ax2_pos = get(ax2, 'Position');
                annotation(fig, 'textbox', [ax2_pos(1), ax2_pos(2) + ax2_pos(4), ax2_pos(3), 0.05], ...
                    'String', shared_opt.subtitle2, ...
                    'FontSize', dynamic_subtitle_fontsize, ...
                    'FontWeight', 'bold', ...
                    'Color', 'k', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom', ...
                    'EdgeColor', 'none', ...
                    'FitBoxToText', 'on');
            end
            
            % Add colorbar for right topoplot if requested
            if strcmpi(shared_opt.showcolorbars, 'on')
                % Create invisible reference image for colorbar with proper gradient
                ax_ref2 = axes('Parent', fig, 'Position', [0.001 0.002 0.001 0.001], 'Visible', 'off');
                imagesc(ax_ref2, linspace(shared_opt.minmax2(1), shared_opt.minmax2(2), 100)');
                set(ax_ref2, 'CLim', shared_opt.minmax2);
                colormap(ax_ref2, shared_colormap2);
                
                % Create colorbar linked to the reference image
                ax2_pos = get(ax2, 'Position');
                cb2_x = ax2_pos(1) + ax2_pos(3) - 0.02;  % Position to the right of topoplot
                cb2_y = ax2_pos(2) + 0.15;  % Vertically centered
                cb2_width = 0.015;
                cb2_height = ax2_pos(4) - 0.3;
                h_cb2 = colorbar(ax_ref2, 'Position', [cb2_x cb2_y cb2_width cb2_height]);
                
                % Set tick marks (5 ticks: min, mid-low, 0, mid-high, max)
                tick_values2 = [shared_opt.minmax2(1), shared_opt.minmax2(1)/2, 0, ...
                               shared_opt.minmax2(2)/2, shared_opt.minmax2(2)];
                tick_labels2 = arrayfun(@(x) sprintf('%.1f', x), tick_values2, 'UniformOutput', false);
                set(h_cb2, 'Ticks', tick_values2, 'TickLabels', tick_labels2, ...
                          'Color', 'k', 'Box', 'on', 'FontSize', dynamic_fontsize * 0.6);
                
                % Hide the reference axes completely
                set(ax_ref2, 'XTick', [], 'YTick', [], 'Box', 'off');
            end
            
            warning(oldWarning);
            
        else
            % SINGLE LAYOUT - Original behavior
            ax = axes('Parent', fig, ...
                      'Units', 'normalized', ...
                      'Position', [0.1 0.1 0.8 0.8], ...
                      'Color', 'white');
            
            % Generate visualization with proper warning suppression
            if strcmpi(shared_opt.mode, '2d')
                % Suppress warning locally
                oldWarning = warning('off', 'all');
                
                limo_topoplot(frame_data, shared_eloc_locs, ...
                    'maplimits', shared_opt.minmax, ...
                    'electrodes', 'off', ...
                    'shading', 'interp', ...
                    'style', 'both', ...
                    'numcontour', 6, ...
                    'whitebk', 'on', ...
                    'contourlinewidth', topo_linewidth, ...
                    shared_opt.topoplotopt{:});

                set(gca, 'Color', 'none');
                set(gcf, 'Color', 'white');
                
                % Set head outline width
                head_lines = findobj(gca, 'Type', 'Line', 'Color', 'k');
                set(head_lines, 'LineWidth', head_linewidth);
                
                % Add circular black outline
                hold(ax, 'on');
                theta = linspace(0, 2*pi, 360);
                
                xlims = get(ax, 'XLim');
                ylims = get(ax, 'YLim');
                radius = min(diff(xlims), diff(ylims))/2 * 0.93;
                
                center_x = mean(xlims);
                center_y = mean(ylims);
                
                circle_x = center_x + radius * cos(theta);
                circle_y = center_y + radius * sin(theta);
                
                patch(ax, 'XData', circle_x, 'YData', circle_y, 'ZData', ones(size(circle_x))*10, ...
                    'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', head_linewidth / 3, 'Clipping', 'off');
                
                hold(ax, 'off');

                warning(oldWarning);
            else
                if isempty(camera_positions) || size(camera_positions, 1) < f
                    error('Camera positions not properly initialized for frame %d', f);
                end
                headplot(frame_data, 'tmp.spl', ...
                    'view', camera_positions(f,:), ...
                    shared_opt.headplotopt{:});
            end
        end
        
        % Add time/frame annotations for BOTH topoplots in dual layout
        if strcmpi(shared_opt.layout, 'dual')
            if strcmpi(shared_opt.framenum, 'on')
                % Add to left topoplot
                text(ax1, 0.02, 0.98, int2str(f), ...
                    'Units', 'normalized', 'FontSize', dynamic_fontsize, ...
                    'VerticalAlignment', 'top', 'Color', 'k');
                % Add to right topoplot
                text(ax2, 0.02, 0.98, int2str(f), ...
                    'Units', 'normalized', 'FontSize', dynamic_fontsize, ...
                    'VerticalAlignment', 'top', 'Color', 'k');
            elseif strcmpi(shared_opt.time, 'on')
                time_str = sprintf('%.3f s', shared_opt.startsec + (original_frame_idx-1)/shared_srate);
                % Add to left topoplot
                text(ax1, 0.02, 0.98, time_str, ...
                    'Units', 'normalized', 'FontSize', dynamic_fontsize, ...
                    'VerticalAlignment', 'top', 'Color', 'k');
                % Add to right topoplot
                text(ax2, 0.02, 0.98, time_str, ...
                    'Units', 'normalized', 'FontSize', dynamic_fontsize, ...
                    'VerticalAlignment', 'top', 'Color', 'k');
            end
        else
            % Original single layout time/frame annotation
            if strcmpi(shared_opt.framenum, 'on')
                text(0.02, 0.98, int2str(f), ...
                    'Units', 'normalized', 'FontSize', dynamic_fontsize, ...
                    'VerticalAlignment', 'top', 'Color', 'k');
            elseif strcmpi(shared_opt.time, 'on')
                text(0.02, 0.98, sprintf('%.3f s', shared_opt.startsec + (original_frame_idx-1)/shared_srate), ...
                    'Units', 'normalized', 'FontSize', dynamic_fontsize, ...
                    'VerticalAlignment', 'top', 'Color', 'k');
            end
        end
        
        % Ensure rendering completes
        drawnow;
        
        % Capture frame
        frame = getframe(fig);
        
        % Only resize if necessary
        if size(frame.cdata, 1) ~= shared_opt.resolution(2) || size(frame.cdata, 2) ~= shared_opt.resolution(1)
            frame.cdata = imresize(frame.cdata, [shared_opt.resolution(2), shared_opt.resolution(1)], 'bicubic');
        end
        
        MovieFrames{f} = frame;
        
        % Clean up
        close(fig);
        
        % Send progress update
        send(dq, f);
        
    catch ME
        warning('Error in frame %d: %s', f, ME.message);
        % Create blank frame with correct resolution
        MovieFrames{f} = struct('cdata', ones(opt.resolution(2), opt.resolution(1), 3, 'uint8')*255, 'colormap', []);
        send(dq, f);
    end
end

fprintf('\nFrame generation complete!\n');

%%%%%%%%%%%%%%%%%%%%%%%%% ASSEMBLY WITH TITLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assembly based on timecourse setting and layout
if strcmpi(opt.timecourse, 'off')
    fprintf('Assembling movie without timecourse...\n');
    
    % Create high-res figure for assembly
    fig_simple = figure('Position', [100 100 opt.resolution(1) opt.resolution(2)], ...
                       'Color', 'white', ...
                       'Units', 'pixels', ...
                       'Resize', 'off');
    
    % Initialize movie
    Movie = struct('cdata', [], 'colormap', []);
    Movie(mframes) = struct('cdata', [], 'colormap', []);
    
    for f = 1:mframes
        % Clear figure
        clf(fig_simple);
        set(fig_simple, 'Color', 'white');
        
        % Add main title at the top - lower position to avoid cutoff
        if ~isempty(opt.title)
            baseline_height = 1000;
            scaling_factor = opt.resolution(2) / baseline_height;
            dynamic_title_fontsize = max(8, round(32 * scaling_factor));

            annotation('textbox', [0 0.92 1 0.05], ...
                      'String', opt.title, ...
                      'HorizontalAlignment', 'center', ...
                      'FontSize', dynamic_title_fontsize, ...
                      'FontWeight', 'bold', ...
                      'EdgeColor', 'none', ...
                      'Color', 'black');
        end
        
        % Create axis for the frame
        ax = axes('Position', [0.01 0.01 0.98 0.98], ...
                  'Color', 'white');
        
        % Display the frame
        if size(MovieFrames{f}.cdata, 1) ~= opt.resolution(2) || size(MovieFrames{f}.cdata, 2) ~= opt.resolution(1)
            resized_frame = imresize(MovieFrames{f}.cdata, [opt.resolution(2), opt.resolution(1)], 'bicubic');
            image(resized_frame);
        else
            image(MovieFrames{f}.cdata);
        end
        axis image;
        axis off;
        
        % Capture with title
        captured_frame = getframe(fig_simple);
        
        % Fix high-DPI display issue
        if size(captured_frame.cdata, 1) ~= opt.resolution(2) || size(captured_frame.cdata, 2) ~= opt.resolution(1)
            captured_frame.cdata = imresize(captured_frame.cdata, [opt.resolution(2), opt.resolution(1)], 'bicubic');
        end
        
        Movie(:,f) = captured_frame;
    end
    close(fig_simple);
    fprintf('Done! Movie assembled.\n');
    
else
    % Complex assembly with timecourse
    fprintf('Assembling movie with timecourse...\n');
    
    % Determine figure size based on layout
    if strcmpi(opt.layout, 'dual')
        % For dual layout with timecourse: vertical stacking
        fig_width = opt.resolution(1);
        fig_height = opt.resolution(2);
    else
        % For single layout with timecourse: original behavior
        if opt.resolution(1) == opt.resolution(2)
            fig_width = opt.resolution(1) + 400;
            fig_height = opt.resolution(2);
        else
            fig_width = opt.resolution(1) + 400;
            fig_height = opt.resolution(2);
        end
    end
    
    % Create high-res figure for final assembly
    fig_final = figure('Position', [100 100 fig_width fig_height], ...
                      'Color', 'white', ...
                      'Units', 'pixels', ...
                      'Resize', 'off');
    
    % Initialize movie
    Movie = struct('cdata', [], 'colormap', []);
    Movie(mframes) = struct('cdata', [], 'colormap', []);
    
    % Determine which channels to display in timecourse
    if ~isempty(opt.timecourse_channels)
        % Use specified channels
        selected_channels = opt.timecourse_channels;
        % Validate indices
        if any(selected_channels < 1) || any(selected_channels > chans)
            error('timecourse_channels contains invalid channel indices');
        end
    elseif ~isempty(opt.timecourse_maxchans)
        % Auto-select evenly spaced channels
        if opt.timecourse_maxchans >= chans
            selected_channels = 1:chans;
        else
            % Evenly space the channels
            step = floor(chans / opt.timecourse_maxchans);
            selected_channels = 1:step:chans;
            % Ensure we have exactly the requested number
            if length(selected_channels) > opt.timecourse_maxchans
                selected_channels = selected_channels(1:opt.timecourse_maxchans);
            end
        end
    else
        % Use all channels (default)
        selected_channels = 1:chans;
    end
    
    num_selected = length(selected_channels);
    selected_data = data(selected_channels, :);
    if isstruct(eloc_locs)
        selected_labels = eloc_locs(selected_channels);
    else
        selected_labels = [];
    end
    
    % Determine layout based on number of channels and topoplot layout
    if strcmpi(opt.layout, 'dual')
        % DUAL LAYOUT WITH TIMECOURSE
        % Adjust positions for dual layout (bottom half)
        base_y = 0.05;
        base_height = 0.40;
    else
        % SINGLE LAYOUT WITH TIMECOURSE
        % Positions for single layout (right side)
        base_y = 0.05;
        base_height = 0.9;
    end
    
    % Create axes based on number of selected channels
    timecourse_axes = [];
    axes_positions = {};
    
    switch num_selected
        case 1
            % Single channel layout - full width
            if strcmpi(opt.layout, 'dual')
                axes_positions = {[0.05, base_y, 0.9, base_height]};
            else
                axes_positions = {[0.75, base_y, 0.2, base_height]};
            end
            
        case 2
            % Side-by-side layout
            if strcmpi(opt.layout, 'dual')
                axes_positions = {[0.05, base_y, 0.43, base_height], ...
                                 [0.52, base_y, 0.43, base_height]};
            else
                % For single layout, stack them vertically
                half_height = base_height * 0.48;
                gap = base_height * 0.04;
                axes_positions = {[0.75, base_y + half_height + gap, 0.2, half_height], ...
                                 [0.75, base_y, 0.2, half_height]};
            end
            
        case 3
            % 2 top, 1 bottom centered
            if strcmpi(opt.layout, 'dual')
                half_height = base_height * 0.48;
                gap = base_height * 0.04;
                axes_positions = {[0.05, base_y + half_height + gap, 0.43, half_height], ...
                                 [0.52, base_y + half_height + gap, 0.43, half_height], ...
                                 [0.26, base_y, 0.43, half_height]};
            else
                % For single layout, stack them vertically
                third_height = base_height * 0.31;
                gap = base_height * 0.035;
                axes_positions = {[0.75, base_y + 2*(third_height + gap), 0.2, third_height], ...
                                 [0.75, base_y + third_height + gap, 0.2, third_height], ...
                                 [0.75, base_y, 0.2, third_height]};
            end
            
        case 4
            % 2x2 grid
            if strcmpi(opt.layout, 'dual')
                half_height = base_height * 0.48;
                gap = base_height * 0.04;
                axes_positions = {[0.05, base_y + half_height + gap, 0.43, half_height], ...
                                 [0.52, base_y + half_height + gap, 0.43, half_height], ...
                                 [0.05, base_y, 0.43, half_height], ...
                                 [0.52, base_y, 0.43, half_height]};
            else
                % For single layout, stack them vertically
                quarter_height = base_height * 0.23;
                gap = base_height * 0.027;
                axes_positions = {[0.75, base_y + 3*(quarter_height + gap), 0.2, quarter_height], ...
                                 [0.75, base_y + 2*(quarter_height + gap), 0.2, quarter_height], ...
                                 [0.75, base_y + quarter_height + gap, 0.2, quarter_height], ...
                                 [0.75, base_y, 0.2, quarter_height]};
            end
            
        otherwise
            % Traditional stacked layout for 5+ channels
            if strcmpi(opt.layout, 'dual')
                axes_positions = {[0.05, base_y, 0.9, base_height]};
            else
                axes_positions = {[0.75, base_y, 0.2, base_height]};
            end
    end
    
    % Create the axes
    for i = 1:length(axes_positions)
        timecourse_axes(i) = axes('Units','Normalized','Position',axes_positions{i}, ...
                                 'Color', 'white');
    end
    
    % Plot data based on number of channels
    initial_xlims = {};
    initial_ylims = {};
    all_limits = [];
    
    if num_selected <= 4
        % Individual plots for each channel (1-4 channels)
        for i = 1:num_selected
            axes(timecourse_axes(i));
            
            % Plot single channel
            channel_data = selected_data(i, :);
            time_vector = (0:length(channel_data)-1) / srate + opt.startsec;
            plot(time_vector, -channel_data, 'k', 'LineWidth', 1);
            
            % Store limits for this axis
            xlim([time_vector(1) time_vector(end)]);
            if strcmpi(opt.timecourse_scaling, 'linked')
                % Calculate common y-limits later
                all_limits = [all_limits; min(-channel_data) max(-channel_data)];
            else
                % Independent scaling
                ylim([min(-channel_data) max(-channel_data)]);
            end
            
            initial_xlims{i} = get(gca, 'XLim');
            initial_ylims{i} = get(gca, 'YLim');
            
            % Customize appearance
            set(gca, 'GridLineStyle', 'none', 'Xgrid', 'off', 'Ygrid', 'on');
            set(gca, 'Color', 'white');
            set(gca, 'XColor', 'black', 'YColor', 'black');
            set(gca, 'TickLabelInterpreter', 'none');
            set(gca, 'XTickLabelRotation', 0);
            
            % Add channel label
            if isstruct(selected_labels)
                channel_label = selected_labels(i).labels;
            else
                channel_label = sprintf('Ch %d', selected_channels(i));
            end
            
            % Position label based on number of channels
            switch num_selected
                case 1
                    ylabel(channel_label, 'FontWeight', 'bold', 'FontSize', 12);
                    xlabel('Time (s)', 'FontSize', 10);
                case {2, 3, 4}
                    title(channel_label, 'FontSize', 10, 'FontWeight', 'bold');
                    if i == num_selected || (num_selected == 3 && i == 3) || ...
                       (num_selected == 4 && (i == 3 || i == 4))
                        xlabel('Time (s)', 'FontSize', 9);
                    end
            end
            
            % Font size adjustment
            if strcmpi(opt.layout, 'dual')
                set(gca, 'FontSize', 10);
            else
                set(gca, 'FontSize', 9);
            end
        end
        
        % Apply linked scaling if requested
        if strcmpi(opt.timecourse_scaling, 'linked') && num_selected > 1
            common_ylim = [min(all_limits(:,1)) max(all_limits(:,2))];
            for i = 1:num_selected
                axes(timecourse_axes(i));
                ylim(common_ylim);
                initial_ylims{i} = common_ylim;
            end
        end
        
    else
        % Traditional multi-channel plot (5+ channels)
        axes(timecourse_axes(1));
        
        % Create temporary file for channel labels
        if isstruct(selected_labels)
            fid = fopen('tmp_file.loc', 'w');
            
            num_channels = length(selected_labels);
            channel_step = ceil(num_channels / 12);
            
            for iChan = 1:num_channels
                if mod(iChan-1, channel_step) == 0
                    fprintf(fid, '0 0 0 %s\n', selected_labels(iChan).labels);
                else
                    fprintf(fid, '0 0 0  \n');
                end
            end
            fclose(fid);
            
            % Plot timecourse
            eegplotold('noui', -selected_data, srate, 0, 'tmp_file.loc', opt.startsec, 'r');
        else
            eegplotold('noui', -selected_data, srate, 0, [], opt.startsec, 'r');
        end
        
        % Get limits and customize appearance
        limits = get(gca, 'Ylim');
        set(gca, 'GridLineStyle', 'none', 'Xgrid', 'off', 'Ygrid', 'on');
        set(gca, 'Color', 'white');
        set(gca, 'XColor', 'black', 'YColor', 'black');
        
        % Ensure tick labels are visible
        set(gca, 'TickLabelInterpreter', 'none');
        set(gca, 'XTickLabelRotation', 0);
        
        % Adjust font size based on layout
        if strcmpi(opt.layout, 'dual')
            set(gca, 'FontSize', 12);
        else
            set(gca, 'FontSize', 10);
        end
        set(gca, 'FontName', 'Helvetica');
        
        % Store initial axis limits
        initial_xlims{1} = get(gca, 'XLim');
        initial_ylims{1} = get(gca, 'YLim');
    end
    
    % Set manual mode for all axes
    for i = 1:length(timecourse_axes)
        axes(timecourse_axes(i));
        set(gca, 'XTickMode', 'manual');
        set(gca, 'XTickLabelMode', 'manual');
    end
    
    % Add title at the top
    if ~isempty(opt.title)
        baseline_height = 1000;
        scaling_factor = opt.resolution(2) / baseline_height;
        dynamic_title_fontsize = max(8, round(32 * scaling_factor));

        if strcmpi(opt.layout, 'dual')
            % For dual layout, title needs to be higher
            annotation('textbox', [0 0.97 1 0.03], ...
                      'String', opt.title, ...
                      'HorizontalAlignment', 'center', ...
                      'FontSize', dynamic_title_fontsize, ...
                      'FontWeight', 'bold', ...
                      'EdgeColor', 'none', ...
                      'Color', 'black');
        else
            annotation('textbox', [0 0.95 1 0.05], ...
                      'String', opt.title, ...
                      'HorizontalAlignment', 'center', ...
                      'FontSize', dynamic_title_fontsize, ...
                      'FontWeight', 'bold', ...
                      'EdgeColor', 'none', ...
                      'Color', 'black');
        end
    end
    
    % Create topoplot axis/axes
    if strcmpi(opt.layout, 'dual')
        % Two topoplots in top half
        axtopoplot1 = axes('Units','Normalized','Position',[0.05 0.52 0.42 0.42], ...
                          'Color', 'white');
        axtopoplot2 = axes('Units','Normalized','Position',[0.53 0.52 0.42 0.42], ...
                          'Color', 'white');
        
        % Add subtitles
        if ~isempty(opt.subtitle1)
            baseline_height = 1000;
            scaling_factor = fig_height / baseline_height;
            dynamic_subtitle_fontsize = max(8, round(subtitle_base_fontsize * scaling_factor));
            
            annotation('textbox', [0.05 0.92 0.42 0.03], ...
                      'String', opt.subtitle1, ...
                      'HorizontalAlignment', 'center', ...
                      'FontSize', dynamic_subtitle_fontsize, ...
                      'FontWeight', 'bold', ...
                      'EdgeColor', 'none', ...
                      'Color', 'black');
        end
        
        if ~isempty(opt.subtitle2)
            annotation('textbox', [0.53 0.92 0.42 0.03], ...
                      'String', opt.subtitle2, ...
                      'HorizontalAlignment', 'center', ...
                      'FontSize', dynamic_subtitle_fontsize, ...
                      'FontWeight', 'bold', ...
                      'EdgeColor', 'none', ...
                      'FitBoxToText', 'on', ...
                      'Color', 'black');
        end
    else
        % Single topoplot (original)
        axtopoplot = axes('Units','Normalized','Position',[0.05 0.05 0.65 0.85], ...
                          'Color', 'white');
    end
    
    % Combine frames with timecourse
    for f = 1:mframes
        indFrame = opt.movieframes(f);
        
        % Calculate time
        x1 = opt.startsec + (indFrame-1)/srate;
        
        % Update timecourse markers for all axes
        for ax_idx = 1:length(timecourse_axes)
            axes(timecourse_axes(ax_idx));
            
            % Remove old red line
            oldLine = findobj(timecourse_axes(ax_idx), 'Type', 'line', 'Color', 'r');
            if ~isempty(oldLine)
                delete(oldLine);
            end
            
            % Get y-limits for this axis
            if num_selected <= 4
                % For individual channel plots
                ylims = initial_ylims{ax_idx};
            else
                % For multi-channel plot
                ylims = limits;
            end
            
            % Add new red marker
            hold on;
            if num_selected <= 4
                % For individual channel plots, use time in seconds
                time_point = x1;
                line([time_point time_point], ylims, 'color', 'r', 'LineWidth', 2);
            else
                % For multi-channel plot, use frame index
                line([indFrame indFrame], ylims, 'color', 'r', 'LineWidth', 2);
            end
            hold off;
            
            % Preserve axis limits
            set(gca, 'XLim', initial_xlims{min(ax_idx, length(initial_xlims))});
            set(gca, 'YLim', initial_ylims{min(ax_idx, length(initial_ylims))});
            
            % Update X-axis label (only for first axis or bottom axes)
            if (num_selected == 1) || ...
               (num_selected == 2 && strcmpi(opt.layout, 'dual')) || ...
               (num_selected == 3 && ax_idx == 3) || ...
               (num_selected == 4 && (ax_idx == 3 || ax_idx == 4)) || ...
               (num_selected > 4 && ax_idx == 1)
                
                set(gca, 'XTickMode', 'manual');
                set(gca, 'XTickLabelMode', 'manual');
                
                if num_selected <= 4
                    % For individual plots, show time value
                    set(gca, 'XTick', time_point);
                    set(gca, 'XTickLabel', sprintf('%.3f s', x1));
                else
                    % For multi-channel plot, use frame index
                    set(gca, 'XTick', indFrame);
                    set(gca, 'XTickLabel', sprintf('%.3f s', x1));
                end
            end
        end
        
        % Force update
        drawnow expose;
        
        % Display topoplot frame(s)
        if strcmpi(opt.layout, 'dual')
            % Process the pre-generated dual frame
            % The MovieFrames{f} already contains both topoplots
            % We need to extract and display them separately
            
            % For dual layout, the frame already has both topoplots
            % We'll display the entire frame in a single axes
            axes(axtopoplot1)
            cla;
            set(gca, 'Color', 'white');
            
            % Extract left half of the frame for first topoplot
            frame_width = size(MovieFrames{f}.cdata, 2);
            frame_height = size(MovieFrames{f}.cdata, 1);
            left_half = MovieFrames{f}.cdata(:, 1:floor(frame_width/2), :);
            
            image(left_half);
            axis image;
            axis off;
            
            axes(axtopoplot2)
            cla;
            set(gca, 'Color', 'white');
            
            % Extract right half of the frame for second topoplot
            right_half = MovieFrames{f}.cdata(:, ceil(frame_width/2):end, :);
            
            image(right_half);
            axis image;
            axis off;
        else
            % Single topoplot (original)
            axes(axtopoplot)
            cla;
            set(gca, 'Color', 'white');
            image(MovieFrames{f}.cdata);
            axis image;
            axis off;
        end
        
        % Capture
        drawnow;
        captured_frame = getframe(fig_final);
        
        % Calculate expected resolution
        if strcmpi(opt.layout, 'dual')
            expected_width = fig_width;
            expected_height = fig_height;
        else
            expected_width = opt.resolution(1) + 400;
            expected_height = opt.resolution(2);
        end
        
        % Ensure frame matches expected resolution
        if size(captured_frame.cdata, 1) ~= expected_height || size(captured_frame.cdata, 2) ~= expected_width
            captured_frame.cdata = imresize(captured_frame.cdata, [expected_height, expected_width], 'bicubic');
        end
        
        Movie(:,f) = captured_frame;
        
        if mod(f, 10) == 0
            fprintf('Assembled %d/%d frames (%.1f%%)\n', f, mframes, (f/mframes)*100);
        end
    end
    close(fig_final);
end

% Re-enable warnings
warning('on', 'MATLAB:handle_graphics:exceptions:SceneNode');

% Cleanup
if exist('tmp_file.loc', 'file')
    delete('tmp_file.loc');
end
if strcmpi(opt.mode, '3d') && strcmpi(opt.layout, 'single') && exist('tmp.spl', 'file')
    delete('tmp.spl');
end

totalTime = toc(startTime);
fprintf('\n=== PERFORMANCE SUMMARY ===\n');
fprintf('Total frames: %d\n', mframes);
fprintf('Workers used: %d\n', pool.NumWorkers);
fprintf('Layout mode: %s\n', opt.layout);
if strcmpi(opt.layout, 'dual')
    fprintf('Dual topoplot subtitles: "%s" | "%s"\n', opt.subtitle1, opt.subtitle2);
    % Determine colormap name for display
    if ~isempty(opt.colormap2)
        colormap2_name = 'custom';
    elseif opt.minmax2(1) < 0 && opt.minmax2(2) > 0
        colormap2_name = 'cool';
    else
        colormap2_name = 'hot';
    end
    fprintf('Colormaps: jet (left) | %s (right)\n', colormap2_name);
    fprintf('Colorbars: %s\n', opt.showcolorbars);
end
fprintf('Total time: %.1f seconds\n', totalTime);
fprintf('Average time per frame: %.2f seconds\n', totalTime/mframes);
fprintf('Frame resolution: %dx%d pixels\n', opt.resolution(1), opt.resolution(2));
fprintf('===========================\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%% Helper function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function camera_positions = calculate_camera_positions(camerapath, movieframes)
    mframes = length(movieframes);
    camera_positions = zeros(mframes, 2);
    
    newpan = mframes + 1;
    posrow = 2;
    if size(camerapath,1) > 1
        newpan = camerapath(posrow,1);
    end
    
    azimuth   = camerapath(1,1);
    az_step   = camerapath(1,2);
    elevation = camerapath(1,3);
    el_step   = camerapath(1,4);
    
    for f = 1:mframes
        indFrame = movieframes(f);
        camera_positions(f,:) = [azimuth, elevation];
        
        if indFrame == newpan
            az_step = camerapath(posrow,2);
            el_step = camerapath(posrow,4);
            posrow = posrow+1;
            if size(camerapath,1)>=posrow
                newpan = camerapath(posrow,1);
            else
                newpan = mframes+1;
            end
        end
        
        azimuth = azimuth+az_step;
        elevation = elevation+el_step;
        elevation = max(-89.99, min(89.99, elevation));
    end
end