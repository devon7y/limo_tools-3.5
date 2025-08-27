function [Movie, Colormap] = eegmovie_parallel(data,srate,eloc_locs,varargin)
% eegmovie_parallel() - Optimized parallel version of eegmovie
%                       Generates movie frames in parallel for faster processing
%
% Additional parameters beyond standard eegmovie:
%   'resolution' = [width height] in pixels. Default: [1200 1000]

set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultAxesColor', 'white');

if nargin<1
    help eegmovie_parallel
    return
end

% Check if Parallel Computing Toolbox is available
if ~license('test', 'Distrib_Computing_Toolbox')
    warning('Parallel Computing Toolbox not available. Falling back to regular eegmovie.');
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
                             'resolution'  'real'    []                [1200 1000];  % [width height]
                             'topoplotopt' 'cell'    {}    {};
                             'headplotopt' 'cell'    {}    {} }, 'eegmovie_parallel');
if ischar(opt), error(opt); end

% Calculate data limits
if opt.minmax ==0,
    datamin = min(min(data));
    datamax = max(max(data));
    absmax  = max([abs(datamin), abs(datamax)]);
    fudge   = 0.05*(datamax-datamin);
    datamin = -absmax-fudge;
    datamax =  absmax+fudge;
    opt.minmax = [datamin datamax];
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

% Pre-calculate camera positions for 3D mode
camera_positions = [];
if strcmpi(opt.mode, '3d')
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

% Create shared variables for workers
selected_data = data(:, opt.movieframes);
shared_data = selected_data;
shared_eloc_locs = eloc_locs;
shared_opt = opt;
shared_srate = srate;

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
parfor f = 1:mframes
    try
        % Create figure at exactly the target resolution
        fig = figure('Visible', 'off', ...
                    'Units', 'pixels', ...
                    'Position', [0 0 shared_opt.resolution(1) shared_opt.resolution(2)], ...
                    'Color', 'white', ...
                    'PaperPositionMode', 'auto', ...
                    'Resize', 'off');
        
        % --- DYNAMIC SCALING ---
        % Define a baseline resolution and calculate scaling factor
        baseline_height = 1000; % Corresponds to default [1200 1000] resolution
        scaling_factor = shared_opt.resolution(2) / baseline_height;
        
        % Scale font size and line widths, with minimum values
        dynamic_fontsize = max(8, round(28 * scaling_factor));
        head_linewidth = max(0.5, 10 * scaling_factor);
        topo_linewidth = max(0.5, scaling_factor); % Thinner than head outline
        
        % Get frame data
        frame_data = shared_data(:, f);
        original_frame_idx = shared_opt.movieframes(f);
        
        % Create axis centered in the figure
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

            set(gca, 'Color', 'none');  % Make axis transparent
            set(gcf, 'Color', 'white');  % Ensure figure is white
            
            % Explicitly set head outline width, as it has a hardcoded width in topoplot
            head_lines = findobj(gca, 'Type', 'Line', 'Color', 'k');
            set(head_lines, 'LineWidth', head_linewidth);
            
            % Add circular black outline around topoplot
            hold(ax, 'on');
            theta = linspace(0, 2*pi, 360);
            
            xlims = get(ax, 'XLim');
            ylims = get(ax, 'YLim');
            radius = min(diff(xlims), diff(ylims))/2 * 0.93; % Finetuning the radius
            
            center_x = mean(xlims);
            center_y = mean(ylims);
            
            circle_x = center_x + radius * cos(theta);
            circle_y = center_y + radius * sin(theta);
            
            % Draw as a patch on top of everything else
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
        
        % Add time/frame annotations with dynamic font size
        if strcmpi(shared_opt.framenum, 'on')
            text(0.02, 0.98, int2str(f), ...
                'Units', 'normalized', 'FontSize', dynamic_fontsize, ...
                'VerticalAlignment', 'top', 'Color', 'k');
        elseif strcmpi(shared_opt.time, 'on')
            text(0.02, 0.98, sprintf('%.3f s', shared_opt.startsec + (original_frame_idx-1)/shared_srate), ...
                'Units', 'normalized', 'FontSize', dynamic_fontsize, ...
                'VerticalAlignment', 'top', 'Color', 'k');
        end
        
        % Ensure rendering completes
        drawnow;
        
        % Capture frame - simple getframe
        frame = getframe(fig);
        
        % Only resize if necessary (for high-DPI displays)
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

% Assembly based on timecourse setting
if strcmpi(opt.timecourse, 'off')
    fprintf('Assembling movie without timecourse...\n');
    
    % Create high-res figure for assembly with WHITE background
    % Use the exact resolution specified by the user
    fig_simple = figure('Position', [100 100 opt.resolution(1) opt.resolution(2)], ...
                       'Color', 'white', ...
                       'Units', 'pixels', ...
                       'Resize', 'off');
    
    % Initialize movie with proper dimensions
    % moviein is deprecated, use struct array instead
    Movie = struct('cdata', [], 'colormap', []);
    Movie(mframes) = struct('cdata', [], 'colormap', []);
    
    for f = 1:mframes
        % Clear figure
        clf(fig_simple);
        set(fig_simple, 'Color', 'white'); % Ensure white background
        
        % Add title at the top with BLACK text, scaled by resolution
        if ~isempty(opt.title)
            baseline_height = 1000;
            scaling_factor = opt.resolution(2) / baseline_height;
            dynamic_title_fontsize = max(8, round(32 * scaling_factor)); % Base size 16

            annotation('textbox', [0 0.95 1 0.05], ...
                      'String', opt.title, ...
                      'HorizontalAlignment', 'center', ...
                      'FontSize', dynamic_title_fontsize, ...
                      'FontWeight', 'bold', ...
                      'EdgeColor', 'none', ...
                      'Color', 'black');  % BLACK text
        end
        
        % Create axis for the frame with WHITE background
        % For square output, use square axes positioning
        ax = axes('Position', [0.0 0.0 1.0 1.0], ...
                      'Color', 'white');
        
        % Display the frame properly
        % Check if frame needs resizing to fit the figure
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
        
        % Fix high-DPI display issue - simple resize to exact dimensions
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
    
    % Create high-res figure for final assembly with WHITE background
    % Note: with timecourse, width needs extra space
    if opt.resolution(1) == opt.resolution(2)
        % For square output, maintain square aspect for the topoplot portion
        fig_final = figure('Position', [100 100 opt.resolution(1)+400 opt.resolution(2)], ...
                          'Color', 'white', ...
                          'Units', 'pixels', ...
                          'Resize', 'off');
    else
        fig_final = figure('Position', [100 100 opt.resolution(1)+400 opt.resolution(2)], ...
                          'Color', 'white');
    end
    
    % Initialize movie with proper dimensions
    % moviein is deprecated, use struct array instead
    Movie = struct('cdata', [], 'colormap', []);
    Movie(mframes) = struct('cdata', [], 'colormap', []);
    
    % Setup timecourse plot ONCE
    axeegplot = axes('Units','Normalized','Position',[.75 .05 .2 .9], ...
                     'Color', 'white');
    
    % Create temporary file for channel labels - FIX THE LABEL ISSUE
    if isstruct(eloc_locs)
        fid = fopen('tmp_file.loc', 'w');
        
        % Write ALL channels but only show labels for subset
        num_channels = length(eloc_locs);
        channel_step = ceil(num_channels / 12); % Show ~12 labels
        
        for iChan = 1:num_channels
            if mod(iChan-1, channel_step) == 0
                % Write the actual label for display
                fprintf(fid, '0 0 0 %s\n', eloc_locs(iChan).labels);
            else
                % Write a space (not empty) to maintain channel count
                fprintf(fid, '0 0 0  \n'); % Single space
            end
        end
        fclose(fid);
        
        % Plot with dark grey traces
        eegplotold('noui', -data, srate, 0, 'tmp_file.loc', opt.startsec, 'r');
    else
        eegplotold('noui', -data, srate, 0, eloc_locs, opt.startsec, 'r');
    end
    
    % Get limits and customize appearance
    limits = get(axeegplot,'Ylim');
    set(axeegplot,'GridLineStyle','none','Xgrid','off','Ygrid','on');
    set(axeegplot,'Color','white');
    set(axeegplot,'XColor','black','YColor','black');
    
    % Ensure tick labels are always visible
    set(axeegplot, 'TickLabelInterpreter', 'none');
    set(axeegplot, 'XTickLabelRotation', 0);
    set(axeegplot, 'FontSize', 10);
    set(axeegplot, 'FontName', 'Helvetica');
    
    % Set manual mode to prevent automatic clearing
    set(axeegplot, 'XTickMode', 'manual');
    set(axeegplot, 'XTickLabelMode', 'manual');
    
    % Store initial axis limits
    initial_xlim = get(axeegplot, 'XLim');
    initial_ylim = get(axeegplot, 'YLim');
    
    % Add title ONCE at the top, scaled by resolution
    if ~isempty(opt.title)
        baseline_height = 1000;
        scaling_factor = opt.resolution(2) / baseline_height;
        dynamic_title_fontsize = max(8, round(32 * scaling_factor)); % Base size 16

        annotation('textbox', [0 0.95 1 0.05], ...
                  'String', opt.title, ...
                  'HorizontalAlignment', 'center', ...
                  'FontSize', dynamic_title_fontsize, ...
                  'FontWeight', 'bold', ...
                  'EdgeColor', 'none', ...
                  'Color', 'black');
    end
    
    % Create topoplot axis ONCE
    axtopoplot = axes('Units','Normalized','Position',[0.05 0.05 0.65 0.85], ...
                      'Color', 'white');
    initial_xlim = get(axeegplot, 'XLim');
    
    % Combine frames with timecourse
    for f = 1:mframes
        indFrame = opt.movieframes(f);
        
        % Update timecourse marker
        axes(axeegplot)
        
        % Remove old red line only
        oldLine = findobj(axeegplot, 'Type', 'line', 'Color', 'r');
        if ~isempty(oldLine)
            delete(oldLine);
        end
        
        % Calculate time
        x1 = opt.startsec + (indFrame-1)/srate;
        
        % Add new red marker
        hold on;
        l1 = line([indFrame indFrame], limits, 'color', 'r', 'LineWidth', 2);
        hold off;
        
        % Preserve axis limits
        set(axeegplot, 'XLim', initial_xlim);
        set(axeegplot, 'YLim', initial_ylim);
        
        % Update X-axis label - with manual mode
        set(axeegplot, 'XTickMode', 'manual');
        set(axeegplot, 'XTickLabelMode', 'manual');
        set(axeegplot, 'XTick', indFrame);
        set(axeegplot, 'XTickLabel', sprintf('%.3f s', x1));
        
        % Force update
        drawnow;
        
        % Display topoplot frame
        axes(axtopoplot)
        cla;
        set(gca, 'Color', 'white');
        image(MovieFrames{f}.cdata);
        axis image;
        axis off;
        
        % Capture
        drawnow;
        captured_frame = getframe(fig_final);
        
        % Calculate expected resolution for timecourse mode (width has extra 400px)
        expected_width = opt.resolution(1) + 400;
        expected_height = opt.resolution(2);
        
        % Ensure frame matches expected resolution (fixes high-DPI display issue)
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
if strcmpi(opt.mode, '3d') && exist('tmp.spl', 'file')
    delete('tmp.spl');
end

totalTime = toc(startTime);
fprintf('\n=== PERFORMANCE SUMMARY ===\n');
fprintf('Total frames: %d\n', mframes);
fprintf('Workers used: %d\n', pool.NumWorkers);
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
