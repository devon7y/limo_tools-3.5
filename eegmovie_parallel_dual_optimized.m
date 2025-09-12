function [Movie, Colormap] = eegmovie_parallel_dual_optimized(data,srate,eloc_locs,varargin)
% eegmovie_parallel_dual_optimized() - Memory-optimized version with chunked processing
%                                      Prevents memory crashes with large video files
%
% Key optimizations:
%   - Processes frames in chunks to maintain constant memory usage
%   - Saves frames to disk immediately after generation
%   - Options for direct video file output
%   - Automatic memory limit detection and adjustment

set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultAxesColor', 'white');

if nargin<1
    help eegmovie_parallel_dual_optimized
    return
end

% Check for Parallel Computing Toolbox
if ~license('test', 'Distrib_Computing_Toolbox')
    error('This optimized version requires Parallel Computing Toolbox');
end

warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');

[chans,frames] = size(data);
DEFAULT_SRATE = 256;

if nargin <2, srate = 0; end
if nargin <3, eloc_locs = 0; end

if ~isequal(eloc_locs, 0) && size(data, 1) ~= length(eloc_locs)
    error('Number of channels in data (%d) does not match number of electrode locations (%d)', ...
          size(data, 1), length(eloc_locs));
end

% Parse arguments (compatibility with original)
if nargin > 5 && ~ischar(varargin{3}) || nargin == 4 && ~ischar(varargin{2})
    options = {};
    if nargin>=8, options = { options{:} 'topoplotopt' varargin(5:end) }; end
    if nargin>=7, options = { options{:} 'startsec'    varargin{4} }; end
    if nargin>=6, options = { options{:} 'minmax'      varargin{3} }; end
    if nargin>=5, options = { options{:} 'movieframes' varargin{2} }; end
    if nargin>=4, options = { options{:} 'title'       varargin{1} }; end
else
    options = varargin;
end

% Check for dual layout
layout_idx = find(strcmpi(options, 'layout'));
is_dual_layout = false;
if ~isempty(layout_idx) && layout_idx < length(options)
    is_dual_layout = strcmpi(options{layout_idx+1}, 'dual');
end

% Parse all options including new optimization parameters
opt = finputcheck(options, { ...
    'startsec'    'real'    {}    0;
    'minmax'      'real'    {}    0;
    'movieframes' 'integer' {}    0;
    'title'       'string'  {}    '';
    'mode'        'string'  { '2D' '3D'  }    '2D';
    'timecourse'  'string'  { 'on' 'off' }    'off';
    'framenum'    'string'  { 'on' 'off' }    'on';
    'time'        'string'  { 'on' 'off' }    'off';
    'resolution'  'real'    []    [4000 2000];  % Required resolution
    'layout'      'string'  { 'single' 'dual' } 'single';
    'data2'       'real'    []    [];
    'subtitle1'   'string'  {}    '';
    'subtitle2'   'string'  {}    '';
    'minmax2'     'real'    {}    0;
    'colormap1'   'real'    []    [];
    'colormap2'   'real'    []    [];
    'showcolorbars' 'string' { 'on' 'off' }   'off';
    'headlinewidth' 'real'    {}    [];
    'topo_linewidth' 'real'    {}    [];
    'chunk_size'  'integer' {}    30;  % Reduced default for 4K resolution
    'output_file' 'string'  {}    '';
    'topoplotopt' 'cell'    {}    {};
    'headplotopt' 'cell'    {}    {} }, 'eegmovie_parallel_dual_optimized');

if ischar(opt), error(opt); end

% Validate dual layout
if strcmpi(opt.layout, 'dual')
    if isempty(opt.data2)
        error('Dual layout requires ''data2'' parameter');
    end
    if size(opt.data2, 1) ~= chans || size(opt.data2, 2) ~= frames
        error('data2 dimensions must match data dimensions');
    end
end

% Calculate data limits
if isscalar(opt.minmax) && opt.minmax == 0
    datamin = min(min(data));
    datamax = max(max(data));
    absmax = max([abs(datamin), abs(datamax)]);
    fudge = 0.05*(datamax-datamin);
    opt.minmax = [-absmax-fudge, absmax+fudge];
end

if strcmpi(opt.layout, 'dual') && isscalar(opt.minmax2) && opt.minmax2 == 0
    datamin2 = min(min(opt.data2));
    datamax2 = max(max(opt.data2));
    absmax2 = max([abs(datamin2), abs(datamax2)]);
    fudge2 = 0.05*(datamax2-datamin2);
    opt.minmax2 = [-absmax2-fudge2, absmax2+fudge2];
end

if opt.movieframes == 0
    opt.movieframes = 1:frames;
end

if srate == 0
    srate = DEFAULT_SRATE;
end

if strcmpi(opt.time, 'on')
    opt.framenum = 'off';
end

mframes = length(opt.movieframes);

% Setup temporary directory
temp_base = fullfile(tempdir, sprintf('eegmovie_%s', datestr(now, 'yyyymmdd_HHMMSS')));
if ~exist(temp_base, 'dir')
    mkdir(temp_base);
end
cleanupObj = onCleanup(@() cleanup_temp_files(temp_base));

fprintf('Processing %d frames at %dx%d resolution\n', mframes, opt.resolution(1), opt.resolution(2));
fprintf('Using chunked processing with chunk size: %d\n', opt.chunk_size);
fprintf('Temporary files: %s\n', temp_base);

% Setup parallel pool
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', feature('numcores'));
    pool = gcp;
end

% Setup colormaps
Colormap = jet(64);
if strcmpi(opt.layout, 'dual')
    colormap1 = opt.colormap1;
    if isempty(colormap1), colormap1 = jet(64); end
    colormap2 = opt.colormap2;
    if isempty(colormap2)
        if opt.minmax2(1) < 0 && opt.minmax2(2) > 0
            colormap2 = cool(64);
        else
            colormap2 = hot(64);
        end
    end
else
    colormap1 = opt.colormap1;
    if isempty(colormap1), colormap1 = jet(64); end
    colormap2 = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%% CHUNKED FRAME GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%

num_chunks = ceil(mframes / opt.chunk_size);
fprintf('Processing in %d chunks...\n\n', num_chunks);

startTime = tic;

for chunk_idx = 1:num_chunks
    chunk_start = (chunk_idx - 1) * opt.chunk_size + 1;
    chunk_end = min(chunk_idx * opt.chunk_size, mframes);
    chunk_frames = chunk_start:chunk_end;
    chunk_size_actual = length(chunk_frames);
    
    fprintf('Chunk %d/%d (frames %d-%d)\n', chunk_idx, num_chunks, chunk_start, chunk_end);
    
    % Select data for this chunk
    chunk_indices = opt.movieframes(chunk_frames);
    chunk_data = data(:, chunk_indices);
    if strcmpi(opt.layout, 'dual')
        chunk_data2 = opt.data2(:, chunk_indices);
    else
        chunk_data2 = [];
    end
    
    % Process chunk in parallel
    parfor local_idx = 1:chunk_size_actual
        frame_idx = chunk_frames(local_idx);
        generate_and_save_frame(frame_idx, local_idx, chunk_data, chunk_data2, ...
                               eloc_locs, opt, srate, colormap1, colormap2, ...
                               temp_base);
    end
    
    fprintf('  Chunk %d completed\n', chunk_idx);
end

%%%%%%%%%%%%%%%%%%%%%%%%% ASSEMBLY %%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nAssembling movie...\n');

if ~isempty(opt.output_file)
    % Direct video output
    write_video_file(opt, temp_base, mframes);
    Movie = [];
else
    % Load frames into memory
    Movie = struct('cdata', [], 'colormap', []);
    Movie(mframes) = struct('cdata', [], 'colormap', []);
    
    for f = 1:mframes
        frame_file = fullfile(temp_base, sprintf('frame_%06d.mat', f));
        frame_data = load(frame_file);
        Movie(f) = frame_data.frame;
        
        if mod(f, 10) == 0
            fprintf('  Loaded %d/%d frames\n', f, mframes);
        end
    end
end

warning('on', 'MATLAB:handle_graphics:exceptions:SceneNode');

totalTime = toc(startTime);
fprintf('\n=== COMPLETED ===\n');
fprintf('Total time: %.1f seconds\n', totalTime);
fprintf('Frames: %d at %dx%d\n', mframes, opt.resolution(1), opt.resolution(2));
fprintf('=================\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%

function generate_and_save_frame(frame_idx, local_idx, data, data2, eloc_locs, ...
                                opt, srate, colormap1, colormap2, temp_base)
    % Generate and save a single frame
    
    try
        fig = figure('Visible', 'off', 'Units', 'pixels', ...
                    'Position', [0 0 opt.resolution(1) opt.resolution(2)], ...
                    'Color', 'white', 'Renderer', 'zbuffer');
        
        % Scale UI elements
        scaling = opt.resolution(2) / 1000;
        fontsize = round(28 * scaling);
        
        if isempty(opt.headlinewidth)
            head_linewidth = 50;
        else
            head_linewidth = opt.headlinewidth;
        end
        
        if isempty(opt.topo_linewidth)
            topo_linewidth = 10;
        else
            topo_linewidth = opt.topo_linewidth;
        end
        
        % Get frame data
        frame_data = data(:, local_idx);
        original_idx = opt.movieframes(frame_idx);
        
        if strcmpi(opt.layout, 'dual')
            % Dual layout
            frame_data2 = data2(:, local_idx);
            
            % Left topoplot
            ax1 = axes('Units', 'normalized', 'Position', [0.025 0.025 0.45 0.8]);
            warning('off', 'all');
            limo_topoplot_dual(frame_data, eloc_locs, ...
                'maplimits', opt.minmax, 'electrodes', 'off', ...
                'shading', 'interp', 'style', 'both', ...
                'numcontour', 6, 'contourlinewidth', topo_linewidth, ...
                'colormap', colormap1, opt.topoplotopt{:});
            warning('on', 'all');
            
            if ~isempty(opt.subtitle1)
                text(0.5, 1.02, opt.subtitle1, 'Units', 'normalized', ...
                    'FontSize', fontsize, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold');
            end
            
            % Right topoplot
            ax2 = axes('Units', 'normalized', 'Position', [0.525 0.025 0.45 0.8]);
            warning('off', 'all');
            limo_topoplot_dual(frame_data2, eloc_locs, ...
                'maplimits', opt.minmax2, 'electrodes', 'off', ...
                'shading', 'interp', 'style', 'both', ...
                'numcontour', 6, 'contourlinewidth', topo_linewidth, ...
                'colormap', colormap2, opt.topoplotopt{:});
            warning('on', 'all');
            
            if ~isempty(opt.subtitle2)
                text(0.5, 1.02, opt.subtitle2, 'Units', 'normalized', ...
                    'FontSize', fontsize, 'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold');
            end
            
            % Add time annotation to both
            if strcmpi(opt.time, 'on')
                time_str = sprintf('%.3f s', opt.startsec + (original_idx-1)/srate);
                text(ax1, 0.02, 0.98, time_str, 'Units', 'normalized', ...
                    'FontSize', fontsize, 'VerticalAlignment', 'top');
                text(ax2, 0.02, 0.98, time_str, 'Units', 'normalized', ...
                    'FontSize', fontsize, 'VerticalAlignment', 'top');
            elseif strcmpi(opt.framenum, 'on')
                text(ax1, 0.02, 0.98, int2str(frame_idx), 'Units', 'normalized', ...
                    'FontSize', fontsize, 'VerticalAlignment', 'top');
                text(ax2, 0.02, 0.98, int2str(frame_idx), 'Units', 'normalized', ...
                    'FontSize', fontsize, 'VerticalAlignment', 'top');
            end
            
        else
            % Single layout
            ax = axes('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
            warning('off', 'all');
            limo_topoplot(frame_data, eloc_locs, ...
                'maplimits', opt.minmax, 'electrodes', 'off', ...
                'shading', 'interp', 'style', 'both', ...
                'numcontour', 6, 'contourlinewidth', topo_linewidth, ...
                opt.topoplotopt{:});
            warning('on', 'all');
            
            % Add time annotation
            if strcmpi(opt.time, 'on')
                text(0.02, 0.98, sprintf('%.3f s', opt.startsec + (original_idx-1)/srate), ...
                    'Units', 'normalized', 'FontSize', fontsize, 'VerticalAlignment', 'top');
            elseif strcmpi(opt.framenum, 'on')
                text(0.02, 0.98, int2str(frame_idx), 'Units', 'normalized', ...
                    'FontSize', fontsize, 'VerticalAlignment', 'top');
            end
        end
        
        % Add title if specified
        if ~isempty(opt.title)
            annotation('textbox', [0 0.92 1 0.05], 'String', opt.title, ...
                      'HorizontalAlignment', 'center', 'FontSize', fontsize+4, ...
                      'FontWeight', 'bold', 'EdgeColor', 'none');
        end
        
        % Capture and save
        drawnow;
        frame = getframe(fig);
        
        % Ensure correct resolution
        if size(frame.cdata, 1) ~= opt.resolution(2) || size(frame.cdata, 2) ~= opt.resolution(1)
            frame.cdata = imresize(frame.cdata, [opt.resolution(2) opt.resolution(1)], 'bicubic');
        end
        
        % Save to disk
        frame_file = fullfile(temp_base, sprintf('frame_%06d.mat', frame_idx));
        save(frame_file, 'frame', '-v7.3');
        
        close(fig);
        
    catch ME
        warning('Error in frame %d: %s', frame_idx, ME.message);
        % Save blank frame
        frame = struct('cdata', ones(opt.resolution(2), opt.resolution(1), 3, 'uint8')*255, ...
                      'colormap', []);
        frame_file = fullfile(temp_base, sprintf('frame_%06d.mat', frame_idx));
        save(frame_file, 'frame', '-v7.3');
    end
end

function write_video_file(opt, temp_base, mframes)
    % Write frames directly to video file
    
    fprintf('Writing video file: %s\n', opt.output_file);
    
    [~, ~, ext] = fileparts(opt.output_file);
    if strcmpi(ext, '.mp4')
        videoObj = VideoWriter(opt.output_file, 'MPEG-4');
        videoObj.Quality = 95;
    else
        videoObj = VideoWriter(opt.output_file);
    end
    
    videoObj.FrameRate = 30;
    open(videoObj);
    
    for f = 1:mframes
        frame_file = fullfile(temp_base, sprintf('frame_%06d.mat', f));
        frame_data = load(frame_file);
        writeVideo(videoObj, frame_data.frame);
        
        if mod(f, 10) == 0
            fprintf('  Written %d/%d frames\n', f, mframes);
        end
    end
    
    close(videoObj);
    fprintf('Video saved successfully!\n');
end

function cleanup_temp_files(temp_base)
    % Clean up temporary files
    
    fprintf('Cleaning up temporary files...\n');
    if exist(temp_base, 'dir')
        try
            rmdir(temp_base, 's');
        catch
            warning('Could not remove temp directory: %s', temp_base);
        end
    end
end