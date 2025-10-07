function [Movie, Colormap] = eegmovie_parallel_dual_optimized(data,srate,eloc_locs,varargin)
% eegmovie_parallel_dual_optimized() - Memory-optimized version with streaming/chunked processing
%                                      Generates and writes frames progressively to prevent memory crashes
%
% Usage:
%   [Movie, Colormap] = eegmovie_parallel_dual_optimized(data, srate, eloc_locs, 'key', value, ...)
%
% Inputs:
%   data       - EEG data matrix (channels x frames)
%   srate      - Sampling rate in Hz
%   eloc_locs  - Electrode locations structure
%
% Optional parameters (same as original plus new ones):
%   'chunk_size'    - Number of frames to process at once (default: 50)
%   'temp_dir'      - Directory for temporary files (default: system temp)
%   'output_file'   - Direct output to video file (skips Movie return)
%   'memory_limit'  - Maximum memory usage in GB (default: auto-detect)
%   [All original parameters still supported]

set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultAxesColor', 'white');

if nargin<1
    help eegmovie_parallel_dual_optimized
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
DEFAULT_SRATE = 256;

if nargin <2, srate = 0; end
if nargin <3, eloc_locs = 0; end

% Verify data and electrode locations compatibility
if ~isequal(eloc_locs, 0) && size(data, 1) ~= length(eloc_locs)
    error('Number of channels in data (%d) does not match number of electrode locations (%d)', size(data, 1), length(eloc_locs));
end

% Parse arguments
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

% Check for layout option early
layout_idx = find(strcmpi(options, 'layout'));
is_dual_layout = false;
if ~isempty(layout_idx) && layout_idx < length(options)
    is_dual_layout = strcmpi(options{layout_idx+1}, 'dual');
end

% Set default resolution - ALWAYS 4000x2000 as required
default_resolution = [4000 2000];

% Parse all options including new optimization parameters
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
                             'layout'      'string'  { 'single' 'dual' } 'single';
                             'data2'       'real'    []                [];
                             'subtitle1'   'string'  {}                '';
                             'subtitle2'   'string'  {}                '';
                             'minmax2'     'real'    {}                0;
                             'colormap1'   'real'    []                [];
                             'colormap2'   'real'    []                [];
                             'showcolorbars' 'string' { 'on' 'off' }   'off';
                             'timecourse_channels' 'real' []    [];
                             'timecourse_maxchans' 'integer' []    [];
                             'timecourse_scaling' 'string' {'linked' 'independent'} 'linked';
                             'headlinewidth' 'real'    {}    [];
                             'topo_linewidth' 'real'    {}    [];
                             'chunk_size'  'integer' {}    50;
                             'temp_dir'    'string'  {}    tempdir;
                             'output_file' 'string'  {}    '';
                             'memory_limit' 'real'   {}    0;
                             'topoplotopt' 'cell'    {}    {};
                             'headplotopt' 'cell'    {}    {} }, 'eegmovie_parallel_dual_optimized');
if ischar(opt), error(opt); end

% CONTINUE IN PART 2