function chunk_info = limo_save_boot_chunks(data, chunk_dir, base_name, chunk_start, chunk_size, var_name)
% LIMO_SAVE_BOOT_CHUNKS - Save bootstrap chunks to disk
%
% INPUTS:
%   data       - Bootstrap data for this chunk
%   chunk_dir  - Directory to save chunks
%   base_name  - Base filename for chunks
%   chunk_start- Starting bootstrap index
%   chunk_size - Size of this chunk
%   var_name   - Variable name to use when saving
%
% OUTPUT:
%   chunk_info - Structure with chunk metadata

if nargin < 6
    var_name = 'H0_data';
end

% Create chunk directory if needed
if ~exist(chunk_dir, 'dir')
    mkdir(chunk_dir);
end

% Generate chunk filename
chunk_num = ceil(chunk_start / chunk_size);
chunk_file = fullfile(chunk_dir, sprintf('%s_chunk_%03d.mat', base_name, chunk_num));

% Store chunk metadata
chunk_info.file = chunk_file;
chunk_info.start_idx = chunk_start;
chunk_info.end_idx = chunk_start + size(data, ndims(data)) - 1;
chunk_info.var_name = var_name;
chunk_info.dimensions = size(data);

% Convert to double to ensure consistency
data = double(data);

% Save with the correct variable name
eval([var_name ' = data;']);
save(chunk_file, var_name, '-v7.3');

% Save metadata
meta_file = fullfile(chunk_dir, sprintf('%s_meta.mat', base_name));
if exist(meta_file, 'file')
    load(meta_file, 'chunk_metadata');
    chunk_metadata(chunk_num) = chunk_info;
else
    chunk_metadata = chunk_info;
end
save(meta_file, 'chunk_metadata', '-v7.3');

fprintf('  Saved chunk %d to %s\n', chunk_num, chunk_file);
end