function [success, final_data] = limo_merge_boot_chunks_incremental(chunk_dir, base_name, output_file, varargin)
% LIMO_MERGE_BOOT_CHUNKS_INCREMENTAL - Merge bootstrap chunks with incremental support
%
% INPUTS:
%   chunk_dir    - Directory containing chunks
%   base_name    - Base filename for chunks
%   output_file  - Output filename for merged data
%   varargin     - Optional: 'var_name', string
%                           'delete_chunks', true/false
%                           'incremental', true/false
%                           'existing_bootstraps', N
%
% OUTPUTS:
%   success     - true if merge successful
%   final_data  - Merged data (optional)

% Parse optional inputs
delete_chunks = false;
var_name = '';
incremental = false;
existing_bootstraps = 0;

for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'delete_chunks'
            delete_chunks = varargin{i+1};
        case 'var_name'
            var_name = varargin{i+1};
        case 'incremental'
            incremental = varargin{i+1};
        case 'existing_bootstraps'
            existing_bootstraps = varargin{i+1};
    end
end

success = false;
final_data = [];

try
    % Load metadata
    meta_file = fullfile(chunk_dir, sprintf('%s_meta.mat', base_name));
    if ~exist(meta_file, 'file')
        error('Metadata file not found: %s', meta_file);
    end
    load(meta_file, 'chunk_metadata');

    % Get variable name from first chunk if not specified
    if isempty(var_name)
        var_name = chunk_metadata(1).var_name;
    end

    % Determine dimensions
    n_chunks = length(chunk_metadata);

    % Calculate total bootstraps from chunks
    total_new_bootstraps = 0;
    for i = 1:n_chunks
        chunk_size = chunk_metadata(i).end_idx - chunk_metadata(i).start_idx + 1;
        total_new_bootstraps = total_new_bootstraps + chunk_size;
    end

    % Load first chunk to get dimensions
    first_chunk_file = chunk_metadata(1).file;
    if ~exist(first_chunk_file, 'file')
        error('First chunk file not found: %s', first_chunk_file);
    end
    load(first_chunk_file, var_name);
    eval(['first_data = ' var_name ';']);

    % Get base dimensions (all except last dimension which is bootstraps)
    base_dims = size(first_data);
    n_dims = ndims(first_data);

    % Determine final dimensions
    if incremental && exist(output_file, 'file')
        % Load existing data
        fprintf('Loading existing bootstrap file for incremental merge...\n');
        existing_data = load(output_file, var_name);
        eval(['old_data = existing_data.' var_name ';']);

        % Create extended array
        final_dims = base_dims;
        final_dims(n_dims) = existing_bootstraps + total_new_bootstraps;
        final_data = NaN(final_dims, 'double');

        % Copy existing data
        fprintf('Copying %d existing bootstraps...\n', existing_bootstraps);
        if n_dims == 3
            final_data(:,:,1:existing_bootstraps) = double(old_data);
        elseif n_dims == 4
            final_data(:,:,:,1:existing_bootstraps) = double(old_data);
        elseif n_dims == 5
            final_data(:,:,:,:,1:existing_bootstraps) = double(old_data);
        else
            error('Unsupported dimensionality: %d', n_dims);
        end
        clear old_data existing_data;

        % Set starting index for new data
        merge_start_idx = existing_bootstraps + 1;
    else
        % Create new array
        final_dims = base_dims;
        final_dims(n_dims) = total_new_bootstraps;
        final_data = NaN(final_dims, 'double');
        merge_start_idx = 1;
    end

    % Merge chunks
    fprintf('Merging %d chunks into %s\n', n_chunks, output_file);
    current_idx = merge_start_idx;

    for i = 1:n_chunks
        chunk_file = chunk_metadata(i).file;

        if ~exist(chunk_file, 'file')
            warning('Chunk file not found: %s', chunk_file);
            continue;
        end

        % Load chunk
        load(chunk_file, var_name);
        eval(['chunk_data = ' var_name ';']);
        chunk_bootstraps = size(chunk_data, n_dims);

        % Insert into final array
        end_idx = current_idx + chunk_bootstraps - 1;

        if n_dims == 3
            final_data(:,:,current_idx:end_idx) = double(chunk_data);
        elseif n_dims == 4
            final_data(:,:,:,current_idx:end_idx) = double(chunk_data);
        elseif n_dims == 5
            final_data(:,:,:,:,current_idx:end_idx) = double(chunk_data);
        end

        current_idx = end_idx + 1;
        fprintf('  Merged chunk %d/%d (bootstraps %d-%d)\n', i, n_chunks, ...
            end_idx - chunk_bootstraps + 1, end_idx);

        clear chunk_data;
    end

    % Save merged file with correct variable name
    eval([var_name ' = final_data;']);
    fprintf('Saving merged data to %s with variable name: %s\n', output_file, var_name);

    % Create output directory if needed
    output_dir = fileparts(output_file);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    save(output_file, var_name, '-v7.3');

    % Delete chunks if requested
    if delete_chunks
        for i = 1:n_chunks
            if exist(chunk_metadata(i).file, 'file')
                delete(chunk_metadata(i).file);
            end
        end
        delete(meta_file);
        fprintf('Deleted %d chunk files\n', n_chunks);
    end

    success = true;
    total_bootstraps = size(final_data, n_dims);
    fprintf('Successfully merged %d chunks (%d total bootstraps: %d existing + %d new)\n', ...
        n_chunks, total_bootstraps, existing_bootstraps, total_new_bootstraps);

catch ME
    fprintf('Error merging chunks: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
    success = false;
end
end
