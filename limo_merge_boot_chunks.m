function [success, final_data] = limo_merge_boot_chunks(chunk_dir, base_name, output_file, varargin)
% LIMO_MERGE_BOOT_CHUNKS - Merge bootstrap chunks into final file
%
% INPUTS:
%   chunk_dir    - Directory containing chunks
%   base_name    - Base filename for chunks
%   output_file  - Output filename for merged data
%   varargin     - Optional: 'delete_chunks', true/false
%                           'verify', true/false
%                           'var_name', string
%
% OUTPUTS:
%   success     - true if merge successful
%   final_data  - Merged data (optional)

% Parse optional inputs
delete_chunks = false;
verify = true;
var_name = '';

for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'delete_chunks'
            delete_chunks = varargin{i+1};
        case 'verify'
            verify = varargin{i+1};
        case 'var_name'
            var_name = varargin{i+1};
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
    
    % Determine full dimensions
    n_chunks = length(chunk_metadata);
    first_chunk_file = chunk_metadata(1).file;
    load(first_chunk_file, var_name);
    eval(['first_data = ' var_name ';']);
    
    % Get base dimensions (all except last dimension which is bootstraps)
    base_dims = size(first_data);
    n_dims = ndims(first_data);
    
    % Calculate total bootstraps
    total_bootstraps = 0;
    for i = 1:n_chunks
        if isfield(chunk_metadata, 'end_idx') && isfield(chunk_metadata, 'start_idx')
            chunk_size = chunk_metadata(i).end_idx - chunk_metadata(i).start_idx + 1;
        else
            % Fallback: load chunk to get size
            load(chunk_metadata(i).file, var_name);
            eval(['tmp_data = ' var_name ';']);
            chunk_size = size(tmp_data, n_dims);
            clear tmp_data;
        end
        total_bootstraps = total_bootstraps + chunk_size;
    end
    
    % Pre-allocate final array
    final_dims = base_dims;
    final_dims(n_dims) = total_bootstraps;
    final_data = NaN(final_dims, 'double');
    
    % Merge chunks
    fprintf('Merging %d chunks into %s\n', n_chunks, output_file);
    current_idx = 1;
    
    for i = 1:n_chunks
        chunk_file = fullfile(chunk_dir, sprintf('%s_chunk_%03d.mat', base_name, i));
        
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
            current_idx - chunk_bootstraps, end_idx);
        
        clear chunk_data;
    end
    
    % Verify dimensions
    if verify
        actual_bootstraps = current_idx - 1;
        if actual_bootstraps ~= total_bootstraps
            warning('Expected %d bootstraps, got %d', total_bootstraps, actual_bootstraps);
        end
        
        % Check for NaN slices that shouldn't be there
        if n_dims == 4
            for b = 1:size(final_data, 4)
                if all(isnan(final_data(:,:,:,b)), 'all')
                    warning('Bootstrap %d is all NaN', b);
                end
            end
        elseif n_dims == 5
            for b = 1:size(final_data, 5)
                if all(isnan(final_data(:,:,:,:,b)), 'all')
                    warning('Bootstrap %d is all NaN', b);
                end
            end
        end
    end
    
    % Save merged file with correct variable name
    eval([var_name ' = final_data;']);
    fprintf('Saving merged data to %s with variable name: %s\n', output_file, var_name);
    save(output_file, var_name, '-v7.3');
    
    % Delete chunks if requested
    if delete_chunks
        for i = 1:n_chunks
            chunk_file = fullfile(chunk_dir, sprintf('%s_chunk_%03d.mat', base_name, i));
            if exist(chunk_file, 'file')
                delete(chunk_file);
            end
        end
        delete(meta_file);
        fprintf('Deleted %d chunk files\n', n_chunks);
    end
    
    success = true;
    fprintf('Successfully merged %d chunks (%d bootstraps total)\n', n_chunks, size(final_data, n_dims));
    
catch ME
    fprintf('Error merging chunks: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
    success = false;
end
end