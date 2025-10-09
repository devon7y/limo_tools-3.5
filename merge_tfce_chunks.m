function output_file = merge_tfce_chunks(chunk_dir, output_file)
%% Merge TFCE Chunks - Indexed Insertion Method
% This function merges individual TFCE chunk files into a single file
% using pre-allocation and indexed insertion, mirroring the logic in
% limo_tfce_handling.m's process_H0_in_batches_fixed function.
%
% USAGE:
%   output_file = merge_tfce_chunks(chunk_dir, output_file)
%
% INPUTS:
%   chunk_dir   - Directory containing tfce_batch_*.mat files
%   output_file - (Optional) Path for merged output file
%                 If not provided, will save to parent directory's H0 folder
%
% EXAMPLE:
%   merge_tfce_chunks('/Volumes/T7/ERP Files/Epoched Files 50/paired_ttest4/tfce/tfce_chunks')

% Set default output file if not provided
if nargin < 2 || isempty(output_file)
    parent_dir = fileparts(fileparts(chunk_dir)); % Go up two levels from tfce/tfce_chunks
    output_file = fullfile(parent_dir, 'H0', 'tfce_H0_merged.mat');
end

fprintf('=== TFCE Chunk Merger (Indexed Insertion Method) ===\n\n');
fprintf('Chunk directory: %s\n', chunk_dir);
fprintf('Output file: %s\n\n', output_file);

% Get list of chunk files
chunk_files = dir(fullfile(chunk_dir, 'tfce_batch_*.mat'));
n_chunks = length(chunk_files);

fprintf('Found %d chunk files\n\n', n_chunks);

if n_chunks == 0
    error('No chunk files found in directory!');
end

% Sort files by number to ensure correct order
[~, sort_idx] = sort({chunk_files.name});
chunk_files = chunk_files(sort_idx);

% --- Determine final dimensions and pre-allocate ---
fprintf('Determining final dimensions from chunks...\n');

% Load first chunk to get base dimensions
first_chunk_info = load(fullfile(chunk_dir, chunk_files(1).name), 'batch_tfce');
if ~isfield(first_chunk_info, 'batch_tfce')
    error('Variable "batch_tfce" not found in first chunk. Chunks may be incompatible.');
end
base_dims = size(first_chunk_info.batch_tfce);
nd = length(base_dims);

% Load last chunk to get final index
last_chunk_info = load(fullfile(chunk_dir, chunk_files(end).name), 'output_end_idx');
if ~isfield(last_chunk_info, 'output_end_idx')
    error('Variable "output_end_idx" not found in last chunk. Chunks may be incompatible.');
end
total_bootstraps = last_chunk_info.output_end_idx;

% Define final output dimensions
output_dims = base_dims;
output_dims(nd) = total_bootstraps; % Replace last dim (chunk size) with total size

fprintf('Final matrix dimensions will be: %s\n', mat2str(output_dims));

% Pre-allocate the final matrix
tfce_H0_score = NaN(output_dims, 'double');
fprintf('Pre-allocated final matrix.\n\n');

% --- Merge all chunks using indexed insertion ---
fprintf('Merging chunks:\n');
successful_merges = 0;
for i = 1:n_chunks
    fprintf('  Merging chunk %d/%d... ', i, n_chunks);
    
    try
        chunk = load(fullfile(chunk_dir, chunk_files(i).name));
        
        % Check for required variables
        if ~isfield(chunk, 'batch_tfce') || ~isfield(chunk, 'output_start_idx') || ~isfield(chunk, 'output_end_idx')
            warning('Skipping chunk %d: missing required variables.', i);
            continue;
        end
        
        % Define the index range for this chunk
        idx_range = chunk.output_start_idx:chunk.output_end_idx;
        
        % Insert the chunk data into the pre-allocated matrix
        if nd == 2
            tfce_H0_score(:, idx_range) = chunk.batch_tfce;
        elseif nd == 3
            tfce_H0_score(:, :, idx_range) = chunk.batch_tfce;
        elseif nd == 4
            tfce_H0_score(:, :, :, idx_range) = chunk.batch_tfce;
        else
            error('Unsupported data dimensionality (%d dimensions)', nd);
        end
        
        successful_merges = successful_merges + 1;
        fprintf('Done\n');
        
    catch ME
        warning('Failed to merge chunk %d: %s', i, ME.message);
    end
end

% --- Finalize and Save ---
n_total_valid = size(tfce_H0_score, nd);
% Note: 'all_valid_bootstraps' cannot be reconstructed from chunks alone,
% as this metadata is created during the main processing loop.
% We save the total count based on the final matrix size.
all_valid_bootstraps = 1:n_total_valid;

fprintf('\n=== Merge Summary ===\n');
fprintf('Successfully merged %d/%d chunks.\n', successful_merges, n_chunks);
fprintf('Final TFCE dimensions: %s\n', mat2str(size(tfce_H0_score)));
fprintf('Total valid bootstraps: %d\n', n_total_valid);

% Create output directory if needed
output_dir = fileparts(output_file);
if ~exist(output_dir, 'dir')
    fprintf('\nCreating output directory: %s\n', output_dir);
    mkdir(output_dir);
end

% Save merged result
fprintf('\nSaving merged TFCE result to:\n  %s\n', output_file);
save(output_file, 'tfce_H0_score', '-v7.3');
save(output_file, 'all_valid_bootstraps', 'n_total_valid', '-append');

fprintf('\n✓ Merge complete!\n');
if exist(output_file, 'file')
    file_info = dir(output_file);
    fprintf('Output file size: %.2f MB\n', file_info.bytes / 1024^2);
end

end
