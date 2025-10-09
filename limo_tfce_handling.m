function [tfce_score,thresholded_maps] = limo_tfce_handling(varargin)

% routine to create tfce files commensurate to boostrapped files
% MEMORY-EFFICIENT VERSION with DIMENSION FIX - processes bootstraps in batches
%
% FORMAT  [tfce_score,thresholded_maps] = limo_tfce_handling(filename,'checkfile','yes')
%         [tfce_score,thresholded_maps] = limo_tfce_handling(filename,'batch_size',100,'temp_dir','/path/to/temp')
%
% INPUTS filename is the stat file that need to be tfced (if H0 exist it is done too)
%        'checkfile' is 'yes' by default - if 'no' and tfce files already exist,
%                     it overwrites without asking otherwise user is prompted
%        'batch_size' number of bootstraps to process at once (default: 50)
%        'temp_dir' directory for temporary files (default: a 'tfce_chunks' subfolder in the 'tfce' directory)
%        'max_workers' maximum parallel workers (default: 4)
%        'incremental' 'yes' (default) or 'no' - enable incremental processing to add clusters to existing distributions
%        'keep_chunks' 'yes' (default) or 'no' - keep temporary chunk files after processing
%
% OUTPUTS tfce_* files are saved on the drive in a tfce folder
%         H0_tfce_* files are saved on the drive in the H0 folder
%         tfce_score is the tfce_score for the data of the file in
%         thresholded_maps is all the thresholded maps i.e. for each dh
%                          value in sum(extent(h)^E*height^H*dh)
%         if a boostrap file exist, thresholded_maps is a 1 * 2 cell array
%         with {1} the maps for the observed data and {2} a cell of cells
%         with the maps of each boostrap
% ------------------------------------------------------------------
% Memory-Efficient Version with Dimension Fix - Cyril R. Pernet & Assistant - 2025

%% Check inputs

file = varargin{1};
[filepath,filename,ext] = fileparts(file);
if isempty(filepath); filepath = pwd; end
if isempty(ext)
    file = dir(fullfile(filepath,[filename '*']));
    filename = file.name;
end

if exist(fullfile(filepath,[filename ext]),'file')
    if ~exist(fullfile(filepath,'LIMO.mat'),'file')
        error('no LIMO.mat found next to %s',filename)
    else
        filename = [filename ext];
        LIMO = load(fullfile(filepath,'LIMO.mat'));
        LIMO = LIMO.(cell2mat(fieldnames(LIMO)));
    end
else
    error('can''t find %s',varargin{1})
end

% Default parameters
checkfile = 'yes';
batch_size = 192;  % Process 192 bootstraps at a time (optimized for 32 workers)
temp_dir = fullfile(LIMO.dir, 'tfce', 'tfce_chunks');  % Default directory for temp files
incremental = 'yes';  % Enable incremental processing by default
keep_chunks = 'yes';  % Keep chunk files after merging by default
N = getenv('NUMBER_OF_PROCESSORS');
if isempty(N)
    N = feature('numcores'); % Fallback to physical cores
else
    N = str2double(N);
end
max_workers = N;  % Set workers to the number of logical cores

% Parse additional inputs
for i=2:2:nargin
    if contains(varargin{i},'checkfile','IgnoreCase',true)
        checkfile = varargin{i+1};
    elseif contains(varargin{i},'batch_size','IgnoreCase',true)
        batch_size = varargin{i+1};
    elseif contains(varargin{i},'temp_dir','IgnoreCase',true)
        temp_dir = varargin{i+1};
    elseif contains(varargin{i},'max_workers','IgnoreCase',true)
        max_workers = varargin{i+1};
    elseif contains(varargin{i},'incremental','IgnoreCase',true)
        incremental = varargin{i+1};
    elseif contains(varargin{i},'keep_chunks','IgnoreCase',true)
        keep_chunks = varargin{i+1};
    end
end

% Create temp directory if it doesn't exist
if ~exist(temp_dir,'dir')
    mkdir(temp_dir);
end

% Display memory-saving settings
fprintf('\n=== MEMORY-EFFICIENT TFCE SETTINGS ===\n');
fprintf('Batch size: %d bootstraps\n', batch_size);
fprintf('Temp directory: %s\n', temp_dir);
fprintf('Max parallel workers: %d\n', max_workers);
fprintf('Incremental processing: %s\n', incremental);
fprintf('Keep chunk files: %s\n', keep_chunks);
fprintf('=====================================\n\n');

%% quick user check
% ----------------

% files to create
tfce_file    = fullfile(LIMO.dir,['tfce' filesep 'tfce_' filename]);
H0_tfce_file = fullfile(LIMO.dir,['H0' filesep 'tfce_H0_' filename]);
% given filename input, we expect H0 to be
H0filename   = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);

if strcmpi(checkfile,'yes')
    if exist(tfce_file,'file')
        % Check if running in headless mode (no display available)
        try
            feature('ShowFigureWindows');
            is_headless = false;
        catch
            is_headless = true;
        end

        if is_headless || ~usejava('desktop')
            % Headless mode: automatically overwrite
            fprintf('TFCE file exists. Running in headless mode - automatically overwriting.\n');
            answer = 'Yes';
        else
            % Interactive mode: ask user
            answer = questdlg('tfce file already exist - overwrite?','data check','Yes','No','Yes');
        end

        if strcmp(answer,'Yes')
            LIMO.design.tfce = 1;
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
        else
            LIMO.design.tfce = 0;
            return
        end
    else
        LIMO.design.tfce = 1;
        if exist(LIMO.dir,'dir')
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
            if ~exist(fullfile(LIMO.dir,'tfce'),'dir')
                mkdir(fullfile(LIMO.dir,'tfce'));
            end
        else
            save(fullfile(pwd,'LIMO.mat'),'LIMO','-v7.3')
            if ~exist(fullfile(pwd,'tfce'),'dir')
                mkdir(fullfile(pwd,'tfce'));
            end
        end
    end
end

if isfield(LIMO.design,'bootstrap')
    nboot = LIMO.design.bootstrap;
end

% check if there is a neighbouring matrix
% (since TFCE integrates over clusters)
if ~isfield(LIMO.data,'neighbouring_matrix')
    warning('no neighbouring matrix found, this is required for TFCE')
    [~, LIMO.data.neighbouring_matrix] = limo_expected_chanlocs;
    if isempty(LIMO.data.neighbouring_matrix)
        return
    else
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
    end
end

% create tfce folder
if ~exist(fullfile(LIMO.dir,'tfce'),'dir')
    mkdir(fullfile(LIMO.dir,'tfce'))
end
fprintf('Thresholding %s using TFCE \n',filename);

% Set up parallel pool with limited workers
current_pool = gcp('nocreate');
if isempty(current_pool)
    parpool('local', max_workers);
elseif current_pool.NumWorkers > max_workers
    delete(current_pool);
    parpool('local', max_workers);
end

% -------------------------------------------------------------------------------
if contains(filename,'R2') || ...
        contains(filename,'semi_partial') % these files last dimension is R2, F, p
    % ----------------------------------------------------------------------------
    R2 = load(fullfile(LIMO.dir,'R2.mat'));
    R2 = R2.(cell2mat(fieldnames(R2)));
    if size(R2,1) == 1
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            [tfce_score(1,:,:),thresholded_maps] = limo_tfce(2, squeeze(R2(:,:,:,2)),[]); % no neighbouring, time-freq cluster
        else
            [tfce_score(1,:),thresholded_maps]   = limo_tfce(1, squeeze(R2(:,:,2)),LIMO.data.neighbouring_matrix);
        end
    else
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            [tfce_score,thresholded_maps] = limo_tfce(3, squeeze(R2(:,:,:,2)),LIMO.data.neighbouring_matrix);
        else
            [tfce_score,thresholded_maps] = limo_tfce(2, squeeze(R2(:,:,2)),LIMO.data.neighbouring_matrix);
        end
    end
    save(tfce_file,'tfce_score','-v7.3'); clear R2 ;
    
    if exist(H0filename,'file')
        fprintf('Applying TFCE to null data using batch processing... \n')
        process_H0_in_batches_fixed(H0filename, H0_tfce_file, LIMO, batch_size, temp_dir, 'R2', incremental, keep_chunks);
        
        % Handle thresholded maps
        tmp                 = thresholded_maps;     clear thresholded_maps;
        thresholded_maps{1} = tmp;                  clear tmp
        thresholded_maps{2} = [];  % Too memory intensive to store all bootstrap thresholded maps
    end
    
    % -------------------------------------------------------------------------------------------------------------
elseif contains(filename,'con') && ~contains(filename,'conditions') || ... % FIX: exclude ANOVA condition files
        (contains(LIMO.design.name,'One sample','IgnoreCase',true) || ...
        contains(LIMO.design.name,'Two samples','IgnoreCase',true) || ...
        contains(LIMO.design.name,'Paired','IgnoreCase',true))  % these file last dimension is mean, se, df, t and p
    % --------------------------------------------------------------------------------------------------------------
    tval = load(filename);
    tval = tval.(cell2mat(fieldnames(tval)));
    
    % FIX: Check if this is actually a t-test file structure
    if ndims(tval) >= 3 && size(tval, ndims(tval)) < 4
        % This looks like an F-value file (F, p), not a t-test file
        % Redirect to F-value handling section
        Fval = tval;
        if size(Fval,1) == 1
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                [tfce_score(1,:,:),thresholded_maps] = limo_tfce(2,squeeze(Fval(:,:,:,1)),[]);
            else
                [tfce_score(1,:),thresholded_maps] = limo_tfce(1,squeeze(Fval(:,:,1)),LIMO.data.neighbouring_matrix);
            end
        else
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                [tfce_score,thresholded_maps] =  limo_tfce(3,squeeze(Fval(:,:,:,1)),LIMO.data.neighbouring_matrix);
            else
                [tfce_score,thresholded_maps] =  limo_tfce(2,squeeze(Fval(:,:,1)),LIMO.data.neighbouring_matrix);
            end
        end
    else
        % Normal t-test file processing
        if size(tval,1) == 1
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                [tfce_score(1,:,:),thresholded_maps] = limo_tfce(2, squeeze(tval(:,:,:,end-1)),[]); % no neighbouring, time-freq cluster
            else
                [tfce_score(1,:),thresholded_maps]   = limo_tfce(1, squeeze(tval(:,:,end-1)),LIMO.data.neighbouring_matrix);
            end
        else
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                [tfce_score,thresholded_maps] = limo_tfce(3, squeeze(tval(:,:,:,end-1)),LIMO.data.neighbouring_matrix);
            else
                [tfce_score,thresholded_maps] = limo_tfce(2, squeeze(tval(:,:,4)),LIMO.data.neighbouring_matrix);
            end
        end
    end
    save(tfce_file,'tfce_score','-v7.3'); clear tval ;
    
    if exist(H0filename,'file')
        fprintf('Applying TFCE to null data using batch processing... \n')
        process_H0_in_batches_fixed(H0filename, H0_tfce_file, LIMO, batch_size, temp_dir, 'ttest', incremental, keep_chunks);
        
        % Handle thresholded maps
        tmp                 = thresholded_maps;     clear thresholded_maps;
        thresholded_maps{1} = tmp;                  clear tmp
        thresholded_maps{2} = [];  % Too memory intensive to store all bootstrap thresholded maps
    end
    
    % ------------------------------------------
else % anything else last dimension is F and p (including ANOVA files)
    % ------------------------------------------
    Fval = load(filename);
    Fval = Fval.(cell2mat(fieldnames(Fval)));
    if contains(filename,'ess')
        Fval = Fval(:,:,end-1:end);
    end

    if size(Fval,1) == 1
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            [tfce_score(1,:,:),thresholded_maps] = limo_tfce(2,squeeze(Fval(:,:,:,1)),[]);
        else
            [tfce_score(1,:),thresholded_maps] = limo_tfce(1,squeeze(Fval(:,:,1)),LIMO.data.neighbouring_matrix);
        end
    else
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            [tfce_score,thresholded_maps] =  limo_tfce(3,squeeze(Fval(:,:,:,1)),LIMO.data.neighbouring_matrix);
        else
            [tfce_score,thresholded_maps] =  limo_tfce(2,squeeze(Fval(:,:,1)),LIMO.data.neighbouring_matrix);
        end
    end
    save(tfce_file,'tfce_score','-v7.3'); clear Fval;
    
    if exist(H0filename,'file')
        fprintf('Applying TFCE to null data using batch processing... \n')
        process_H0_in_batches_fixed(H0filename, H0_tfce_file, LIMO, batch_size, temp_dir, 'F', incremental, keep_chunks);
        
        % Handle thresholded maps
        tmp                 = thresholded_maps;     clear thresholded_maps;
        thresholded_maps{1} = tmp;                  clear tmp
        thresholded_maps{2} = [];  % Too memory intensive to store all bootstrap thresholded maps
    end
end

% Clean up temp directory (optional)
if strcmpi(keep_chunks, 'no')
    fprintf('Cleaning up temporary chunk files...\n');
    temp_pattern = fullfile(temp_dir, 'tfce_batch_*.mat');
    temp_files = dir(temp_pattern);
    for i = 1:length(temp_files)
        delete(fullfile(temp_files(i).folder, temp_files(i).name));
    end
else
    fprintf('Keeping temporary chunk files in: %s\n', temp_dir);
end

end

%% Helper function for batch processing of H0 data with DIMENSION FIX

function process_H0_in_batches_fixed(H0filename, H0_tfce_file, LIMO, batch_size, temp_dir, test_type, incremental, keep_chunks)
    
    % Use matfile to access H0 data without loading it all
    fprintf('Opening H0 file using memory-mapped access...\n');
    
    % First check if file exists and is readable
    if ~exist(H0filename, 'file')
        error('H0 file does not exist: %s', H0filename);
    end
    
    try
        m = matfile(H0filename);
    catch ME
        error('Cannot open H0 file: %s\nError: %s', H0filename, ME.message);
    end
    
    % Get variable name and dimensions
    vars = who(m);
    if isempty(vars)
        error('H0 file appears to be empty or corrupted: %s', H0filename);
    end
    
    H0_varname = vars{1};
    fprintf('H0 variable name: %s\n', H0_varname);
    
    % Get dimensions safely
    try
        H0_dims = size(m, H0_varname);
        fprintf('H0 dimensions: %s\n', mat2str(H0_dims));
    catch ME
        error('Cannot read dimensions from H0 file: %s\nError: %s', H0filename, ME.message);
    end
    
    % CRITICAL FIX: Determine the correct dimension structure
    % The H0 files should have structure: [channels, time, 2, bootstraps]
    % where dimension containing size 2 is the F/p statistics
    
    % Find which dimension has size 2 (statistics dimension)
    stats_dim_pos = find(H0_dims == 2);
    if isempty(stats_dim_pos)
        error('No dimension with size 2 found. H0 file may be missing statistics dimension.');
    end
    
    % For most LIMO outputs, stats should be in position ndims-1
    expected_stats_pos = length(H0_dims) - 1;
    
    if length(stats_dim_pos) > 1
        % Multiple dimensions with size 2, use the expected position
        if any(stats_dim_pos == expected_stats_pos)
            stats_dim_pos = expected_stats_pos;
        else
            warning('Multiple dimensions with size 2. Using position %d', stats_dim_pos(end));
            stats_dim_pos = stats_dim_pos(end);
        end
    end
    
    fprintf('Statistics dimension identified at position %d\n', stats_dim_pos);
    
    % Determine which statistic to extract based on test type
    if strcmpi(test_type, 'R2') || strcmpi(test_type, 'F')
        stat_index = 1;  % F values
        fprintf('Extracting F values (index 1 in statistics dimension)\n');
    elseif strcmpi(test_type, 'ttest')
        stat_index = 1;  % For t-tests stored as F values
        fprintf('Extracting F values from t-test (index 1 in statistics dimension)\n');
    else
        stat_index = 1;  % Default to F values
    end
    
    % Total bootstraps is the last dimension
    total_bootstraps = H0_dims(end);
    actual_bootstraps = min(total_bootstraps, LIMO.design.bootstrap);
    
    fprintf('Total bootstraps in file: %d\n', total_bootstraps);
    
    % ========= INCREMENTAL PROCESSING LOGIC =========
    existing_bootstraps = 0;
    existing_valid_bootstraps = [];
    bootstraps_to_process = 1:actual_bootstraps;
    
    if strcmpi(incremental, 'yes') && exist(H0_tfce_file, 'file')
        fprintf('\n--- Checking existing TFCE file for incremental processing ---\n');
        
        try
            % Load existing file to check dimensions and bootstrap count
            existing_data = matfile(H0_tfce_file);
            existing_vars = who(existing_data);
            
            if any(strcmp(existing_vars, 'tfce_H0_score'))
                existing_dims = size(existing_data, 'tfce_H0_score');
                existing_bootstraps = existing_dims(end);
                
                fprintf('Found existing TFCE file with %d bootstraps\n', existing_bootstraps);
                
                % Check for metadata about valid bootstraps
                if any(strcmp(existing_vars, 'valid_bootstraps'))
                    existing_valid_bootstraps = existing_data.valid_bootstraps;
                    fprintf('Found metadata: %d valid bootstraps in existing file\n', length(existing_valid_bootstraps));
                else
                    % Assume all existing bootstraps are valid
                    existing_valid_bootstraps = 1:existing_bootstraps;
                    fprintf('No metadata found, assuming all %d existing bootstraps are valid\n', existing_bootstraps);
                end
                
                % Determine which bootstraps still need processing
                if actual_bootstraps > existing_bootstraps
                    bootstraps_to_process = (existing_bootstraps + 1):actual_bootstraps;
                    fprintf('Will process additional bootstraps: %d to %d\n', existing_bootstraps + 1, actual_bootstraps);
                else
                    fprintf('All bootstraps already processed. No additional processing needed.\n');
                    return;
                end
            else
                fprintf('Existing file found but no tfce_H0_score variable. Starting fresh.\n');
            end
        catch ME
            fprintf('Could not read existing TFCE file: %s\n', ME.message);
            fprintf('Starting fresh processing...\n');
        end
    else
        if strcmpi(incremental, 'yes')
            fprintf('No existing TFCE file found. Starting fresh processing.\n');
        else
            fprintf('Incremental processing disabled. Starting fresh processing.\n');
        end
    end
    
    fprintf('Processing bootstraps: %d\n', length(bootstraps_to_process));
    
    % ========= CORRUPTION DETECTION =========
    fprintf('\n--- Checking for corrupted regions ---\n');
    
    % Test accessibility
    test_chunk_size = 100;
    accessible_bootstraps = true(actual_bootstraps, 1);
    corrupted_regions = [];
    
    for test_start = 1:test_chunk_size:actual_bootstraps
        test_end = min(test_start + test_chunk_size - 1, actual_bootstraps);
        
        try
            % Test read based on dimensions
            if length(H0_dims) == 4
                if stats_dim_pos == 3
                    test_read = m.(H0_varname)(1,1,stat_index,test_start);
                else
                    test_read = m.(H0_varname)(1,1,1,test_start);
                end
            elseif length(H0_dims) == 5
                if stats_dim_pos == 4
                    test_read = m.(H0_varname)(1,1,1,stat_index,test_start);
                else
                    test_read = m.(H0_varname)(1,1,1,1,test_start);
                end
            end
        catch
            fprintf('  Bootstraps %d-%d: CORRUPTED - will skip\n', test_start, test_end);
            accessible_bootstraps(test_start:test_end) = false;
            if isempty(corrupted_regions)
                corrupted_regions = [test_start, test_end];
            else
                corrupted_regions(end+1,:) = [test_start, test_end];
            end
        end
    end
    
    % Filter accessibility check to only bootstraps we need to process
    bootstraps_accessible = accessible_bootstraps(bootstraps_to_process);
    valid_bootstraps_to_process = bootstraps_to_process(bootstraps_accessible);
    n_valid_new = length(valid_bootstraps_to_process);
    n_corrupted_new = length(bootstraps_to_process) - n_valid_new;
    
    fprintf('\nCorruption summary for new bootstraps:\n');
    fprintf('  Valid new bootstraps: %d (%.1f%%)\n', n_valid_new, 100*n_valid_new/length(bootstraps_to_process));
    if n_corrupted_new > 0
        fprintf('  Corrupted new bootstraps: %d (%.1f%%)\n', n_corrupted_new, 100*n_corrupted_new/length(bootstraps_to_process));
        fprintf('  Processing only valid new bootstraps\n');
    end
    
    % Combine existing and new valid bootstraps
    if existing_bootstraps > 0
        all_valid_bootstraps = [existing_valid_bootstraps(:); valid_bootstraps_to_process(:)];
    else
        all_valid_bootstraps = valid_bootstraps_to_process;
    end
    n_total_valid = length(all_valid_bootstraps);
    
    % ========= DETERMINE OUTPUT DIMENSIONS =========
    % CRITICAL: Output should match expected LIMO structure
    % Remove the statistics dimension and use total valid bootstrap count
    
    output_dims = H0_dims;
    output_dims(stats_dim_pos) = [];  % Remove statistics dimension
    output_dims(end) = n_total_valid; % Use total valid bootstrap count
    
    fprintf('\nOutput TFCE dimensions: %s\n', mat2str(output_dims));
    
    % Process batches (only for new bootstraps)
    if n_valid_new == 0
        fprintf('\nNo new bootstraps to process.\n');
        return;
    end
    
    n_batches = ceil(n_valid_new / batch_size);
    fprintf('\nProcessing %d new valid bootstraps in %d batches\n', n_valid_new, n_batches);
    
    % Process each batch
    for batch = 1:n_batches
        fprintf('\nBatch %d/%d: ', batch, n_batches);
        
        % Calculate indices for new bootstraps only
        valid_start_idx = (batch - 1) * batch_size + 1;
        valid_end_idx = min(batch * batch_size, n_valid_new);
        batch_indices = valid_bootstraps_to_process(valid_start_idx:valid_end_idx);
        current_batch_size = length(batch_indices);
        
        % Adjust indices for final output (account for existing bootstraps)
        output_start_idx = existing_bootstraps + valid_start_idx;
        output_end_idx = existing_bootstraps + valid_end_idx;
        
        fprintf('Loading %d bootstraps...', current_batch_size);
        
        % Extract batch data with correct indexing
        try
            % Pre-allocate based on data structure (without stats dimension)
            if length(H0_dims) == 4  % [channels, time, stats, boots]
                batch_data = NaN(H0_dims(1), H0_dims(2), current_batch_size);
                
                % Read data
                for idx = 1:current_batch_size
                    boot_num = batch_indices(idx);
                    if stats_dim_pos == 3
                        batch_data(:,:,idx) = m.(H0_varname)(:,:,stat_index,boot_num);
                    else
                        % Handle unexpected dimension order
                        temp = m.(H0_varname)(:,:,:,boot_num);
                        batch_data(:,:,idx) = temp(:,:,stat_index);
                    end
                end
                
            elseif length(H0_dims) == 5  % [channels, freq, time, stats, boots] for TF
                batch_data = NaN(H0_dims(1), H0_dims(2), H0_dims(3), current_batch_size);
                
                % Read data
                for idx = 1:current_batch_size
                    boot_num = batch_indices(idx);
                    if stats_dim_pos == 4
                        batch_data(:,:,:,idx) = m.(H0_varname)(:,:,:,stat_index,boot_num);
                    else
                        % Handle unexpected dimension order
                        temp = m.(H0_varname)(:,:,:,:,boot_num);
                        batch_data(:,:,:,idx) = squeeze(temp(:,:,:,stat_index));
                    end
                end
            else
                error('Unexpected H0 dimensions: %s', mat2str(H0_dims));
            end
            
        catch ME
            error('Failed to read batch %d: %s', batch, ME.message);
        end
        
        % Process batch with TFCE
        fprintf(' Computing TFCE...\n');
        
        % Pre-allocate results matching batch_data dimensions
        batch_tfce = NaN(size(batch_data), 'double');
        
        % Process each bootstrap
        par_tfce = cell(1, current_batch_size);
        parfor b = 1:current_batch_size
            if H0_dims(1) == 1  % Single channel
                if length(size(batch_data)) == 3  % Time-frequency
                    par_tfce{b} = limo_tfce(2, squeeze(batch_data(:,:,b)), [], 0);
                else  % Time only
                    par_tfce{b} = limo_tfce(1, squeeze(batch_data(:,b)), LIMO.data.neighbouring_matrix, 0);
                end
            else  % Multiple channels
                if length(size(batch_data)) == 4  % Time-frequency
                    par_tfce{b} = limo_tfce(3, squeeze(batch_data(:,:,:,b)), LIMO.data.neighbouring_matrix, 0);
                else  % Time only
                    par_tfce{b} = limo_tfce(2, squeeze(batch_data(:,:,b)), LIMO.data.neighbouring_matrix, 0);
                end
            end
        end
        
        % Store results
        for b = 1:current_batch_size
            if length(size(batch_tfce)) == 3
                batch_tfce(:,:,b) = double(par_tfce{b});
            elseif length(size(batch_tfce)) == 4
                batch_tfce(:,:,:,b) = double(par_tfce{b});
            else
                batch_tfce(:,b) = double(par_tfce{b});
            end
        end
        
        % Save batch
        temp_file = fullfile(temp_dir, sprintf('tfce_batch_%03d.mat', batch));
        save(temp_file, 'batch_tfce', 'output_start_idx', 'output_end_idx', '-v7.3');
        
        fprintf('  Saved batch to %s\n', temp_file);
        clear batch_data batch_tfce par_tfce
    end
    
    % ========= MERGE BATCHES =========
    fprintf('\n--- Merging batch results ---\n');
    
    % Handle output file creation or extension
    if existing_bootstraps == 0 || ~strcmpi(incremental, 'yes')
        % Create new file
        if exist(H0_tfce_file, 'file')
            delete(H0_tfce_file);
        end
        fprintf('Creating new output file with dimensions: %s\n', mat2str(output_dims));
        tfce_H0_score = NaN(output_dims, 'double');
        
        % Copy existing data if we're in incremental mode
        if existing_bootstraps > 0 && strcmpi(incremental, 'yes')
            fprintf('  Copying %d existing bootstrap results...\n', existing_bootstraps);
            % Note: This path should not be taken in current logic, but kept for safety
            existing_data = matfile(H0_tfce_file);
            if length(output_dims) == 3
                tfce_H0_score(:,:,1:existing_bootstraps) = existing_data.tfce_H0_score;
            elseif length(output_dims) == 4
                tfce_H0_score(:,:,:,1:existing_bootstraps) = existing_data.tfce_H0_score;
            elseif length(output_dims) == 2
                tfce_H0_score(:,1:existing_bootstraps) = existing_data.tfce_H0_score;
            end
        end
    else
        % Extend existing file
        fprintf('Extending existing output file to dimensions: %s\n', mat2str(output_dims));
        
        % Load existing data
        existing_data = load(H0_tfce_file, 'tfce_H0_score');
        old_tfce = existing_data.tfce_H0_score;
        
        % Create extended array
        tfce_H0_score = NaN(output_dims, 'double');
        
        % Copy existing data
        if length(output_dims) == 3
            tfce_H0_score(:,:,1:existing_bootstraps) = old_tfce;
        elseif length(output_dims) == 4
            tfce_H0_score(:,:,:,1:existing_bootstraps) = old_tfce;
        elseif length(output_dims) == 2
            tfce_H0_score(:,1:existing_bootstraps) = old_tfce;
        end
        
        clear old_tfce existing_data;
    end
    
    % Merge batches
    successful_merges = 0;
    for batch = 1:n_batches
        temp_file = fullfile(temp_dir, sprintf('tfce_batch_%03d.mat', batch));
        
        if exist(temp_file, 'file')
            fprintf('  Merging batch %d/%d...', batch, n_batches);
            
            try
                batch_result = load(temp_file);
                
                % Insert into output array
                idx_range = batch_result.output_start_idx:batch_result.output_end_idx;
                
                if length(output_dims) == 3
                    tfce_H0_score(:,:,idx_range) = batch_result.batch_tfce;
                elseif length(output_dims) == 4
                    tfce_H0_score(:,:,:,idx_range) = batch_result.batch_tfce;
                elseif length(output_dims) == 2
                    tfce_H0_score(:,idx_range) = batch_result.batch_tfce;
                end
                
                successful_merges = successful_merges + 1;
                fprintf(' Done\n');
                
                % Delete temp file (optional)
                if strcmpi(keep_chunks, 'no')
                    delete(temp_file);
                end
            catch ME
                warning('Failed to merge batch %d: %s', batch, ME.message);
            end
        end
    end
    
    % Save final result
    fprintf('\nSaving final TFCE H0 result to: %s\n', H0_tfce_file);
    save(H0_tfce_file, 'tfce_H0_score', '-v7.3');

    % Save metadata about valid bootstraps
    fprintf('Saving metadata about valid bootstraps...\n');
    save(H0_tfce_file, 'all_valid_bootstraps', 'n_total_valid', '-append');
    
    % Save corruption information if there were any corrupted bootstraps
    if n_corrupted_new > 0
        corrupted_new_bootstraps = bootstraps_to_process(~bootstraps_accessible);
        save(H0_tfce_file, 'corrupted_new_bootstraps', 'n_corrupted_new', '-append');
    end
    
    fprintf('\nTFCE batch processing complete!\n');
    fprintf('Successfully merged %d/%d batches\n', successful_merges, n_batches);
    fprintf('Total output contains %d TFCE scores (%d existing + %d new)\n', n_total_valid, existing_bootstraps, n_valid_new);
    
    if n_corrupted_new > 0
        fprintf('\nNOTE: %d new corrupted bootstraps were excluded from analysis.\n', n_corrupted_new);
    end
end