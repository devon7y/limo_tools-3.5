function [tfce_score,thresholded_maps] = limo_tfce_handling(varargin)

% routine to create tfce files commensurate to boostrapped files
% MEMORY-EFFICIENT VERSION - processes bootstraps in batches
%
% FORMAT  [tfce_score,thresholded_maps] = limo_tfce_handling(filename,'checkfile','yes')
%         [tfce_score,thresholded_maps] = limo_tfce_handling(filename,'batch_size',100,'temp_dir','/path/to/temp')
%
% INPUTS filename is the stat file that need to be tfced (if H0 exist it is done too)
%        'checkfile' is 'yes' by default - if 'no' and tfce files already exist,
%                     it overwrites without asking otherwise user is prompted
%        'batch_size' number of bootstraps to process at once (default: 50)
%        'temp_dir' directory for temporary files (default: system temp)
%        'max_workers' maximum parallel workers (default: 4)
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
% Memory-Efficient Version - Cyril R. Pernet & Assistant - 2025

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
batch_size = 50;  % Process 50 bootstraps at a time
temp_dir = tempdir;  % Default system temp directory
max_workers = 4;  % Limit parallel workers

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
        answer = questdlg('tfce file already exist - overwrite?','data check','Yes','No','Yes');
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
        process_H0_in_batches(H0filename, H0_tfce_file, LIMO, batch_size, temp_dir, 'R2');
        
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
        process_H0_in_batches(H0filename, H0_tfce_file, LIMO, batch_size, temp_dir, 'ttest');
        
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
        process_H0_in_batches(H0filename, H0_tfce_file, LIMO, batch_size, temp_dir, 'F');
        
        % Handle thresholded maps
        tmp                 = thresholded_maps;     clear thresholded_maps;
        thresholded_maps{1} = tmp;                  clear tmp
        thresholded_maps{2} = [];  % Too memory intensive to store all bootstrap thresholded maps
    end
end

% Clean up temp directory
temp_pattern = fullfile(temp_dir, 'tfce_batch_*.mat');
delete(temp_pattern);

end

%% Helper function for batch processing of H0 data
function process_H0_in_batches(H0filename, H0_tfce_file, LIMO, batch_size, temp_dir, test_type)
    
    % Use matfile to access H0 data without loading it all
    fprintf('Opening H0 file using memory-mapped access...\n');
    m = matfile(H0filename);
    
    % Get variable name and dimensions
    vars = who(m);
    H0_varname = vars{1};
    H0_dims = size(m, H0_varname);
    
    % Determine data structure based on test type
    if strcmpi(test_type, 'R2')
        stat_dim = 2;  % F values are in dimension 2
        total_bootstraps = H0_dims(end);
    elseif strcmpi(test_type, 'ttest')
        if length(H0_dims) >= 4 && H0_dims(end-1) < 4
            stat_dim = 1;  % F-value structure
        else
            stat_dim = H0_dims(end-1) - 1;  % t values are in end-1 dimension
        end
        total_bootstraps = H0_dims(end);
    else % F-value files
        stat_dim = 1;  % F values are in dimension 1
        total_bootstraps = H0_dims(end);
    end
    
    % Calculate number of batches
    n_batches = ceil(total_bootstraps / batch_size);
    
    fprintf('Processing %d bootstraps in %d batches of %d\n', total_bootstraps, n_batches, batch_size);
    
    % Initialize result dimensions based on data type
    if H0_dims(1) == 1  % Single channel
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            tfce_dims = [1, H0_dims(2), H0_dims(3), total_bootstraps];
        else
            tfce_dims = [1, H0_dims(2), total_bootstraps];
        end
    else  % Multiple channels
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            tfce_dims = [H0_dims(1), H0_dims(2), H0_dims(3), total_bootstraps];
        else
            tfce_dims = [H0_dims(1), H0_dims(2), total_bootstraps];
        end
    end
    
    % Process batches
    for batch = 1:n_batches
        fprintf('\nProcessing batch %d/%d...', batch, n_batches);
        
        % Calculate batch indices
        start_idx = (batch - 1) * batch_size + 1;
        end_idx = min(batch * batch_size, total_bootstraps);
        batch_indices = start_idx:end_idx;
        current_batch_size = length(batch_indices);
        
        % Extract batch data using matfile
        fprintf(' Loading bootstraps %d-%d...', start_idx, end_idx);
        
        if strcmpi(test_type, 'R2')
            if length(H0_dims) == 5  % Time-frequency
                batch_data = m.(H0_varname)(:,:,:,stat_dim,batch_indices);
            elseif length(H0_dims) == 4
                batch_data = m.(H0_varname)(:,:,stat_dim,batch_indices);
            else
                error('Unexpected dimensions for R2 test: %s', mat2str(H0_dims));
            end
        
        elseif strcmpi(test_type, 'ttest')
            if length(H0_dims) >= 4 && H0_dims(end-1) < 4  % F-value structure
                if length(H0_dims) == 5  % Time-frequency
                    batch_data = m.(H0_varname)(:,:,:,stat_dim,batch_indices);
                elseif length(H0_dims) == 4
                    batch_data = m.(H0_varname)(:,:,stat_dim,batch_indices);
                else
                    error('Unexpected dimensions for t-test: %s', mat2str(H0_dims));
                end
            else
                error('Unsupported H0_dims structure for t-test: %s', mat2str(H0_dims));
            end
        
        else  % F-value files
            if length(H0_dims) == 5  % Hypothetical higher-D TF case
                batch_data = m.(H0_varname)(:,:,:,stat_dim,batch_indices);
            elseif length(H0_dims) == 4  % This is your case
                batch_data = m.(H0_varname)(:,:,stat_dim,batch_indices);
            else
                error('Unexpected dimensions for F-value test: %s', mat2str(H0_dims));
            end
        end
        
        % Process batch with parallel computation
        fprintf(' Computing TFCE...\n');
        
        % Pre-allocate batch results
        if H0_dims(1) == 1  % Single channel
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                batch_tfce = NaN(1, H0_dims(2), H0_dims(3), current_batch_size);
            else
                batch_tfce = NaN(1, H0_dims(2), current_batch_size);
            end
        else  % Multiple channels
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                batch_tfce = NaN(H0_dims(1), H0_dims(2), H0_dims(3), current_batch_size);
            else
                batch_tfce = NaN(H0_dims(1), H0_dims(2), current_batch_size);
            end
        end
        
        % Process each bootstrap in the batch
        par_tfce = cell(1, current_batch_size);
        parfor b = 1:current_batch_size
            if H0_dims(1) == 1  % Single channel
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    par_tfce{b} = limo_tfce(2, squeeze(batch_data(:,:,:,b)), [], 0);
                else
                    par_tfce{b} = limo_tfce(1, squeeze(batch_data(:,:,b)), LIMO.data.neighbouring_matrix, 0);
                end
            else  % Multiple channels
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    par_tfce{b} = limo_tfce(3, squeeze(batch_data(:,:,:,b)), LIMO.data.neighbouring_matrix, 0);
                else
                    par_tfce{b} = limo_tfce(2, squeeze(batch_data(:,:,b)), LIMO.data.neighbouring_matrix, 0);
                end
            end
        end
        
        % Convert cell array back to numeric array
        for b = 1:current_batch_size
            if H0_dims(1) == 1
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    batch_tfce(1,:,:,b) = par_tfce{b};
                else
                    batch_tfce(1,:,b) = par_tfce{b};
                end
            else
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    batch_tfce(:,:,:,b) = par_tfce{b};
                else
                    batch_tfce(:,:,b) = par_tfce{b};
                end
            end
        end
        
        % Save batch results to temporary file
        temp_file = fullfile(temp_dir, sprintf('tfce_batch_%03d.mat', batch));
        fprintf('  Saving batch to %s\n', temp_file);
        save(temp_file, 'batch_tfce', 'batch_indices', '-v7.3');
        
        % Clear variables to free memory
        clear batch_data batch_tfce
        
    end
    
    % Merge batch results
    fprintf('\nMerging batch results...\n');
    
    % Create output file with matfile for memory efficiency
    if exist(H0_tfce_file, 'file')
        delete(H0_tfce_file);
    end
    m_out = matfile(H0_tfce_file, 'Writable', true);
    m_out.tfce_H0_score = NaN(tfce_dims);
    
    % Load and merge each batch
    for batch = 1:n_batches
        fprintf('  Merging batch %d/%d...', batch, n_batches);
        
        temp_file = fullfile(temp_dir, sprintf('tfce_batch_%03d.mat', batch));
        batch_result = load(temp_file);
        
        % Write to output file
        if length(tfce_dims) == 4  % Time-frequency or multi-channel time
            m_out.tfce_H0_score(:,:,:,batch_result.batch_indices) = batch_result.batch_tfce;
        else  % Single channel time
            m_out.tfce_H0_score(:,:,batch_result.batch_indices) = batch_result.batch_tfce;
        end
        
        % Delete temp file
        delete(temp_file);
        fprintf(' Done\n');
    end
    
    fprintf('Batch processing complete!\n');
end