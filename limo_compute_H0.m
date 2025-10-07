function varargout=limo_compute_H0(varargin)

% this function allows comnputing the null distribution of any of the
% LIMO EEG 2nd level statistical tests using bootstrap
%
% FORMAT [H0_data,boot_table] = limo_compute_H0(0,data,nboot,test);
%        status  = limo_compute_H0(type,data,label,nboot)
%        status  = limo_compute_H0(type,data,label,nboot,'incremental','yes')
%        status  = limo_compute_H0(type,data,label,nboot,'existing_bootstraps',N)
%
% -------------------------------------------------------------------------
% [H0_data,boot_table] = limo_compute_H0(0,data,nboot,test)
% the generic case to compute H0 using limo stat functions
%
% INPUTS: 0 stands for type = 0 and indicates the generic case
%        data this is a 4D or 3D matrix of data
%        nboot is the number of boostrap to use (recommended at least 800)
%        test is one of limo statistical test: 'limo_trimci' 'limo_yuen'
%                                              'limo_yuend' 'limo_glm'
%                                              'limo_'
%
% OUTPUTS: H0_data is the matrix of statistical values (t/F and p) needed
%          to do a correction for multiple comparisons
%          boot_table is the resampling table used
%
% -------------------------------------------------------------------------
% status = limo_compute_H0(type,data,label,boot_table)
% this is the format used by LIMO EEG, and in this case it creates data on
% the drive using dedicated names and format - see limo_random_robust for
% details of inputs and outputs
%
% ADDITIONAL PARAMETERS (for incremental bootstrapping):
%        'incremental' - 'yes' or 'no' (default 'no') - enable incremental processing
%        'existing_bootstraps' - number of existing bootstraps to extend from
%        'keep_chunks' - 'yes' or 'no' (default 'no') - keep temporary chunk files
%
% see also limo_create_boot_table
%
% ------------------------------
%  Copyright (C) LIMO Team 2019

type = varargin{1};

% Parse additional parameters for incremental processing
incremental = 'no';
existing_bootstraps = 0;
keep_chunks = 'no';

% Look for optional parameters
for i = 1:2:length(varargin)-1
    if i+1 <= length(varargin)
        if strcmpi(varargin{i}, 'incremental')
            incremental = varargin{i+1};
        elseif strcmpi(varargin{i}, 'existing_bootstraps')
            existing_bootstraps = varargin{i+1};
        elseif strcmpi(varargin{i}, 'keep_chunks')
            keep_chunks = varargin{i+1};
        end
    end
end

% get the data and boot_table (fix typos in variable names)
if type == 0
    test = varargin{end};
    if strcmpi(test,'limo_yuen_ttest') || strcmpi(test,'limo_yuend_ttest')
        data1 = varargin{2};
        data2 = varargin{3};
        nboot = varargin{4};
    else
        data = varargin{2};
        nboot = varargin{3};
    end
elseif type == 1 || type > 3
    data = varargin{2};
    nboot = varargin{3};
else
    data1 = varargin{2};
    data2 = varargin{3};
    nboot = varargin{4};
end

fprintf('\n=== INCREMENTAL BOOTSTRAP H0 COMPUTATION ===\n');
fprintf('Incremental mode: %s\n', incremental);
if strcmpi(incremental, 'yes')
    fprintf('Existing bootstraps: %d\n', existing_bootstraps);
    fprintf('New bootstraps to compute: %d\n', nboot - existing_bootstraps);
end
fprintf('Keep chunk files: %s\n', keep_chunks);
fprintf('===============================================\n\n');

clear varargin

% -------------------------------------------------------------------------

switch type
    case {1}
        % one sample t-test
        % -----------------
        bootex = 1;
        boot_name = sprintf('H0_one_sample_ttest_parameter_%g',parameter);
        boot_file_path = ['H0', filesep, boot_name, '.mat'];
        
        % Handle incremental processing
        if exist(boot_file_path, 'file') && strcmpi(incremental, 'yes')
            % Check existing bootstrap count
            try
                existing_info = whos('-file', boot_file_path);
                for fi = 1:length(existing_info)
                    if startsWith(existing_info(fi).name, 'H0_')
                        existing_size = existing_info(fi).size;
                        current_existing = existing_size(end);  % Last dimension
                        fprintf('Found existing bootstrap file with %d bootstraps\n', current_existing);
                        
                        if current_existing >= nboot
                            fprintf('Requested bootstraps (%d) already exist. No additional processing needed.\n', nboot);
                            bootex = 0;
                        else
                            fprintf('Will add %d new bootstraps to existing %d\n', nboot - current_existing, current_existing);
                            existing_bootstraps = current_existing;
                        end
                        break;
                    end
                end
            catch
                fprintf('Could not read existing bootstrap file. Starting fresh.\n');
                existing_bootstraps = 0;
            end
        elseif exist(boot_file_path, 'file') && strcmpi(incremental, 'no')
            % Traditional behavior
            answer = questdlg('a boostrap file already exist - overwrite?','data check','Yes','No','Yes');
            if strcmp(answer,'Yes')
                bootex = 1;
            else
                bootex = 0;
            end
        end
        
        if bootex == 1;
            if ~exist('H0', 'dir'), mkdir H0; end
            
            % Handle incremental vs fresh processing
            if strcmpi(incremental, 'yes') && existing_bootstraps > 0
                % Load existing data for extension
                fprintf('Loading existing bootstrap data for extension...\n');
                existing_data = load(boot_file_path);
                field_names = fieldnames(existing_data);
                bootstrap_var = [];
                for fn = 1:length(field_names)
                    if startsWith(field_names{fn}, 'H0_')
                        bootstrap_var = field_names{fn};
                        break;
                    end
                end
                
                if ~isempty(bootstrap_var)
                    H0_one_sample_existing = existing_data.(bootstrap_var);
                    % Create extended array
                    H0_one_sample = NaN(size(data,1), size(data,2), 2, nboot);
                    H0_one_sample(:,:,:,1:existing_bootstraps) = H0_one_sample_existing;
                    clear existing_data H0_one_sample_existing;
                    
                    % Only process new bootstraps
                    bootstrap_start = existing_bootstraps + 1;
                    bootstrap_end = nboot;
                    fprintf('Processing bootstraps %d to %d\n', bootstrap_start, bootstrap_end);
                else
                    error('Could not find bootstrap variable in existing file');
                end
            else
                % Fresh processing
                H0_one_sample = NaN(size(data,1), size(data,2), 2, nboot);
                bootstrap_start = 1;
                bootstrap_end = nboot;
                fprintf('Processing all bootstraps %d to %d\n', bootstrap_start, bootstrap_end);
            end
            
            % create centered data to estimate H0
            centered_data = data - repmat(limo_trimmed_mean(data),[1 1 size(data,3)]);
            % get boot table
            disp('making boot table ...')
            boot_table = limo_create_boot_table(data, nboot);
            save(['H0', filesep, 'boot_table'], 'boot_table')
            
            % get results under H0 - only for new bootstraps
            for electrode = 1:size(data,1)
                fprintf('bootstrap: electrode %g parameter %g (bootstraps %d-%d)\n', electrode, parameter, bootstrap_start, bootstrap_end);
                tmp = centered_data(electrode,:,:); Y = tmp(1,:,find(~isnan(tmp(1,1,:))));
                
                if exist('parfor','file') ~=0
                    % Use parallel processing for new bootstraps only
                    t = cell(1, bootstrap_end - bootstrap_start + 1);
                    p = cell(1, bootstrap_end - bootstrap_start + 1);
                    
                    parfor b_idx = 1:(bootstrap_end - bootstrap_start + 1)
                        actual_b = bootstrap_start + b_idx - 1;
                        [t{b_idx},~,~,~,p{b_idx},~,~] = limo_trimci(Y(1,:,boot_table{electrode}(:,actual_b)));
                    end
                    
                    % Store results for new bootstraps
                    for b_idx = 1:(bootstrap_end - bootstrap_start + 1)
                        actual_b = bootstrap_start + b_idx - 1;
                        H0_one_sample(electrode,:,1,actual_b) = t{b_idx};
                        H0_one_sample(electrode,:,2,actual_b) = p{b_idx};
                    end
                    clear t p;
                else
                    % Sequential processing for new bootstraps only
                    for actual_b = bootstrap_start:bootstrap_end
                        [H0_one_sample(electrode,:,1,actual_b),~,~,~,H0_one_sample(electrode,:,2,actual_b),~,~] = ...
                            limo_trimci(Y(1,:,boot_table{electrode}(:,actual_b)));
                    end
                end
                clear tmp Y
            end % closes for electrode
            
            if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                H0_one_sample = limo_tf_5d_reshape(H0_one_sample);
            end
            
            % Save with metadata for incremental processing
            save(['H0', filesep, boot_name], 'H0_one_sample', '-v7.3');
            
            % Save bootstrap metadata
            bootstrap_info.total_bootstraps = nboot;
            bootstrap_info.existing_bootstraps = existing_bootstraps;
            bootstrap_info.new_bootstraps = bootstrap_end - bootstrap_start + 1;
            bootstrap_info.incremental_mode = incremental;
            bootstrap_info.creation_time = datestr(now);
            save(['H0', filesep, boot_name], 'bootstrap_info', '-append');
            
            fprintf('Bootstrap H0 computation complete. Total bootstraps: %d\n', nboot);
        end
        
        % ----------------------------------------------------------------------
        % ----------------------------------------------------------------------
        %                    THIS IS TYPE 0 - GENRERIC CASE
        % ----------------------------------------------------------------------
        % ----------------------------------------------------------------------
    case {0}
        % one sample t-test
        % -----------------
        if strcmpi(test,'limo_trimci')
            
            reshape = 0;
            if numel(size(data)) == 4
                reshape = 1; [elect,freq,time,obs]=size(data);
                data = limo_tf_4d_reshape(data,[elect,freq*time,obs]);
            end
            
            % Handle incremental processing for generic case
            if strcmpi(incremental, 'yes') && existing_bootstraps > 0
                H0 = NaN(size(data,1), size(data,2), 2, nboot);
                % Note: For generic case, we assume fresh computation
                % User should save/load existing H0 data externally
                bootstrap_start = existing_bootstraps + 1;
                bootstrap_end = nboot;
                fprintf('Generic case: Processing new bootstraps %d to %d\n', bootstrap_start, bootstrap_end);
            else
                H0 = NaN(size(data,1), size(data,2), 2, nboot); % stores T and p values for each boot under H0
                bootstrap_start = 1;
                bootstrap_end = nboot;
            end
            
            centered_data = data - repmat(limo_trimmed_mean(data),[1 1 size(data,3)]);
            disp('making boot table ...'); boot_table = limo_create_boot_table(data,nboot);
            
            for electrode = 1:size(data,1)
                fprintf('bootstrap: electrode %g (bootstraps %d-%d)\n', electrode, bootstrap_start, bootstrap_end);
                tmp = centered_data(electrode,:,:); Y = tmp(1,:,find(~isnan(tmp(1,1,:))));
                
                % Process only new bootstraps
                t = cell(1, bootstrap_end - bootstrap_start + 1);
                p = cell(1, bootstrap_end - bootstrap_start + 1);
                
                parfor b_idx = 1:(bootstrap_end - bootstrap_start + 1)
                    actual_b = bootstrap_start + b_idx - 1;
                    [t{b_idx},~,~,~,p{b_idx},~,~] = limo_trimci(Y(1,:,boot_table{electrode}(:,actual_b)));
                end
                
                % Store results for new bootstraps
                for b_idx = 1:(bootstrap_end - bootstrap_start + 1)
                    actual_b = bootstrap_start + b_idx - 1;
                    H0(electrode,:,1,actual_b) = t{b_idx};
                    H0(electrode,:,2,actual_b) = p{b_idx};
                end
                clear tmp Y t p
            end % closes for electrode
            
            if reshape == 1
                H0 = limo_tf_5d_reshape(H0,[elect,freq,time,nboot]);
            end
            
            % Add metadata for incremental processing
            if nargout > 2
                bootstrap_info.total_bootstraps = nboot;
                bootstrap_info.existing_bootstraps = existing_bootstraps;
                bootstrap_info.new_bootstraps = bootstrap_end - bootstrap_start + 1;
                bootstrap_info.incremental_mode = incremental;
                bootstrap_info.creation_time = datestr(now);
                varargout{3} = bootstrap_info;
            end
        end
        clear data; varargout{2} = boot_table; clear boot_table  % Fixed typo
        varargout{1} = H0; clear H0;
end


