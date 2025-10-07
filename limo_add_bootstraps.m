function limo_add_bootstraps(LIMO_path, additional_bootstraps, chunk_size)
% LIMO_ADD_BOOTSTRAPS Add more bootstraps to existing analysis
%
% This function allows you to incrementally add more bootstraps to an
% existing analysis without regenerating all previous bootstraps.
%
% INPUTS:
%   LIMO_path            - Path to LIMO.mat file
%   additional_bootstraps - Number of bootstraps to add
%   chunk_size           - Bootstraps per chunk (default: 100)
%
% Example:
%   % Add 1000 more bootstraps to existing analysis
%   limo_add_bootstraps('LIMO.mat', 1000);
%
%   % Add 5000 bootstraps in chunks of 50
%   limo_add_bootstraps('LIMO.mat', 5000, 50);
%
% Devon Yanitski & Claude 2025

if nargin < 3
    chunk_size = 100;
end

% Load LIMO structure
fprintf('Loading LIMO structure...\n');
load(LIMO_path, 'LIMO');

% Check current bootstrap status
chunk_dir = fullfile(LIMO.dir, 'H0', 'chunks');
current_bootstraps = 0;

% Count existing bootstraps from chunks
if exist(chunk_dir, 'dir')
    chunk_files = dir(fullfile(chunk_dir, '*_chunk_*.mat'));
    if ~isempty(chunk_files)
        % Find highest bootstrap number
        for i = 1:length(chunk_files)
            tokens = regexp(chunk_files(i).name, '_chunk_(\d+)_(\d+)\.mat', 'tokens');
            if ~isempty(tokens)
                end_boot = str2double(tokens{1}{2});
                current_bootstraps = max(current_bootstraps, end_boot);
            end
        end
    end
end

% Check merged files
if current_bootstraps == 0
    % Check if merged H0 files exist
    H0_files = dir(fullfile(LIMO.dir, 'H0', 'H0_*.mat'));
    for i = 1:length(H0_files)
        if ~contains(H0_files(i).name, 'chunk') && ~contains(H0_files(i).name, 'info')
            try
                m = matfile(fullfile(LIMO.dir, 'H0', H0_files(i).name));
                vars = who(m);
                if ~isempty(vars)
                    dims = size(m, vars{1});
                    current_bootstraps = max(current_bootstraps, dims(end));
                end
            catch
                % Skip corrupted files
            end
            break; % Only need to check one file
        end
    end
end

fprintf('\nCurrent bootstrap status:\n');
fprintf('  Existing bootstraps: %d\n', current_bootstraps);
fprintf('  Additional bootstraps requested: %d\n', additional_bootstraps);
fprintf('  Total after addition: %d\n', current_bootstraps + additional_bootstraps);



% Update LIMO structure
old_bootstrap = LIMO.design.bootstrap;
LIMO.design.bootstrap = current_bootstraps + additional_bootstraps;
save(LIMO_path, 'LIMO');

% Check for parallel processing capabilities
fprintf('\n=== CHECKING PARALLEL PROCESSING ===\n');
limo_check_ppool;

% Generate additional bootstraps
fprintf('\n=== GENERATING ADDITIONAL BOOTSTRAPS ===\n');
fprintf('Starting from bootstrap %d\n', current_bootstraps + 1);

% Create chunk directory if needed
if ~exist(chunk_dir, 'dir')
    mkdir(chunk_dir);
end

% Load or create boot table
boot_table_file = fullfile(LIMO.dir, 'H0', 'boot_table.mat');
if exist(boot_table_file, 'file')
    fprintf('Loading existing boot table...\n');
    boot_table_struct = load(boot_table_file);
    boot_table = boot_table_struct.boot_table;
else
    fprintf('Creating new boot table...\n');
    % We'll create it after loading the data to know the dimensions
    boot_table = [];
end

% Determine analysis type by checking existing files
analysis_type = [];
parameter = [];

% Check for different analysis types based on existing files
if exist(fullfile(LIMO.dir, 'one_sample_ttest_parameter_1.mat'), 'file')
    analysis_type = 'one_sample';
    % Find parameter number from filename
    files = dir(fullfile(LIMO.dir, 'one_sample_ttest_parameter_*.mat'));
    if ~isempty(files)
        tokens = regexp(files(1).name, 'parameter_(\d+)', 'tokens');
        if ~isempty(tokens)
            parameter = str2double(tokens{1}{1});
        end
    end
elseif exist(fullfile(LIMO.dir, 'two_samples_ttest_parameter_1.mat'), 'file')
    analysis_type = 'two_samples';
    files = dir(fullfile(LIMO.dir, 'two_samples_ttest_parameter_*.mat'));
    if ~isempty(files)
        tokens = regexp(files(1).name, 'parameter_(\d+)', 'tokens');
        if ~isempty(tokens)
            parameter = str2double(tokens{1}{1});
        end
    end
elseif ~isempty(dir(fullfile(LIMO.dir, 'paired_samples_ttest_parameter_*.mat')))
    analysis_type = 'paired_samples';
    files = dir(fullfile(LIMO.dir, 'paired_samples_ttest_parameter_*.mat'));
    if ~isempty(files)
        % Extract parameter from filename - handle both single and multiple parameters
        tokens = regexp(files(1).name, 'parameter_(.+)\.mat', 'tokens');
        if ~isempty(tokens)
            param_str = tokens{1}{1};
            % Handle cases like "1" or "1_2" etc.
            parameter = str2num(param_str); %#ok<ST2NM>
            if isempty(parameter)
                parameter = 1; % Default fallback
            end
        end
    end
elseif exist(fullfile(LIMO.dir, 'R2.mat'), 'file')
    analysis_type = 'regression';
elseif exist(fullfile(LIMO.dir, 'Condition_effect_1.mat'), 'file')
    analysis_type = 'anova';
elseif ~isempty(dir(fullfile(LIMO.dir, 'Rep_ANOVA_*.mat')))
    analysis_type = 'repeated_anova';

    % Handle repeated measures ANOVA
    if isfield(LIMO.design, 'repeated_measure')
        if length(unique(LIMO.data.Cat)) == 1
            if length(LIMO.design.repeated_measure) == 1
                rep_type = 1; % Simple repeated measure
            else
                rep_type = 2; % Multiple factors
            end
        else
            if length(LIMO.design.repeated_measure) == 1
                rep_type = 3; % Group × repeated measure
            else
                rep_type = 4; % Group × multiple factors
            end
        end

        % Load additional required variables
        factor_levels = LIMO.design.repeated_measure;
        gp_vector = LIMO.data.Cat;
        C = LIMO.design.C;

        if rep_type == 3 || rep_type == 4
            gp_values = unique(gp_vector);
            k = length(gp_values);
            X = NaN(size(gp_vector,1),k+1);
            for g = 1:k
                X(:,g) = gp_vector == gp_values(g);
            end
            X(:,end) = 1;
        else
            X = [];
        end

        % Get filenames
        Rep_filenames = {};
        IRep_filenames = {};
        nb_effects = length(LIMO.design.C);

        for i = 1:nb_effects
            if contains(LIMO.design.effects{i},'Main effect')
                if isfield(LIMO.design,'factor_names')
                    Rep_filenames{i} = sprintf('Rep_ANOVA_Main_effect_%g_%s.mat',i,LIMO.design.factor_names{i});
                else
                    Rep_filenames{i} = sprintf('Rep_ANOVA_Main_effect_%g.mat',i);
                end
            elseif contains(LIMO.design.effects{i},'Interaction')
                Interaction = LIMO.design.effects{i}(length('Interaction')+1:end);
                Interaction(isspace(Interaction)) = [];
                Rep_filenames{i} = sprintf('Rep_ANOVA_Interaction_Factors_%s.mat',Interaction);
            end
        end

        if rep_type == 3 || rep_type == 4
            for i = 1:size(nb_effects,1)
                if contains(LIMO.design.effects{i},'Main effect')
                    if isfield(LIMO.design,'factor_names')
                        IRep_filenames{i} = sprintf('Rep_ANOVA_Interaction_gp_Factor_%g_%s.mat',i,LIMO.design.factor_names{i});
                    else
                        IRep_filenames{i} = sprintf('Rep_ANOVA_Interaction_gp_Factor_%g.mat',i);
                    end
                else
                    Interaction = LIMO.design.effects{i}(length('Interaction')+1:end);
                    Interaction(isspace(Interaction)) = [];
                    IRep_filenames{i} = sprintf('Rep_ANOVA_Interaction_gp_Factors_%s.mat',Interaction);
                end
            end
        end
    end
else
    error('Cannot determine analysis type from existing files. Supported types: one-sample t-test, two-samples t-test, paired t-test, regression, ANOVA, repeated measures ANOVA');
end

fprintf('Detected analysis type: %s\n', analysis_type);
if ~isempty(parameter)
    fprintf('Parameter: %s\n', num2str(parameter));
end

% Prepare data for bootstrap processing based on analysis type
switch analysis_type
    case 'one_sample'
        % Load data for one-sample t-test with memory optimization
        data_files = find_limo_data_files(LIMO.dir, 'Y');
        if isempty(data_files)
            error('Cannot find data files for one-sample t-test. Tried: Y*.mat, Yr.mat');
        end

        data_file_path = fullfile(LIMO.dir, data_files(1).name);
        file_info = dir(data_file_path);
        file_size_mb = file_info.bytes / 1024^2;

        if file_size_mb > 100  % Use memory mapping for files > 100MB
            fprintf('Large data file detected (%.1f MB) - using memory-mapped access\n', file_size_mb);
            m = matfile(data_file_path);
            vars = who(m);
            if ~isempty(vars)
                data = m.(vars{1});
            else
                error('No variables found in data file');
            end
        else
            original_data = load(data_file_path);
            data = original_data.(cell2mat(fieldnames(original_data)));
            clear original_data;
        end

        % Center data under H0 (mean = 0)
        if strcmpi(LIMO.design.method, 'Trimmed Mean')
            centered_data = data - repmat(limo_trimmed_mean(data), [1 1 size(data,3)]);
        else
            centered_data = data - repmat(nanmean(data,3), [1 1 size(data,3)]);
        end

    case 'two_samples'
        % Load data for two-samples t-test with memory optimization
        data1_files = find_limo_data_files(LIMO.dir, 'Y1');
        data2_files = find_limo_data_files(LIMO.dir, 'Y2');

        if isempty(data1_files) || isempty(data2_files)
            error('Cannot find data files for two-samples t-test. Tried: Y1*.mat, Y1r.mat, Yr1.mat and Y2*.mat, Y2r.mat, Yr2.mat');
        end

        % Load data1 with memory optimization
        data1_file_path = fullfile(LIMO.dir, data1_files(1).name);
        file1_info = dir(data1_file_path);
        file1_size_mb = file1_info.bytes / 1024^2;

        if file1_size_mb > 100
            fprintf('Large data file 1 detected (%.1f MB) - using memory-mapped access\n', file1_size_mb);
            m1 = matfile(data1_file_path);
            vars1 = who(m1);
            data1 = m1.(vars1{1});
        else
            original_data1 = load(data1_file_path);
            data1 = original_data1.(cell2mat(fieldnames(original_data1)));
            clear original_data1;
        end

        % Load data2 with memory optimization
        data2_file_path = fullfile(LIMO.dir, data2_files(1).name);
        file2_info = dir(data2_file_path);
        file2_size_mb = file2_info.bytes / 1024^2;

        if file2_size_mb > 100
            fprintf('Large data file 2 detected (%.1f MB) - using memory-mapped access\n', file2_size_mb);
            m2 = matfile(data2_file_path);
            vars2 = who(m2);
            data2 = m2.(vars2{1});
        else
            original_data2 = load(data2_file_path);
            data2 = original_data2.(cell2mat(fieldnames(original_data2)));
            clear original_data2;
        end

        % Center both datasets under H0
        if contains(LIMO.design.method, 'Trimmed Mean', 'IgnoreCase', true)
            data1_centered = data1 - repmat(limo_trimmed_mean(data1), [1 1 size(data1,3)]);
            data2_centered = data2 - repmat(limo_trimmed_mean(data2), [1 1 size(data2,3)]);
        else
            data1_centered = data1 - repmat(nanmean(data1,3), [1 1 size(data1,3)]);
            data2_centered = data2 - repmat(nanmean(data2,3), [1 1 size(data2,3)]);
        end
        centered_data = {data1_centered, data2_centered};

    case 'paired_samples'
        % Load data for paired t-test with memory optimization
        data1_files = find_limo_data_files(LIMO.dir, 'Y1');
        data2_files = find_limo_data_files(LIMO.dir, 'Y2');

        if isempty(data1_files) || isempty(data2_files)
            error('Cannot find data files for paired t-test. Tried: Y1*.mat, Y1r.mat, Yr1.mat and Y2*.mat, Y2r.mat, Yr2.mat');
        end

        % Load data1 with memory optimization
        data1_file_path = fullfile(LIMO.dir, data1_files(1).name);
        file1_info = dir(data1_file_path);
        file1_size_mb = file1_info.bytes / 1024^2;

        if file1_size_mb > 100
            fprintf('Large data file 1 detected (%.1f MB) - using memory-mapped access\n', file1_size_mb);
            m1 = matfile(data1_file_path);
            vars1 = who(m1);
            data1 = m1.(vars1{1});
        else
            original_data1 = load(data1_file_path);
            data1 = original_data1.(cell2mat(fieldnames(original_data1)));
            clear original_data1;
        end

        % Load data2 with memory optimization
        data2_file_path = fullfile(LIMO.dir, data2_files(1).name);
        file2_info = dir(data2_file_path);
        file2_size_mb = file2_info.bytes / 1024^2;

        if file2_size_mb > 100
            fprintf('Large data file 2 detected (%.1f MB) - using memory-mapped access\n', file2_size_mb);
            m2 = matfile(data2_file_path);
            vars2 = who(m2);
            data2 = m2.(vars2{1});
        else
            original_data2 = load(data2_file_path);
            data2 = original_data2.(cell2mat(fieldnames(original_data2)));
            clear original_data2;
        end

        % Center both datasets under H0 (difference = 0)
        if contains(LIMO.design.method, 'Trimmed Mean', 'IgnoreCase', true)
            data1_centered = data1 - repmat(limo_trimmed_mean(data1), [1 1 size(data1,3)]);
            data2_centered = data2 - repmat(limo_trimmed_mean(data2), [1 1 size(data2,3)]);
        else
            data1_centered = data1 - repmat(nanmean(data1,3), [1 1 size(data1,3)]);
            data2_centered = data2 - repmat(nanmean(data2,3), [1 1 size(data2,3)]);
        end
        centered_data = {data1_centered, data2_centered};

    case 'repeated_anova'
        % Load data for repeated measures ANOVA with memory optimization
        data_files = dir(fullfile(LIMO.dir, 'Yr.mat'));
        if isempty(data_files)
            error('Cannot find original data file Yr.mat for bootstrap processing');
        end

        data_file_path = fullfile(LIMO.dir, 'Yr.mat');
        file_info = dir(data_file_path);
        file_size_mb = file_info.bytes / 1024^2;

        if file_size_mb > 100
            fprintf('Large repeated ANOVA data file detected (%.1f MB) - using memory-mapped access\n', file_size_mb);
            m = matfile(data_file_path);
            data = m.Yr;
        else
            original_data = load(data_file_path);
            data = original_data.Yr;
            clear original_data;
        end

        if rep_type == 1 || rep_type == 2
            % Center data for each condition
            centered_data = NaN(size(data));
            nb_conditions = prod(factor_levels);

            for condition = 1:nb_conditions
                if contains(LIMO.design.method, 'Trimmed Mean', 'IgnoreCase', true)
                    avg = repmat(limo_trimmed_mean(data(:,:,:,condition), 3), [1 1 size(data,3)]);
                else
                    avg = repmat(nanmean(data(:,:,:,condition), 3), [1 1 size(data,3)]);
                end
                centered_data(:,:,:,condition) = data(:,:,:,condition) - avg;
            end

        elseif rep_type == 3 || rep_type == 4
            % Center data within each group separately
            centered_data = NaN(size(data));
            nb_conditions = prod(factor_levels);

            for gp = 1:LIMO.design.nb_conditions
                gp_index = find(gp_vector == gp);
                for condition = 1:nb_conditions
                    if contains(LIMO.design.method, 'Trimmed Mean', 'IgnoreCase', true)
                        avg = repmat(limo_trimmed_mean(data(:,:,gp_index,condition), 3), [1 1 length(gp_index)]);
                    else
                        avg = repmat(nanmean(data(:,:,gp_index,condition), 3), [1 1 length(gp_index)]);
                    end
                    centered_data(:,:,gp_index,condition) = data(:,:,gp_index,condition) - avg;
                end
            end
        end

    otherwise
        error('Unsupported analysis type: %s. Currently supported: one_sample, two_samples, paired_samples, repeated_anova', analysis_type);
end

% Save centered data for future use
save(fullfile(LIMO.dir, 'H0', 'centered_data.mat'), 'centered_data', '-v7.3');

% Create or extend boot table based on analysis type
if isempty(boot_table)
    fprintf('Creating boot table for %d bootstraps...\n', LIMO.design.bootstrap);
    switch analysis_type
        case {'one_sample', 'paired_samples'}
            if strcmp(analysis_type, 'paired_samples')
                boot_table = limo_create_boot_table(data1, LIMO.design.bootstrap);
            else
                boot_table = limo_create_boot_table(data, LIMO.design.bootstrap);
            end
        case 'two_samples'
            boot_table1 = limo_create_boot_table(data1, LIMO.design.bootstrap);
            boot_table2 = limo_create_boot_table(data2, LIMO.design.bootstrap);
            boot_table = boot_table1; % Primary boot table
            save(fullfile(LIMO.dir, 'H0', 'boot_table2.mat'), 'boot_table2', '-v7.3');
        case 'repeated_anova'
            boot_table = limo_create_boot_table(data(:,:,:,1), LIMO.design.bootstrap);
    end
elseif size(boot_table{1}, 2) < LIMO.design.bootstrap
    fprintf('Extending boot table from %d to %d bootstraps...\n', size(boot_table{1}, 2), LIMO.design.bootstrap);
    n_subjects = size(boot_table{1}, 1);

    for channel = 1:length(boot_table)
        current_boots = size(boot_table{channel}, 2);
        additional_boots = LIMO.design.bootstrap - current_boots;

        % Generate additional random indices
        new_boots = zeros(n_subjects, additional_boots);
        for b = 1:additional_boots
            new_boots(:,b) = randi(n_subjects, n_subjects, 1);
        end

        boot_table{channel} = [boot_table{channel}, new_boots];
    end

    % For two-samples, also extend boot_table2
    if strcmp(analysis_type, 'two_samples')
        boot_table2_file = fullfile(LIMO.dir, 'H0', 'boot_table2.mat');
        if exist(boot_table2_file, 'file')
            boot_table2_struct = load(boot_table2_file);
            boot_table2 = boot_table2_struct.boot_table2;

            for channel = 1:length(boot_table2)
                current_boots = size(boot_table2{channel}, 2);
                additional_boots = LIMO.design.bootstrap - current_boots;

                % Generate additional random indices
                new_boots = zeros(n_subjects, additional_boots);
                for b = 1:additional_boots
                    new_boots(:,b) = randi(n_subjects, n_subjects, 1);
                end

                boot_table2{channel} = [boot_table2{channel}, new_boots];
            end
            save(boot_table2_file, 'boot_table2', '-v7.3');
        end
    end
end

% Save boot table
save(fullfile(LIMO.dir, 'H0', 'boot_table.mat'), 'boot_table', '-v7.3');

% Process additional bootstraps in chunks
n_chunks = ceil(additional_bootstraps / chunk_size);
fprintf('Processing %d additional bootstraps in %d chunks...\n', additional_bootstraps, n_chunks);

% Initialize bootstrap result arrays based on analysis type
switch analysis_type
    case 'one_sample'
        H0_one_sample = NaN(size(centered_data,1), size(centered_data,2), 2, additional_bootstraps);
    case 'two_samples'
        H0_two_samples = NaN(size(centered_data{1},1), size(centered_data{1},2), 2, additional_bootstraps);
    case 'paired_samples'
        H0_paired_samples = NaN(size(centered_data{1},1), size(centered_data{1},2), 2, additional_bootstraps);
    case 'repeated_anova'
        if rep_type == 1
            tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1), size(centered_data,2), 1, 2, additional_bootstraps);
        elseif rep_type == 2
            tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1), size(centered_data,2), nb_effects, 2, additional_bootstraps);
        elseif rep_type == 3
            tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1), size(centered_data,2), 1, 2, additional_bootstraps);
            H0_Rep_ANOVA_Gp_effect = NaN(size(centered_data,1), size(centered_data,2), 2, additional_bootstraps);
            tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(centered_data,1), size(centered_data,2), 1, 2, additional_bootstraps);
        elseif rep_type == 4
            tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1), size(centered_data,2), nb_effects, 2, additional_bootstraps);
            H0_Rep_ANOVA_Gp_effect = NaN(size(centered_data,1), size(centered_data,2), 2, additional_bootstraps);
            tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(centered_data,1), size(centered_data,2), nb_effects, 2, additional_bootstraps);
        end
end

% Process bootstraps in chunks
for chunk = 1:n_chunks
    chunk_start_local = (chunk - 1) * chunk_size + 1;
    chunk_end_local = min(chunk * chunk_size, additional_bootstraps);

    % Global bootstrap indices (offset by existing bootstraps)
    chunk_start_global = current_bootstraps + chunk_start_local;
    chunk_end_global = current_bootstraps + chunk_end_local;

    fprintf('\nProcessing chunk %d/%d (local bootstraps %d-%d, global %d-%d)\n', ...
        chunk, n_chunks, chunk_start_local, chunk_end_local, chunk_start_global, chunk_end_global);

    % Debug: Check boot table dimensions
    fprintf('Boot table dimensions: %d channels, %d bootstraps\n', length(boot_table), size(boot_table{1}, 2));
    fprintf('B_global range will be: %d to %d\n', current_bootstraps + chunk_start_local, current_bootstraps + chunk_end_local);

    % Prepare variables for parallel processing
    local_chunk_size = chunk_end_local - chunk_start_local + 1;
    local_analysis_type = analysis_type;
    local_centered_data = centered_data;
    local_boot_table = boot_table;
    local_LIMO_method = LIMO.design.method;
    if exist('boot_table2', 'var')
        local_boot_table2 = boot_table2;
    else
        local_boot_table2 = [];
    end

    % Initialize variables for parallel access (all analysis types need these defined)
    local_gp_vector = [];
    local_C = [];
    local_X = [];
    local_factor_levels = [];
    local_rep_type = [];

    % Copy repeated anova variables for parallel access if needed
    if strcmp(analysis_type, 'repeated_anova')
        if exist('gp_vector', 'var')
            local_gp_vector = gp_vector;
        end
        if exist('C', 'var')
            local_C = C;
        end
        if exist('X', 'var')
            local_X = X;
        end
        if exist('factor_levels', 'var')
            local_factor_levels = factor_levels;
        end
        if exist('rep_type', 'var')
            local_rep_type = rep_type;
        end
    end

    % Process each bootstrap in this chunk using parallel processing
    % Use cell arrays to collect results from parfor
    bootstrap_results = cell(local_chunk_size, 1);

    fprintf('About to start parfor loop with %d iterations...\n', local_chunk_size);
    tic; % Start timing

    % Process bootstraps using the same logic as limo_random_robust.m
    for B_idx = 1:local_chunk_size
        B_local = chunk_start_local + B_idx - 1;
        B_global = current_bootstraps + B_local;

        if mod(B_idx, 5) == 0 || B_idx == 1
            fprintf('  Processing bootstrap %d/%d (global: %d)\n', B_idx, local_chunk_size, B_global);
        end

        local_result = struct();
        switch local_analysis_type
            case 'one_sample'
                local_result.t_vals = NaN(size(local_centered_data,1), size(local_centered_data,2));
                local_result.p_vals = NaN(size(local_centered_data,1), size(local_centered_data,2));
                % Process each channel for this bootstrap (matching limo_process_bootstrap_chunk)
                for channel = 1:size(local_centered_data,1)
                    tmp = local_centered_data(channel,:,:);
                    Y = tmp(1,:,find(~isnan(tmp(1,1,:))));

                    if ~isempty(Y)
                        % Get bootstrap sample using boot table indices
                        boot_sample = Y(1,:,local_boot_table{channel}(:,B_global));

                        % Compute one-sample t-test
                        if strcmpi(local_LIMO_method, 'Trimmed Mean')
                            [t_val, ~, ~, ~, p_val] = limo_trimci(boot_sample);
                        else
                            [~, ~, ~, ~, ~, t_val, p_val] = limo_ttest(1, boot_sample, 0, 0.05);
                        end

                        local_result.t_vals(channel,:) = t_val;
                        local_result.p_vals(channel,:) = p_val;
                    end
                end

            case 'two_samples'
                array = intersect(find(~isnan(local_centered_data{1}(:,1,1))), find(~isnan(local_centered_data{2}(:,1,1))));

                local_result.t_vals = NaN(size(local_centered_data{1},1), size(local_centered_data{1},2));
                local_result.p_vals = NaN(size(local_centered_data{1},1), size(local_centered_data{1},2));

                for e = 1:length(array)
                    channel = array(e);

                    % Get data for first group
                    tmp1 = local_centered_data{1}(channel,:,:);
                    Y1_full = tmp1(1,:,find(~isnan(tmp1(1,1,:))));
                    if ~isempty(Y1_full)
                        Y1 = Y1_full(1,:,local_boot_table{channel}(:,B_global));
                    else
                        continue;
                    end

                    % Get data for second group
                    tmp2 = local_centered_data{2}(channel,:,:);
                    Y2_full = tmp2(1,:,find(~isnan(tmp2(1,1,:))));
                    if ~isempty(Y2_full)
                        if ~isempty(local_boot_table2)
                            Y2 = Y2_full(1,:,local_boot_table2{channel}(:,B_global));
                        else
                            Y2 = Y2_full(1,:,local_boot_table{channel}(:,B_global));
                        end
                    else
                        continue;
                    end

                    % Compute two-samples t-test
                    if contains(local_LIMO_method, 'Trimmed Mean', 'IgnoreCase', true)
                        [t_val, ~, ~, ~, p_val] = limo_yuen_ttest(Y1, Y2);
                    else
                        [~, ~, ~, ~, ~, t_val, p_val] = limo_ttest(2, Y1, Y2, 0.05);
                    end

                    local_result.t_vals(channel,:) = t_val;
                    local_result.p_vals(channel,:) = p_val;
                end

            case 'paired_samples'
                array = intersect(find(~isnan(local_centered_data{1}(:,1,1))), find(~isnan(local_centered_data{2}(:,1,1))));

                local_result.t_vals = NaN(size(local_centered_data{1},1), size(local_centered_data{1},2));
                local_result.p_vals = NaN(size(local_centered_data{1},1), size(local_centered_data{1},2));

                for e = 1:length(array)
                    channel = array(e);

                    % Get data for both conditions (same bootstrap indices for paired data)
                    tmp1 = local_centered_data{1}(channel,:,:);
                    Y1_full = tmp1(1,:,find(~isnan(tmp1(1,1,:))));
                    tmp2 = local_centered_data{2}(channel,:,:);
                    Y2_full = tmp2(1,:,find(~isnan(tmp2(1,1,:))));

                    if ~isempty(Y1_full) && ~isempty(Y2_full)
                        Y1 = Y1_full(1,:,local_boot_table{channel}(:,B_global));
                        Y2 = Y2_full(1,:,local_boot_table{channel}(:,B_global));

                        % Compute paired t-test
                        if contains(local_LIMO_method, 'Trimmed Mean', 'IgnoreCase', true)
                            [t_val, ~, ~, ~, p_val] = limo_yuend_ttest(Y1, Y2);
                        else
                            [~, ~, ~, ~, ~, t_val, p_val] = limo_ttest(1, Y1, Y2, 0.05);
                        end

                        local_result.t_vals(channel,:) = t_val;
                        local_result.p_vals(channel,:) = p_val;
                    end
                end

            case 'repeated_anova'
                array = find(~isnan(local_centered_data(:,1,1,1)));

                if local_rep_type == 1
                    local_result.rep_F = NaN(size(local_centered_data,1), size(local_centered_data,2), 1);
                    local_result.rep_p = NaN(size(local_centered_data,1), size(local_centered_data,2), 1);
                elseif local_rep_type == 2
                    nb_effects_local = length(local_C);
                    local_result.rep_F = NaN(size(local_centered_data,1), size(local_centered_data,2), nb_effects_local);
                    local_result.rep_p = NaN(size(local_centered_data,1), size(local_centered_data,2), nb_effects_local);
                elseif local_rep_type == 3
                    local_result.rep_F = NaN(size(local_centered_data,1), size(local_centered_data,2), 1);
                    local_result.rep_p = NaN(size(local_centered_data,1), size(local_centered_data,2), 1);
                    local_result.gp_F = NaN(size(local_centered_data,1), size(local_centered_data,2));
                    local_result.gp_p = NaN(size(local_centered_data,1), size(local_centered_data,2));
                    local_result.int_F = NaN(size(local_centered_data,1), size(local_centered_data,2), 1);
                    local_result.int_p = NaN(size(local_centered_data,1), size(local_centered_data,2), 1);
                elseif local_rep_type == 4
                    nb_effects_local = length(local_C);
                    local_result.rep_F = NaN(size(local_centered_data,1), size(local_centered_data,2), nb_effects_local);
                    local_result.rep_p = NaN(size(local_centered_data,1), size(local_centered_data,2), nb_effects_local);
                    local_result.gp_F = NaN(size(local_centered_data,1), size(local_centered_data,2));
                    local_result.gp_p = NaN(size(local_centered_data,1), size(local_centered_data,2));
                    local_result.int_F = NaN(size(local_centered_data,1), size(local_centered_data,2), nb_effects_local);
                    local_result.int_p = NaN(size(local_centered_data,1), size(local_centered_data,2), nb_effects_local);
                end

                for e = 1:length(array)
                    channel = array(e);

                    % Get bootstrapped data for this channel
                    tmp = squeeze(local_centered_data(channel,:,local_boot_table{channel}(:,B_global),:));
                    if size(local_centered_data,2) == 1
                        Y = ones(1,size(tmp,1),size(tmp,2)); Y(1,:,:) = tmp;
                        gp = local_gp_vector(find(~isnan(Y(1,:,1))),:);
                        Y = Y(:,find(~isnan(Y(1,:,1))),:);
                    else
                        Y = tmp(:,find(~isnan(tmp(1,:,1))),:);
                        gp = local_gp_vector(find(~isnan(tmp(1,:,1))));
                    end

                    if local_rep_type == 3 || local_rep_type == 4
                        XB = local_X(find(~isnan(tmp(1,:,1))),:);
                    else
                        XB = [];
                    end

                    % Compute analysis for this bootstrap sample
                    if local_rep_type == 1
                        if contains(local_LIMO_method, 'Trimmed Mean', 'IgnoreCase', true)
                            result = limo_robust_rep_anova(Y, gp, local_factor_levels, local_C);
                        else
                            result = limo_rep_anova(Y, gp, local_factor_levels, local_C);
                        end
                        local_result.rep_F(channel,:,1) = result.F;
                        local_result.rep_p(channel,:,1) = result.p;
                    elseif local_rep_type == 2
                        if contains(local_LIMO_method, 'Trimmed Mean', 'IgnoreCase', true)
                            result = limo_robust_rep_anova(Y, gp, local_factor_levels, local_C);
                        else
                            result = limo_rep_anova(Y, gp, local_factor_levels, local_C);
                        end
                        local_result.rep_F(channel,:,:) = result.F';
                        local_result.rep_p(channel,:,:) = result.p';
                    elseif local_rep_type == 3
                        if contains(local_LIMO_method, 'Trimmed Mean', 'IgnoreCase', true)
                            result = limo_robust_rep_anova(Y, gp, local_factor_levels, local_C, XB);
                        else
                            result = limo_rep_anova(Y, gp, local_factor_levels, local_C, XB);
                        end
                        local_result.rep_F(channel,:,1) = result.repeated_measure.F;
                        local_result.rep_p(channel,:,1) = result.repeated_measure.p;
                        local_result.gp_F(channel,:) = result.gp.F;
                        local_result.gp_p(channel,:) = result.gp.p;
                        local_result.int_F(channel,:,1) = result.interaction.F;
                        local_result.int_p(channel,:,1) = result.interaction.p;
                    elseif local_rep_type == 4
                        if contains(local_LIMO_method, 'Trimmed Mean', 'IgnoreCase', true)
                            result = limo_robust_rep_anova(Y, gp, local_factor_levels, local_C, XB);
                        else
                            result = limo_rep_anova(Y, gp, local_factor_levels, local_C, XB);
                        end
                        local_result.rep_F(channel,:,:) = result.repeated_measure.F';
                        local_result.rep_p(channel,:,:) = result.repeated_measure.p';
                        local_result.gp_F(channel,:) = result.gp.F;
                        local_result.gp_p(channel,:) = result.gp.p;
                        local_result.int_F(channel,:,:) = result.interaction.F';
                        local_result.int_p(channel,:,:) = result.interaction.p';
                    end
                end
        end

        bootstrap_results{B_idx} = local_result;
    end

    parfor_time = toc; % End timing
    fprintf('Parfor loop completed in %.2f seconds\n', parfor_time);

    % Copy chunk results back to global arrays from parfor results
    switch analysis_type
        case 'one_sample'
            for B_idx = 1:local_chunk_size
                result = bootstrap_results{B_idx};
                H0_one_sample(:,:,1,chunk_start_local+B_idx-1) = result.t_vals;
                H0_one_sample(:,:,2,chunk_start_local+B_idx-1) = result.p_vals;
            end
        case 'two_samples'
            for B_idx = 1:local_chunk_size
                result = bootstrap_results{B_idx};
                H0_two_samples(:,:,1,chunk_start_local+B_idx-1) = result.t_vals;
                H0_two_samples(:,:,2,chunk_start_local+B_idx-1) = result.p_vals;
            end
        case 'paired_samples'
            for B_idx = 1:local_chunk_size
                result = bootstrap_results{B_idx};
                H0_paired_samples(:,:,1,chunk_start_local+B_idx-1) = result.t_vals;
                H0_paired_samples(:,:,2,chunk_start_local+B_idx-1) = result.p_vals;
            end
        case 'repeated_anova'
            for B_idx = 1:local_chunk_size
                result = bootstrap_results{B_idx};
                if isfield(result, 'rep_F')
                    tmp_boot_H0_Rep_ANOVA(:,:,:,1,chunk_start_local+B_idx-1) = result.rep_F;
                    tmp_boot_H0_Rep_ANOVA(:,:,:,2,chunk_start_local+B_idx-1) = result.rep_p;
                end
                if rep_type == 3 || rep_type == 4
                    if isfield(result, 'gp_F')
                        H0_Rep_ANOVA_Gp_effect(:,:,1,chunk_start_local+B_idx-1) = result.gp_F;
                        H0_Rep_ANOVA_Gp_effect(:,:,2,chunk_start_local+B_idx-1) = result.gp_p;
                    end
                    if isfield(result, 'int_F')
                        tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(:,:,:,1,chunk_start_local+B_idx-1) = result.int_F;
                        tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(:,:,:,2,chunk_start_local+B_idx-1) = result.int_p;
                    end
                end
            end
    end

    % Save chunk results with new filename format
    switch analysis_type
        case 'one_sample'
            base_name = sprintf('H0_one_sample_ttest_parameter_%s', num2str(parameter));
            chunk_data = H0_one_sample(:,:,:,chunk_start_local:chunk_end_local);
            filename = sprintf('%s_chunk_%d-%d.mat', base_name, chunk_start_global, chunk_end_global);
            S = struct('H0_one_sample', chunk_data);
            save(fullfile(chunk_dir, filename), '-struct', 'S', '-v7.3');
            fprintf('  Saved chunk to %s\n', filename);

        case 'two_samples'
            base_name = sprintf('H0_two_samples_ttest_parameter_%s', num2str(parameter));
            chunk_data = H0_two_samples(:,:,:,chunk_start_local:chunk_end_local);
            filename = sprintf('%s_chunk_%d-%d.mat', base_name, chunk_start_global, chunk_end_global);
            S = struct('H0_two_samples', chunk_data);
            save(fullfile(chunk_dir, filename), '-struct', 'S', '-v7.3');
            fprintf('  Saved chunk to %s\n', filename);

        case 'paired_samples'
            base_name = sprintf('H0_paired_samples_ttest_parameter_%s', num2str(parameter));
            chunk_data = H0_paired_samples(:,:,:,chunk_start_local:chunk_end_local);
            filename = sprintf('%s_chunk_%d-%d.mat', base_name, chunk_start_global, chunk_end_global);
            S = struct('H0_paired_samples', chunk_data);
            save(fullfile(chunk_dir, filename), '-struct', 'S', '-v7.3');
            fprintf('  Saved chunk to %s\n', filename);

        case 'repeated_anova'
            % Save main repeated ANOVA chunk
            base_name = 'H0_Rep_ANOVA';
            chunk_data = tmp_boot_H0_Rep_ANOVA(:,:,:,:,chunk_start_local:chunk_end_local);
            filename = sprintf('%s_chunk_%d-%d.mat', base_name, chunk_start_global, chunk_end_global);
            S = struct('H0_Rep_ANOVA', chunk_data);
            save(fullfile(chunk_dir, filename), '-struct', 'S', '-v7.3');
            fprintf('  Saved chunk to %s\n', filename);

            if rep_type == 3 || rep_type == 4
                % Save group effect chunk
                gp_base_name = 'H0_Rep_ANOVA_Gp_effect';
                gp_chunk_data = H0_Rep_ANOVA_Gp_effect(:,:,:,chunk_start_local:chunk_end_local);
                filename = sprintf('%s_chunk_%d-%d.mat', gp_base_name, chunk_start_global, chunk_end_global);
                S = struct('H0_Rep_ANOVA_Gp_effect', gp_chunk_data);
                save(fullfile(chunk_dir, filename), '-struct', 'S', '-v7.3');
                fprintf('  Saved chunk to %s\n', filename);

                % Save interaction chunk
                int_base_name = 'H0_Rep_ANOVA_Interaction_with_gp';
                int_chunk_data = tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(:,:,:,:,chunk_start_local:chunk_end_local);
                filename = sprintf('%s_chunk_%d-%d.mat', int_base_name, chunk_start_global, chunk_end_global);
                S = struct('H0_Rep_ANOVA_Interaction_with_gp', int_chunk_data);
                save(fullfile(chunk_dir, filename), '-struct', 'S', '-v7.3');
                fprintf('  Saved chunk to %s\n', filename);
            end
    end

    fprintf('  Chunk %d/%d completed\n', chunk, n_chunks);
end

fprintf('\n=== ADDITIONAL BOOTSTRAPS COMPLETE ===\n');
fprintf('Total bootstraps now available: %d\n', current_bootstraps + additional_bootstraps);

% Always merge chunks into final H0 files
fprintf('\n=== MERGING ALL CHUNKS ===\n');

switch analysis_type
    case 'one_sample'
        fprintf('Merging one-sample t-test chunks...\n');
        existing_file = fullfile(LIMO.dir, 'H0', sprintf('H0_one_sample_ttest_parameter_%s.mat', num2str(parameter)));
        chunk_pattern = sprintf('H0_one_sample_ttest_parameter_%s', num2str(parameter));

        merged_data = limo_merge_additional_chunks(chunk_dir, chunk_pattern, ...
            existing_file, current_bootstraps, additional_bootstraps, 1);

        if ~isempty(merged_data)
            H0_one_sample = merged_data;
            save(existing_file, 'H0_one_sample', '-v7.3');
            fprintf('  Saved merged data to %s\n', existing_file);
        end

    case 'two_samples'
        fprintf('Merging two-samples t-test chunks...\n');
        existing_file = fullfile(LIMO.dir, 'H0', sprintf('H0_two_samples_ttest_parameter_%s.mat', num2str(parameter)));
        chunk_pattern = sprintf('H0_two_samples_ttest_parameter_%s', num2str(parameter));

        merged_data = limo_merge_additional_chunks(chunk_dir, chunk_pattern, ...
            existing_file, current_bootstraps, additional_bootstraps, 1);

        if ~isempty(merged_data)
            H0_two_samples = merged_data;
            save(existing_file, 'H0_two_samples', '-v7.3');
            fprintf('  Saved merged data to %s\n', existing_file);
        end

    case 'paired_samples'
        fprintf('Merging paired t-test chunks...\n');
        existing_file = fullfile(LIMO.dir, 'H0', sprintf('H0_paired_samples_ttest_parameter_%s.mat', num2str(parameter)));
        chunk_pattern = sprintf('H0_paired_samples_ttest_parameter_%s', num2str(parameter));

        merged_data = limo_merge_additional_chunks(chunk_dir, chunk_pattern, ...
            existing_file, current_bootstraps, additional_bootstraps, 1);

        if ~isempty(merged_data)
            H0_paired_samples = merged_data;
            save(existing_file, 'H0_paired_samples', '-v7.3');
            fprintf('  Saved merged data to %s\n', existing_file);
        end

    case 'repeated_anova'
        % Merge main effect files
        for i = 1:nb_effects
            fprintf('Merging effect %d/%d: %s...\n', i, nb_effects, Rep_filenames{i});

            % Find existing bootstrap file
            existing_file = fullfile(LIMO.dir, 'H0', sprintf('H0_%s', Rep_filenames{i}));

            % Merge chunks for this effect
            merged_data = limo_merge_additional_chunks(chunk_dir, 'H0_Rep_ANOVA', ...
                existing_file, current_bootstraps, additional_bootstraps, i);

            if ~isempty(merged_data)
                H0_Rep_ANOVA = merged_data;
                save(existing_file, 'H0_Rep_ANOVA', '-v7.3');
                fprintf('  Saved merged data to %s\n', existing_file);
            end
        end

        if rep_type == 3 || rep_type == 4
            % Merge group effect
            fprintf('Merging group effect...\n');
            existing_gp_file = fullfile(LIMO.dir, 'H0', 'H0_Rep_ANOVA_Gp_effect.mat');
            merged_gp_data = limo_merge_additional_chunks(chunk_dir, 'H0_Rep_ANOVA_Gp_effect', ...
                existing_gp_file, current_bootstraps, additional_bootstraps, 1);

            if ~isempty(merged_gp_data)
                H0_Rep_ANOVA_Gp_effect = merged_gp_data;
                save(existing_gp_file, 'H0_Rep_ANOVA_Gp_effect', '-v7.3');
                fprintf('  Saved merged group effect data\n');
            end

            % Merge interaction effects
            for i = 1:nb_effects
                fprintf('Merging interaction effect %d/%d: %s...\n', i, nb_effects, IRep_filenames{i});
                existing_int_file = fullfile(LIMO.dir, 'H0', sprintf('H0_%s', IRep_filenames{i}));
                merged_int_data = limo_merge_additional_chunks(chunk_dir, 'H0_Rep_ANOVA_Interaction_with_gp', ...
                    existing_int_file, current_bootstraps, additional_bootstraps, i);

                if ~isempty(merged_int_data)
                    H0_Rep_ANOVA_Interaction_with_gp = merged_int_data;
                    save(existing_int_file, 'H0_Rep_ANOVA_Interaction_with_gp', '-v7.3');
                    fprintf('  Saved merged interaction data\n');
                end
            end
        end
end

fprintf('Merging complete!\n');
fprintf('Chunks preserved in: %s\n', chunk_dir);

end

% =========================================================================
% HELPER FUNCTION FOR MERGING ADDITIONAL BOOTSTRAP CHUNKS
% =========================================================================

function merged_data = limo_merge_additional_chunks(chunk_dir, chunk_pattern, existing_file, existing_bootstraps, additional_bootstraps, effect_index)
    % Merge additional bootstrap chunks with existing bootstrap data
    %
    % INPUTS:
    %   chunk_dir - directory containing chunk files
    %   chunk_pattern - pattern for chunk file names (e.g., 'H0_Rep_ANOVA')
    %   existing_file - path to existing bootstrap file
    %   existing_bootstraps - number of existing bootstraps
    %   additional_bootstraps - number of new bootstraps to add
    %   effect_index - which effect to extract from multi-effect chunks
    %
    % OUTPUT:
    %   merged_data - combined existing + new bootstrap data

    merged_data = [];

    % Find chunk files for this pattern
    chunk_files = dir(fullfile(chunk_dir, sprintf('%s_chunk_*.mat', chunk_pattern)));
    if isempty(chunk_files)
        fprintf('  Warning: No chunk files found for pattern %s\n', chunk_pattern);
        return;
    end

    fprintf('  Found %d chunk files to merge\n', length(chunk_files));

    % Load existing data if file exists
    existing_data = [];
    if exist(existing_file, 'file')
        fprintf('  Loading existing bootstrap file...\n');
        existing_struct = load(existing_file);
        field_names = fieldnames(existing_struct);

        % Find the bootstrap variable
        bootstrap_var = '';
        for i = 1:length(field_names)
            if contains(field_names{i}, 'H0_')
                bootstrap_var = field_names{i};
                break;
            end
        end

        if ~isempty(bootstrap_var)
            existing_data = existing_struct.(bootstrap_var);
            fprintf('  Loaded existing data with %d bootstraps\n', size(existing_data, ndims(existing_data)));
        end
    end

    % Load and concatenate chunk data
    % Pre-scan chunks to determine total size for accumulator
    total_bootstraps_to_merge = 0;
    if ~isempty(chunk_files)
        for i = 1:length(chunk_files)
            try
                chunk_info = whos('-file', fullfile(chunk_files(i).folder, chunk_files(i).name));
                dims = chunk_info(1).size; % Assuming one variable per chunk file
                total_bootstraps_to_merge = total_bootstraps_to_merge + dims(end);
            catch ME
                 fprintf('  Warning: Could not read chunk %s for sizing: %s\n', chunk_files(i).name, ME.message);
            end
        end
    end

    new_chunks_data = [];
    total_new_bootstraps = 0;
    current_bootstrap_idx = 1;

    for i = 1:length(chunk_files)
        chunk_file = fullfile(chunk_files(i).folder, chunk_files(i).name);
        fprintf('  Loading chunk %d/%d: %s\n', i, length(chunk_files), chunk_files(i).name);

        try
            chunk_struct = load(chunk_file);

            % Determine the variable name in the chunk
            % The limo_save_boot_chunks function saves data using the var_name directly
            field_names = fieldnames(chunk_struct);

            % Find the data variable (should be the H0_* variable)
            chunk_data = [];
            for j = 1:length(field_names)
                if contains(field_names{j}, 'H0_')
                    chunk_data = chunk_struct.(field_names{j});
                    break;
                end
            end

            % Fallback to old naming if no H0_ variable found
            if isempty(chunk_data)
                if contains(chunk_pattern, 'Gp_effect')
                    chunk_data = chunk_struct.gp_chunk_data;
                elseif contains(chunk_pattern, 'Interaction')
                    chunk_data = chunk_struct.int_chunk_data;
                else
                    chunk_data = chunk_struct.chunk_data;
                end
            end

            % Extract the specific effect if multi-effect data
            if ndims(chunk_data) >= 4 && size(chunk_data, 3) > 1 && exist('effect_index', 'var') && ~contains(chunk_pattern, 'ttest')
                if ndims(chunk_data) == 5
                    % Format: (channels, time, effects, stats, bootstraps)
                    chunk_data = squeeze(chunk_data(:, :, effect_index, :, :));
                elseif ndims(chunk_data) == 4 && ~contains(chunk_pattern, 'Gp_effect')
                    % Format: (channels, time, effects, bootstraps) for stats
                    chunk_data = squeeze(chunk_data(:, :, effect_index, :));
                end
            end

            % Debug: Check chunk data dimensions
            fprintf('    Loaded chunk with dimensions: [%s]\n', num2str(size(chunk_data)));

            % Accumulate chunk data
            if isempty(new_chunks_data)
                chunk_dims = size(chunk_data);
                chunk_dims(end) = total_bootstraps_to_merge; % Set total expected bootstraps
                new_chunks_data = NaN(chunk_dims);
                fprintf('    Initialized accumulator with dimensions: [%s]\n', num2str(size(new_chunks_data)));
            end

            chunk_size = size(chunk_data, ndims(chunk_data));
            end_idx = current_bootstrap_idx + chunk_size - 1;

            % Debug: Check dimension compatibility
            if ndims(chunk_data) ~= ndims(new_chunks_data)
                fprintf('    Warning: Dimension mismatch - chunk: %dD, accumulator: %dD\n', ...
                    ndims(chunk_data), ndims(new_chunks_data));
            end

            % Store chunk data with proper dimension handling
            if ndims(new_chunks_data) == 4 && ndims(chunk_data) == 4
                new_chunks_data(:, :, :, current_bootstrap_idx:end_idx) = chunk_data;
            elseif ndims(new_chunks_data) == 3 && ndims(chunk_data) == 3
                new_chunks_data(:, :, current_bootstrap_idx:end_idx) = chunk_data;
            elseif ndims(new_chunks_data) == 2 && ndims(chunk_data) == 2
                new_chunks_data(:, current_bootstrap_idx:end_idx) = chunk_data;
            else
                % Handle dimension mismatches
                if ndims(new_chunks_data) == 4 && ndims(chunk_data) == 3
                    % Try to reshape chunk_data to 4D by assuming missing dimension is 2 (T and p)
                    if size(chunk_data, 3) == chunk_size
                        % Data is [channels, frames, bootstraps] - need to make it [channels, frames, 2, bootstraps/2]
                        fprintf('    Attempting to reshape 3D chunk to 4D format\n');
                        if mod(chunk_size, 2) == 0
                            reshaped_chunk = reshape(chunk_data, [size(chunk_data,1), size(chunk_data,2), 2, chunk_size/2]);
                            new_chunk_size = chunk_size/2;
                            new_end_idx = current_bootstrap_idx + new_chunk_size - 1;
                            new_chunks_data(:, :, :, current_bootstrap_idx:new_end_idx) = reshaped_chunk;
                            end_idx = new_end_idx; % Update end_idx
                        else
                            error('Cannot reshape 3D chunk to 4D - bootstrap dimension not divisible by 2');
                        end
                    else
                        error('Incompatible dimensions between chunk (%dD) and accumulator (%dD)', ...
                            ndims(chunk_data), ndims(new_chunks_data));
                    end
                else
                    error('Incompatible dimensions between chunk (%dD) and accumulator (%dD)', ...
                        ndims(chunk_data), ndims(new_chunks_data));
                end
            end

            current_bootstrap_idx = end_idx + 1;
            total_new_bootstraps = total_new_bootstraps + chunk_size;

        catch ME
            fprintf('  Warning: Failed to load chunk %s: %s\n', chunk_files(i).name, ME.message);
        end
    end

    % Combine existing and new data
    if ~isempty(existing_data) && ~isempty(new_chunks_data)
        fprintf('  Combining existing (%d) and new (%d) bootstraps\n', ...
            size(existing_data, ndims(existing_data)), total_new_bootstraps);

        % Determine final dimensions
        if ndims(existing_data) == ndims(new_chunks_data)
            final_dims = size(existing_data);
            final_dims(end) = existing_bootstraps + total_new_bootstraps;
            merged_data = NaN(final_dims);

            % Copy existing data
            if ndims(merged_data) == 4
                merged_data(:, :, :, 1:existing_bootstraps) = existing_data;
                merged_data(:, :, :, existing_bootstraps+1:end) = new_chunks_data(:, :, :, 1:total_new_bootstraps);
            elseif ndims(merged_data) == 3
                merged_data(:, :, 1:existing_bootstraps) = existing_data;
                merged_data(:, :, existing_bootstraps+1:end) = new_chunks_data(:, :, 1:total_new_bootstraps);
            else
                merged_data(:, 1:existing_bootstraps) = existing_data;
                merged_data(:, existing_bootstraps+1:end) = new_chunks_data(:, 1:total_new_bootstraps);
            end
        else
            fprintf('  Warning: Dimension mismatch between existing and new data\n');
            fprintf('  Existing data dimensions: [%s], New data dimensions: [%s]\n', ...
                num2str(size(existing_data)), num2str(size(new_chunks_data)));
            % Handle variable dimensions
            if ndims(new_chunks_data) == 4
                actual_bootstraps = min(total_new_bootstraps, size(new_chunks_data, 4));
                merged_data = new_chunks_data(:, :, :, 1:actual_bootstraps);
            elseif ndims(new_chunks_data) == 3
                actual_bootstraps = min(total_new_bootstraps, size(new_chunks_data, 3));
                merged_data = new_chunks_data(:, :, 1:actual_bootstraps);
            else
                actual_bootstraps = min(total_new_bootstraps, size(new_chunks_data, 2));
                merged_data = new_chunks_data(:, 1:actual_bootstraps);
            end
            fprintf('  Using %d bootstraps from new data\n', actual_bootstraps);
        end

    elseif ~isempty(new_chunks_data)
        fprintf('  Using only new bootstrap data (%d bootstraps)\n', total_new_bootstraps);
        fprintf('  New data dimensions: [%s]\n', num2str(size(new_chunks_data)));
        % Handle variable dimensions
        if ndims(new_chunks_data) == 4
            actual_bootstraps = min(total_new_bootstraps, size(new_chunks_data, 4));
            merged_data = new_chunks_data(:, :, :, 1:actual_bootstraps);
        elseif ndims(new_chunks_data) == 3
            actual_bootstraps = min(total_new_bootstraps, size(new_chunks_data, 3));
            merged_data = new_chunks_data(:, :, 1:actual_bootstraps);
        else
            actual_bootstraps = min(total_new_bootstraps, size(new_chunks_data, 2));
            merged_data = new_chunks_data(:, 1:actual_bootstraps);
        end
        fprintf('  Final merged data dimensions: [%s]\n', num2str(size(merged_data)));
        fprintf('  Using %d bootstraps from new data\n', actual_bootstraps);
    else
        fprintf('  Warning: No valid chunk data found\n');
        return;
    end

    fprintf('  Successfully merged data: final size = [%s]\n', num2str(size(merged_data)));
end

% =========================================================================
% HELPER FUNCTION FOR FINDING LIMO DATA FILES
% =========================================================================

function data_files = find_limo_data_files(limo_dir, base_name)
    % FIND_LIMO_DATA_FILES - Find LIMO data files with various naming conventions
    %
    % INPUTS:
    %   limo_dir  - Directory to search in
    %   base_name - Base name pattern (e.g., 'Y', 'Y1', 'Y2')
    %
    % OUTPUT:
    %   data_files - Directory listing of found files

    data_files = [];

    % Try different naming conventions
    if strcmp(base_name, 'Y')
        % For one-sample: Y*.mat, Yr.mat
        patterns = {'Y*.mat', 'Yr.mat'};
    elseif strcmp(base_name, 'Y1')
        % For first group: Y1*.mat, Y1r.mat, Yr1.mat
        patterns = {'Y1*.mat', 'Y1r.mat', 'Yr1.mat'};
    elseif strcmp(base_name, 'Y2')
        % For second group: Y2*.mat, Y2r.mat, Yr2.mat
        patterns = {'Y2*.mat', 'Y2r.mat', 'Yr2.mat'};
    else
        % Generic pattern
        patterns = {[base_name '*.mat'], [base_name 'r.mat']};
    end

    % Try each pattern until we find files
    for i = 1:length(patterns)
        data_files = dir(fullfile(limo_dir, patterns{i}));
        if ~isempty(data_files)
            fprintf('Found data files using pattern: %s\n', patterns{i});
            break;
        end
    end
end

