function chunk_results = limo_process_bootstrap_chunk(analysis_type, centered_data, ...
    boot_indices, chunk_start, chunk_end, LIMO, options)
% LIMO_PROCESS_BOOTSTRAP_CHUNK - Process a single chunk of bootstrap iterations
%
% INPUTS:
%   analysis_type - Integer 1-6 corresponding to analysis types in limo_random_robust
%   centered_data - Data centered under H0
%   boot_indices  - Bootstrap indices table  
%   chunk_start   - Starting bootstrap index
%   chunk_end     - Ending bootstrap index
%   LIMO          - LIMO structure with design info
%   options       - Struct with analysis-specific parameters
%                   .parameter (for t-tests)
%                   .C (contrast matrix for ANOVA)
%                   .factor_levels (for repeated measures)
%                   .gp_vector (group vector)
%                   .X (design matrix for between-subjects effects)
%
% OUTPUT:
%   chunk_results - Struct containing bootstrap results for this chunk

% Initialize output structure
chunk_results = struct();
n_bootstraps = chunk_end - chunk_start + 1;

% Get data dimensions
if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
    is_tf = true;
else
    is_tf = false;
end

switch analysis_type
    case 1 % One-sample t-test
        % Pre-allocate
        H0_one_sample = NaN(size(centered_data,1), size(centered_data,2), 2, n_bootstraps);
        
        % Process each bootstrap
        for b = 1:n_bootstraps
            actual_boot = chunk_start + b - 1;
            if mod(b, 10) == 0
                fprintf('  Bootstrap %d/%d (global: %d)\n', b, n_bootstraps, actual_boot);
            end
            
            % Process each channel
            for channel = 1:size(centered_data,1)
                tmp = centered_data(channel,:,:);
                Y = tmp(1,:,find(~isnan(tmp(1,1,:))));
                
                if ~isempty(Y)
                    boot_sample = Y(1,:,boot_indices{channel}(:,actual_boot));
                    
                    if strcmpi(LIMO.design.method,'Trimmed Mean')
                        [t_val,~,~,~,p_val,~,~] = limo_trimci(boot_sample);
                    elseif strcmpi(LIMO.design.method,'Mean')
                        [~,~,~,~,~,t_val,p_val] = limo_ttest(1,boot_sample,0,5/100);
                    end
                    
                    H0_one_sample(channel,:,1,b) = t_val;
                    H0_one_sample(channel,:,2,b) = p_val;
                end
            end
        end
        
        chunk_results.H0_one_sample = H0_one_sample;
        chunk_results.var_name = sprintf('H0_one_sample_ttest_parameter_%g', options.parameter);
        
    case 2 % Two-samples t-test
        % Pre-allocate
        H0_two_samples = NaN(size(centered_data{1},1), size(centered_data{1},2), 2, n_bootstraps);
        data1_centered = centered_data{1};
        data2_centered = centered_data{2};
        boot_table1 = options.boot_table1;
        boot_table2 = options.boot_table2;
        
        % Process each bootstrap
        for b = 1:n_bootstraps
            actual_boot = chunk_start + b - 1;
            if mod(b, 10) == 0
                fprintf('  Bootstrap %d/%d (global: %d)\n', b, n_bootstraps, actual_boot);
            end
            
            array = intersect(find(~isnan(data1_centered(:,1,1))),find(~isnan(data2_centered(:,1,1))));
            for e = 1:size(array,1)
                channel = array(e);
                tmp = data1_centered(channel,:,:); 
                Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); 
                tmp = data2_centered(channel,:,:); 
                Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); 
                
                if ~isempty(Y1) && ~isempty(Y2)
                    boot_Y1 = Y1(1,:,boot_table1{channel}(:,actual_boot));
                    boot_Y2 = Y2(1,:,boot_table2{channel}(:,actual_boot));
                    
                    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true) || ...
                            contains(LIMO.design.method,'Welch','IgnoreCase',true)
                        [t_val,~,~,~,p_val,~,~] = limo_yuen_ttest(boot_Y1,boot_Y2);
                    else
                        [~,~,~,~,~,t_val,p_val] = limo_ttest(2,boot_Y1,boot_Y2,.05);
                    end
                    
                    H0_two_samples(channel,:,1,b) = t_val;
                    H0_two_samples(channel,:,2,b) = p_val;
                end
            end
        end
        
        chunk_results.H0_two_samples = H0_two_samples;
        chunk_results.var_name = sprintf('H0_two_samples_ttest_parameter_%g_%g', options.parameter);
        
    case 3 % Paired t-test
        % Pre-allocate
        H0_paired_samples = NaN(size(centered_data{1},1), size(centered_data{1},2), 2, n_bootstraps);
        data1_centered = centered_data{1};
        data2_centered = centered_data{2};
        
        % Process each bootstrap
        for b = 1:n_bootstraps
            actual_boot = chunk_start + b - 1;
            if mod(b, 10) == 0
                fprintf('  Bootstrap %d/%d (global: %d)\n', b, n_bootstraps, actual_boot);
            end
            
            array = intersect(find(~isnan(data1_centered(:,1,1))),find(~isnan(data2_centered(:,1,1))));
            for e = 1:size(array,1)
                channel = array(e);
                tmp = data1_centered(channel,:,:); 
                Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); 
                tmp = data2_centered(channel,:,:); 
                Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); 
                
                if ~isempty(Y1) && ~isempty(Y2)
                    boot_Y1 = Y1(1,:,boot_indices{channel}(:,actual_boot));
                    boot_Y2 = Y2(1,:,boot_indices{channel}(:,actual_boot));
                    
                    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                        [t_val,~,~,~,p_val,~,~] = limo_yuend_ttest(boot_Y1,boot_Y2);
                    else
                        [~,~,~,~,~,t_val,p_val] = limo_ttest(1,boot_Y1,boot_Y2);
                    end
                    
                    H0_paired_samples(channel,:,1,b) = t_val;
                    H0_paired_samples(channel,:,2,b) = p_val;
                end
            end
        end
        
        chunk_results.H0_paired_samples = H0_paired_samples;
        chunk_results.var_name = sprintf('H0_paired_samples_ttest_parameter_%s', num2str(options.parameter')');
        
    case 4 % Regression
        % Call limo_glm_boot for regression bootstrapping
        % This requires the Yr file and X matrix
        nboot = n_bootstraps;
        boot_indices_chunk = boot_indices(:, chunk_start:chunk_end);
        
        % The regression bootstrap is handled by limo_glm_boot
        % We need to adapt it to work with chunks
        fprintf('Processing regression bootstrap chunk (bootstraps %d-%d)\n', chunk_start, chunk_end);
        
        % Note: This requires modification of limo_glm_boot to accept chunk parameters
        % For now, we'll use a simplified approach
        chunk_results.H0_regression = [];
        chunk_results.var_name = 'H0_regression';
        warning('Regression chunking requires modification of limo_glm_boot');
        
    case 5 % N-way ANOVA/ANCOVA
        if LIMO.design.fullfactorial == 0 && LIMO.design.nb_continuous == 0
            % Robust 1-way ANOVA
            H0_Condition_effect = NaN(size(centered_data,1), size(centered_data,2), 2, n_bootstraps);
            
            array = find(~isnan(centered_data(:,1,1)));
            for b = 1:n_bootstraps
                actual_boot = chunk_start + b - 1;
                if mod(b, 10) == 0
                    fprintf('  Bootstrap %d/%d (global: %d)\n', b, n_bootstraps, actual_boot);
                end
                
                for channel = 1:size(array,1)
                    e = array(channel);
                    index = find(~isnan(squeeze(centered_data(e,1,:))));
                    X = LIMO.design.X(index,1:end-1);
                    if sum(sum(X) == 0) == 0
                        boot_data = squeeze(centered_data(e,:,boot_indices{e}(:,actual_boot)));
                        [H0_Condition_effect(e,:,1,b), H0_Condition_effect(e,:,2,b)] = ...
                            limo_robust_1way_anova(boot_data,X,20);
                    end
                end
            end
            
            chunk_results.H0_Condition_effect = H0_Condition_effect;
            chunk_results.var_name = 'H0_Condition_effect_1';
        else
            % Standard ANOVA/ANCOVA - requires limo_glm_boot
            chunk_results.H0_anova = [];
            chunk_results.var_name = 'H0_anova';
            warning('Standard ANOVA/ANCOVA chunking requires modification of limo_glm_boot');
        end
        
    case 6 % Repeated measures ANOVA
        % Extract options
        C = options.C;
        factor_levels = options.factor_levels;
        gp_vector = options.gp_vector;
        type = options.type;
        
        % Determine array dimensions based on type
        if type == 1
            tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1), size(centered_data,2), 1, 2, n_bootstraps);
        elseif type == 2
            tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1), size(centered_data,2), length(C), 2, n_bootstraps);
        elseif type == 3
            tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1), size(centered_data,2), 1, 2, n_bootstraps);
            H0_Rep_ANOVA_Gp_effect = NaN(size(centered_data,1), size(centered_data,2), 2, n_bootstraps);
            tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(centered_data,1), size(centered_data,2), 1, 2, n_bootstraps);
        else % type == 4
            tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1), size(centered_data,2), length(C), 2, n_bootstraps);
            H0_Rep_ANOVA_Gp_effect = NaN(size(centered_data,1), size(centered_data,2), 2, n_bootstraps);
            tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(centered_data,1), size(centered_data,2), length(C), 2, n_bootstraps);
        end
        
        % Process bootstraps
        array = find(~isnan(centered_data(:,1,1,1)));
        for b = 1:n_bootstraps
            actual_boot = chunk_start + b - 1;
            if mod(b, 10) == 0
                fprintf('  Bootstrap %d/%d (global: %d)\n', b, n_bootstraps, actual_boot);
            end
            
            for e = 1:length(array)
                channel = array(e);
                
                % Get bootstrapped data
                tmp = squeeze(centered_data(channel,:,boot_indices{channel}(:,actual_boot),:));
                
                if size(centered_data,2) == 1
                    Y = ones(1,size(tmp,1),size(tmp,2));
                    Y(1,:,:) = tmp;
                    gp = gp_vector(find(~isnan(Y(1,:,1))),:);
                    Y = Y(:,find(~isnan(Y(1,:,1))),:);
                else
                    Y = tmp(:,find(~isnan(tmp(1,:,1))),:);
                    gp = gp_vector(find(~isnan(tmp(1,:,1))));
                end
                
                if type == 3 || type == 4
                    X = options.X;
                    XB = X(find(~isnan(tmp(1,:,1))),:);
                else
                    XB = [];
                end
                
                % Run analysis based on type
                if type == 1
                    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C);
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C);
                    end
                    tmp_boot_H0_Rep_ANOVA(channel,:,1,1,b) = result.F;
                    tmp_boot_H0_Rep_ANOVA(channel,:,1,2,b) = result.p;
                    
                elseif type == 2
                    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C);
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C);
                    end
                    tmp_boot_H0_Rep_ANOVA(channel,:,:,1,b) = result.F';
                    tmp_boot_H0_Rep_ANOVA(channel,:,:,2,b) = result.p';
                    
                elseif type == 3
                    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB);
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C,XB);
                    end
                    tmp_boot_H0_Rep_ANOVA(channel,:,1,1,b) = result.repeated_measure.F;
                    tmp_boot_H0_Rep_ANOVA(channel,:,1,2,b) = result.repeated_measure.p;
                    H0_Rep_ANOVA_Gp_effect(channel,:,1,b) = result.gp.F;
                    H0_Rep_ANOVA_Gp_effect(channel,:,2,b) = result.gp.p;
                    tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(channel,:,1,1,b) = result.interaction.F;
                    tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(channel,:,1,2,b) = result.interaction.p;
                    
                elseif type == 4
                    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB);
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C,XB);
                    end
                    tmp_boot_H0_Rep_ANOVA(channel,:,:,1,b) = result.repeated_measure.F';
                    tmp_boot_H0_Rep_ANOVA(channel,:,:,2,b) = result.repeated_measure.p';
                    H0_Rep_ANOVA_Gp_effect(channel,:,1,b) = result.gp.F;
                    H0_Rep_ANOVA_Gp_effect(channel,:,2,b) = result.gp.p;
                    tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(channel,:,:,1,b) = result.interaction.F';
                    tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(channel,:,:,2,b) = result.interaction.p';
                end
            end
        end
        
        % Store results
        chunk_results.H0_Rep_ANOVA = tmp_boot_H0_Rep_ANOVA;
        if type == 3 || type == 4
            chunk_results.H0_Rep_ANOVA_Gp_effect = H0_Rep_ANOVA_Gp_effect;
            chunk_results.H0_Rep_ANOVA_Interaction_with_gp = tmp_boot_H0_Rep_ANOVA_Interaction_with_gp;
        end
        chunk_results.type = type;
        chunk_results.var_names = options.var_names;
        
    otherwise
        error('Unknown analysis type: %d', analysis_type);
end

% Apply Time-Frequency reshape if needed
if is_tf
    fields = fieldnames(chunk_results);
    for f = 1:length(fields)
        if ~isempty(strfind(fields{f}, 'H0_')) && ~isempty(strfind(fields{f}, 'var_name'))
            chunk_results.(fields{f}) = limo_tf_5d_reshape(chunk_results.(fields{f}));
        end
    end
end

end
