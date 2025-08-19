function LIMOPath = limo_random_robust(varargin)

% This function makes the result files for the random effects of various tests
% as well as organizes and makes files for boostrap. It is interfaced with
% limo_random_effect which itself interfaces with the user to select and pass
% the data in the appropriate format. Limo_random_robust calls low level
% functions to perform the actual computation and add the trimmed mean or mean
% of parameters for the correponding test (helps for vizualizing effects)
%
% FORMAT LIMOPath = limo_random_robust(test,data,label,LIMO)
%
% INPUTS
%
% limo_random_robust(1,y,parameter number,LIMO)
%                    1 = a one-sample t-test
%                    y = data (dim channels, time or freq, subjects)
%                      = data (dim channels, freq, time, subjects)
%                      = the name of the Yr file
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%
% limo_random_robust(2,y1,y2,parameter number,LIMO)
%                    2 = two samples t-test
%                    y1 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y1r file
%                    y2 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y2r file
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%
% limo_random_robust(3,y1,y2,parameter number,LIMO)
%                    3 = paired t-test
%                    y1 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y1r file
%                    y2 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y2r file
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%
% limo_random_robust(4,y,X,parameter number,LIMO)
%                    4 = regression analysis
%                    y = data (dim channels, time or freq, subjects)
%                      = data (dim channels, freq, time, subjects)
%                    X = continuous regressor(s)
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%
% limo_random_robust(5,y,cat,cont,LIMO,'go',option)
%                    5 = N-way ANOVA/ANCOVA
%                    y = data (dim channels, time or freq, subjects)
%                      = data (dim channels, freq, time, subjects)
%                      = the name of the Yr file
%                    cat = categorical variable(s)
%                    cont = continuous regressors (covariates)
%                    LIMO the basic structure with data, design and channel info
%                    'go' is optional and prompt or not the design
%                         options are 'yes' (prompt, default),
%                         or 'no' (no prompt, usuful for scripting)
%
% limo_random_robust(6,y,gp,factor_levels,LIMO,'go',option)
%                    6 = Repeated measures ANOVA/ANCOVA using multivariate approach
%                    y = data (dim channels, time or freq, subjects, measures)
%                      = data (dim channels, freq, time, subjects, measures)
%                      = the name of the Yr file
%                    gp = a vector defining gps
%                    factor_levels = a vector specifying the levels of each repeated measure factor
%                    LIMO the basic structure with data, design and channel info
%                         or the full name of the LIMO file
%                    'go' is optional and prompt or not the pseudo-design
%                         options are 'yes' (prompt, default),
%                         or 'no' (no prompt, usuful for scripting)
%
% OUTPUT
% write on the disk matrices correponding to the test (Yr and LIMO.mat are generated in limo_random_select,
% and for Regression, ANOVA, the LIMO.mat structure is updated)
%
% 1 one_sample_parameter_X (channels, frames [time, freq or freq-time], [mean value, se, df, t, p])
%   H0_one_sample_ttest_parameter_X (channels, frames, [T values under H0, p values under H0], LIMO.design.bootstrap)
%
% 2 two_samples_parameter_X (channels, frames [time, freq or freq-time], [mean value, se, df, t, p])
%   H0_two_samples_ttest_parameter_X (channels, frames, [T values under H0, p values under H0], LIMO.design.bootstrap)
%
% 3 paired_samples_parameter_X (channels, frames [time, freq or freq-time], [mean value, se, df, t, p])
%   H0_paired_samples_ttest_parameter_X (channels, frames, [T values under H0, p values under H0], LIMO.design.bootstrap)
%
% 4 R2 (channels, frames [time, freq or freq-time], [F p values])
%   H0_R2 (channels, frames, [F p], LIMO.design.bootstrap)
%   Covariate_effect_X (channels, frames [time, freq or freq-time], [F p values])
%   H0_Covariate_effect_X (channels, frames, [F p], LIMO.design.bootstrap)
%
% 5 Condition_effect_X (channels, frames [time, freq or freq-time], [F p values])
%   H0_Condition_effect_X (channels, frames, [F p], LIMO.design.bootstrap)
%   Covariate_effect_X (channels, frames [time, freq or freq-time], [F p values])
%   H0_Covariate_effect_X (channels, frames, [F p], LIMO.design.bootstrap)
%
% 6 Rep_ANOVA_Factor_X (channels, frames [time, freq or freq-time], [F p values])
%   Rep_ANOVA_Gp_effect (channels, frames [time, freq or freq-time], [F p values])
%   Rep_ANOVA_Interaction_gp_Factor_X (channels, frames [time, freq or freq-time], [F p values])
%   H0_XXXXX same as above, including LIMO.design.bootstrap on the last dimension
%
% LIMOPath = LIMO.dir or [] if failed
%
% See also LIMO_TRIMCI LIMO_YUEN_TTEST LIMO_YUEND_TTEST LIMO_ROBUST_1WAY_ANOVA
% LIMO_GLM1 LIMO_EEG(4) LIMO_EEG_TF(4) LIMO_REP_ANOVA LIMO_CREATE_BOOT_TABLE
% ------------------------------
%  Copyright (C) LIMO Team 2020

warning off
%% inputs checks
LIMOPath = [];
if nargin == 0
    help limo_random_robust
    return
else
    type  = varargin{1};
    if type <1 || type > 6
        error('type argument must be between 1 and 6')
    end
end

%% Variables to persist for bootstrap
persistent_data = struct();
persistent_parameter = [];

%% start

switch type
    %--------------------------------------------------------------------------
    % One Sample t-test // bootstrap-t method
    %--------------------------------------------------------------------------
    case {1}
        
        if ischar(varargin{2})
            data      = load(varargin{2});
            data      = data.(cell2mat(fieldnames(data)));
        else
            data  = varargin{2};
        end
        parameter = varargin{3};
        if ischar(varargin{4})
            LIMO = load(varargin{4});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{4})
            LIMO = varargin{4};
        end
        clear varargin
        cd(LIMO.dir);
        
        if strcmp(LIMO.Analysis,'Time-Frequency') 
            data = limo_tf_4d_reshape(data);
        end
        
        % ------------------------------------------------
        % check the data structure
        for e=1:size(data,1)
            tmp = isnan(data(e,1,:));
            if length(tmp) == sum(isnan(tmp))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty - analysis aborded']);
                return
            elseif (length(tmp) - sum(isnan(tmp))) < 3
                errordlg([LIMO.Type ' ' num2str(e) ' has less than 3 subjects - analysis aborded']);
                return
            end
        end
        clear tmp
        
        if ~isfield(LIMO.design,'method')
            LIMO.design.method = 'Trimmed Mean';
        end
        
        % ------------------------------------------------
        % make a one_sample file per parameter (channels, frames, [mean value, se, df, t, p])
        one_sample = NaN(size(data,1), size(data,2), 5);
        name       = sprintf('one_sample_ttest_parameter_%g',parameter);
        
        for channel = 1:size(data,1) % run per channel because we have to remove NaNs
            fprintf('analyse parameter %g channel %g \n',parameter, channel);
            tmp = data(channel,:,:);
            if nansum(tmp(1,1,:)) == 0
                warning('empty channel - filling with NaNs')
                one_sample(channel,:,:) = NaN;
            else
                Y = tmp(1,:,find(~isnan(tmp(1,1,:))));
                if strcmpi(LIMO.design.method,'Trimmed Mean')
                    [one_sample(channel,:,4),one_sample(channel,:,1),~,one_sample(channel,:,2), ...
                        one_sample(channel,:,5),~,one_sample(channel,:,3)] = limo_trimci(Y);
                elseif strcmpi(LIMO.design.method,'Mean')
                    [one_sample(channel,:,1),one_sample(channel,:,3),~,sd,n, ...
                        one_sample(channel,:,4),one_sample(channel,:,5)] = limo_ttest(1,Y,0,5/100);
                    one_sample(channel,:,2) = sd./sqrt(n);
                else
                    error('unrecognized LIMO.design.method: %s',LIMO.design.method)
                end
                clear tmp Y
            end
        end
        
        if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
            one_sample = limo_tf_4d_reshape(one_sample);
        end
        save (name,'one_sample', '-v7.3')
        LIMOPath = LIMO.dir;
        
        % Store for bootstrap
        persistent_data.data = data;
        persistent_parameter = parameter;
        
    %--------------------------------------------------------------------------
    % Two Samples t-test // percentile bootstrap technique
    %--------------------------------------------------------------------------
    case {2}
        
        if ischar(varargin{2})
            data1  = load(varargin{2});
            data1  = data1.(cell2mat(fieldnames(data1)));
        else
            data1  = varargin{2};
        end
        if ischar(varargin{3})
            data2  = load(varargin{3});
            data2  = data2.(cell2mat(fieldnames(data2)));
        else
            data2  = varargin{3};
        end
        parameter = varargin{4};
        if ischar(varargin{5})
            LIMO = load(varargin{5});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{5})
            LIMO = varargin{5};
        end
        cd(LIMO.dir);
        clear varargin
       
        if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
            data1 = limo_tf_4d_reshape(data1);
            data2 = limo_tf_4d_reshape(data2);
        end
        
        % ------------------------------------------------
        % check the data structure
        if size(data1,1) ~= size(data2,1)
            error(['groups have a different number of ' LIMO.Type])
        end
        
        for e=1:size(data1,1)
            tmp = isnan(data1(e,1,:));
            tmp2 = isnan(data2(e,1,:));
            if length(tmp) == sum(isnan(tmp))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty in group 1 - analysis aborded']);
                return
            elseif (length(tmp) - sum(isnan(tmp))) < 3
                errordlg([LIMO.Type ' ' num2str(e) ' has less than 3 subjects in group 1 - analysis aborded']);
                return
            elseif length(tmp2) == sum(isnan(tmp2))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty in group 2 - analysis aborded']);
                return
            elseif (length(tmp2) - sum(isnan(tmp2))) < 3
                errordlg([LIMO.Type ' ' num2str(e) ' has less than 3 subjects in group 2 - analysis aborded']);
                return
            end
        end
        clear tmp tmp2
        
        if ~isfield(LIMO.design,'method')
            LIMO.design.method = 'Trimmed Mean';
        end
        
        % ------------------------------------------------
        % make a two_samples file per parameter (channels, frames, [mean value, se, df, t, p])
        two_samples = NaN(size(data1,1), size(data1,2),5);
        name = sprintf('two_samples_ttest_parameter_%g_%g',parameter);
        
        array = intersect(find(~isnan(data1(:,1,1))),find(~isnan(data2(:,1,1))));
        for e = 1:size(array,1)
            channel = array(e);
            fprintf('analyse parameter %g channel %g',parameter, channel); disp(' ');
            tmp = data1(channel,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            tmp = data2(channel,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                [two_samples(channel,:,4),two_samples(channel,:,1),two_samples(channel,:,2),...
                    ~,two_samples(channel,:,5),~,two_samples(channel,:,3)]=limo_yuen_ttest(Y1,Y2);
            else % if strcmpi(LIMO.design.method,'Mean')
                [two_samples(channel,:,1),two_samples(channel,:,3),~,sd,~,two_samples(channel,:,4),...
                    two_samples(channel,:,5)]=limo_ttest(2,Y1,Y2,.05);
                sd = sd.^2; a = sd(1,:)./size(Y1,3); b = sd(1,:)./size(Y2,3);
                two_samples(channel,:,2) = sqrt(a + b);
            end
            clear Y1 Y2
        end
        
        if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
            two_samples = limo_tf_4d_reshape(two_samples);
        end
        save (name,'two_samples', '-v7.3')
        LIMOPath = LIMO.dir;
        
        % Store for bootstrap
        persistent_data.data1 = data1;
        persistent_data.data2 = data2;
        persistent_parameter = parameter;
        
    %--------------------------------------------------------------------------
    % Paired t-test // percentile bootstrap technique
    %--------------------------------------------------------------------------
    case {3}
        
        if ischar(varargin{2})
            data1  = load(varargin{2});
            data1  = data1.(cell2mat(fieldnames(data1)));
        else
            data1  = varargin{2};
        end
        if ischar(varargin{3})
            data2  = load(varargin{3});
            data2  = data2.(cell2mat(fieldnames(data2)));
        else
            data2  = varargin{3};
        end
        parameter = varargin{4};
        if ischar(varargin{5})
            LIMO = load(varargin{5});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{5})
            LIMO = varargin{5};
        end
        cd(LIMO.dir);
        clear varargin
        
        if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
            data1 = limo_tf_4d_reshape(data1);
            data2 = limo_tf_4d_reshape(data2);
        end
        
        % ------------------------------------------------
        % check the data structure
        if size(data1,1) ~= size(data2,1)
            error(['samples have a different number of ' LIMO.Type ',not a paired t-tests'])
        end
        
        for e=1:size(data1,1)
            tmp = isnan(data1(e,1,:));
            tmp2 = isnan(data2(e,1,:));
            if length(tmp) ~= length(isnan(tmp2))
                errordlg([LIMO.Type ' ' num2str(e) ' has unpaired data - analysis aborded, not a paired t-test']);
                return
            elseif length(tmp) == sum(isnan(tmp))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty - analysis aborded']);
                return
            elseif (length(tmp) - sum(isnan(tmp))) < 3
                errordlg([LIMO.Type ' ' num2str(e) ' has less than 3 subjects - analysis aborded']);
                return
            end
        end
        clear tmp tmp2
        
        if ~isfield(LIMO.design,'method')
            LIMO.design.method = 'Trimmed Mean';
        end
        
        % ------------------------------------------------
        % make a paired_samples file per parameter (channels, frames, [mean value, se, df, t, p])
        paired_samples = NaN(size(data1,1), size(data1,2),5);
        name = sprintf('paired_samples_ttest_parameter_%s',num2str(parameter')');
        
        array = intersect(find(~isnan(data1(:,1,1))),find(~isnan(data2(:,1,1))));
        for e = 1:size(array,1)
            channel = array(e);
            fprintf('analyse parameter %s channel %g',num2str(parameter')', channel); disp(' ');
            tmp = data1(channel,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            tmp = data2(channel,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                [paired_samples(channel,:,4),paired_samples(channel,:,1),paired_samples(channel,:,2),...
                    ~,paired_samples(channel,:,5),~,paired_samples(channel,:,3)]=limo_yuend_ttest(Y1,Y2); 
            else % strcmpi(LIMO.design.method,'Mean')
                [paired_samples(channel,:,1),paired_samples(channel,:,3),~,sd,n,paired_samples(channel,:,4),...
                    paired_samples(channel,:,5)]=limo_ttest(1,Y1,Y2,.05);
                paired_samples(channel,:,2) = sd./sqrt(n);
            end
            clear Y1 Y2
        end
        
        if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
            paired_samples = limo_tf_4d_reshape(paired_samples);
        end
        save (name,'paired_samples', '-v7.3')
        LIMOPath = LIMO.dir; 
        
        % Store for bootstrap
        persistent_data.data1 = data1;
        persistent_data.data2 = data2;
        persistent_parameter = parameter;
        
    %------------------------------------------------------------------
    % Regression // percentile bootstrap under H0
    %------------------------------------------------------------------
    case {4}
        
        if ischar(varargin{2})
            data = load(varargin{2});
            data = data.cell2mat(fieldnames(data));
        else
            data = varargin{2};
        end
        regressors = varargin{3}; % the predictors across subjects like e.g. age
        parameter  = varargin{4}; % the parameters from 1st level matrices the regression is computed on (just for name)
        if ischar(varargin{5})
            LIMO = load(varargin{5});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{5})
            LIMO = varargin{5};
        end
        cd(LIMO.dir);
        
        if nargin >5
            for in = 6:2:nargin
                if strcmpi(varargin{in},'zscore')
                    answer = varargin{in+1};
                elseif strcmpi(varargin{in},'go')
                    go = varargin{in+1};
                end
            end
        else
            go = 'No';
        end
        clear varargin
        cd(LIMO.dir);
        
        % ------------------------------------------------
        % check the data structure
        for e=1:size(data,1)
            if strcmp(LIMO.Analysis,'Time-Frequency')
                tmp = isnan(data(e,1,1,:));
            else
                tmp = isnan(data(e,1,:));
            end
            
            if length(tmp) == sum(isnan(tmp))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty - analysis aborded']);
                return
            elseif (length(tmp) - sum(isnan(tmp))) < 3
                errordlg([LIMO.Type ' ' num2str(e) ' has less than 3 subjects - analysis aborded']);
                return
            elseif (length(tmp) - sum(isnan(tmp))) < 6
                warndlg([LIMO.Type ' ' num2str(e) ' has less than 6 subjects - regression results will likely be biased']);
            end
        end
        
        % ------------------------------------------------
        % update the LIMO structure
        LIMO.data.Cat                = 0;
        LIMO.data.Cont               = regressors;
        LIMO.data.data_dir           = pwd;
        LIMO.design.type_of_analysis = 'Mass-univariate';
        LIMO.design.method           = 'IRLS'; 
        LIMO.design.fullfactorial    = 0;
        LIMO.design.status           = 'to do';
        
        if ~exist('answer','var') 
            answer = questdlg('zscore regressor(s)?','Regression option','Yes','No','Yes');
        end
        
        if isempty(answer)
            return
        elseif strcmpi(answer,'Yes')
            LIMO.design.zscore = 1;
        else
            LIMO.design.zscore = 0;
        end
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');

        % make design matrix and files
        if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
            [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
                LIMO.design.nb_continuous] = limo_design_matrix_tf(data, LIMO,1);
        else
            [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
                LIMO.design.nb_continuous] = limo_design_matrix(data, LIMO,1);
        end
        
        % ------------------------------------------------
        % do the analysis
        if strcmpi(go,'no')
            go = questdlg('run the analysis?','Start GLM analysis','Yes','No','Yes');
        end
        close('LIMO design');        

        if strcmpi(go,'Yes')
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
            if nargout ~= 0, LIMOPath = fullfile(pwd,'LIMO.mat'); end
            clear data regressors files; limo_eeg(4); disp('regression analysis done');
        else
            return
        end
                
    %--------------------------------------------------------------------------
    % N-ways ANOVA / ANCOVA
    %--------------------------------------------------------------------------
    case {5}
        
        data  = varargin{2};
        if ischar(data)
            data = load(data);
            data = data.(cell2mat(fieldnames(data)));
        end
        cat   = varargin{3}; 
        cont  = varargin{4}; 
        if ischar(varargin{5})
            LIMO = load(varargin{5});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{5})
            LIMO = varargin{5};
        end
        cd(LIMO.dir);
        
        if nargin ==7
            if strcmpi(varargin{6},'go')
                design_check = varargin{7};
            end
        end
        clear varargin
        cd(LIMO.dir);        
        
        % ------------------------------------------------
        % check the data structure
        for e=1:size(data,1)
            if strcmp(LIMO.Analysis,'Time-Frequency')
                tmp = isnan(data(e,1,1,:));
            else
                tmp = isnan(data(e,1,:));
            end
            
            if length(tmp) == sum(isnan(tmp))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty - analysis aborded']);
                return
            end
        end
        
        % ------------------------------------------------
        % update the LIMO structure
        LIMO.design.type_of_analysis  = 'Mass-univariate';
        LIMO.data.Cat                 = cat;
        LIMO.data.Cont                = cont;
        LIMO.data.data_dir            = pwd;
        LIMO.design.zscore            = 1;
        if size(cat,2) > 1
            LIMO.design.fullfactorial = 1;
        else
            LIMO.design.fullfactorial = 0;
        end
        
        % make design matrix and files
        if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
            [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
                LIMO.design.nb_continuous] = limo_design_matrix_tf(data, LIMO,1);
        else
            [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
                LIMO.design.nb_continuous] = limo_design_matrix(data, LIMO,1);
        end
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');

        % ------------------------------------------------
        % do the analysis
        if ~exist('design_check','var')
            go = questdlg('run the analysis?','Start GLM analysis','Yes','No','Yes');
        else
            if strcmpi(design_check,'yes')
                go = 'yes';
            else
                go = questdlg('run the analysis?','Start GLM analysis','Yes','No','Yes');
            end
        end
        close('LIMO design');
        
        if strcmpi(go,'Yes')
            if LIMO.design.fullfactorial == 0 && LIMO.design.nb_continuous == 0
                if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                    data = limo_tf_4d_reshape(data);
                end
                Yhat             = NaN(size(data));
                Condition_effect = NaN(size(data,1),size(data,2),2);
                LIMO.design.df    = NaN(size(data,1),1);
                LIMO.design.dfe   = NaN(size(data,1),size(data,2));
                array            = find(~isnan(data(:,1,1)));
                for e=1:size(array,1)
                    channel = array(e); fprintf('processing channel %g \n',channel);
                    [Condition_effect(channel,:,1), Condition_effect(channel,:,2),...
                        Yhat(channel,:,:), LIMO.design.df(channel),LIMO.design.dfe(channel,:)] = ...
                        limo_robust_1way_anova(squeeze(data(channel,:,:)),LIMO.design.X(:,1:end-1),20); % no intercept in this model
                end
                delete(fullfile(LIMO.dir,'Betas.mat')); % no betas here
                delete(fullfile(LIMO.dir,'R2.mat'));    % no R2
                
                if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                    Condition_effect = limo_tf_4d_reshape(Condition_effect);
                    Yhat = limo_tf_4d_reshape(Yhat); % these are the trimmed mean
                    data = limo_tf_4d_reshape(data);
                end
                save('Condition_effect_1.mat','Condition_effect', '-v7.3');
                clear  Condition_effect_1
                save('Yhat.mat','Yhat', '-v7.3');                                          
                Res = data - Yhat;
                clear Yhat
                save('Res.mat','Res', '-v7.3');                                            
                clear Res
                LIMO.design.status = 'done';
                LIMO.design.method = 'Generalized Welch''s method';
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
                
                % Store for bootstrap
                persistent_data.data = data;
            else
                LIMO.design.method = 'IRLS'; % will switch to OLS if N<50
                LIMO.design.status = 'to do';
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
                LIMOPath = LIMO.dir; 
                clear data LIMO
                limo_eeg(4); return
            end
        else
            disp('Analysis aborted')
            return
        end
        
    %----------------------------------------------------------------------------------------------
    % Repeated Measure ANOVA (multivariate approach) - bootstrap centering data
    %----------------------------------------------------------------------------------------------
    case {6}
        
        if ischar(varargin{2})
            data = load(varargin{2});
            data = data.(cell2mat((fieldnames(data))));
        else
            data = varargin{2}; % e,f,subjects,measures
        end
        
        factor_levels     = varargin{4}; % vector eg [2 3] --> LIMO.design.repeated_measure
        gp_vector         = varargin{3}; % length of data, indices groups --> LIMO.data.Cat  
                                         % size(LIMO.design.X,2) = % N*prod(factor_levels)+1
        if ischar(varargin{5})
            LIMO = load(varargin{5});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{5})
            LIMO = varargin{5};
        end
        
        if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
            tmp = NaN(size(data,1), size(data,2)*size(data,3),size(data,4),size(data,5));
            for measure = 1:size(data,5)
                if size(data,1) == 1
                    tmp(:,:,:,measure) = limo_tf_4d_reshape(data(1,:,:,:,measure));
                else
                    tmp(:,:,:,measure) = limo_tf_4d_reshape(squeeze(data(:,:,:,:,measure)));
                end
            end
            clear data; data=tmp; clear tmp;
            if size(data,3) ~= length(gp_vector)
                error('gp vector is not commensurate to the data (dimension 4)')
            end
            if size(data,4) ~= prod(factor_levels)
                error('factor_levels are not commensurate to the data (dimension 5)')
            end
        else
            if size(data,3) ~= length(gp_vector)
                error('gp vector is not commensurate to the data (dimension 3)')
            end
            if size(data,4) ~= prod(factor_levels)
                error('factor_levels are not commensurate to the data (dimension 4)')
            end
        end
        
        if nargin > 5
            if strcmpi(varargin{6},'go')
                go = varargin{7};
            end
        else
            go = 'no';
        end
        clear varargin
        
        % ------------------------------------------------
        % update the LIMO structure
        if isfield(LIMO,'dir')
            cd(LIMO.dir)
        end
        LIMO.data.Cat                = gp_vector;
        LIMO.data.Cont               = 0;
        LIMO.data.data_dir           = pwd;
        LIMO.design.type_of_analysis = 'Mass-univariate';
        LIMO.design.nb_conditions    = length(unique(gp_vector));
        LIMO.design.nb_interactions  = 0;
        LIMO.design.nb_continuous    = 0;
        LIMO.design.fullfactorial    = 0;
        LIMO.design.zscore           = 0;
        LIMO.design.repeated_measure = factor_levels;
        
        % specific stuff for repeated measures
        % from the input we know which case to handle
        if unique(gp_vector) == 1
            % one sample
            if length(factor_levels) ==1
                rep_type = 1; % simple repeated measure ANOVA (1 effect)
            elseif length(factor_levels) >1
                rep_type = 2; % multiple factors (effects and interactions)
            end
        else
            % k samples
            if length(factor_levels) ==1
                rep_type = 3; % gp * repeated measure (2 effects + interactions)
            elseif length(factor_levels) >1
                rep_type = 4; % gp * multiple factors (effects and interactions)
            end
        end
        
        % ------------------------------------------------
        % Complete case 6 analysis
        % ------------------------------------------------
        
        % make files to be stored
        % -----------------------
        if rep_type == 1 % one factor
            LIMO.design.method     = 'Mean'; % change to Trimmed Mean for robust ANOVA
            C                      = [eye(size(data,4)-1) ones(size(data,4)-1,1).*-1]; % contrast
            tmp_Rep_ANOVA          = NaN(size(data,1),size(data,2),1,2); % store F and p
            LIMO.design.effects{1} = 'Main effect';
            LIMO.design.C{1}       = C;
            x                      = kron(eye(prod(factor_levels)),ones(size(data,3),1));
            LIMO.design.X          = [x ones(size(x,1),1)];
            
        elseif rep_type == 2 % many factors
            LIMO.design.method = 'Mean';
            C                  = limo_OrthogContrasts(factor_levels);
            tmp_Rep_ANOVA      = NaN(size(data,1),size(data,2),length(C),2); % store F and p for each within factor and interactions
            LIMO.design.C      = C;
            index = length(factor_levels)+1;
            
            for i= 1:length(factor_levels)
                LIMO.design.effects{i} = ['Main effect ' num2str(i)];
            end
            
            for i= 2:length(factor_levels)
                n = nchoosek([1:length(factor_levels)],i);
                for j=1:size(n,1)
                    LIMO.design.effects{index} = ['Interaction ' num2str(n(j,:))]; index = index+1;
                end
            end
            x             = kron(eye(prod(factor_levels)),ones(size(data,3),1));
            LIMO.design.X = [x ones(size(x,1),1)];
            
        elseif rep_type == 3 % one factor within and one factor between
            LIMO.design.method                = 'Mean';
            gp_values                         = unique(gp_vector);
            k                                 = length(gp_values);
            X                                 = NaN(size(gp_vector,1),k+1);
            for g =1:k
                X(:,g) = gp_vector == gp_values(g);
            end
            X(:,end) = 1; % design matrix for gp effects
            C                                 = [eye(size(data,4)-1) ones(size(data,4)-1,1).*-1]; % contrast
            tmp_Rep_ANOVA                     = NaN(size(data,1),size(data,2),1,2);
            Rep_ANOVA_Gp_effect               = NaN(size(data,1),size(data,2),2);
            tmp_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),1,2);
            LIMO.design.C{1}                  = C;
            LIMO.design.effects{1}            = 'Main effect';
            x                                 = kron(X(:,1:k),eye(prod(factor_levels)));
            LIMO.design.X                     = [x sum(x,2)]; % just for display
            LIMO.design.nb_interactions       = length(LIMO.design.C);
           
        elseif rep_type == 4 % many factors within and one factor between
            LIMO.design.method                = 'Mean';
            gp_values                         = unique(gp_vector);
            k                                 = length(gp_values);
            X                                 = NaN(size(gp_vector,1),k+1);
            for g =1:k
                X(:,g) = gp_vector == gp_values(g);
            end
            X(:,end)                          = 1; % design matrix for gp effects
            C                                 = limo_OrthogContrasts(factor_levels);
            tmp_Rep_ANOVA                     = NaN(size(data,1),size(data,2),length(C),2);
            Rep_ANOVA_Gp_effect               = NaN(size(data,1),size(data,2),2);
            tmp_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),length(C),2);
            LIMO.design.C                     = C;
            for i= 1:length(factor_levels)
                LIMO.design.effects{i} = ['Main effect ' num2str(i)];
            end
            index                             = length(factor_levels)+1;
            for i= 2:length(factor_levels)
                n = nchoosek([1:length(factor_levels)],i);
                for j=1:size(n,1)
                    LIMO.design.effects{index} = ['Interaction ' num2str(n(j,:))];
                    index = index+1;
                end
            end
            x                                  = kron(X(:,1:k),eye(prod(factor_levels)));
            LIMO.design.X                      = [x sum(x,2)]; % just for display
            LIMO.design.nb_interactions        = length(LIMO.design.C);
       end
        
        if isempty(dir('Rep_ANOVA_Factor*.mat'))
            
            % check the design with user
            % --------------------------
            if ~strcmpi(go,'Yes')
                figure('Name','LIMO design'); set(gcf,'Color','w');
                imagesc(LIMO.design.X); colormap('gray');
                title('ANOVA model','FontSize',16);xlabel('regressors');
                ylabel('subjects'); drawnow;
                go = questdlg('start the analysis?');
                close('LIMO design')
                if ~strcmpi(go,'Yes')
                    return
                end
            end
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
            
            % do the analysis
            % ---------------
            array = find(~isnan(data(:,1,1,1)));
            for e = 1:length(array)
                channel = array(e);
                fprintf('analyse channel %g/%g\n ...', channel,size(data,1));
                tmp = squeeze(data(channel,:,:,:));
                if size(data,2) == 1
                    Y  = ones(1,size(tmp,1),size(tmp,2)); Y(1,:,:) = tmp;
                    gp = gp_vector(find(~isnan(Y(1,:,1))),:);
                    Y  = Y(:,find(~isnan(Y(1,:,1))),:);
                else
                    Y  = tmp(:,find(~isnan(tmp(1,:,1))),:);
                    gp = gp_vector(find(~isnan(tmp(1,:,1))),:);
                end
                
                if rep_type == 3 || rep_type == 4
                    XB = X(find(~isnan(tmp(1,:,1))),:);
                end
                
                if rep_type == 1
                    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C); % trimmed means
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C); % usual means
                    end
                    tmp_Rep_ANOVA(channel,:,1,1) = result.F;
                    tmp_Rep_ANOVA(channel,:,1,2) = result.p;
                    LIMO.design.df(channel)      = result.df;
                    LIMO.design.dfe(channel)     = result.dfe;
                elseif rep_type == 2
                    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C); % trimmed means
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C); % usual means
                    end
                    tmp_Rep_ANOVA(channel,:,:,1) = result.F';
                    tmp_Rep_ANOVA(channel,:,:,2) = result.p';
                    LIMO.design.df(channel,:)    = result.df;
                    LIMO.design.dfe(channel,:)   = result.dfe;
                elseif rep_type == 3
                    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB); % trimmed means
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C,XB); % usual means
                    end
                    tmp_Rep_ANOVA(channel,:,1,1)                   = result.repeated_measure.F;
                    tmp_Rep_ANOVA(channel,:,1,2)                   = result.repeated_measure.p;
                    LIMO.design.df(channel)                        = result.repeated_measure.df;
                    LIMO.design.dfe(channel)                       = result.repeated_measure.dfe;
                    Rep_ANOVA_Gp_effect(channel,:,1)               = result.gp.F;
                    Rep_ANOVA_Gp_effect(channel,:,2)               = result.gp.p;
                    LIMO.design.group.df(channel)                  = result.gp.df;
                    LIMO.design.group.dfe(channel)                 = result.gp.dfe;
                    tmp_Rep_ANOVA_Interaction_with_gp(channel,:,1,1) = result.interaction.F;
                    tmp_Rep_ANOVA_Interaction_with_gp(channel,:,1,2) = result.interaction.p;
                    LIMO.design.interaction.df(channel)            = result.interaction.df;
                    LIMO.design.interaction.dfe(channel)           = result.interaction.dfe;
                elseif rep_type == 4
                    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB); % trimmed means
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C,XB); % usual means
                    end
                    tmp_Rep_ANOVA(channel,:,:,1)                     = result.repeated_measure.F';
                    tmp_Rep_ANOVA(channel,:,:,2)                     = result.repeated_measure.p';
                    LIMO.design.df(channel,:)                        = result.repeated_measure.df;
                    LIMO.design.dfe(channel,:)                       = result.repeated_measure.dfe;
                    Rep_ANOVA_Gp_effect(channel,:,1)                 = result.gp.F;
                    Rep_ANOVA_Gp_effect(channel,:,2)                 = result.gp.p;
                    LIMO.design.group.df(channel)                    = result.gp.df;
                    LIMO.design.group.dfe(channel)                   = result.gp.dfe;
                    tmp_Rep_ANOVA_Interaction_with_gp(channel,:,:,1) = result.interaction.F';
                    tmp_Rep_ANOVA_Interaction_with_gp(channel,:,:,2) = result.interaction.p';
                    LIMO.design.interaction.df(channel,:)            = result.interaction.df;
                    LIMO.design.interaction.dfe(channel,:)           = result.interaction.dfe;
                end
                
                nb_effects = size(tmp_Rep_ANOVA,3);
                clear tmp Y gp result
            end
            
            % save stuff
            % ---------
            Rep_filenames = cell(1,nb_effects);
            for i=1:nb_effects
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
                
                % save each factor effect as F/p values
                % use reshape instead of squeeze in case there is only 1 channel
                Rep_ANOVA = reshape(tmp_Rep_ANOVA(:,:,i,:),...
                    [size(tmp_Rep_ANOVA,1) size(tmp_Rep_ANOVA,2) size(tmp_Rep_ANOVA,4)]);
                if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                    Rep_ANOVA = limo_tf_4d_reshape(Rep_ANOVA);
                end
                save(Rep_filenames{i},'Rep_ANOVA', '-v7.3');
                if nargout ~= 0, LIMOPath{i} = [fullfile(pwd,Rep_filenames{i}),'.mat']; end
            end
            
            if rep_type == 3 || rep_type ==4
                IRep_filenames = cell(1,nb_effects);
                for i=1:size(tmp_Rep_ANOVA_Interaction_with_gp,3)
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
                    
                    % save each interaction effect as F/p values
                    Rep_ANOVA_Interaction_with_gp = reshape(tmp_Rep_ANOVA_Interaction_with_gp(:,:,i,:),...
                        [size(tmp_Rep_ANOVA_Interaction_with_gp,1) size(tmp_Rep_ANOVA_Interaction_with_gp,2) size(tmp_Rep_ANOVA_Interaction_with_gp,4)]);
                    if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                        Rep_ANOVA_Interaction_with_gp = limo_tf_4d_reshape(Rep_ANOVA_Interaction_with_gp);
                    end
                    save(IRep_filenames{i},'Rep_ANOVA_Interaction_with_gp', '-v7.3');
                    clear Rep_ANOVA_Interaction_with_gp;
                    if nargout ~= 0, LIMOPath = [fullfile(pwd,IRep_filenames{i}),'.mat']; end
                end
                
                % Main group effet
                save('Rep_ANOVA_Gp_effect.mat','Rep_ANOVA_Gp_effect','-v7.3');
                if nargout ~= 0, LIMOPath = fullfile(pwd,'Rep_ANOVA_Gp_effect.mat'); end
            end
        end
        
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
        
        % if skipping the above
        if ~exist('nb_effects','var')
            nb_effects = length(LIMO.design.C);
        end
        
        % clear up all tmp files
        clear tmp_Rep_ANOVA
        if rep_type == 3 || rep_type == 4
            clear Rep_ANOVA_Gp_effect tmp_Rep_ANOVA_Interaction_with_gp
        end
        
        % Store for bootstrap - we'll need the original data
        persistent_data.data = data;
        persistent_data.gp_vector = gp_vector;
        persistent_data.factor_levels = factor_levels;
        persistent_data.rep_type = rep_type;
        if exist('Rep_filenames', 'var')
            persistent_data.Rep_filenames = Rep_filenames;
        end
        if exist('IRep_filenames', 'var')
            persistent_data.IRep_filenames = IRep_filenames;
        end
        if exist('X', 'var')
            persistent_data.X = X;
        end
        if exist('C', 'var')
            persistent_data.C = C;
        end
        
        LIMOPath = LIMO.dir;
        disp('Repeated Measures ANOVA done')
        
end % switch type

%% ========== CENTRALIZED CHUNKED BOOTSTRAP PROCESSING ==========
% This section handles bootstrapping for all analysis types (1-6)

if exist('LIMO', 'var') && isfield(LIMO, 'design') && isfield(LIMO.design, 'bootstrap') && LIMO.design.bootstrap > 0
    fprintf('\n=== CENTRALIZED BOOTSTRAP PROCESSING ===\n');
    
    % Determine if we should run bootstrap based on the analysis type
    run_bootstrap = true;
    
    % Check if bootstrap files already exist
    if usejava('desktop')
        boot_files = dir(fullfile(LIMO.dir, 'H0', 'H0_*.mat'));
        if ~isempty(boot_files)
            answer = questdlg('Bootstrap files already exist - overwrite?', 'Data check', 'Yes', 'No', 'Yes');
            if ~strcmp(answer, 'Yes')
                run_bootstrap = false;
            end
        end
    end
    
    if run_bootstrap && type >= 1 && type <= 6
        % Check for parallel pool
        limo_check_ppool;
        
        % Create H0 directory
        if ~exist(fullfile(LIMO.dir, 'H0'), 'dir')
            mkdir(fullfile(LIMO.dir, 'H0'));
        end
        
        % Chunking parameters
        chunk_size = 100;  % Adjust based on memory constraints
        n_chunks = ceil(LIMO.design.bootstrap / chunk_size);
        chunk_dir = fullfile(LIMO.dir, 'H0', 'chunks');
        
        fprintf('Total bootstraps: %d\n', LIMO.design.bootstrap);
        fprintf('Chunk size: %d\n', chunk_size);
        fprintf('Number of chunks: %d\n', n_chunks);
        
        % Create chunk directory
        if ~exist(chunk_dir, 'dir')
            mkdir(chunk_dir);
        end
        
        % Prepare centered data and boot tables based on analysis type
        options = struct();
        parameter = persistent_parameter;
        
        if type == 1  % One-sample t-test
            data = persistent_data.data;
            % Center data under H0
            if strcmpi(LIMO.design.method, 'Trimmed Mean')
                centered_data = data - repmat(limo_trimmed_mean(data), [1 1 size(data,3)]);
            else
                centered_data = data - repmat(nanmean(data,3), [1 1 size(data,3)]);
            end
            boot_table = limo_create_boot_table(data, LIMO.design.bootstrap);
            save(fullfile(LIMO.dir, 'H0', 'boot_table'), 'boot_table');
            options.parameter = parameter;
            
        elseif type == 2  % Two-samples t-test
            data1 = persistent_data.data1;
            data2 = persistent_data.data2;
            % Center both datasets
            if contains(LIMO.design.method, 'Trimmed Mean', 'IgnoreCase', true)
                data1_centered = data1 - repmat(limo_trimmed_mean(data1), [1 1 size(data1,3)]);
                data2_centered = data2 - repmat(limo_trimmed_mean(data2), [1 1 size(data2,3)]);
            else
                data1_centered = data1 - repmat(nanmean(data1,3), [1 1 size(data1,3)]);
                data2_centered = data2 - repmat(nanmean(data2,3), [1 1 size(data2,3)]);
            end
            centered_data = {data1_centered, data2_centered};
            boot_table1 = limo_create_boot_table(data1, LIMO.design.bootstrap);
            boot_table2 = limo_create_boot_table(data2, LIMO.design.bootstrap);
            save(fullfile(LIMO.dir, 'H0', 'boot_table1'), 'boot_table1');
            save(fullfile(LIMO.dir, 'H0', 'boot_table2'), 'boot_table2');
            boot_table = boot_table1;  % Primary boot table
            options.boot_table1 = boot_table1;
            options.boot_table2 = boot_table2;
            options.parameter = parameter;
            
        elseif type == 3  % Paired t-test
            data1 = persistent_data.data1;
            data2 = persistent_data.data2;
            % Center both datasets
            if contains(LIMO.design.method, 'Trimmed Mean', 'IgnoreCase', true)
                data1_centered = data1 - repmat(limo_trimmed_mean(data1), [1 1 size(data1,3)]);
                data2_centered = data2 - repmat(limo_trimmed_mean(data2), [1 1 size(data2,3)]);
            else
                data1_centered = data1 - repmat(nanmean(data1,3), [1 1 size(data1,3)]);
                data2_centered = data2 - repmat(nanmean(data2,3), [1 1 size(data2,3)]);
            end
            centered_data = {data1_centered, data2_centered};
            boot_table = limo_create_boot_table(data1, LIMO.design.bootstrap);
            save(fullfile(LIMO.dir, 'H0', 'boot_table'), 'boot_table');
            options.parameter = parameter;
            
        elseif type == 4  % Regression
            fprintf('Note: Regression bootstrapping handled by limo_eeg(4)\n');
            run_bootstrap = false;
            
        elseif type == 5  % N-way ANOVA/ANCOVA
            if LIMO.design.fullfactorial == 0 && LIMO.design.nb_continuous == 0
                data = persistent_data.data;
                % Robust 1-way ANOVA
                if strcmp(LIMO.Analysis, 'Time-Frequency') || strcmp(LIMO.Analysis, 'ITC')
                    data = limo_tf_4d_reshape(data);
                end
                % Center data for each condition
                for c = 1:(size(LIMO.design.X,2)-1)
                    index = find(LIMO.design.X(:,c));
                    data(:,:,index) = data(:,:,index) - repmat(limo_trimmed_mean(data(:,:,index)), [1 1 length(index)]);
                end
                centered_data = data;
                boot_table = limo_create_boot_table(data, LIMO.design.bootstrap);
                save(fullfile(LIMO.dir, 'H0', 'boot_table'), 'boot_table');
            else
                fprintf('Note: Standard ANOVA/ANCOVA bootstrapping handled by limo_eeg(4)\n');
                run_bootstrap = false;
            end
            
        elseif type == 6  % Repeated Measures ANOVA
            data = persistent_data.data;
            gp_vector = persistent_data.gp_vector;
            factor_levels = persistent_data.factor_levels;
            rep_type = persistent_data.rep_type;
            
            % Center data for each cell based on rep_type
            centered_data = NaN(size(data,1),size(data,2),size(data,3),size(data,4));
            nb_conditions = prod(factor_levels);
            
            if rep_type == 1 || rep_type == 2
                % Within-subjects only: center each condition across subjects
                for condition=1:nb_conditions
                    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                        avg = repmat(limo_trimmed_mean(data(:,:,:,condition),3),[1 1 size(data,3)]);
                    else
                        avg = repmat(nanmean(data(:,:,:,condition),3),[1 1 size(data,3)]);
                    end
                    centered_data(:,:,:,condition) = data(:,:,:,condition) - avg;
                end
            else
                % Mixed design: center within each group separately
                for gp=1:LIMO.design.nb_conditions
                    gp_index = find(LIMO.data.Cat == gp);
                    for condition=1:nb_conditions
                        if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                            avg = repmat(limo_trimmed_mean(data(:,:,gp_index,condition),3),[1 1 length(gp_index)]);
                        else
                            avg = repmat(nanmean(data(:,:,gp_index,condition),3),[1 1 length(gp_index)]);
                        end
                        centered_data(:,:,gp_index,condition) = data(:,:,gp_index,condition) - avg;
                    end
                end
            end
            
            % Create boot table (same for all conditions)
            boot_table = limo_create_boot_table(data(:,:,:,1), LIMO.design.bootstrap);
            save(fullfile(LIMO.dir, 'H0', 'boot_table'), 'boot_table');
            save(fullfile(LIMO.dir, 'H0', 'centered_data'), 'centered_data', '-v7.3');
            
            % Set additional options for repeated measures
            options.rep_type = rep_type;
            options.factor_levels = factor_levels;
            options.gp_vector = gp_vector;
            if isfield(persistent_data, 'C')
                options.C = persistent_data.C;
            end
            if isfield(persistent_data, 'X')
                options.X = persistent_data.X;
            end
            if isfield(persistent_data, 'Rep_filenames')
                options.Rep_filenames = persistent_data.Rep_filenames;
            end
            if isfield(persistent_data, 'IRep_filenames')
                options.IRep_filenames = persistent_data.IRep_filenames;
            end
            
        end  % switch type
        
        % Process chunks if bootstrap should run (excluding type 4 regression and type 6 repeated measures)
        if run_bootstrap && type ~= 4 && type ~= 6
            fprintf('\nProcessing bootstrap in %d chunks...\n', n_chunks);
            
            for chunk = 1:n_chunks
                chunk_start = (chunk - 1) * chunk_size + 1;
                chunk_end = min(chunk * chunk_size, LIMO.design.bootstrap);
                
                fprintf('\n--- Processing chunk %d/%d (bootstraps %d-%d) ---\n', ...
                    chunk, n_chunks, chunk_start, chunk_end);
                
                % Process this chunk
                chunk_results = limo_process_bootstrap_chunk(type, centered_data, ...
                    boot_table, chunk_start, chunk_end, LIMO, options);
                
                % Save chunk results
                if ~isempty(chunk_results.var_name)
                    fields = fieldnames(chunk_results);
                    for f = 1:length(fields)
                        field_name = fields{f};
                        if startsWith(field_name, 'H0_') && ~strcmp(field_name, 'var_name') && ~strcmp(field_name, 'var_names')
                            var_name = chunk_results.var_name;
                            base_name = var_name;
                            limo_save_boot_chunks(chunk_results.(field_name), chunk_dir, ...
                                base_name, chunk_start, chunk_size, field_name);
                        end
                    end
                end
            end
            
            % Merge chunks
            fprintf('\n=== MERGING BOOTSTRAP CHUNKS ===\n');
            
            % Determine which files to merge based on analysis type
            if type == 1
                var_name = sprintf('H0_one_sample_ttest_parameter_%g', parameter);
                output_file = fullfile(LIMO.dir, 'H0', [var_name '.mat']);
                limo_merge_boot_chunks(chunk_dir, var_name, output_file, ...
                    'var_name', 'H0_one_sample', 'delete_chunks', false);
                
            elseif type == 2
                var_name = sprintf('H0_two_samples_ttest_parameter_%g_%g', parameter);
                output_file = fullfile(LIMO.dir, 'H0', [var_name '.mat']);
                limo_merge_boot_chunks(chunk_dir, var_name, output_file, ...
                    'var_name', 'H0_two_samples', 'delete_chunks', false);
                
            elseif type == 3
                var_name = sprintf('H0_paired_samples_ttest_parameter_%s', num2str(parameter')');
                output_file = fullfile(LIMO.dir, 'H0', [var_name '.mat']);
                limo_merge_boot_chunks(chunk_dir, var_name, output_file, ...
                    'var_name', 'H0_paired_samples', 'delete_chunks', false);
                
            elseif type == 5
                if LIMO.design.fullfactorial == 0 && LIMO.design.nb_continuous == 0
                    output_file = fullfile(LIMO.dir, 'H0', 'H0_Condition_effect_1.mat');
                    limo_merge_boot_chunks(chunk_dir, 'H0_Condition_effect_1', output_file, ...
                        'var_name', 'H0_Condition_effect', 'delete_chunks', false);
                end
            end
            
            % Optional cleanup
            if usejava('desktop')
                answer = questdlg('Delete chunk files?', 'Cleanup', 'Yes', 'No', 'No');
                if strcmp(answer, 'Yes')
                    rmdir(chunk_dir, 's');
                    fprintf('Chunk directory deleted\n');
                end
            end
            
            fprintf('\n=== BOOTSTRAP PROCESSING COMPLETE ===\n');
        end  % if run_bootstrap
        
        % Special handling for type 6 (Repeated Measures ANOVA) bootstrap
        if run_bootstrap && type == 6
            fprintf('\n=== REPEATED MEASURES ANOVA BOOTSTRAP ===\n');
            
            % Use original bootstrap method for repeated measures
            % Load centered data
            centered_data = load(fullfile(LIMO.dir, 'H0', 'centered_data'));
            centered_data = centered_data.centered_data;
            boot_table = load(fullfile(LIMO.dir, 'H0', 'boot_table'));
            boot_table = boot_table.boot_table;
            
            % Get parameters from options
            rep_type = options.rep_type;
            factor_levels = options.factor_levels;
            gp_vector = options.gp_vector;
            C = options.C;
            
            % Reconstruct X for rep_type 3 and 4 if needed
            if rep_type == 3 || rep_type == 4
                if isfield(options, 'X') && ~isempty(options.X)
                    X = options.X;
                else
                    % Recreate X from gp_vector
                    gp_values = unique(gp_vector);
                    k = length(gp_values);
                    X = NaN(size(gp_vector,1),k+1);
                    for g = 1:k
                        X(:,g) = gp_vector == gp_values(g);
                    end
                    X(:,end) = 1; % intercept
                end
            else
                X = []; % Not needed for rep_type 1 and 2
            end
            
            % Prepare bootstrap arrays based on rep_type
            if rep_type == 1
                tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1),size(centered_data,2),1,2,LIMO.design.bootstrap);
            elseif rep_type == 2
                tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1),size(centered_data,2),length(C),2,LIMO.design.bootstrap);
            elseif rep_type == 3
                tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1),size(centered_data,2),1,2,LIMO.design.bootstrap);
                H0_Rep_ANOVA_Gp_effect = NaN(size(centered_data,1),size(centered_data,2),2,LIMO.design.bootstrap);
                tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(centered_data,1),size(centered_data,2),1,2,LIMO.design.bootstrap);
            else
                tmp_boot_H0_Rep_ANOVA = NaN(size(centered_data,1),size(centered_data,2),length(C),2,LIMO.design.bootstrap);
                H0_Rep_ANOVA_Gp_effect = NaN(size(centered_data,1),size(centered_data,2),2,LIMO.design.bootstrap);
                tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(centered_data,1),size(centered_data,2),length(C),2,LIMO.design.bootstrap);
            end
            
            % Bootstrap computation using chunked approach to avoid memory issues
            fprintf('Computing bootstrap for Repeated Measures ANOVA...\n');
            chunk_size = 50; % Smaller chunks for repeated measures due to 4D data
            n_boot_chunks = ceil(LIMO.design.bootstrap / chunk_size);
            
            % Make sure all variables are available for parfor
            local_rep_type = rep_type;
            local_factor_levels = factor_levels;
            local_gp_vector = gp_vector;
            local_C = C;
            local_X = X;
            local_method = LIMO.design.method;
            
            for boot_chunk = 1:n_boot_chunks
                boot_start = (boot_chunk - 1) * chunk_size + 1;
                boot_end = min(boot_chunk * chunk_size, LIMO.design.bootstrap);
                fprintf('Bootstrap chunk %d/%d (bootstraps %d-%d)\n', boot_chunk, n_boot_chunks, boot_start, boot_end);
                
                % Process this chunk of bootstraps
                parfor B = boot_start:boot_end
                    local_B = B - boot_start + 1;
                    array = find(~isnan(centered_data(:,1,1,1)));
                    
                    % Initialize local arrays for this bootstrap (always initialize all variables for parfor)
                    if local_rep_type == 1 || local_rep_type == 3
                        local_tmp_boot = NaN(size(centered_data,1),size(centered_data,2),1,2);
                    else % rep_type == 2 or 4
                        local_tmp_boot = NaN(size(centered_data,1),size(centered_data,2),length(local_C),2);
                    end
                    local_gp_boot = NaN(size(centered_data,1),size(centered_data,2),2); % Always initialize
                    if local_rep_type == 1 || local_rep_type == 3
                        local_int_boot = NaN(size(centered_data,1),size(centered_data,2),1,2); % Always initialize
                    else % rep_type == 2 or 4
                        local_int_boot = NaN(size(centered_data,1),size(centered_data,2),length(local_C),2); % Always initialize
                    end
                    
                    for e = 1:length(array)
                        channel = array(e);
                        % Get bootstrapped data for this channel
                        tmp = squeeze(centered_data(channel,:,boot_table{channel}(:,B),:));
                        if size(centered_data,2) == 1
                            Y = ones(1,size(tmp,1),size(tmp,2)); Y(1,:,:) = tmp;
                            gp = local_gp_vector(find(~isnan(Y(1,:,1))),:);
                            Y = Y(:,find(~isnan(Y(1,:,1))),:);
                        else
                            Y = tmp(:,find(~isnan(tmp(1,:,1))),:);
                            gp = local_gp_vector(find(~isnan(tmp(1,:,1))));
                        end
                        
                        if local_rep_type == 3 || local_rep_type == 4
                            if ~isempty(local_X)
                                XB = local_X(find(~isnan(tmp(1,:,1))),:);
                            else
                                XB = [];
                            end
                        else
                            XB = [];
                        end
                        
                        % Compute analysis for this bootstrap sample
                        if local_rep_type == 1
                            if contains(local_method,'Trimmed Mean','IgnoreCase',true)
                                result = limo_robust_rep_anova(Y,gp,local_factor_levels,local_C);
                            else
                                result = limo_rep_anova(Y,gp,local_factor_levels,local_C);
                            end
                            local_tmp_boot(channel,:,1,1) = result.F;
                            local_tmp_boot(channel,:,1,2) = result.p;
                        elseif local_rep_type == 2
                            if contains(local_method,'Trimmed Mean','IgnoreCase',true)
                                result = limo_robust_rep_anova(Y,gp,local_factor_levels,local_C);
                            else
                                result = limo_rep_anova(Y,gp,local_factor_levels,local_C);
                            end
                            local_tmp_boot(channel,:,:,1) = result.F';
                            local_tmp_boot(channel,:,:,2) = result.p';
                        elseif local_rep_type == 3
                            if contains(local_method,'Trimmed Mean','IgnoreCase',true)
                                result = limo_robust_rep_anova(Y,gp,local_factor_levels,local_C,XB);
                            else
                                result = limo_rep_anova(Y,gp,local_factor_levels,local_C,XB);
                            end
                            local_tmp_boot(channel,:,1,1) = result.repeated_measure.F;
                            local_tmp_boot(channel,:,1,2) = result.repeated_measure.p;
                            local_gp_boot(channel,:,1) = result.gp.F;
                            local_gp_boot(channel,:,2) = result.gp.p;
                            local_int_boot(channel,:,1,1) = result.interaction.F;
                            local_int_boot(channel,:,1,2) = result.interaction.p;
                        elseif local_rep_type == 4
                            if contains(local_method,'Trimmed Mean','IgnoreCase',true)
                                result = limo_robust_rep_anova(Y,gp,local_factor_levels,local_C,XB);
                            else
                                result = limo_rep_anova(Y,gp,local_factor_levels,local_C,XB);
                            end
                            local_tmp_boot(channel,:,:,1) = result.repeated_measure.F';
                            local_tmp_boot(channel,:,:,2) = result.repeated_measure.p';
                            local_gp_boot(channel,:,1) = result.gp.F;
                            local_gp_boot(channel,:,2) = result.gp.p;
                            local_int_boot(channel,:,:,1) = result.interaction.F';
                            local_int_boot(channel,:,:,2) = result.interaction.p';
                        end
                    end
                    
                    % Store results for this bootstrap
                    tmp_boot_H0_Rep_ANOVA(:,:,:,:,B) = local_tmp_boot;
                    if local_rep_type == 3 || local_rep_type == 4
                        H0_Rep_ANOVA_Gp_effect(:,:,:,B) = local_gp_boot;
                        tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(:,:,:,:,B) = local_int_boot;
                    end
                end
            end
            
            % Save bootstrap results using chunked saving to avoid memory issues
            fprintf('Saving repeated measures bootstrap results...\n');
            if isfield(options, 'Rep_filenames')
                Rep_filenames = options.Rep_filenames;
                for i=1:size(tmp_boot_H0_Rep_ANOVA,3)
                    name = sprintf('H0_%s',Rep_filenames{i});
                    H0_Rep_ANOVA = squeeze(tmp_boot_H0_Rep_ANOVA(:,:,i,:,:));
                    if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                        H0_Rep_ANOVA = limo_tf_5d_reshape(H0_Rep_ANOVA);
                    end
                    % Use chunked saving for large files
                    limo_save_rep_anova_bootstrap(H0_Rep_ANOVA, fullfile(LIMO.dir, 'H0', name));
                end
                
                if rep_type == 3 || rep_type == 4
                    % Save group effects
                    if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                        H0_Rep_ANOVA_Gp_effect = limo_tf_5d_reshape(H0_Rep_ANOVA_Gp_effect);
                    end
                    limo_save_rep_anova_bootstrap(H0_Rep_ANOVA_Gp_effect, fullfile(LIMO.dir, 'H0', 'H0_Rep_ANOVA_Gp_effect.mat'));
                    
                    % Save interactions
                    if isfield(options, 'IRep_filenames')
                        IRep_filenames = options.IRep_filenames;
                        for i=1:size(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp,3)
                            name = sprintf('H0_%s',IRep_filenames{i});
                            H0_Rep_ANOVA_Interaction_with_gp = squeeze(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(:,:,i,:,:));
                            if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                                H0_Rep_ANOVA_Interaction_with_gp = limo_tf_5d_reshape(H0_Rep_ANOVA_Interaction_with_gp);
                            end
                            limo_save_rep_anova_bootstrap(H0_Rep_ANOVA_Interaction_with_gp, fullfile(LIMO.dir, 'H0', name));
                        end
                    end
                end
            end
            
            fprintf('=== REPEATED MEASURES ANOVA BOOTSTRAP COMPLETE ===\n');
        end  % if type == 6 bootstrap
        
    end  % if run_bootstrap && type >= 1 && type <= 6
end  % if LIMO.design.bootstrap > 0

%% ========== END CENTRALIZED BOOTSTRAP ==========

% Final TFCE handling if needed
if exist('LIMO', 'var') && isfield(LIMO, 'design') && LIMO.design.tfce ~= 0
    if exist('name', 'var')
        limo_tfce_handling(fullfile(LIMO.dir, name));
        LIMO.design.tfce = 1;
    elseif type == 6 && exist('Rep_filenames', 'var')
        % Handle TFCE for repeated measures ANOVA files
        fprintf('Thresholding bootstrapped Rep ANOVA using TFCE \n');
        for i=1:length(Rep_filenames)
            limo_tfce_handling(fullfile(LIMO.dir, Rep_filenames{i}));
            if exist('IRep_filenames','var')
                if i == 1
                    limo_tfce_handling(fullfile(LIMO.dir, 'Rep_ANOVA_Gp_effect.mat'));
                end
                limo_tfce_handling(fullfile(LIMO.dir, IRep_filenames{i}));
            end
        end
        LIMO.design.tfce = 1;
    end
end

% Save final LIMO structure
if exist('LIMO', 'var')
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
end

warning on
end % main function

% =========================================================================
% HELPER FUNCTION FOR REPEATED MEASURES BOOTSTRAP SAVING
% =========================================================================

function limo_save_rep_anova_bootstrap(data, filepath)
    % Save repeated measures ANOVA bootstrap data with memory optimization
    % For very large files, use matfile for streaming writes
    
    dims = size(data);
    
    % For smaller files, use regular save
    if numel(data) < 1e8  % Less than ~100MB for double precision
        if contains(filepath, 'Rep_ANOVA')
            H0_Rep_ANOVA = data;
            save(filepath, 'H0_Rep_ANOVA', '-v7.3');
        elseif contains(filepath, 'Gp_effect')
            H0_Rep_ANOVA_Gp_effect = data;
            save(filepath, 'H0_Rep_ANOVA_Gp_effect', '-v7.3');
        else
            H0_Rep_ANOVA_Interaction_with_gp = data;
            save(filepath, 'H0_Rep_ANOVA_Interaction_with_gp', '-v7.3');
        end
        return;
    end
    
    % For larger files, use streaming approach
    fprintf('    Large bootstrap file detected - using streaming save\n');
    
    % Create matfile
    if exist(filepath, 'file')
        delete(filepath);
    end
    
    m = matfile(filepath, 'Writable', true);
    
    % Set up the variable based on filename
    if contains(filepath, 'Rep_ANOVA') && ~contains(filepath, 'Gp_effect') && ~contains(filepath, 'Interaction')
        m.H0_Rep_ANOVA = data;
    elseif contains(filepath, 'Gp_effect')
        m.H0_Rep_ANOVA_Gp_effect = data;
    else
        m.H0_Rep_ANOVA_Interaction_with_gp = data;
    end
    
    fprintf('    Bootstrap file saved successfully\n');
end
