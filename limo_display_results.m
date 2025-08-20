function limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,varargin)

% This function displays various results
% The arguments specify cases for the
% different kind of figures, thresholds etc ..
%
% FORMAT:
% limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,options)
%
% INPUTS:
%   Type      = type of images/plot to do
%               1 - 2D images with a intensity plotted as function of time (x) and electrodes (y)
%               2 - topographic plot a la eeglab
%               3 - plot the ERP data (original or modeled)
%   Filename  = Name of the file to image
%   PathName  = Path of the file to image
%   p         = threshold p value e.g. 0.05
%   MCC       = Multiple Comparison technique
%               1=None, 2= Cluster, 3=TFCE, 4=T max
%   LIMO      = LIMO structure
%   flag      = indicates to allow surfing the figure (1) or not (0)
%
% OPTIONAL INPUTS  (Usage: {''key'', value, ... })
% 'channels' : Provide the index of the channel to be used.
% 'regressor': Provide the index of the regressor to be used.
% 'plot3type': Type of plots to show when 'Type' is 3. Select between {'Original', 'Modeled', 'Adjusted'}
% 'sumstats' : Course plot summary statistics 'Mean' or 'Trimmed'
% 'restrict' : for time-frequency data, plot restrict plot to 'Time' or 'Frequency'
% 'dimvalue' : for time-frequency data, what value to resctrict on (e.g. restrict to 'Time' with dimvalue 5Hz)
%
% Although the function is mainly intented to be used via the GUI, some figures
% can be generated automatically, for instance limo_display_results(1,'R2.mat',pwd,0.05,5,LIMO,0);
% would load the R2.mat file from the current directory, and plot all
% electrodes/time frames F values thresholded using tfce at alpha 0.05
% topoplot and ERP like figures can't be automated since they require user
% input
%
% Cyril Pernet, Guillaume Rousselet, Carl Gaspar,
% Nicolas Chauveau, Andrew Stewart, Ramon Martinez-Cancino
%
% see also limo_stat_values limo_display_image topoplot limo_course_plot
% ----------------------------------------------------------------------
%  Copyright (C) LIMO Team 2019

try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else
        g = [];
    end
catch
    disp('limo_display_results() error: calling convention {''key'', value, ... } error'); return;
end

try g.channels;  catch, g.channels  = [];  end % No default values
try g.regressor; catch, g.regressor = [];  end % No default values
try g.plot3type; catch, g.plot3type = [];  end % No default values
try g.sumstats;  catch, g.sumstats  = [];  end % No default values
try g.restrict;  catch, g.restrict  = [];  end % No default values
try g.dimvalue;  catch, g.dimvalue  = [];  end % No default values

toplot = load(fullfile(PathName,FileName));
toplot = toplot.(cell2mat(fieldnames(toplot)));
if nargin <= 6
    flag = 1;
end

choice = 'use theoretical p values'; % threshold based on what is computed since H0 is used for clustering
% see limo_stat_values - discontinuated empirical threshold (misleading)

% Load LIMO structure if a path was provided
if ischar(LIMO)
    load(LIMO, 'LIMO');
end

if LIMO.design.bootstrap == 0
    if MCC == 2
        errordlg2('Clustering thresholding necessitates boostrap - invalid choice');
    elseif MCC == 3
        errordlg2('TFCE thresholding necessitates boostrap - invalid choice');
    elseif MCC == 4
        errordlg2('Maximum stat thresholding necessitates bootstrap - invalid choice');
    end
    MCC = 1;
end

if LIMO.design.bootstrap == 1 && LIMO.design.tfce == 0 && MCC == 3
    errordlg2('TFCE thresholding hasn''t been computed - invalid choice');
    MCC =1;
end

% -------------------------------------------------------------------------
% -------------------      LEVEL 1     ------------------------------------
% -------------------  SINGLE SUBJECT  ------------------------------------
% -------------------------------------------------------------------------
if LIMO.Level == 1 
    
    switch Type
        
        case{1}
            
            %--------------------------
            % imagesc of the results
            %--------------------------
            
            if  strcmpi(LIMO.design.type_of_analysis,'Mass-univariate')
                
                % univariate results from 1st level analysis
                % ------------------------------------------
                
                % if previously plotted recover data from the cache
                data_cached = 0;
                if isfield(LIMO,'cache')
                    try
                        if strcmpi(LIMO.cache.fig.name, FileName) && ...
                                LIMO.cache.fig.MCC == MCC && ...
                                LIMO.cache.fig.threshold == p
                            
                            disp('using cached data');
                            mask = LIMO.cache.fig.mask;
                            if isempty(mask)
                                data_cached = 0;
                            elseif sum(mask(:)) == 0
                                warndlg('  no values under threshold  ','no significant effect','modal');
                                return
                            else
                                M           = LIMO.cache.fig.pval;
                                mytitle     = LIMO.cache.fig.title;
                                toplot      = LIMO.cache.fig.stats;
                                data_cached = 1;
                                assignin('base','p_values',M)
                                assignin('base','mask',mask)
                            end
                        end
                    catch no_cache
                        fprintf('could not load cached data %s',no_cache.message)
                        data_cached = 0;
                    end
                end
                
                % ------------------
                % compute the plot
                % ------------------
                if data_cached == 0
                    
                    [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO,choice);
                    
                    if isempty(mask)
                        disp('no values computed'); return
                    elseif sum(mask(:)) == 0
                        warndlg('  no values under threshold  ','no significant effect','modal');
                        LIMO.cache.fig.name       = FileName;
                        LIMO.cache.fig.MCC        = MCC;
                        LIMO.cache.fig.stats      = [];
                        LIMO.cache.fig.threshold  = p;
                        LIMO.cache.fig.pval       = M;
                        LIMO.cache.fig.mask       = mask;
                        LIMO.cache.fig.title      = mytitle;
                        % do an exception for designs with just the constant
                        if strcmpi(FileName,'R2.mat') && size(LIMO.design.X,2)==1
                            mask = ones(size(mask)); LIMO.cache.fig.mask = mask;
                            mytitle = 'R^2 Coef unthresholded'; LIMO.cache.fig.title = mytitle;
                            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
                        else
                            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
                            return
                        end
                    else
                        assignin('base','p_values',M)
                        assignin('base','mask',mask)
                    end
                    
                    if contains(FileName,'R2','IgnoreCase',true)
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,1)); % plot R2 values instead of F
                        else
                            toplot = squeeze(toplot(:,:,1));
                        end
                        assignin('base','R2_values',toplot)
                        
                    elseif contains(FileName,'Condition_effect','IgnoreCase',true) || ...
                            contains(FileName,'Covariate_effect','IgnoreCase',true) || ...
                            contains(FileName,'Interaction_effect','IgnoreCase',true) || ...
                            contains(FileName,'semi_partial_coef','IgnoreCase',true)
                        
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,1)); % plot F values
                        else
                            toplot = squeeze(toplot(:,:,1));
                        end
                        
                        if contains(FileName,'semi_partial_coef','IgnoreCase',true)
                            assignin('base','semi_partial_coef',toplot)
                        else
                            assignin('base','F_values',toplot)
                        end
                        
                    elseif strcmpi(FileName(1:4),'con_')
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,4)); % plot T values
                        else
                            toplot = squeeze(toplot(:,:,4));
                        end
                        assignin('base','T_values',toplot)
                        
                    elseif strcmpi(FileName(1:4),'ess_')
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,end-1)); % plot F values
                        else
                            toplot = squeeze(toplot(:,:,end-1));
                        end
                        assignin('base','F_values',toplot)
                        
                    else
                        errordlg2('file not supported');
                        return
                    end
                end
                
                % -------------------------------------------------------------------------
                %              Actual plot takes place here
                % -------------------------------------------------------------------------
                if ~isempty(toplot)
                    
                    % cache the results for next time
                    if data_cached == 0 && ~all(mask(:)==1)
                        LIMO.cache.fig.name       = FileName;
                        LIMO.cache.fig.MCC        = MCC;
                        LIMO.cache.fig.stats      = toplot;
                        LIMO.cache.fig.threshold  = p;
                        LIMO.cache.fig.pval       = M;
                        LIMO.cache.fig.mask       = mask;
                        LIMO.cache.fig.title      = mytitle;
                        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
                    end
                    
                    if ndims(toplot)==3
                        limo_display_image_tf(LIMO,toplot,mask,mytitle,flag);
                    else
                        limo_display_image(LIMO,toplot,mask,mytitle,flag)
                    end
                end
                
            else
                
                % mutivariate results from 1st level analysis
                % ------------------------------------------
                if strncmp(FileName,'R2',2) || strncmp(FileName,'Condition_effect',16) || strncmp(FileName,'Covariate_effect',16)   % MANOVA PLOTTING
                    if strncmp(FileName,'R2',2)
                        
                        R2_EV     = load(fullfile(LIMO.dir,'R2_EV.mat'));
                        R2_EV     = R2_EV.R2_EV;
                        EV        = R2_EV(1:size(R2_EV,1),:); % no point plotting 0, just pick 5 1st Eigen values
                        R2_EV_var = load(fullfile(LIMO.dir,'R2_EV_var.mat'));
                        R2_EV_var = R2_EV_var.R2_EV_var;
                        test      =  sum(R2_EV_var(1,:) > 95) / size(R2_EV_var,2); % If more than 50% of the time-frames have a
                        % first eigenvalue with a proportion higher than 90%, the results of Roy's test are displayed,
                        if test > .50
                            choice = 'Roy';
                        else
                            choice = 'Pillai';
                        end
                        clear R2_EV;
                        
                        F_values(:,1) = squeeze(toplot(:,2));
                        F_values(:,2) = squeeze(toplot(:,4));
                        [M, mask, mytitle] = limo_mstat_values(Type,FileName,p,MCC,LIMO,choice);
                        if isempty(mask)
                            return
                        elseif sum(mask(:)) == 0
                            warndlg('  no values under threshold  ','no significant effect','modal');
                            return
                        else
                            toplot = squeeze(toplot(:,1)); % plot R2 values instead of F
                            assignin('base','F_values',F_values)
                            assignin('base','p_values',M)
                            assignin('base','mask',mask)
                            clear R2
                        end
                        
                    else
                        
                        if strcmpi(FileName(end-6:end),'_EV.mat')
                            FileName = [FileName(1:end-7) '.mat'];
                            toplot   = load(fullfile(PathName,FileName));
                            toplot   = toplot.(cell2mat(fieldnames(toplot)));
                        end
                        name   = sprintf('%s_%g_EV',FileName(1:end-4),str2double(FileName(max(strfind(FileName,'_')):end-4)));
                        EV     = load(fullfile(LIMO.dir,name));
                        EV     = EV.(cell2mat(fieldnames(EV)));
                        EV     = EV(1:size(Condition_effect_EV,1),:); % no point plotting 0, just pick 5 1st Eigen values
                        name   = sprintf('%s_%g_EV_var',FileName(1:end-4),str2double(FileName(max(strfind(FileName,'_')):end-4)));
                        EV_var = load(fullfile(LIMO.dir,name));
                        EV_var = EV_var.(cell2mat(fieldnames(EV_var)));
                        EV_var = EV_var(1:size(EV_var,1),:);
                        test =  sum(EV_var(1,:) > 95) / size(EV_var,2); % If more than 50% of the time-frames have a
                        %first eigenvalue with a proportion higher than 90%, the results of Roy's test are displayed,
                        if test > .50
                            choice = 'Roy';
                        else
                            choice = 'Pillai';
                        end
                        
                        F_values(:,1) = squeeze(toplot(:,1));
                        F_values(:,2) = squeeze(toplot(:,3));
                        [M, mask, mytitle] = limo_mstat_values(Type,FileName,p,MCC,LIMO,choice);
                        if isempty(mask)
                            return
                        elseif sum(mask(:)) == 0
                            warndlg('  no values under threshold  ','no significant effect','modal');
                            return
                        else
                            if strcmpi(choice,'Roy')
                                toplot = F_values(:,1);
                            else
                                toplot = F_values(:,2);
                            end
                            assignin('base','F_values',F_values)
                            assignin('base','p_values',M)
                            assignin('base','mask',mask)
                            clear R2
                        end
                    end
                    
                    figure; set(gcf,'Color','w');
                    % imagesc eigen values
                    subplot(3,3,[4 5 7 8]);
                    timevect = linspace(LIMO.data.start,LIMO.data.end,size(EV,2));
                    scale = EV; scale(scale==0)=NaN;
                    imagesc(timevect,1:size(EV,1),scale);
                    color_images_(scale,LIMO);  colorbar
                    ylabel('Eigen Values','Fontsize',14)
                    set(gca,'YTickLabel',{'1','2','3','4','5'});
                    title('non-zero Eigen values','Fontsize',14)
                    
                    % imagesc effect values
                    subplot(3,3,[1 2]);
                    scale = toplot'.*mask; scale(scale==0)=NaN;
                    imagesc(timevect,1,scale);
                    caxis([min(scale(:)), max(scale(:))]);
                    color_images_(scale,LIMO); xlabel(' ')
                    title(mytitle,'Fontsize',18); colorbar
                    ylabel(' '); set(gca,'YTickLabel',{''});
                    
                    % ERP plot1 - Roy -
                    subplot(3,3,6);
                    plot(timevect, F_values(:,1),'LineWidth',3); grid on; axis tight
                    mytitle2 = sprintf('F values - Roy');
                    title(mytitle2,'FontSize',14)
                    
                    % ERP plot2 - Pillai -
                    subplot(3,3,9);
                    plot(timevect, F_values(:,2),'LineWidth',3); grid on; axis tight
                    mytitle2 = sprintf('F value - Pillai');
                    title(mytitle2,'FontSize',14)
                end % end of MANOVA PLOTTING
                
                if strncmp(FileName,'Discriminant_coeff',18) || strncmp(FileName,'Discriminant_scores',19)
                    Discriminant_coeff      = load(fullfile(LIMO.dir,'Discriminant_coeff'));
                    Discriminant_coeff      = Discriminant_coeff.Discriminant_coeff;
                    Discriminant_scores     = load(fullfile(LIMO.dir,'Discriminant_scores'));
                    Discriminant_scores     = Discriminant_scores.Discriminant_scores;
                    Condition_effect_EV_var = load(fullfile(LIMO.dir,'Condition_effect_1_EV_var.mat'));
                    Condition_effect_EV_var = Condition_effect_EV_var.Condition_effect_EV_var;
                    
                    time = linspace(LIMO.data.start,LIMO.data.end, size(Discriminant_coeff,2));
                    input_title = sprintf('which time-frame to plot (in ms)?: ');
                    timepoint = inputdlg(input_title,'Plotting option');
                    t = dsearchn(time', str2double(timepoint{1}));
                    groupcolors = 'rgbcwmryk';
                    groupsymbols = 'xo*+.sdv<>';
                    [class,~] = find(LIMO.design.X(:,1:LIMO.design.nb_conditions)');
                    k = LIMO.design.nb_conditions;
                    
                    if k>2
                        figure;set(gcf,'Color','w');
                        subplot(2,2,[1 2]); % 2D plot of two discriminant functions
                        gscatter(squeeze(Discriminant_scores(1,t,:)), squeeze(Discriminant_scores(2,t,:)), class, groupcolors(1:k), groupsymbols(1:k));
                        grid on; axis tight;
                        xlabel(['Z1, var: ' num2str(round(Condition_effect_EV_var(1,t)),2) '%'],'Fontsize',14);
                        ylabel(['Z2, var: ' num2str(round(Condition_effect_EV_var(2,t)),2) '%'],'Fontsize',14);
                        title(['Results of the discriminant analysis at ' num2str(time(t)) 'ms'], 'Fontsize', 18);
                        z1 = subplot(2,2,3); % First discriminant coeff
                        cc = limo_color_images(Discriminant_coeff(:,t,1)); % get a color map commensurate to that
                        topoplot(Discriminant_coeff(:,t,1),LIMO.data.chanlocs, 'electrodes','off','style','map','whitebk', 'on','colormap',cc);colorbar;
                        title('Z1','Fontsize',14); colormap(z1, 'hot');
                        z2 = subplot(2,2,4); % Second discriminant coeff
                        cc = limo_color_images(Discriminant_coeff(:,t,2)); % get a color map commensurate to that
                        topoplot(Discriminant_coeff(:,t,2),LIMO.data.chanlocs, 'electrodes','off','style','map','whitebk', 'on','colormap',cc);colorbar;
                        title('Z2','Fontsize',14); colormap(z2, 'hot');
                    elseif k==2
                        figure;set(gcf,'Color','w');
                        subplot(2,2,[1 2]); % 1D plot of two discriminant functions
                        data = squeeze(Discriminant_scores(1,t,:));
                        class1 = data(class == 1);
                        class2 = data(class == 2);
                        histogram(class1, 'BinWidth', 0.1);
                        hold on
                        histogram(class2,'BinWidth',0.1);
                        hold off
                        legend show
                        grid on; axis tight;
                        xlabel(['Z1, var: ' num2str(round(Condition_effect_EV_var(1,t)),2) '%'],'Fontsize',14);
                        title(['Results of the discriminant analysis at ' num2str(time(t)) 'ms'], 'Fontsize', 18);
                        z1 = subplot(2,2,[3,4]); % First discriminant coeff
                        cc = limo_color_images(Discriminant_coeff(:,t,1)); % get a color map commensurate to that
                        topoplot(Discriminant_coeff(:,t,1),LIMO.data.chanlocs, 'electrodes','off','style','map','whitebk', 'on','colormap',cc);colorbar;
                        title('Z1','Fontsize',14); colormap(z1, 'hot');
                    end
                    % Create tabbed figure with 2D and 3D views for discriminant analysis
                    if ~strcmpi(LIMO.Analysis,'Time-Frequency') % Only for 2D data
                        % Create main tabbed figure
                        main_fig = figure('Color','w','Name','Discriminant coefficients Z1 - Results');
                        tab_group = uitabgroup(main_fig);
                        
                        % Size figure to match MATLAB desktop
                        try
                            desktop_pos = get(0, 'ScreenSize');  % [left, bottom, width, height] of screen
                            % Use 90% of screen size with some margin
                            fig_width = desktop_pos(3) * 0.9;
                            fig_height = desktop_pos(4) * 0.8;
                            fig_left = desktop_pos(3) * 0.05;  % 5% margin from left
                            fig_bottom = desktop_pos(4) * 0.1;  % 10% margin from bottom
                            set(main_fig, 'Position', [fig_left, fig_bottom, fig_width, fig_height]);
                        catch
                            % Fallback to default if sizing fails
                            set(main_fig, 'Position', [100, 100, 1200, 800]);
                        end
                        
                        % Create 2D view tab first
                        tab_2d = uitab(tab_group, 'Title', '2D View');
                        limo_display_image(LIMO,abs(Discriminant_coeff(:,:,1)),abs(Discriminant_coeff(:,:,1)),'Discriminant coefficients Z1',flag,tab_2d);
                        
                        % Create 3D view tab  
                        tab_3d = uitab(tab_group, 'Title', '3D View');
                        axes_3d = axes('Parent', tab_3d);
                        
                        % Prepare data
                        if strcmpi(LIMO.Analysis,'Time')
                            if exist('timevect','var')
                                xvect = timevect;
                            else
                                xvect = linspace(LIMO.data.start, LIMO.data.end, size(toplot,2));
                            end
                            xlabel_text = 'Time (ms)';
                        else % Frequency
                            if exist('freqvect','var')
                                xvect = freqvect;
                            else
                                xvect = linspace(LIMO.data.freqlist(1), LIMO.data.freqlist(end), size(toplot,2));
                            end
                            xlabel_text = 'Frequency (Hz)';
                        end
                        
                        % Create channel vector
                        num_channels = size(toplot,1);
                        channel_vector = 1:num_channels;
                        
                        % Create meshgrid for 3D plotting
                        [X, Z] = meshgrid(xvect, channel_vector);
                        
                        % Plot 3D surface
                        surf(axes_3d, X, toplot, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
                        
                        % Customize the plot
                        % Set Y-axis limits to encompass full range of statistical values
                        y_min = min(toplot(:));
                        y_max = max(toplot(:));
                        y_range = y_max - y_min;
                        y_margin = y_range * 0.05; % Add 5% margin
                        ylim(axes_3d, [y_min - y_margin, y_max + y_margin]);
                        xlabel(axes_3d, xlabel_text, 'FontSize', 14);
                        ylabel(axes_3d, 'F values', 'FontSize', 14);
                        zlabel(axes_3d, 'Channels', 'FontSize', 14);
                        title(axes_3d, 'Discriminant coefficients Z1 - 3D View', 'FontSize', 16);
                        
                        % Set colormap
                        colormap(axes_3d, limo_color_images(toplot));
                        cb = colorbar(axes_3d);
                set(cb, 'Color', 'k'); % Make colorbar text black
                        
                        % Adjust view angle
                        view(axes_3d, [-45 30]);
                        
                        % Set z-axis ticks and labels
                        if num_channels <= 20
                            set(axes_3d, 'ZTick', channel_vector);
                            if exist('label_electrodes','var') && ~isempty(label_electrodes)
                                set(axes_3d, 'ZTickLabel', flipud(label_electrodes(:)));
                            end
                        else
                            % For many channels, show fewer labels
                            tick_indices = round(linspace(1, num_channels, min(10, num_channels)));
                            set(axes_3d, 'ZTick', tick_indices);
                            if exist('label_electrodes','var') && ~isempty(label_electrodes)
                                set(axes_3d, 'ZTickLabel', label_electrodes(tick_indices));
                            end
                        end
                        
                        % Add grid and lighting
                        grid(axes_3d, 'on');
                        light('Position', [-1 -1 2], 'Style', 'local');
                        lighting gouraud;
                        material dull;
                        
                        % Apply cluster visualization if available
                        if exist('mask','var') && ~isempty(mask) && MCC == 2 && max(mask(:)) > 1
                            hold(axes_3d, 'on');
                            % Create cluster-specific visualization  
                            n_cluster = max(mask(:));
                            % Use distinguishable colors, avoiding light colors
                            cluster_colors = [
                                0.0000 0.4470 0.7410;  % Blue
                                0.8500 0.3250 0.0980;  % Orange  
                                0.9290 0.6940 0.1250;  % Yellow
                                0.4940 0.1840 0.5560;  % Purple
                                0.4660 0.6740 0.1880;  % Green
                                0.3010 0.7450 0.9330;  % Cyan
                                0.6350 0.0780 0.1840   % Dark Red
                            ];
                            % Repeat colors if more clusters than predefined colors
                            if n_cluster > size(cluster_colors, 1)
                                cluster_colors = repmat(cluster_colors, ceil(n_cluster/size(cluster_colors,1)), 1);
                            end
                            
                            % Store handles for legend
                            h_clusters = [];
                            
                            for cluster_id = 1:n_cluster
                                cluster_data = toplot;
                                cluster_data(mask ~= cluster_id) = NaN;
                                
                                % Plot each cluster with a different color
                                h_cluster = surf(X, cluster_data, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                                set(h_cluster, 'FaceColor', cluster_colors(cluster_id, :));
                                h_clusters(cluster_id) = h_cluster;
                            end
                            
                            % Update title for 3D view
                            title([mytitle ' - 3D View'], 'FontSize', 16);
                            
                            % Add threshold plane for reference
                            if p > 0 && p <= 1  % Valid p-value threshold
                                % Compute statistical threshold
                                % Try to get degrees of freedom from LIMO structure
                                if isfield(LIMO, 'model') && isfield(LIMO.model, 'model_df')
                                    df_vals = LIMO.model.model_df;
                                    if length(df_vals) >= 2
                                        stat_threshold = finv(1-p, df_vals(1), df_vals(2));  % F-distribution
                                        is_ttest = false;
                                    else
                                        stat_threshold = tinv(1-p/2, df_vals(1));       % t-distribution (two-tailed)
                                        is_ttest = true;
                                    end
                                else
                                    % Use a default threshold and detect t-test from title
                                    stat_threshold = -log10(p) * 2; % Rough approximation
                                    is_ttest = contains(lower(mytitle), {'ttest', 't values', 'paired', 'one sample'});
                                end
                                
                                % Create threshold plane as solid surface
                                [X_thresh, Z_thresh] = meshgrid(xvect, channel_vector);
                                threshold_handles = [];
                                threshold_labels = {};
                                
                                % Always add positive threshold
                                Y_thresh_pos = ones(size(X_thresh)) * stat_threshold;
                                h_thresh_pos = surf(X_thresh, Y_thresh_pos, Z_thresh, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'FaceColor', 'red');
                                threshold_handles(end+1) = h_thresh_pos;
                                
                                % Only add negative threshold for explicit t-tests with significant negative clusters
                                has_negative_clusters = is_ttest && any(toplot(:) < -stat_threshold) && any(mask(:) > 0 & toplot(:) < 0);
                                if has_negative_clusters
                                    Y_thresh_neg = ones(size(X_thresh)) * (-stat_threshold);
                                    h_thresh_neg = surf(X_thresh, Y_thresh_neg, Z_thresh, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'FaceColor', 'blue');
                                    threshold_handles(end+1) = h_thresh_neg;
                                    threshold_labels = {sprintf('+Threshold (p=%.3f)', p), sprintf('-Threshold (p=%.3f)', p)};
                                else
                                    threshold_labels = {sprintf('Threshold (p=%.3f)', p)};
                                end
                                
                                % Create legend with proper handles and light mode styling
                                legend_entries = cell(n_cluster + length(threshold_labels), 1);
                                legend_handles = [h_clusters, threshold_handles];
                                for i = 1:n_cluster
                                    legend_entries{i} = sprintf('Cluster %d', i);
                                end
                                for i = 1:length(threshold_labels)
                                    legend_entries{n_cluster + i} = threshold_labels{i};
                                end
                                h_legend = legend(legend_handles, legend_entries, 'Location', 'best');
                                set(h_legend, 'Color', 'white', 'TextColor', 'black', 'EdgeColor', 'black');
                            end
                        end
                        
                        % Store the 3D plot data
                        assignin('base', 'Plot3D_X', X);
                        assignin('base', 'Plot3D_Y', toplot);
                        assignin('base', 'Plot3D_Z', Z);
                    else
                        % If no 3D plot created for discriminant analysis, use normal display
                        limo_display_image(LIMO,abs(Discriminant_coeff(:,:,1)),abs(Discriminant_coeff(:,:,1)),'Discriminant coefficients Z1',flag)
                    end
                    
                    %                     figure;set(gcf,'Color','w');
                    %                     for t=1:size(Discriminant_coeff,2)
                    %                     topoplot(Discriminant_coeff(:,t,1),LIMO.data.chanlocs, 'electrodes','numbers','style','map');
                    %                     title(['Discriminant values first discriminant at timepoint ' num2str(t) ' corresponding to ' num2str(time(t)) ' ms']);
                    %                     pause(.01)
                    %                     end;
                end
                
                if strncmp(FileName,'Linear_Classification',21)
                    Linear_Classification = load(fullfile(LIMO.dir,'Linear_Classification'));
                    Linear_Classification = Linear_Classification.Linear_Classification;
                    [~, mask, mytitle] = limo_mstat_values(Type,FileName,p,MCC,LIMO,choice);
                    timevect = linspace(LIMO.data.start,LIMO.data.end,size(Linear_Classification,1));
                    figure;set(gcf,'Color','w');
                    subplot(3,1,[1 2]); % lineplot
                    plot(timevect,Linear_Classification(:,2),'LineWidth',3);title(mytitle, 'Fontsize', 18);
                    ylabel('decoding accuracies', 'Fontsize', 14);grid on; axis tight; hold on;
                    plot(timevect, Linear_Classification(:,2) + 2*Linear_Classification(:,3), 'k-','LineWidth',1); hold on;
                    plot(timevect, Linear_Classification(:,2) - 2*Linear_Classification(:,3), 'k-','LineWidth',1)
                    line([0,0],[0,1], 'color', 'black')
                    subplot(3, 1, 3); % imagesc accuracies
                    toplot = Linear_Classification(:,2); scale = toplot'.*mask;scale(scale==0)=NaN;
                    imagesc(timevect,1,scale);xlabel('Time in ms');
                    color_images_(scale,LIMO);
                    ylabel(' '); set(gca,'YTickLabel',{''});
                end
                
                if strncmp(FileName,'Quadratic_Classification',24)
                    Quadratic_Classification = load(fullfile(LIMO.dir,'Quadratic_Classification'));
                    Quadratic_Classification = Quadratic_Classification.Quadratic_Classification;
                    timevect = linspace(LIMO.data.start,LIMO.data.end,size(Quadratic_Classification,1));
                    figure;set(gcf,'Color','w');
                    subplot(3,1,[1 2]); % lineplot
                    plot(timevect,Quadratic_Classification(:,2),'LineWidth',3);title('CV quadratic decoding accuracies +/- 2SD', 'Fontsize', 18);
                    ylabel('decoding accuracies', 'Fontsize', 14);grid on; axis tight; hold on
                    plot(timevect, Quadratic_Classification(:,2) + 2*Quadratic_Classification(:,3), 'k-','LineWidth',1); hold on;
                    plot(timevect, Quadratic_Classification(:,2) - 2*Quadratic_Classification(:,3), 'k-','LineWidth',1)
                    line([0,0],[0,1], 'color', 'black')
                    subplot(3,1, 3); % imagesc plot training accuracies
                    scale = Quadratic_Classification(:,2)'; scale(scale==0)=NaN;
                    imagesc(timevect,1,scale);xlabel('Time in ms');
                    color_images_(scale,LIMO);
                    ylabel(' '); set(gca,'YTickLabel',{''});
                end
            end
            
            
        case{2}
            
            %--------------------------
            % topoplot
            %--------------------------
            
            % univariate results from 1st level analysis
            % ------------------------------------------
            
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                warndlg('topoplot not supported for 3D data')
            else
                EEG.trials = 1;
                EEG.chanlocs = LIMO.data.chanlocs;
                if strcmpi(LIMO.Analysis,'Time')
                    EEG.xmin  = LIMO.data.start / 1000;% in msec
                    EEG.xmax  = LIMO.data.end / 1000;  % in msec
                    EEG.times = LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000; % in sec;
                    if length(EEG.times) > 2
                        EEG.times = [EEG.times(1) EEG.times(end)];
                    end
                elseif strcmpi(LIMO.Analysis,'Frequency')
                    EEG.xmin = LIMO.data.freqlist(1);
                    EEG.xmax = LIMO.data.freqlist(end);
                    freqlist = inputdlg('specify frequency range e.g. [5:2:40]','Choose Frequencies to plot');
                    if isempty(freqlist)
                        return
                    else
                        if contains(cell2mat(freqlist),':')
                            EEG.freq = eval(cell2mat(freqlist));
                        else
                            EEG.freq = str2double(cell2mat(freqlist));
                        end
                        
                        if min(EEG.freq)<EEG.xmin || max(EEG.freq)>EEG.xmax
                            errordlg('slected frequency out of bound'); return
                        end
                    end
                end
                
                if contains(FileName,'R2','IgnoreCase',true)
                    if size(LIMO.design.X,2)==1
                        EEG.data    = squeeze(toplot(:,:,1));
                        EEG.setname = 'R2 values for the mean';
                    else
                        EEG.data    = squeeze(toplot(:,:,2));
                        EEG.setname = 'R2 - F values';
                    end
                    call_topolot(EEG,FileName,LIMO.Analysis)
                elseif contains(FileName,'Condition_effect','IgnoreCase',true) || ...
                        contains(FileName,'Covariate_effect','IgnoreCase',true) || ...
                        contains(FileName,'Interaction_effect','IgnoreCase',true)  || ...
                        contains(FileName,'Condition_effect','IgnoreCase',true)
                    EEG.data    = squeeze(toplot(:,:,1));
                    call_topolot(EEG,FileName,LIMO.Analysis)
                elseif contains(FileName,'con','IgnoreCase',true) || contains(FileName,'ess','IgnoreCase',true)
                    EEG.data    = squeeze(toplot(:,:,end-1));
                    call_topolot(EEG,FileName,LIMO.Analysis)
                elseif contains(FileName,'semi partial_coef.mat','IgnoreCase',true)
                    regressor = str2double(cell2mat(inputdlg('which regressor(s) to plot (e.g. 1:3)','Plotting option')));
                    if max(regressor) > size(toplot,3); errordlg('error in regressor number'); return; end
                    for b = regressor
                        EEG.data     = squeeze(toplot(:,:,b,1));
                        call_topolot(EEG,FileName,LIMO.Analysis)
                    end
                else
                    disp('file not supported');
                    return
                end
                
                if contains(FileName,'con','IgnoreCase',true)
                    assignin('base','T_values',EEG.data);
                else
                    assignin('base','F_values',EEG.data);
                end
            end
            
        case{3}
            
            %--------------------------
            % Time course / Power
            %--------------------------
            
            % which variable(s) to plot
            % ----------------------
            if isempty(g.regressor)
                input_title = sprintf('which regressor to plot?: 1 to %g ',size(LIMO.design.X,2));
                regressor   = inputdlg(input_title,'Plotting option');
            else
                regressor = g.regressor;
            end
            
            if isempty(regressor); disp('selection aborded'); return; end
            regressor = cell2mat(regressor);
            if isempty(regressor); disp('selection aborded'); return; end
            if ~contains(regressor,'['); regressor=['[' regressor ']']; end
            if ischar(regressor); regressor=str2num(regressor); end %#ok<ST2NM>
            regressor = sort(regressor);
            
            if max(regressor) > size(LIMO.design.X,2)
                errordlg('invalid regressor number'); 
            end
            
            categorical = sum(LIMO.design.nb_conditions) + sum(LIMO.design.nb_interactions);
            if max(regressor) == size(LIMO.design.X,2)
                tmp = regressor(1:end-1);
            else
                tmp = regressor;
            end
            
            cat = sum(tmp<=categorical); cont = sum(tmp>categorical);
            if cat >=1 && cont >=1
                errordlg2('you can''t plot categorical and continuous regressors together'); return
            end
            
            % which data type to make
            % ------------------------
            if isempty(g.plot3type) && ~any(strcmpi(g.plot3type,{'Original','Modelled','Adjusted'}))
                extra = questdlg('Which data type to plot?','Options','Original','Modelled','Adjusted','Adjusted');
            else
                extra = g.plot3type;
            end
            if isempty(extra)
                return
            elseif strcmpi(extra,'Original')
                if regressor == size(LIMO.design.X,2)
                    errordlg('you can''t plot adjusted mean for original data'); return
                end
            end
            
            % timing /frequency info
            % -----------------------
            if strcmpi(LIMO.Analysis,'Time')
                timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
            elseif strcmpi(LIMO.Analysis,'Frequency')
                freqvect=LIMO.data.freqlist';
            elseif strcmpi(LIMO.Analysis,'Time-Frequency')
                timevect = linspace(LIMO.data.start,LIMO.data.end,LIMO.data.size4D(3));
                freqvect = linspace(LIMO.data.lowf,LIMO.data.highf,LIMO.data.size4D(2));
            end
            
            % which channel/frequency to plot
            % --------------------------------
            if isempty(g.channels)
                channel = inputdlg('which channel to plot','Plotting option');
            else
                channel = g.channels;
            end
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                disp('loading the 4D data ...')
                frequency = inputdlg('which Frequency to plot','Plotting option');
            else
                frequency = [];
            end
            
            if strcmpi(channel,'') || strcmpi(frequency,'')
                disp('looking for max');
                R2 = load(fullfile(LIMO.dir,'R2.mat'));
                R2 = R2.R2;
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tmp = squeeze(R2(:,:,:,1)); clear R2
                    [e,f,~] = ind2sub(size(tmp),find(tmp==max(tmp(:))));
                    if length(e) ~= 1; e = e(1); f = f(1); end
                    
                    if strcmpi(channel,'')
                        channel = e;
                    else
                        channel = eval(cell2mat(channel));
                    end
                    if size(channel) > 1
                        errordlg('invalid channel choice'); return
                    elseif channel > size(LIMO.data.chanlocs,2) || channel < 1
                        errordlg('invalid channel number'); return
                    end
                    
                    if strcmpi(frequency,'')
                        freq_index = f;
                        frequency = freqvect(freq_index);
                    else
                        frequency = eval(cell2mat(frequency));
                    end
                    if size(frequency) > 1
                        errordlg('invalid frequency choice'); return
                    elseif frequency > LIMO.data.tf_freqs(end) || frequency < LIMO.data.tf_freqs(1)
                        errordlg('invalid frequency number'); return
                    end
                else
                    tmp = squeeze(R2(:,:,1)); clear R2
                    [channel,~] = ind2sub(size(tmp),find(tmp==max(tmp(:))));
                end
                clear tmp
            else
                channel = eval(cell2mat(channel));
                if size(channel) > 1
                    errordlg('invalid channel choice'); return
                elseif channel > size(LIMO.data.chanlocs,2) || channel < 1
                    errordlg('invalid channel number'); return
                end
                
                if ~isempty(frequency)
                    frequency = eval(cell2mat(frequency));
                    if size(frequency) > 1
                        errordlg('invalid frequency choice');
                    elseif frequency > freqvect(end) || frequency < freqvect(1)
                        errordlg('invalid frequency number');
                    end
                    % pick the nearest frequency index
                    [~, freq_index] = min(abs(freqvect-frequency ));
                    frequency = freqvect(freq_index);
                end
            end
            
            % down to business
            % ----------------------
            data_cached = 0;
            if isfield(LIMO,'cache')
                if strcmpi(LIMO.Analysis,'Time-Frequency') && isfield(LIMO.cache,'ERPplot')
                    
                    if mean([LIMO.cache.Courseplot.channel == channel ...
                            LIMO.cache.Courseplot.regressor == regressor ...
                            LIMO.cache.Courseplot.frequency == frequency]) == 1 ...
                            && strcmpi('LIMO.cache.Courseplot.extra',extra)
                        
                        if sum(regressor <= categorical) == length(regressor)
                            average = LIMO.cache.Courseplot.average;
                            ci = LIMO.cache.Courseplot.ci;
                            mytitle = LIMO.cache.Courseplot.title;
                            disp('using cached data');
                            data_cached = 1;
                        else
                            continuous = LIMO.cache.Courseplot.continuous;
                            mytitle = LIMO.cache.Courseplot.title;
                            disp('using cached data');
                            data_cached = 1;
                        end
                    end
                    
                elseif strcmpi(LIMO.Analysis,'Time') && isfield(LIMO.cache,'Courseplot') || ...
                        strcmpi(LIMO.Analysis,'Frequency') && isfield(LIMO.cache,'Courseplot')
                    
                    if length(LIMO.cache.Courseplot.regressor) == length(channel)
                        if mean([LIMO.cache.Courseplot.channel == channel ...
                                LIMO.cache.Courseplot.regressor == regressor]) == 1  ...
                                && strcmpi('LIMO.cache.Courseplot.extra',extra)
                            
                            if sum(regressor <= categorical) == length(regressor)
                                average = LIMO.cache.Courseplot.average;
                                ci = LIMO.cache.Courseplot.ci;
                                mytitle = LIMO.cache.Courseplot.title;
                                disp('using cached data');
                                data_cached = 1;
                            else
                                continuous = LIMO.cache.Courseplot.continuous;
                                mytitle = LIMO.cache.Courseplot.title;
                                disp('using cached data');
                                data_cached = 1;
                            end
                        end
                    end
                end
            end
            
            % no cache = compute
            if data_cached == 0
                
                probs = [p/2; 1-p/2];
                z = norminv(probs);
                
                if strcmpi(extra,'Original')
                    Yr = load(fullfile(LIMO.dir,'Yr.mat')); Yr = Yr.Yr;
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        for i=length(regressor):-1:1
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                data     = squeeze(Yr(channel,freq_index,:,index{i}));
                                mytitle  = sprintf('Original ERSP at \n channel %s (%g) at %g Hz', LIMO.data.chanlocs(channel).labels, channel, frequency);
                            else
                                data     = squeeze(Yr(channel,:,index{i}));
                            end
                            average(i,:) = nanmean(data,2);
                            se           = nanstd(data,0,2) ./ sqrt(numel(index{i}));
                            ci(i,:,:)    = repmat(average(i,:),2,1) + repmat(se',2,1).*repmat(z,1,size(Yr,2));
                        end
                        
                        if strcmpi(LIMO.Analysis,'Time')
                            mytitle = sprintf('Original ERP at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        elseif strcmpi(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Original Power Spectrum at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        end
                    else % continuous variable
                        for i=length(regressor):-1:1
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                continuous(i,:,:) = Yr(channel,freq_index,:,sorting_values);
                                mytitle{i}        = sprintf('Original single trials \n sorted by regressor %g channel %s (%g)', regressor(i), LIMO.data.chanlocs(channel).labels, channel);
                            else
                                continuous(i,:,:) = Yr(channel,:,sorting_values);
                                mytitle{i}        = sprintf('Original single trials \n sorted by regressor %g \n channel %s (%g) at %s Hz', regressor(i), LIMO.data.chanlocs(channel).labels, channel, frequency);
                            end
                        end
                    end
                    clear Yr
                elseif strcmpi(extra,'Modelled')
                    Betas = load(fullfile(LIMO.dir,'Betas.mat'));
                    Betas = Betas.Betas;
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        Betas = squeeze(Betas(channel,freq_index,:,:));
                    else
                        Betas = squeeze(Betas(channel,:,:));
                    end
                    Yh = (LIMO.design.X*Betas')'; % modelled data
                    
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        Yr = load(fullfile(LIMO.dir,'Yr.mat')); Yr = Yr.Yr;
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            Yr = squeeze(Yr(:,freq_index,:,:));
                        end
                        R = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
                        
                        for i=length(regressor):-1:1
                            index{i}     = find(LIMO.design.X(:,regressor(i)));
                            data         = squeeze(Yh(:,index{i}));
                            average(i,:) = mean(data,2);
                            var          = diag(((R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')'*(R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')) / LIMO.model.model_df(2));
                            CI           = sqrt(var/size(index{i},1))*z';
                            ci(i,:,:)    = (repmat(mean(data,2),1,2)+CI)';
                        end
                        
                        if strcmpi(LIMO.Analysis,'Time')
                            mytitle = sprintf('Modelled ERP at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        elseif strcmpi(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Modelled Power Spectrum at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle = sprintf('Modelled ERSP \n channel %s (%g) at %g Hz', LIMO.data.chanlocs(channel).labels, channel, frequency);
                        end
                    else % continuous variable
                        for i=length(regressor):-1:1
                            index{i}                         = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values] = sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            continuous(i,:,:)                = Yh(:,sorting_values);
                            
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                mytitle{i} = sprintf('Modelled single trials \n sorted by regressor %g \n channel %s (%g) at %g Hz', regressor(i), LIMO.data.chanlocs(channel).labels, channel, frequency);
                            else
                                mytitle{i} = sprintf('Modelled single trials \n sorted by regressor %g channel %s (%g)', regressor(i), LIMO.data.chanlocs(channel).labels, channel);
                            end
                        end
                    end
                else % Adjusted
                    allvar = 1:size(LIMO.design.X,2)-1;
                    allvar(regressor)=[];
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        Yr    = load(fullfile(LIMO.dir,'Yr.mat'));
                        Yr    = squeeze(Yr.Yr(channel,freq_index,:,:));
                        Betas = load(fullfile(LIMO.dir,'Betas.mat'));
                        Betas = squeeze(Betas.Betas(channel,freq_index,:,:));
                    else
                        Yr    = load(fullfile(LIMO.dir,'Yr.mat'));
                        Yr    = squeeze(Yr.Yr(channel,:,:));
                        Betas = load(fullfile(LIMO.dir,'Betas.mat'));
                        Betas = squeeze(Betas.Betas(channel,:,:));
                    end
                    confounds = (LIMO.design.X(:,allvar)*Betas(:,allvar)')';
                    Ya        = Yr - confounds; clear Yr Betas confounds;
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        for i=length(regressor):-1:1
                            index{i}     = find(LIMO.design.X(:,regressor(i)));
                            data         = squeeze(Ya(:,index{i}));
                            average(i,:) = nanmean(data,2);
                            se           = nanstd(data,0,2) ./ sqrt(numel(index{i}));
                            ci(i,:,:)    = repmat(average(i,:),2,1) + repmat(se',2,1).*repmat(z,1,size(Ya,1));
                        end
                        if strcmpi(LIMO.Analysis,'Time')
                            mytitle = sprintf('Adjusted ERP at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        elseif strcmpi(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Adjusted Power Spectrum at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle = sprintf('Adjusted ERSP channel %s (%g) at %g Hz', LIMO.data.chanlocs(channel).labels, channel, frequency);
                        end
                    else % continuous variable
                        for i=length(regressor):-1:1
                            index{i}                         = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values] = sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            continuous(i,:,:)                = Ya(:,sorting_values);
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                mytitle{i} = sprintf('Adjusted single trials \n sorted by regressor \n %g channel %s (%g) at %g Hz', regressor(i), LIMO.data.chanlocs(channel).labels, channel, frequency);
                            else
                                mytitle{i} = sprintf('Adjusted single trials \n sorted by regressor %g channel %s (%g)', regressor(i), LIMO.data.chanlocs(channel).labels, channel);
                            end
                        end
                    end
                end
            end
            
            % make the figure(s)
            % ------------------
            figure;set(gcf,'Color','w')
            if sum(regressor <= categorical) == length(regressor)
                for i=1:size(average,1)
                    if i==1
                        colorOrder = get(gca, 'ColorOrder');
                        colorOrder = repmat(colorOrder,ceil(size(average,1)/size(colorOrder,1)),1);
                    end
                    
                    if strcmpi(LIMO.Analysis,'Frequency')
                        try
                            plot(freqvect,average(i,:),'LineWidth',1.5,'Color',colorOrder(i,:)); hold on
                        catch
                            freqvect = linspace(LIMO.data.start,LIMO.data.end,size(average,2));
                            plot(freqvect,average(i,:),'LineWidth',1.5,'Color',colorOrder(i,:)); hold on
                        end
                    else
                        plot(timevect,average(i,:),'LineWidth',1.5,'Color',colorOrder(i,:)); hold on
                    end
                    
                    x = squeeze(ci(i,1,:)); y = squeeze(ci(i,2,:));
                    if strcmpi(LIMO.Analysis,'Frequency')
                        fillhandle = patch([reshape(freqvect, 1, numel(freqvect)) fliplr(reshape(freqvect, 1, numel(freqvect)))], [x' fliplr(y')], colorOrder(i,:));
                    else
                        fillhandle = patch([reshape(timevect, 1, numel(timevect)) fliplr(reshape(timevect, 1, numel(timevect)))], [x',fliplr(y')], colorOrder(i,:));
                    end
                    set(fillhandle,'EdgeColor',colorOrder(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);
                end
                
                % if regressor spans columns of an effect, plot significant time frames
                index = 1; index2 = LIMO.design.nb_conditions(1);
                for i=1:length(LIMO.design.nb_conditions)
                    effect = index:index2;
                    if length(regressor) == length(effect)
                        if mean(regressor == effect) == 1
                            name = sprintf('Condition_effect_%g.mat',i);
                            % load(name);
                            if isfield(LIMO,'cache') && isfield(LIMO,'fig')
                                if strcmpi(LIMO.cache.fig.name,name) && ...
                                        LIMO.cache.fig.MCC == MCC && ...
                                        LIMO.cache.fig.threshold == p
                                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                                        sig = single(LIMO.cache.fig.mask(channel,freq_index,:)); sig(sig==0)=NaN;
                                    else
                                        sig = single(LIMO.cache.fig.mask(channel,:)); sig(sig==0)=NaN;
                                    end
                                end
                            else
                                if strcmpi(LIMO.Analysis,'Time-Frequency')
                                    [~, mask] = limo_stat_values(name,p,MCC,LIMO,choice);
                                    sig = single(squeeze(mask(channel,freq_index,:))); sig(sig==0)=NaN;
                                else
                                    [~, mask] = limo_stat_values(name,p,MCC,LIMO,choice);
                                    sig = single(mask(channel,:)); sig(sig==0)=NaN;
                                end
                            end
                            h = axis;
                            if strcmpi(LIMO.Analysis,'Frequency')
                                plot(freqvect,sig.*h(3),'r*','LineWidth',2)
                            else
                                plot(timevect,sig.*h(3),'r*','LineWidth',2)
                            end
                            break
                        end
                    else
                        index = index+LIMO.design.nb_conditions(i);
                        if i<length(LIMO.design.nb_conditions)
                            index2 = LIMO.design.nb_conditions(i)+LIMO.design.nb_conditions(i+1);
                        end
                    end
                end
                
                if LIMO.design.nb_interactions ~= 0
                    index = sum(LIMO.design.nb_conditions)+1; index2 = sum(LIMO.design.nb_conditions)+LIMO.design.nb_interactions(1);
                    for i=1:length(LIMO.design.nb_interactions)
                        effect = index:index2;
                        if length(regressor) == length(effect)
                            if mean(regressor == effect) == 1
                                name = sprintf('Interaction_effect_%g.mat',i);
                                % load(name);
                                if isfield(LIMO,'cache') && isfield(LIMO,'fig')
                                    if strcmpi(LIMO.cache.fig.name,name) && ...
                                            LIMO.cache.fig.MCC == MCC && ...
                                            LIMO.cache.fig.threshold == p
                                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                                            sig = single(squeeze(LIMO.cache.fig.mask(channel,freq_index,:))); sig(sig==0)=NaN;
                                        else
                                            sig = single(LIMO.cache.fig.mask(channel,:)); sig(sig==0)=NaN;
                                        end
                                    end
                                else
                                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                                        [~, mask] = limo_stat_values(name,p,MCC,LIMO,choice);
                                        sig = single(squeeze(mask(channel,freq_index,:))); sig(sig==0)=NaN;
                                    else
                                        [~, mask] = limo_stat_values(name,p,MCC,LIMO,choice);
                                        sig = single(mask(channel,:)); sig(sig==0)=NaN;
                                    end
                                end
                                h = axis;
                                if strcmpi(LIMO.Analysis,'Frequency')
                                    plot(freqvect,sig.*h(3),'r*','LineWidth',2)
                                else
                                    plot(timevect,sig.*h(3),'r*','LineWidth',2)
                                end
                                break
                            end
                        else
                            index = index+LIMO.design.nb_interactions(i);
                            if i<length(LIMO.design.nb_interactions)
                                index2 = LIMO.design.nb_interactions(i)+LIMO.design.nb_interactions(i+1);
                            end
                        end
                    end
                end
                
                % --
                axis tight; grid on; box on
                title(mytitle,'FontSize',19); drawnow;
                assignin('base','Plotted_data', average)
                set(gca,'FontSize',14,'Layer','Top')
                if strcmpi(LIMO.Analysis,'Frequency')
                    xlabel('Freq in Hz','FontSize',16)
                    ylabel('Power Spectrum in {\mu}V^2/Hz','FontSize',16);
                else
                    xlabel('Time in ms','FontSize',16)
                    ylabel('Amplitude in {\mu}V','FontSize',16)
                end
                
                LIMO.cache.Courseplot.extra     = extra;
                LIMO.cache.Courseplot.average   = average;
                LIMO.cache.Courseplot.channel   = channel;
                LIMO.cache.Courseplot.regressor = regressor;
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    LIMO.cache.Courseplot.frequency = frequency;
                end
                LIMO.cache.Courseplot.ci        = ci;
                LIMO.cache.Courseplot.title     = mytitle;
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
                
            else
                for i=1:size(continuous,1)
                    if i > 1; figure;set(gcf,'Color','w'); end
                    index = find(~isnan(squeeze(continuous(i,1,:))));
                    
                    if strcmpi(LIMO.Analysis,'Frequency')
                        try
                            surf(index,freqvect,squeeze(continuous(i,:,index)));shading interp
                        catch
                            freqvect = linspace(LIMO.data.start,LIMO.data.end,size(continuous,2));
                            surf(index,freqvect,squeeze(continuous(i,:,index)));shading interp
                        end
                        ylabel('Frequency in Hz','FontSize',16)
                        zlabel('Power Spectrum in {\mu}V^2/Hz','FontSize',16)
                    else
                        surf(index,timevect,squeeze(continuous(i,:,index)));shading interp
                        ylabel('Time in ms','FontSize',16)
                        zlabel('Amplitude in {\mu}V','FontSize',16)
                    end
                    % --
                    axis tight; title(mytitle{i},'FontSize',14); drawnow;
                    xlabel('Sorted trials','FontSize',16)
                    try %#ok<TRYNC>
                        set(gca,'XTick',index, 'XTickLabels', reg_values(index));
                    end
                end
                
                LIMO.cache.Courseplot.continuous = continuous;
                LIMO.cache.Courseplot.channel    = channel;
                LIMO.cache.Courseplot.regressor  = regressor;
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    LIMO.cache.Courseplot.frequency = frequency;
                end
                LIMO.cache.Courseplot.title      = mytitle;
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
            end
            
    end % closes switch
    
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % -------------------                               LEVEL 2
    % -------------------                            GROUP EFFECTS
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    
elseif LIMO.Level == 2
    
    if contains(FileName,'LIMO') || contains(FileName,'Y')
        disp('select a statitical result file - plot aborded')
        return
    end
    
    % if previously plotted, recover from the cache
    data_cached = 0;
    if isfield(LIMO,'cache')
        try
            if strcmpi(LIMO.cache.fig.name, FileName) && ...
                    LIMO.cache.fig.MCC == MCC && ...
                    LIMO.cache.fig.threshold == p
                
                disp('using cached data');
                mask = LIMO.cache.fig.mask;
                if isempty(mask)
                    data_cached = 0;
                elseif sum(mask(:)) == 0
                    warndlg('  no values under threshold  ','no significant effect','modal');
                    return
                else
                    toplot      = LIMO.cache.fig.stats;
                    M           = LIMO.cache.fig.pval;
                    mask        = LIMO.cache.fig.mask;
                    mytitle     = LIMO.cache.fig.title;
                    data_cached = 1;
                    assignin('base','stat_values',toplot)
                    assignin('base','p_values',M)
                    assignin('base','mask',mask)
                end
            end
        catch no_cache
            data_cached = 0;
            warning(no_cache,'failed to chache data %s',no_cache.message)
        end
    end
    
    % if there is no cached data, compute and plot
    % -------------------------------------------
    if data_cached == 0
        
        [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO);
        
        if isempty(mask)
            return
        elseif sum(mask(:)) == 0
            warndlg('  no values under threshold  ','no significant effect','modal');
            return
        else
            assignin('base','p_values',squeeze(M))
            assignin('base','mask',squeeze(mask))
        end
        
        if strcmpi(LIMO.Analysis,'Time-Frequency') || strcmpi(LIMO.Analysis,'ITC')
            if contains(FileName,'R2') || ...
                    contains(FileName,'semi_partial')
                toplot = squeeze(toplot(:,:,:,1));
            elseif contains(FileName,'ttest','IgnoreCase',true) || ...
                    contains(FileName,'LI_Map','IgnoreCase',true)
                toplot = squeeze(toplot(:,:,:,4));
            elseif strncmp(FileName,'con_',4)
                toplot = squeeze(toplot(:,:,:,4));
            elseif strncmp(FileName,'ess_',4)
                if ~exist('ess','var')
                    effect_nb = eval(FileName(5:end-4)); %#ok<NASGU>
                end
                toplot = squeeze(toplot(:,:,:,end-1));
            elseif contains(FileName,'Condition') || ...
                    contains(FileName,'Covariate') || ...
                    contains(FileName,'Rep_ANOVA')
                toplot = squeeze(toplot(:,:,:,1));
            else
                disp('file no supported'); return
            end
        else
            if contains(FileName,'R2') || ...
                    contains(FileName,'semi_partial')
                toplot = squeeze(toplot(:,:,1));
            elseif contains(FileName,'ttest','IgnoreCase',true) || ...
                    contains(FileName,'LI_Map','IgnoreCase',true)
                toplot = squeeze(toplot(:,:,4));
            elseif strncmp(FileName,'con_',4)
                toplot = squeeze(toplot(:,:,4));
            elseif strncmp(FileName,'ess_',4)
                 if ~exist('ess','var')
                    effect_nb = eval(FileName(5:end-4)); %#ok<NASGU>
                end
                toplot = squeeze(toplot(:,:,end-1));
            elseif contains(FileName,'Condition') || ...
                    contains(FileName,'Covariate') || ...
                    contains(FileName,'Rep_ANOVA')
                toplot = squeeze(toplot(:,:,1));
            else
                disp('file no supported'); return
            end
        end
        assignin('base','stat_values',toplot)
        data_cached = 0;
    end
    
    % ------------------------------
    %      Image and topoplot
    % ----------------------------
    if Type == 1 || Type == 2
        
        % cache the results for next time
        % ------------------------------
        if data_cached == 0
            LIMO.cache.fig.name       = FileName;
            LIMO.cache.fig.MCC        = MCC;
            LIMO.cache.fig.stats      = toplot;
            LIMO.cache.fig.threshold  = p;
            LIMO.cache.fig.pval       = squeeze(M);
            LIMO.cache.fig.mask       = squeeze(mask);
            LIMO.cache.fig.title      = mytitle;
            if exist(LIMO.dir,'dir')
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
            else
                disp('Cached data in LIMO.mat cannot be updated - LIMO dir doesn''t exist (likely moved files)')
            end
        end
        
        % image all results
        % ------------------
        if Type == 1 && ~strcmpi(LIMO.Analysis,'Time-Frequency') && ~strcmpi(LIMO.Analysis,'ITC')
            % Create 3D surface plot for Level 2
            if Type == 1 && ~strcmpi(LIMO.Analysis,'Time-Frequency') && ~strcmpi(LIMO.Analysis,'ITC')
                % Create 3D figure window
                main_fig = figure('Color','w','Name',[mytitle ' - 3D View']);
                
                % Size figure to match MATLAB desktop
                try
                    desktop_pos = get(0, 'ScreenSize');  % [left, bottom, width, height] of screen
                    % Use 90% of screen size with some margin
                    fig_width = desktop_pos(3) * 0.9;
                    fig_height = desktop_pos(4) * 0.8;
                    fig_left = desktop_pos(3) * 0.05;  % 5% margin from left
                    fig_bottom = desktop_pos(4) * 0.1;  % 10% margin from bottom
                    set(main_fig, 'Position', [fig_left, fig_bottom, fig_width, fig_height]);
                catch
                    % Fallback to default if sizing fails
                    set(main_fig, 'Position', [100, 100, 1200, 800]);
                end
                
                % Prepare data
                if strcmpi(LIMO.Analysis,'Time')
                    if isfield(LIMO.data,'timevect')
                        time_vect = LIMO.data.timevect;
                    else
                        time_vect = linspace(LIMO.data.start, LIMO.data.end, size(toplot,2));
                    end
                    time_label = 'Time (ms)';
                else % Frequency
                    if isfield(LIMO.data,'freqlist')
                        time_vect = LIMO.data.freqlist;
                    else
                        time_vect = linspace(LIMO.data.start, LIMO.data.end, size(toplot,2));
                    end
                    time_label = 'Frequency (Hz)';
                end
                
                % Create channel vector
                num_channels = size(toplot,1);
                channel_vector = 1:num_channels;
                
                % Create meshgrid for 3D plotting
                [T, C] = meshgrid(time_vect, channel_vector);
                
                % Create axes in the 3D figure
                axes_3d = axes('Parent', main_fig);
                
                % Plot 3D surface
                surf(axes_3d, C, T, toplot, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
                
                % Customize the plot
                % Set Z-axis limits to encompass full range of statistical values
                z_min = min(toplot(:));
                z_max = max(toplot(:));
                z_range = z_max - z_min;
                z_margin = z_range * 0.05; % Add 5% margin
                zlim(axes_3d, [z_min - z_margin, z_max + z_margin]);
                set(axes_3d, 'Color', [1 1 1]); % White background for the axes
                xlabel(axes_3d, 'Channels', 'FontSize', 14, 'Color', 'k');
                ylabel(axes_3d, time_label, 'FontSize', 14, 'Color', 'k');
                zlabel(axes_3d, 'Statistical values', 'FontSize', 14, 'Color', 'k');
                title_handle = title(axes_3d, [mytitle ' - 3D View'], 'FontSize', 16);
                set(title_handle, 'Color', 'k');
                set(axes_3d, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k'); % Black axes lines and ticks

                % Set colormap
                colormap(axes_3d, limo_color_images(toplot));
                cb = colorbar(axes_3d);
                set(cb, 'Color', 'k'); % Make colorbar text black
                
                % Adjust view angle
                view(axes_3d, [-45 30]);
                
                % Set x-axis ticks and labels for channels
                if num_channels <= 20
                    set(axes_3d, 'XTick', channel_vector);
                    % Get channel labels
                    if isempty(LIMO.design.electrode)
                        label_electrodes = cell(length(LIMO.data.chanlocs), 1);
                        for i = 1:length(LIMO.data.chanlocs)
                            label_electrodes{i} = LIMO.data.chanlocs(i).labels;
                        end
                        set(axes_3d, 'XTickLabel', label_electrodes);
                    else
                        if length(LIMO.design.electrode) == 1
                            set(axes_3d, 'XTickLabel', {num2str(LIMO.design.electrode)});
                        end
                    end
                else
                    % For many channels, show fewer labels
                    tick_indices = round(linspace(1, num_channels, min(10, num_channels)));
                    set(axes_3d, 'XTick', tick_indices);
                    if isempty(LIMO.design.electrode)
                        tick_labels = cell(1, length(tick_indices));
                        for i_tick = 1:length(tick_indices)
                            tick_labels{i_tick} = LIMO.data.chanlocs(tick_indices(i_tick)).labels;
                        end
                        set(axes_3d, 'XTickLabel', tick_labels);
                    end
                end
                
                % Add grid and lighting
                grid(axes_3d, 'on');
                set(axes_3d, 'GridColor', 'k', 'GridAlpha', 0.5); % Black grid
                light('Position', [-1 -1 2], 'Style', 'local');
                lighting gouraud;
                material dull;
                
                % Apply cluster visualization if available
                if exist('mask','var') && ~isempty(mask) && MCC == 2 && max(mask(:)) > 1
                    hold(axes_3d, 'on');
                    % Create cluster-specific visualization  
                    n_cluster = max(mask(:));
                    % Use distinguishable colors, avoiding light colors
                    cluster_colors = [
                        0.0000 0.4470 0.7410;  % Blue
                        0.8500 0.3250 0.0980;  % Orange  
                        0.9290 0.6940 0.1250;  % Yellow
                        0.4940 0.1840 0.5560;  % Purple
                        0.4660 0.6740 0.1880;  % Green
                        0.3010 0.7450 0.9330;  % Cyan
                        0.6350 0.0780 0.1840   % Dark Red
                    ];
                    
                    h_clusters = gobjects(n_cluster, 1);
                    legend_handles = [];
                    legend_entries = {};
                    
                    for cluster_id = 1:n_cluster
                        cluster_mask = (mask == cluster_id);
                        if sum(cluster_mask(:)) > 0
                            cluster_data = nan(size(toplot));
                            cluster_data(cluster_mask) = toplot(cluster_mask);
                            
                            % Create surface for this cluster
                            if cluster_id <= size(cluster_colors, 1)
                                color_idx = cluster_id;
                            else
                                color_idx = mod(cluster_id - 1, size(cluster_colors, 1)) + 1;
                            end
                            
                            h_cluster = surf(axes_3d, C, T, cluster_data, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                            set(h_cluster, 'FaceColor', cluster_colors(color_idx, :));
                            h_clusters(cluster_id) = h_cluster;
                            
                            % Add to legend
                            legend_handles(end+1) = h_cluster;
                            legend_entries{end+1} = sprintf('Cluster %d', cluster_id);
                        end
                    end
                    
                end
                
                % Add significance threshold planes for reference
                threshold_handles = [];
                threshold_entries = {};
                
                if p > 0 && p <= 1  % Valid p-value threshold
                    hold(axes_3d, 'on');
                    % Compute statistical threshold
                    % Try to get degrees of freedom from LIMO structure
                    if isfield(LIMO, 'model') && isfield(LIMO.model, 'model_df')
                        df = LIMO.model.model_df(2);  % Error degrees of freedom
                    elseif isfield(LIMO, 'design') && isfield(LIMO.design, 'nb_subjects')
                        df = LIMO.design.nb_subjects - 2;  % Simple approximation
                    else
                        df = 30;  % Default fallback
                    end
                    
                    % Determine test type and compute threshold
                    % Check for F-test/ANOVA first (takes precedence)
                    is_ftest = contains(lower(FileName), 'f') || contains(lower(mytitle), 'f') || ...
                               contains(lower(FileName), 'anova') || contains(lower(mytitle), 'anova') || ...
                               contains(lower(FileName), 'interaction') || contains(lower(mytitle), 'interaction');
                    is_ttest = ~is_ftest && (contains(lower(FileName), 't') || contains(lower(mytitle), 't'));
                    
                    if is_ftest
                        % For F-statistics or ANOVA
                        stat_threshold = finv(1 - p, 1, df);  % F-test
                    elseif is_ttest
                        stat_threshold = tinv(1 - p/2, df);  % Two-tailed t-test
                    else
                        % Default to F-test if uncertain
                        stat_threshold = finv(1 - p, 1, df);  % F-test
                        is_ftest = true;
                    end
                    
                    % Create threshold planes
                    [C_thresh, T_thresh] = meshgrid(channel_vector, time_vect);
                    Z_thresh_pos = ones(size(C_thresh)) * stat_threshold;
                    
                    % Plot positive threshold plane (red, top)
                    h_pos_thresh = surf(axes_3d, C_thresh, T_thresh, Z_thresh_pos, 'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                    threshold_handles(end+1) = h_pos_thresh;
                    
                    if is_ttest && ~is_ftest
                        threshold_entries{end+1} = sprintf('Significance threshold (±%.2f)', stat_threshold);
                        % Plot negative threshold plane (blue, bottom) for t-tests only
                        Z_thresh_neg = ones(size(C_thresh)) * (-stat_threshold);
                        h_neg_thresh = surf(axes_3d, C_thresh, T_thresh, Z_thresh_neg, 'FaceColor', 'blue', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                        threshold_handles(end+1) = h_neg_thresh;
                        threshold_entries{end+1} = 'Negative threshold';
                    else
                        % For F-tests/ANOVA, only positive threshold makes sense
                        threshold_entries{end+1} = sprintf('Significance threshold (%.2f)', stat_threshold);
                    end
                end
                
                % Combine legend entries
                if exist('legend_handles', 'var') && ~isempty(legend_handles) && ~isempty(threshold_handles)
                    % Combine cluster and threshold legend entries
                    all_handles = [legend_handles, threshold_handles];
                    all_entries = [legend_entries, threshold_entries];
                    h_legend = legend(axes_3d, all_handles, all_entries, 'Location', 'best');
                    set(h_legend, 'Color', 'white', 'TextColor', 'black', 'EdgeColor', 'black');
                elseif exist('legend_handles', 'var') && ~isempty(legend_handles)
                    % Only cluster legend
                    h_legend = legend(axes_3d, legend_handles, legend_entries, 'Location', 'best');
                    set(h_legend, 'Color', 'white', 'TextColor', 'black', 'EdgeColor', 'black');
                elseif ~isempty(threshold_handles)
                    % Only threshold legend
                    h_legend = legend(axes_3d, threshold_handles, threshold_entries, 'Location', 'best');
                    set(h_legend, 'Color', 'white', 'TextColor', 'black', 'EdgeColor', 'black');
                end
                
                % Store the 3D plot data
                assignin('base', 'Plot3D_X', C);
                assignin('base', 'Plot3D_Y', T);
                assignin('base', 'Plot3D_Z', toplot);
                
                % Use the 2D view tab - but call limo_display_image WITHOUT parent_tab to avoid topoplot issues
                limo_display_image(LIMO,toplot,mask,mytitle,flag)
            else
                % If no 3D plot created, use normal display
                limo_display_image(LIMO,toplot,mask,mytitle,flag)
            end
            
        elseif Type == 1 && strcmpi(LIMO.Analysis,'Time-Frequency') || ...
                Type == 1 && strcmpi(LIMO.Analysis,'ITC')
            if ndims(toplot)==3
                limo_display_image_tf(LIMO,toplot,mask,mytitle,flag);
            else
                limo_display_image(LIMO,squeeze(toplot),squeeze(mask),mytitle,flag)
            end
            
            
        elseif Type == 2
            %--------------------------
            % topoplot
            %--------------------------
            
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                errordlg('topoplot not supported for time-frequency analyses')
            else
                if isfield(LIMO.design,'channel')  % not full scalp
                    if ~isempty(LIMO.design.electrode)
                        msgbox('Only one channel found','No topoplot')
                        return
                    end
                end
            end
            
            if sum(mask(:)) == 0
                warndlg('no values under threshold','no significant effect');
            else
                EEG.data     = toplot;
                EEG.setname  = mytitle;
                EEG.nbchan   = size(EEG.data,1);
                EEG.pnts     = size(EEG.data,2);
                EEG.trials   = 1;
                EEG.chanlocs = LIMO.data.chanlocs;
                EEG.xmin     = LIMO.data.start/1000;
                EEG.xmax     = LIMO.data.end/1000;
                EEG.times    = linspace(EEG.xmin,EEG.xmax,EEG.pnts);
                pop_topoplot(EEG);
                assignin('base','Plotted_data',EEG.data)
            end
            
    elseif Type == 3
        
        %--------------------------
        % Course plot
        %--------------------------       
        
        if contains(FileName,'one_sample','IgnoreCase',true) || contains(FileName,'two_samples','IgnoreCase',true) || ...
                contains(FileName,'paired_samples','IgnoreCase',true) || contains(FileName,'con_','IgnoreCase',true) || ...
                contains(FileName,'ess_','IgnoreCase',true)
            % ------------------------------------------------------------------------------------------------------------
            % stat file dim = (electrodes, frames, [mean value, se, df, t, p])
            % H0 file dim = (electrodes,frames,[t, p],nboot)
            
            data = load(fullfile(PathName,FileName));
            data = data.(cell2mat(fieldnames(data)));
                
                % Customize the plot
                % Set Z-axis limits to encompass full range of statistical values
                z_min = min(toplot(:));
                z_max = max(toplot(:));
                z_range = z_max - z_min;
                z_margin = z_range * 0.05; % Add 5% margin
                zlim(axes_3d, [z_min - z_margin, z_max + z_margin]);
                set(axes_3d, 'Color', [1 1 1]); % White background for the axes
                xlabel(axes_3d, 'Channels', 'FontSize', 14, 'Color', 'k');
                ylabel(axes_3d, time_label, 'FontSize', 14, 'Color', 'k');
                zlabel(axes_3d, 'Statistical values', 'FontSize', 14, 'Color', 'k');
                title_handle = title(axes_3d, [mytitle ' - 3D View'], 'FontSize', 16);
                set(title_handle, 'Color', 'k');
                set(axes_3d, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k'); % Black axes lines and ticks

                % Set colormap
                colormap(axes_3d, limo_color_images(toplot));
                cb = colorbar(axes_3d);
                set(cb, 'Color', 'k'); % Make colorbar text black
                
                % Adjust view angle
                view(axes_3d, [-45 30]);
                
                % Set x-axis ticks and labels for channels
                if num_channels <= 20
                    set(axes_3d, 'XTick', channel_vector);
                    % Get channel labels
                    if isempty(LIMO.design.electrode)
                        label_electrodes = cell(length(LIMO.data.chanlocs), 1);
                        for i = 1:length(LIMO.data.chanlocs)
                            label_electrodes{i} = LIMO.data.chanlocs(i).labels;
                        end
                        set(axes_3d, 'XTickLabel', label_electrodes);
                    else
                        if length(LIMO.design.electrode) == 1
                            set(axes_3d, 'XTickLabel', {num2str(LIMO.design.electrode)});
                        end
                    end
                else
                    % For many channels, show fewer labels
                    tick_indices = round(linspace(1, num_channels, min(10, num_channels)));
                    set(axes_3d, 'XTick', tick_indices);
                    if isempty(LIMO.design.electrode)
                        tick_labels = cell(1, length(tick_indices));
                        for i_tick = 1:length(tick_indices)
                            tick_labels{i_tick} = LIMO.data.chanlocs(tick_indices(i_tick)).labels;
                        end
                        set(axes_3d, 'XTickLabel', tick_labels);
                    end
                end
                
                % Add grid and lighting
                grid(axes_3d, 'on');
                set(axes_3d, 'GridColor', 'k', 'GridAlpha', 0.5); % Black grid
                light('Position', [-1 -1 2], 'Style', 'local');
                lighting gouraud;
                material dull;
                
                % Apply mask if available - make non-significant values transparent
                if exist('mask','var') && ~isempty(mask)
                    masked_data = toplot;
                    masked_data(~mask) = NaN;
                    hold on;
                    
                    % Check if we have cluster correction (mask contains cluster labels > 1)
                    if MCC == 2 && max(mask(:)) > 1  % Cluster correction with multiple clusters
                        % Create cluster-specific visualization
                        n_cluster = max(mask(:));
                        % Use distinguishable colors, avoiding light colors
                        cluster_colors = [
                            0.0000 0.4470 0.7410;  % Blue
                            0.8500 0.3250 0.0980;  % Orange  
                            0.9290 0.6940 0.1250;  % Yellow
                            0.4940 0.1840 0.5560;  % Purple
                            0.4660 0.6740 0.1880;  % Green
                            0.3010 0.7450 0.9330;  % Cyan
                            0.6350 0.0780 0.1840   % Dark Red
                        ];
                        % Repeat colors if more clusters than predefined colors
                        if n_cluster > size(cluster_colors, 1)
                            cluster_colors = repmat(cluster_colors, ceil(n_cluster/size(cluster_colors,1)), 1);
                        end
                        
                        % Store handles for legend
                        h_clusters = [];
                        
                        for cluster_id = 1:n_cluster
                            cluster_data = toplot;
                            cluster_data(mask ~= cluster_id) = NaN;
                            
                            % Plot each cluster with a different color
                            h_cluster = surf(C, T, cluster_data, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                            set(h_cluster, 'FaceColor', cluster_colors(cluster_id, :));
                            h_clusters(cluster_id) = h_cluster;
                        end
                        
                        % Add 3D view title
                        title_handle = title([mytitle ' - 3D View'], 'FontSize', 16);
                        
                        % Add threshold plane for reference
                        if p > 0 && p <= 1  % Valid p-value threshold
                            % Compute statistical threshold
                            % Try to get degrees of freedom from LIMO structure
                            if isfield(LIMO, 'model') && isfield(LIMO.model, 'model_df')
                                df_vals = LIMO.model.model_df;
                                if length(df_vals) >= 2
                                    stat_threshold = finv(1-p, df_vals(1), df_vals(2));  % F-distribution
                                    is_ttest = false;
                                else
                                    stat_threshold = tinv(1-p/2, df_vals(1));       % t-distribution (two-tailed)
                                    is_ttest = true;
                                end
                            else
                                % Use a default threshold and detect t-test from title
                                stat_threshold = -log10(p) * 2; % Rough approximation
                                is_ttest = contains(lower(mytitle), {'ttest', 't values', 'paired', 'one sample'});
                            end
                            
                            % Create threshold plane as solid surface
                            [C_thresh, T_thresh] = meshgrid(channel_vector, time_vect);
                            threshold_handles = [];
                            threshold_labels = {};
                            
                            % Always add positive threshold
                            Z_thresh_pos = ones(size(C_thresh)) * stat_threshold;
                            h_thresh_pos = surf(C_thresh, T_thresh, Z_thresh_pos, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'FaceColor', 'red');
                            threshold_handles(end+1) = h_thresh_pos;
                            
                            % Only add negative threshold for explicit t-tests with significant negative clusters
                            has_negative_clusters = is_ttest && any(toplot(:) < -stat_threshold) && any(mask(:) > 0 & toplot(:) < 0);
                            if has_negative_clusters
                                Z_thresh_neg = ones(size(C_thresh)) * (-stat_threshold);
                                h_thresh_neg = surf(C_thresh, T_thresh, Z_thresh_neg, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'FaceColor', 'blue');
                                threshold_handles(end+1) = h_thresh_neg;
                                threshold_labels = {sprintf('+Threshold (p=%.3f)', p), sprintf('-Threshold (p=%.3f)', p)};
                            else
                                threshold_labels = {sprintf('Threshold (p=%.3f)', p)};
                            end
                            
                            % Create legend with proper handles and light mode styling
                            legend_entries = cell(n_cluster + length(threshold_labels), 1);
                            legend_handles = [h_clusters, threshold_handles];
                            for i = 1:n_cluster
                                legend_entries{i} = sprintf('Cluster %d', i);
                            end
                            for i = 1:length(threshold_labels)
                                legend_entries{n_cluster + i} = threshold_labels{i};
                            end
                            h_legend = legend(legend_handles, legend_entries, 'Location', 'best');
                            set(h_legend, 'Color', 'white', 'TextColor', 'black', 'EdgeColor', 'black');
                        end
                    else
                        % Standard binary mask overlay
                        surf(C, T, masked_data, 'EdgeColor', 'none', 'FaceAlpha', 1);
                        
                        % Add threshold plane for reference (even for binary mask)
                        if p > 0 && p <= 1  % Valid p-value threshold
                            % Compute statistical threshold based on degrees of freedom
                            if exist('df','var') && ~isempty(df)
                                stat_threshold = finv(1-p, df(1), df(2));
                            else
                                % Use a default threshold based on p-value (approximate)
                                stat_threshold = -log10(p) * 2; % Rough approximation
                            end
                            
                            % Create threshold plane
                            [C_thresh, T_thresh] = meshgrid(channel_vector, time_vect);
                            Z_thresh = ones(size(C_thresh)) * stat_threshold;
                            mesh(C_thresh, T_thresh, Z_thresh, 'FaceAlpha', 0.2, 'EdgeColor', 'red', 'LineWidth', 1);
                        end
                    end
                end
                
                % Store the 3D plot data
                assignin('base', 'Plot3D_X', C);
                assignin('base', 'Plot3D_Y', T);
                assignin('base', 'Plot3D_Z', toplot);
                
                % Use the 2D view tab that was created earlier
                limo_display_image(LIMO,toplot,mask,mytitle,flag,tab_2d)
            else
                % If no 3D plot created, use normal display
                limo_display_image(LIMO,toplot,mask,mytitle,flag)
            end
            
        elseif Type == 1 && strcmpi(LIMO.Analysis,'Time-Frequency') || ...
                Type == 1 && strcmpi(LIMO.Analysis,'ITC')
            if ndims(toplot)==3
                limo_display_image_tf(LIMO,toplot,mask,mytitle,flag);
            else
                limo_display_image(LIMO,squeeze(toplot),squeeze(mask),mytitle,flag)
            end
            
        elseif Type == 2
            %--------------------------
            % topoplot
            %--------------------------
            
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                errordlg('topoplot not supported for time-frequency analyses')
            else
                if isfield(LIMO.design,'channel')  % not full scalp
                    if ~isempty(LIMO.design.electrode)
                        msgbox('Only one channel found','No topoplot')
                        return
                    end
                end
            end
            
            if sum(mask(:)) == 0
                warndlg('no values under threshold','no significant effect');
            else
                EEG.data     = toplot;
                EEG.setname  = mytitle;
                EEG.chanlocs = LIMO.data.chanlocs;
                
                if size(toplot,2) == 1
                    opt = {'maplimits','maxmin','verbose','off'};
                    if isfield(LIMO,'Type')
                        if strcmpi(LIMO.Type,'Components')
                            opt = {'maplimits','absmax','electrodes','off','verbose','off'};
                        end
                    end
                    figure; set(gcf,'Color','w','InvertHardCopy','off');
                    topoplot(toplot(:,1),EEG.chanlocs,opt{:});
                    title('Topoplot','FontSize',12)
                else
                    if strcmpi(LIMO.Analysis,'Time')
                        EEG.xmin  = LIMO.data.start/1000; % in sec
                        EEG.xmax  = LIMO.data.end/1000;   % in sec
                        EEG.times = (LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000); % in sec;
                        call_topolot(EEG,FileName,LIMO.Analysis)
                        % pop_topoplot(EEG);
                    elseif strcmpi(LIMO.Analysis,'Frequency')
                        EEG.xmin  = LIMO.data.freqlist(1);
                        EEG.xmax  = LIMO.data.freqlist(end);
                        freqlist  = inputdlg('specify frequency range e.g. [5:2:40]','Choose Frequencies to plot');
                        if isempty(freqlist)
                            return
                        else
                            EEG.freq = str2double(cell2mat(freqlist));
                            if min(EEG.freq)<EEG.xmin || max(EEG.freq)>EEG.xmax
                                errordlg('selected frequency out of bound'); return
                            end
                        end
                        call_topolot(EEG,FileName,LIMO.Analysis)
                    end
                    assignin('base','Plotted_data',EEG.data)
                end
            end
        end
        
    elseif Type == 3
        
        %--------------------------
        % Course plot
        %--------------------------       
        
        if contains(FileName,'one_sample','IgnoreCase',true) || contains(FileName,'two_samples','IgnoreCase',true) || ...
                contains(FileName,'paired_samples','IgnoreCase',true) || contains(FileName,'con_','IgnoreCase',true) || ...
                contains(FileName,'ess_','IgnoreCase',true)
            % ------------------------------------------------------------------------------------------------------------
            % stat file dim = (electrodes, frames, [mean value, se, df, t, p])
            % H0 file dim = (electrodes,frames,[t, p],nboot)
            
            data = load(fullfile(PathName,FileName));
            data = data.(cell2mat(fieldnames(data)));
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                [~,channel,freq,time] = limo_display_reducedim(data(:,:,:,[4 5]),LIMO,g.channels,g.restrict,g.dimvalue);
                data                  = squeeze(data(channel,freq,time,:,:)); % 2D
                sig                   = squeeze(single(mask(channel,freq,time))); %1D
            else
                [~,channel,freq,time] = limo_display_reducedim(data(:,:,[4 5]),LIMO,g.channels);
                data                  = squeeze(data(channel,time,:));
                sig                   = single(mask(channel,:));
            end
            sig(sig==0)=NaN;
            
            % compute
            trimci      = NaN(size(data,1),3);
            trimci(:,2) = data(:,1); % mean values
            if contains(FileName,'ess','IgnoreCase',true)
                start_at = max(strfind(FileName,'_'))+1;
                C = LIMO.contrast{eval(FileName(start_at:end-4))}.C;
                df = rank(C); % rank of the relevant contrast
                trimci(:,1) = squeeze(trimci(:,2))-(finv(1-p./2*size(C,1),df,data(:,3)).*data(:,2));
                trimci(:,3) = squeeze(trimci(:,2))+(finv(1-p./2*size(C,1),df,data(:,3)).*data(:,2));
            else
                trimci(:,1) = squeeze(trimci(:,2))-(tinv(1-p./2,data(:,3)).*data(:,2));
                trimci(:,3) = squeeze(trimci(:,2))+(tinv(1-p./2,data(:,3)).*data(:,2));
            end
            
            % plot
            if strcmpi(LIMO.Analysis,'Time')
                if isfield(LIMO.data,'timevect')
                    xvect = LIMO.data.timevect;
                else
                    xvect = [];
                end
                
                if size(xvect,2) ~= size(toplot,2)
                    xvect              = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
                    LIMO.data.timevect = xvect;
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                end
                
            elseif strcmpi(LIMO.Analysis,'Frequency')
                if isfield(LIMO.data,'freqlist')
                    xvect=LIMO.data.freqlist;
                else
                    xvect = [];
                end
                
                if size(xvect,2) ~= size(toplot,2)
                    xvect              = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
                    LIMO.data.freqlist = xvect;
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                end
                
            elseif strcmpi(LIMO.Analysis,'Time-Frequency')
                if length(time) > 1 && isfield(LIMO.data,'tf_times')
                    xvect = LIMO.data.tf_times;
                elseif length(freq) > 1 && isfield(LIMO.data,'tf_freqs')
                    xvect = LIMO.data.tf_freqs;
                else
                    xvect = [];
                end
                
                if length(time) > 1&& size(xvect,2) ~= size(data,1)
                    xvect              = linspace(LIMO.data.start,LIMO.data.end,size(data,1));
                    LIMO.data.tf_times =  xvect;
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                elseif length(freq) > 1 && size(xvect,2) ~= size(data,1)
                    xvect              = linspace(LIMO.data.lowf,LIMO.data.highf,size(data,1));
                    LIMO.data.tf_freqs =  xvect;
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                end
            end
            
            figure;
            set(gcf,'Color','w')
            plot(xvect,squeeze(trimci(:,2)),'LineWidth',3);
            fillhandle = patch([xvect,fliplr(xvect)], [trimci(:,1)' fliplr(trimci(:,3)')], [1 0 0]);
            set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);% set edge color
            grid on; box on; axis tight
            h = axis;  hold on;
            plot(xvect,(sig./10+1).*h(3),'k.','MarkerSize',20)
            if strcmpi(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(time) > 1
                xlabel('Time in ms','FontSize',14)
                ylabel('Amplitude (A.U.)','FontSize',14)
            elseif strcmpi(LIMO.Analysis,'Frequency') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(freq) > 1
                xlabel('Frequency in Hz','FontSize',14)
                ylabel('Spectral Power (A.U.)','FontSize',14)
            end
            if isempty(LIMO.design.electrode)
                title(sprintf('%s \n%s %s %s (%g)',mytitle,'Mean values',LIMO.Type(1:end-1),LIMO.data.chanlocs(channel).labels,channel),'FontSize',16); drawnow;
            else
                title(sprintf('%s \n%s virtual %s',mytitle,'Mean values',LIMO.Type(1:end-1)),'FontSize',16); drawnow;
            end
            assignin('base','Plotted_data',trimci);
            
            
        elseif contains(LIMO.design.name,'regression','IgnoreCase',true) && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true) || ...
                contains(LIMO.design.name,'ANOVA') && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true) || ...
                contains(LIMO.design.name,'ANCOVA') && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true)
            % --------------------------------------------------------------------------------
            
            % which variable(s) to plot
            % ----------------------
            if size(LIMO.design.X,2) >= 2
                if contains(FileName,'Condition_effect_') || ...
                        contains(FileName,'Covariate_effect_')
                    regressor = eval(FileName(18:end-4));
                    if contains(FileName,'Covariate_effect_')
                        regressor = regressor+sum(LIMO.design.nb_conditions);
                    end
                else
                    input_title = sprintf('which regressor to plot?: 1 to %g ',size(LIMO.design.X,2)-1);
                    regressor = inputdlg(input_title,'Plotting option');
                end
                
                if isempty(regressor)
                    warning on
                    warning('couldn''t figure out the regressor number/column, plot aborded')
                    return
                end
                
                try
                    if iscell(regressor)
                        regressor = sort(eval(cell2mat(regressor)));
                    end
                    if max(regressor) > size(LIMO.design.X,2)
                        errordlg('invalid regressor number');
                    end
                catch reginput_error
                    fprintf('error with regressor numbers/columns line 1373:\n %s',reginput_error.message)
                    return
                end
            else
                regressor = 1;
            end
            
            categorical = sum(LIMO.design.nb_conditions) + sum(LIMO.design.nb_interactions);
            if max(regressor) == size(LIMO.design.X,2)
                tmp = regressor(1:end-1);
            else
                tmp = regressor;
            end
            
            if sum(tmp<=categorical) >=1 && sum(tmp>categorical) >=1
                errordlg('you can''t plot categorical and continuous regressors together'); return
            end
            
            % load the effect
            % --------------
            data = load(FileName);
            data = data.(cell2mat(fieldnames(data)));
            if numel(size(data)) == 3 && size(data,2) == 1
                errordlg2('single time point detected, plot aborded'); return
            elseif numel(size(data)) == 4 && size(data,3) == 1
                errordlg2('single time point detected, plot aborded'); return
            end
            
            % which course plot to make
            % -------------------------
            if isempty(g.plot3type)
                extra = questdlg('Plotting data','Options','Original data','Modelled data','Adjusted data','Modelled data');
            else
                extra = g.plot3type;
                % allow typos
                if contains(extra,'orig','Ignorecase',true)
                    extra = 'Original';
                elseif contains(extra,'Model','Ignorecase',true)
                    extra = 'Modelled';
                elseif contains(extra,'Adj','Ignorecase',true)
                    extra = 'Adjusted';
                else
                    if exist(errodlg2,'file')
                        errordlg2(sprintf('input option ''%s'' invalid',extra)); return
                    else
                        errordlg(sprintf('input option ''%s'' invalid',extra)); return
                    end
                end
            end
            
            if isempty(extra)
                return
            elseif strcmpi(extra,'Original data')
                if regressor == size(LIMO.design.X,2)
                    errordlg('you can''t plot adjusted mean for original data'); return
                end
            end
            
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                [~,channel,freq,time] = limo_display_reducedim(data,LIMO,g.channels,g.restrict,g.dimvalue);
                if length(freq) == 1; g.restrict = 'time';
                else; g.restrict = 'frequency'; end
                sig         = squeeze(single(mask(channel,freq,time))); %1D
            else
                [~,channel] = limo_display_reducedim(data,LIMO,g.channels);
                sig         = single(mask(channel,:));
            end
            sig(sig==0)=NaN;
            clear data
            
            % down to business
            % ----------------------
            probs = [p/2; 1-p/2];
            z     = norminv(probs);
            Yr    = load(fullfile(LIMO.dir,'Yr.mat'));
            Yr    = Yr.Yr;
            
            if contains(extra,'Original','Ignorecase',true) 
                if regressor <= length(LIMO.design.nb_conditions) && ...
                        LIMO.design.nb_conditions ~= 0 % for categorical variables
                    if length(LIMO.design.nb_conditions) == 1
                        start = 1;
                    else
                        start = sum(LIMO.design.nb_conditions(1:regressor-1));
                    end
                    
                    for i = (start+LIMO.design.nb_conditions(regressor)-1):-1:start
                        index{i} = find(LIMO.design.X(:,i));
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            data = squeeze(Yr(channel,freq,:,index{i}));
                        else
                            data = squeeze(Yr(channel,:,index{i}));
                        end
                        average(i,:) = nanmean(data,2);
                        se           = (nanstd(data,0,2) ./ sqrt(numel(index{i})));
                        ci(i,:,:)    = repmat(average(i,:),2,1) + repmat(se',2,1).*repmat(z,1,size(Yr,2));
                    end
                    clear mytitle
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Original subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                    else
                        mytitle = sprintf('Original subjects'' parameters at optimized channel');
                    end
                else % continuous variable
                    for i=max(regressor):-1:min(regressor)
                        index{i}           = find(LIMO.design.X(:,i));
                        [~,sorting_values] = sort(LIMO.design.X(index{i},i));  % continuous variable 3D plot
                        reg_values(i,:)    = LIMO.design.X(sorting_values,i);
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            if strcmpi(g.restrict,'time')
                                continuous(i,:,:) = Yr(channel,freq,:,sorting_values);
                            else
                                continuous(i,:,:) = Yr(channel,:,time,sorting_values);
                            end
                        else
                            continuous(i,:,:) = Yr(channel,:,sorting_values);
                        end
                        clear mytitle
                        if isempty(LIMO.design.electrode)
                            mytitle{i} = sprintf('Original subjects'' parameters \n sorted by regressor %g channel %s (%g)', i, LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle{i} = sprintf('Original subjects'' parameters, \n sorted by regressor %g at optimized channel', i);
                        end
                    end
                    remove                 = find(sum(reg_values == 0,2) == size(continuous,3));
                    continuous(remove,:,:) = [];
                    reg_values(remove,:)   = [];
                end
                
            elseif contains(extra,{'Modelled','Modeled'},'Ignorecase',true)
                if exist('Betas.mat','file') % OLS & IRLS GLM
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        Betas = load('Betas.mat');
                        if strcmpi(g.restrict,'time')
                            Betas = squeeze(Betas.Betas(channel,freq,:,:));
                        else
                            Betas = squeeze(Betas.Betas(channel,:,time,:));
                        end
                        R     = eye(size(Yr,4)) - (LIMO.design.X*pinv(LIMO.design.X));
                    else
                        Betas = load('Betas.mat');
                        Betas = squeeze(Betas.Betas(channel,:,:));
                        R     = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
                    end
                    Yh = (LIMO.design.X*Betas')'; % modelled data
                else % strcmpi(LIMO.design.method,'Generalized Welch's method')
                    [~,~,Yh,~,dfe] = limo_robust_1way_anova(squeeze(Yr(channel,:,:)),LIMO.design.X);
                    Res            = squeeze(Yr(channel,:,:))-Yh;
                end
                
                if regressor <= length(LIMO.design.nb_conditions) && ...
                        LIMO.design.nb_conditions ~= 0 % for categorical variables
                    if length(LIMO.design.nb_conditions) == 1
                        start = 1;
                    else
                        start = sum(LIMO.design.nb_conditions(1:regressor-1));
                    end
                    
                    for i = (start+LIMO.design.nb_conditions(regressor)-1):-1:start
                        index{i}     = find(LIMO.design.X(:,i));
                        data         = squeeze(Yh(:,index{i}));
                        average(i,:) = nanmean(data,2);
                        index{i}     = index{i}(find(~isnan(squeeze(Yr(channel,1,index{i}))))); %#ok<FNDSB>
                        if exist('R','var')
                            var      = diag(((R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')'*(R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')) / LIMO.model.model_df(2));
                        else
                            var      = diag(Res*Res')./dfe;
                        end
                        CI           = sqrt(var/size(index{i},1))*z';
                        ci(i,:,:)    = (repmat(nanmean(data,2),1,2)+CI)';
                    end
                    clear mytitle
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Modelled subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                    else
                        mytitle = sprintf('Modelled subjects'' parameters at optimized channel');
                    end
                    clear Yr
                else % continuous variable
                    for i=max(regressor):-1:min(regressor)
                        index{i}           = find(LIMO.design.X(:,i));
                        [~,sorting_values] = sort(LIMO.design.X(index{i},i));  % continuous variable 3D plot
                        reg_values(i,:)    = LIMO.design.X(sorting_values,i);
                        continuous(i,:,:)  = Yh(:,sorting_values);
                        clear mytitle
                        if isempty(LIMO.design.electrode)
                            mytitle = sprintf('Modelled subjects'' parameters \n sorted by regressor %g channel %s (%g)', ...
                                i, LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle = sprintf('Modelled subjects'' parameters \n sorted by regressor %g at optimized channel', i);
                        end
                    end
                    remove                 = find(sum(reg_values == 0,2) == size(continuous,3));
                    continuous(remove,:,:) = [];
                    reg_values(remove,:)   = [];
                end
                
            elseif contains(extra,'Adjusted','Ignorecase',true)
                if length(LIMO.design.nb_conditions) == 1 && LIMO.design.nb_continuous == 0
                    warning on;
                    if exist('warndlg2','file')
                        warndlg2('Only one condition detected, no adjusted data possible');return
                    else
                        warndlg('Only one condition detected, no adjusted data possible');return
                    end                       
                end
                
                allvar = 1:size(LIMO.design.X,2)-1;
                allvar(regressor)=[]; % all but constant and columns of interest
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    Betas     = load('Betas.mat');
                    if strcmpi(g.restrict,'time')
                        Yr    = squeeze(Yr(channel,freq,:,:));
                        Betas = squeeze(Betas.Betas(channel,freq,:,:));
                    else
                        Yr    = squeeze(Yr(channel,:,time,:));
                        Betas = squeeze(Betas.Betas(channel,:,time,:));
                    end
                    confounds = (LIMO.design.X(:,allvar)*Betas(:,allvar)')';
                else
                    Yr        = squeeze(Yr(channel,:,:));
                    Betas     = load('Betas.mat');
                    Betas     = squeeze(Betas.Betas(channel,:,:));
                    confounds = (LIMO.design.X(:,allvar)*Betas(:,allvar)')';
                end
                Ya = Yr - confounds;
                clear Yr Betas confounds;
                
                if regressor <= length(LIMO.design.nb_conditions) && ...
                        LIMO.design.nb_conditions ~= 0 % for categorical variables
                    if length(LIMO.design.nb_conditions) == 1
                        start = 1;
                    else
                        start = sum(LIMO.design.nb_conditions(1:regressor-1));
                    end
                    
                    for i = (start+LIMO.design.nb_conditions(regressor)-1):-1:start
                        index{i}     = find(LIMO.design.X(:,i));
                        data         = squeeze(Ya(:,index{i})); % use adjusted data
                        average(i,:) = nanmean(data,2);
                        se           = nanstd(data,0,2) ./ sqrt(numel(index{i}));
                        ci(i,:,:)    = repmat(average(i,:),2,1) + repmat(se',2,1).*repmat(z,1,size(Ya,1));
                    end
                    clear mytitle
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Adjusted subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                    else
                        mytitle = sprintf('Adjusted subjects'' parameters at  at optimized channel');
                    end
                    clear Yr
                    
                else % continuous variable ; regressor value already + sum(LIMO.design.nb_conditions)
                    for i=max(regressor):-1:min(regressor)
                        index{i}           = find(LIMO.design.X(:,i));
                        [~,sorting_values] = sort(LIMO.design.X(index{i},i));  % continuous variable 3D plot
                        reg_values(i,:)    = LIMO.design.X(sorting_values,i);
                        continuous(i,:,:)  = Ya(:,sorting_values);
                        clear mytitle
                        if isempty(LIMO.design.electrode)
                            mytitle = sprintf('Adjusted subjects'' parameters \n sorted by regressor %g channel %s (%g)', i, LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle = sprintf('Adjusted subjects'' parameters \n sorted by regressor %g at optimized channel', i);
                        end
                    end
                    remove                 = find(sum(reg_values == 0,2) == size(continuous,3));
                    continuous(remove,:,:) = [];
                    reg_values(remove,:)   = [];
                end
            else
                error('unspecified data type to plot ''Original'',''Modelled'' or ''Adjusted'' ')
            end
            
            % make the figure(s)
            % ------------------
            figure;set(gcf,'Color','w')
            if regressor <= length(LIMO.design.nb_conditions) && ...
                    LIMO.design.nb_conditions ~= 0 % for categorical variables
                for i=1:size(average,1)
                    if strcmpi(LIMO.Analysis,'Time') || strcmpi(g.restrict,'time')
                        timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
                        plot(timevect,squeeze(average(i,:)),'LineWidth',1.5); hold on
                        xlabel('Time in ms','FontSize',14)
                        ylabel('Amplitude (A.U.)','FontSize',14)
                    elseif strcmpi(LIMO.Analysis,'Frequency') || strcmpi(g.restrict,'frequency')
                        freqvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                        plot(freqvect,squeeze(average(i,:)),'LineWidth',1.5); hold on
                        xlabel('Frequency in Hz','FontSize',14)
                        ylabel('Spectral Power (A.U.)','FontSize',14)
                    else
                        error('couldn''t figure out what dimension to plot')
                    end
                    
                    if i==1
                        colorOrder = get(gca, 'ColorOrder');
                        colorOrder = repmat(colorOrder,ceil(size(average,1)/size(colorOrder,1)),1);
                    end
                    x = squeeze(ci(i,1,:)); y = squeeze(ci(i,2,:));
                    fillhandle = patch([timevect fliplr(timevect)], [x',fliplr(y')], colorOrder(i,:));
                    set(fillhandle,'EdgeColor',colorOrder(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                end
                
                h = axis;
                if strcmpi(LIMO.Analysis,'Time')
                    plot(timevect,(sig./10+1).*h(3),'r*','LineWidth',2)
                else
                    plot(freqvect,(sig./10+1).*h(3),'r*','LineWidth',2)
                end
                
                axis tight; grid on; box on
                title(mytitle,'FontSize',16); drawnow;
                assignin('base','Plotted_data', average)
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                set(gca,'FontSize',14);
                if strcmpi(LIMO.Analysis,'Time')
                    ylabel('Amplitude (A.U.)','FontSize',16)
                    xlabel('Time in ms','FontSize',16)
                else
                    ylabel('Spectral Power (A.U.)','FontSize',16)
                    xlabel('Frequency in Hz','FontSize',16)
                end
                
            else % 3D plots
                for i=1:size(continuous,1)
                    if i > 1; figure;set(gcf,'Color','w'); end
                    index = find(~isnan(squeeze(continuous(i,1,:))));
                    if strcmpi(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency')
                        if strcmpi(LIMO.Analysis,'Time')
                            timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in msec
                        else
                            timevect = LIMO.data.tf_times;
                        end
                        surf(index,timevect,squeeze(continuous(i,:,index)));shading interp
                        ylabel('Time in ms','FontSize',16)
                        zlabel('Amplitude (A.U.)','FontSize',16)
                    else
                        surf(index,LIMO.data.freqlist,squeeze(continuous(i,:,index)));shading interp
                        ylabel('Frequency in Hz','FontSize',16)
                        zlabel('Spectral Power (A.U.)','FontSize',16)
                    end
                    % --
                    axis tight; title(mytitle,'FontSize',14); drawnow;
                    xlabel('Sorted variable','FontSize',14)
                    try
                        set(gca,'XTick',index, 'XTickLabels', reg_values(index));
                    catch label_err
                        warning on; warning('could not set X-labels:\n%s',label_err)
                    end
                end
            end
            
            
        elseif contains(LIMO.design.name,'Repeated','IgnoreCase',true)   % All stuffs for repeated measures ANOVA
            % -----------------------------------------------------------------------------
            
            if contains(FileName,'LIMO')
                error('Select summary stat file, nothing to infer from LIMO file')
            end
            
            % which summary stat
            % -------------------
            if ~isempty(g.sumstats) && any(strcmpi(g.sumstats,{'Mean','Trimmed'}))
                extra = g.sumstats;
            else
                if contains(LIMO.design.name,'robust','Ignorecase',true)
                    extra = 'Trimmed Mean';
                else
                    extra = 'Mean';
                end
                % let's not give GUI option and follow the design
                % extra = questdlg('Summarize data using:','Data plot option','Mean','Trimmed Mean','Mean');
                % if isempty(extra)
                %     return
                % end
            end

            if ~contains(FileName,'Rep_ANOVA_Interaction') && ...
                ~contains(FileName,'Rep_ANOVA_Gp') 
                % contains(FileName,{'Rep_ANOVA_Main','Rep_ANOVA'}) 
                % --------------------------------------------------
                                
                % the data to plot are the difference in Yr given LIMO.design.C (see limo_rep_anova)
                if contains(FileName,'Main_effect','IgnoreCase',true)
                    index1     = strfind(FileName,'Main_effect')+length('Main_effect')+1;
                    index2     = max(strfind(FileName,'_'))-1;
                    effect_nb  = eval(FileName(index1:index2));
                elseif contains(FileName,'Interaction','IgnoreCase',true)
                    index1     = strfind(FileName,'Interaction')+length('Interaction')+1;
                    index2     = max(strfind(FileName,'_'))-1;
                    effect_nb  = eval(FileName(index1:index2));
                else
                    index1     = strfind(FileName,'Factor')+length('Factor')+1;
                    effect_nb  = eval(FileName(index1:end));
                end
                C              = LIMO.design.C{effect_nb};
                Data           = load(fullfile(LIMO.dir,'Yr.mat'));
                Data           = Data.(cell2mat(fieldnames(Data)));
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    [~,channel,freq,time] = limo_display_reducedim(Data,LIMO,g.channels,g.restrict,g.dimvalue);
                    Data              = squeeze(Data(channel,freq,time,:,:));
                    sig               = squeeze(single(mask(channel,freq,time)));
                else
                    [~,channel,freq,time] = limo_display_reducedim(Data,LIMO,g.channels);
                    Data              = squeeze(Data(channel,time,:,:)); % note freq/time variables have the same values
                    sig               = single(mask(channel,:));
                end
                sig(sig==0)=NaN;
                
                % compute differences between pairs using C and Cov
                n = size(Data,2);
                if strcmpi(extra,'Mean')
                    for time_or_freq = size(Data,1):-1:1
                        avg(time_or_freq,:) = nanmean(C*squeeze(Data(time_or_freq,:,:))',2);
                        S(time_or_freq,:,:) = nancov(squeeze(Data(time_or_freq,:,:)));
                    end
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Original %s \n %s %s (%g)',mytitle,LIMO.Type(1:end-1),LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Original %s \n virtual %s',mytitle,LIMO.Type(1:end-1));
                    end
                else
                    g=floor((20/100)*n); %% compute for 20% trimmed mean
                    for time_or_freq = size(Data,1):-1:1
                        [v,indices]          = sort(squeeze(Data(time_or_freq,:,:))); % sorted data
                        TD(time_or_freq,:,:) = v((g+1):(n-g),:); % trimmed data - doesn't matter if relationship was kept since we only compute means
                        avg(time_or_freq,:)  = nanmean(C*squeeze(TD(time_or_freq,:,:))',2);
                        v(1:g+1,:)           = repmat(v(g+1,:),g+1,1);
                        v(n-g:end,:)         = repmat(v(n-g,:),g+1,1); % winsorized data
                        [~,reorder]          = sort(indices);
                        for j = size(Data,3):-1:1
                            SD(:,j)          = v(reorder(:,j),j); 
                        end % restore the order of original data
                        S(time_or_freq,:,:)  = cov(SD); % winsorized covariance
                    end
                    clear mytitle
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Trimmed %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Trimmed %s \n optimized channel',mytitle);
                    end
                end
                
                % CI
                dfe = size(Data,2)-size(Data,3)+1;
                % c = avg + 2*(finv(p./(2*size(C,1)),df,dfe).*(sqrt(C*squeeze(S(time_or_freq,:,:))*C'))); % uses Bonferoni inequality
                % b = avg - 2*(finv(p./(2*size(C,1)),df,dfe).*(sqrt(C*squeeze(S(time_or_freq,:,:))*C')));
                bound = (abs(tinv(p./(2*size(C,1)),dfe)).*diag((sqrt(C*squeeze(S(time_or_freq,:,:))*C'))));
                c = avg + repmat(bound', [length(avg),1]);
                b = avg - repmat(bound', [length(avg),1]);
                
                % do the figure
                colours = limo_color_images(size(avg,2));
                if strcmpi(LIMO.Analysis,'Time')
                    if isfield(LIMO.data,'timevect')
                        xvect = LIMO.data.timevect;
                    else
                        LIMO.data.timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
                        save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.timevect;
                    end
                elseif strcmpi(LIMO.Analysis,'Frequency')
                    if isfield(LIMO.data,'timevect')
                        xvect = LIMO.data.freqlist;
                    else
                        LIMO.data.freqlist = linspace(LIMO.data.lowf,LIMO.data.highf,size(toplot,2));
                        save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.freqlist;
                    end
                elseif strcmpi(LIMO.Analysis,'Time-Frequency')
                    if length(time) > 1
                        if isfield(LIMO.data,'tf_times')
                            xvect = LIMO.data.tf_times;
                        else
                            LIMO.data.tf_times = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
                            save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.tf_times;
                        end
                    elseif length(freq) > 1
                        if isfield(LIMO.data,'tf_freqs')
                            xvect = LIMO.data.tf_freqs;
                        else
                            LIMO.data.tf_freqs = linspace(LIMO.data.lowf,LIMO.data.highf,size(toplot,2));
                            save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.tf_freqs;
                        end
                    end
                end
                
                figure;set(gcf,'Color','w')
                for cond = 1:size(c,2)
                    plot(xvect,avg(:,cond)','LineWidth',3,'color',colours(cond,:));
                    fillhandle = patch([xvect fliplr(xvect)], [c(:,cond)',fliplr(b(:,cond)')], colours(cond,:));
                    set(fillhandle,'EdgeColor',colours(cond,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                    hold on
                end
                grid on; box on; axis tight; hold on;
                h = axis;  plot(xvect,(sig./10+1).*h(3),'k.','MarkerSize',20)
                if strcmpi(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(time) > 1
                    xlabel('Time in ms','FontSize',14)
                    ylabel('Amplitude (A.U.)','FontSize',14)
                elseif strcmpi(LIMO.Analysis,'Frequency') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(freq) > 1
                    xlabel('Frequency in Hz','FontSize',14)
                    ylabel('Spectral Power (A.U.)','FontSize',14)
                end
                set(gca,'FontSize',14,'layer','top');
                title(mytitle,'FontSize',16); drawnow;
                trimci = [c ; avg ; b];
                assignin('base','Plotted_data',trimci);
                
                % ----------------------
            elseif contains(FileName,'Rep_ANOVA_Gp')  %% plot pairs of gp differences
                
                % -------------------
                % which ERP to make
                % ------------------
                extra = questdlg('Plotting ERP','ERP Options','Original','Modelled','Original');
                if isempty(extra)
                    return;
                end
                % -----------------------
                Rep_ANOVA_Gp_effect = load(FileName);
                Rep_ANOVA_Gp_effect = Rep_ANOVA_Gp_effect.(cell2mat(fieldnames(Rep_ANOVA_Gp_effect)));
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    [~,channel,freq,time] = limo_display_reducedim(Rep_ANOVA_Gp_effect,LIMO,g.channels,g.restrict,g.dimvalue);
                    Rep_ANOVA_Gp_effect   = squeeze(Rep_ANOVA_Gp_effect(channel,freq,time,:,:));
                    sig                   = squeeze(single(mask(channel,freq,time)));
                else
                    [~,channel,freq,time] = limo_display_reducedim(Rep_ANOVA_Gp_effect,LIMO,g.channels); %#ok<ASGLU>
                    Rep_ANOVA_Gp_effect   = squeeze(Rep_ANOVA_Gp_effect(channel,time,:,:)); % note freq/time variables have the same values
                    sig                   = single(mask(channel,:));
                end
                
                % check channel to plot
                if size(Rep_ANOVA_Gp_effect,1) > 1
                    if isempty(channel)
                        [v,e] = max(Rep_ANOVA_Gp_effect(:,:,1));
                        [~,c] = max(v); channel = e(c);
                    else
                        if ischar(channel)
                            channel = str2double(channel);
                        end
                        
                        if length(channel) > 1
                            error('1 channel only can be plotted')
                        elseif channel > size(Rep_ANOVA_Gp_effect,1)
                            error('channel number invalid')
                        end
                    end
                else
                    if length(LIMO.design.electrode) == 1
                        channel = LIMO.design.electrode;
                    else
                        channel = 1;  % accomodates the fact that all matrices have the channel dim (even = 1)
                    end
                end
                clear Rep_ANOVA_Gp_effect
                
                % compute the pair-wise differences and plot
                Yr = load(fullfile(LIMO.dir,'Yr.mat'));
                Yr = Yr.Yr;
                if strcmpi(extra,'Original')
                    Data = mean(squeeze(Yr(channel,:,:,:)),3);
                    combinations = nchoosek(1:LIMO.design.nb_conditions,2);
                    for d = size(combinations,1):-1:1
                        Effect(:,d) = nanmean(Data(:,find(LIMO.data.Cat == combinations(d,1))),2) - nanmean(Data(:,find(LIMO.data.Cat == combinations(d,2))),2); %#ok<FNDSB>
                    end
                end
                
                % design matrix
                X = zeros(size(Yr,3),LIMO.design.nb_conditions+1);
                X(:,end) = 1;
                for i=1:LIMO.design.nb_conditions
                    X(find(LIMO.data.Cat == i),i) = 1; %#ok<FNDSB>
                end
                
                % data again
                Y = nanmean(squeeze(Yr(channel,:,:,:)),3);
                X = X(find(~isnan(Y(1,:))),:); %#ok<FNDSB>
                Y = Y(:,find(~isnan(Y(1,:))))'; %#ok<FNDSB>
                if strcmpi(extra,'Modelled')
                    beta = pinv(X)*Y; Yhat = X*beta;
                    combinations = nchoosek(1:LIMO.design.nb_conditions,2);
                    for d = 1:size(combinations,1)
                        Effect(:,d) = nanmean(Yhat(find(X(:,combinations(d,1))),:),1)' - nanmean(Yhat(find(X(:,combinations(d,2))),:),1)'; %#ok<FNDSB>
                    end
                end
                
                Res    = (Y'*(eye(size(Y,1)) - (X*pinv(X)))*Y);
                df     = size(Y,1)-rank(X); t = tcdf(1-p,df);
                sigma2 = sum((Res.^2./df),2);
                v      = t.*sqrt(sigma2 ./ norm(X(:,1:end-1)).^2);
                b      = Effect - v(channel,:);
                c      = Effect + v(channel,:);
                if channel == 1
                    if length(LIMO.design.electrode) == 1
                        if strcmpi(extra,'Original')
                            mytitle = sprintf('Mean parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        else
                            mytitle = sprintf('Modelled parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        end
                    else
                        if strcmpi(extra,'Original')
                            mytitle = sprintf('Mean parameter difference between groups \n optimized channel');
                        else
                            mytitle = sprintf('Modelled parameter difference between groups \n optimized channel');
                        end
                    end
                else
                    if strcmpi(extra,'Original')
                        mytitle = sprintf('Mean parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Modelled parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(channel).labels,channel);
                    end
                end
                
                figure;set(gcf,'Color','w'); hold on
                RGB = limo_color_images(size(Effect,2));
                if strcmpi(LIMO.Analysis,'Time')
                    xvect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in sec
                else
                    xvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                end
                
                for d = size(combinations,1):-1:1
                    plot(xvect,Effect(:,d),'LineWidth',3,'Color',RGB(d,:));
                    fillhandle = patch([xvect fliplr(xvect)], [c(:,d)',fliplr(b(:,d)')], RGB(d,:));
                    set(fillhandle,'EdgeColor',RGB(d,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                    Gp_difference{d} = [c(:,d)' ; Effect(:,d)' ; b(:,d)'];
                end
                grid on; box on; axis tight; hold on; 
                h = axis; sig(sig==0)=NaN;
                plot(xvect,(sig./10+1).*h(3),'r.','MarkerSize',20)
                set(gca,'FontSize',14,'layer','top');
                title(mytitle,'FontSize',16); drawnow;
                assignin('base','Plotted_data',Gp_difference);
                
                % -------------------------
            elseif contains(FileName,'Rep_ANOVA_Interaction') % Gp * Repeated measures - plot differences btween condition per gp
                % ------------------------
                
                Rep_ANOVA_Interaction_with_gp = load(FileName);
                Rep_ANOVA_Interaction_with_gp = Rep_ANOVA_Interaction_with_gp.(cell2mat(fieldnames(Rep_ANOVA_Interaction_with_gp)));
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    [~,channel,freq,time]         = limo_display_reducedim(Rep_ANOVA_Interaction_with_gp,LIMO,g.channels,g.restrict,g.dimvalue);
                    Rep_ANOVA_Interaction_with_gp = squeeze(Rep_ANOVA_Interaction_with_gp(channel,freq,time,:,:));
                    sig                           = squeeze(single(mask(channel,freq,time)));
                else
                    [~,channel,freq,time]         = limo_display_reducedim(Rep_ANOVA_Interaction_with_gp,LIMO,g.channels); %#ok<ASGLU>
                    Rep_ANOVA_Interaction_with_gp = squeeze(Rep_ANOVA_Interaction_with_gp(channel,time,:,:)); % note freq/time variables have the same values
                    sig                           = single(mask(channel,:));
                end
                sig(sig==0)=NaN;

                if size(Rep_ANOVA_Interaction_with_gp,1) > 1
                    if isempty(channel) 
                        [v,e] = max(Rep_ANOVA_Interaction_with_gp(:,:,1));
                        [~,c] = max(v); channel = e(c);
                    else
                        if ischar(channel)
                            channel = str2double(channel);
                        end
                        
                        if length(channel) > 1
                            error('1 channel only can be plotted')
                        elseif channel > size(Rep_ANOVA_Interaction_with_gp,1)
                            error('channel number invalid')
                        end
                    end
                else
                    channel = 1;
                end
                
                % the data to plot are the difference in Yr given LIMO.design.C (see limo_rep_anova)
                effect_nb = str2double(FileName(strfind(FileName,'Factor_')+7:strfind(FileName,'Factor_')+6+strfind(FileName(strfind(FileName,'Factor_')+6:end),'_')));
                C         = LIMO.design.C{effect_nb};
                Yr        = load(fullfile(LIMO.dir,'Yr.mat'));
                Yr        = Yr.Yr;
                Data      = squeeze(Yr(channel,:,:,:));
                
                % compute differences between pairs using C and Cov
                if strcmpi(extra,'Mean')
                    for gp = 1:LIMO.design.nb_conditions
                        index = find(LIMO.data.Cat==gp);
                        for time=size(Data,1):-1:1
                            avg(gp,time,:) = nanmean(C*squeeze(Data(time,index,:))',2);
                            S(gp,time,:,:) = cov(squeeze(Data(time,index,:)));
                        end
                    end
                    if ~isempty(LIMO.design.electrode)
                        mytitle = sprintf('Original %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Original %s \n optimized channel',mytitle);
                    end
                else
                    for gp = 1:LIMO.design.nb_conditions
                        index = find(LIMO.data.Cat==gp);
                        n     = length(index);
                        g     = floor((20/100)*n); %% compute for 20% trimmed mean
                        for time_or_freq=size(Data,1):-1:1
                            [v,indices]             = sort(squeeze(Data(time_or_freq,index,:))); % sorted data
                            TD(gp,time_or_freq,:,:) = v((g+1):(n-g),:); % trimmed data - doesn't matter if relationship was kept since we only compute means
                            avg(gp,time_or_freq)    = nanmean(C*squeeze(TD(gp,time_or_freq,:,:))',2);
                            v(1:g+1,:)              = repmat(v(g+1,:),g+1,1);
                            v(n-g:end,:)            = repmat(v(n-g,:),g+1,1); % winsorized data
                            [~,reorder]             = sort(indices);
                            for j = 1:size(Data,3)
                                SD(:,j)             = v(reorder(:,j),j); % restore the order of original data
                            end
                            S(gp,time_or_freq,:,:)  = cov(SD); % winsorized covariance
                        end
                        clear SD
                    end
                    if ~isempty(LIMO.design.electrode)
                        mytitle = sprintf('Trimmed %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Trimmed %s \n optimized channel',mytitle);
                    end
                end
                
                figure; set(gcf,'Color','w'); hold on
                if strcmpi(LIMO.Analysis,'Time')
                    xvect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in sec
                else
                    xvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                end
                
                colorindex = 1;
                trimci     = cell(size(avg,1)*size(avg,3),1);
                RGB        = limo_color_images(size(avg,1)*size(avg,3)); 
                for gp = 1:LIMO.design.nb_conditions
                    % get the variance per comparison
                    for frame = size(avg,2):-1:1
                        varc(:,frame) = diag(sqrt(C*squeeze(S(gp,frame,:,:))*C'));
                    end
                    % plot each comparison
                    for c = 1:size(avg,3)
                        plot(xvect,squeeze(avg(gp,:,c)),'Color',RGB(colorindex,:),'LineWidth',3);
                        % there is an error in the following 2 formulas I cannot fix -- @disbeat
                        index      = find(LIMO.data.Cat==gp);
                        dfe        = size(Data,2)-length(index)+1;
                        up         = avg(gp,:,c) + tinv(p./(2*size(C,1)),dfe).* varc(c,:);
                        down       = avg(gp,:,c) - tinv(p./(2*size(C,1)),dfe).* varc(c,:);
                        fillhandle = patch([xvect fliplr(xvect)], [up,fliplr(down)], RGB(colorindex,:));
                        set(fillhandle,'EdgeColor',RGB(colorindex,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                        trimci{colorindex} = [up ; squeeze(avg(gp,:,c)); down];
                        colorindex = colorindex + 1;
                    end
                end
                
                grid on; box on; axis tight
                h = axis;  hold on;
                plot(xvect,(sig./10+1).*h(3),'r.','MarkerSize',20)
                ylabel('Amplitude (A.U.)','FontSize',14)
                if strcmpi(LIMO.Analysis,'Time')
                    xlabel('Time in ms','FontSize',14)
                else
                    ylabel('Frequency in Hz','FontSize',14)
                end
                set(gca,'FontSize',14,'layer','top');
                title(mytitle,'FontSize',16); drawnow;
                assignin('base','Plotted_data',trimci);
            end
        else
            errordlg('this file is not supported for this kind of plot','Nothing plotted')
        end
    end % closes type
    
elseif strcmpi(LIMO.Level,'LI')
    
    [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO,[],[]);
    
    if Type == 1
        %--------------------------
        % imagesc of the results
        %--------------------------
        if sum(mask(:)) == 0
            warndlg('no values under threshold parameter','no significant effect');
        else
            scale = M.*mask;
            if min(scale(:))<0
                scale(scale==0)=min(scale(:))+(min(scale(:))/10);
            else
                scale(scale==0)=NaN;
            end
            
            figure; set(gcf,'Color','w');
            timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(M,2));
            imagesc(timevect,1:size(M,1),scale);
            title(mytitle,'FontSize',18);
            color_images_(scale,LIMO);
            assignin('base','Plotted_data',scale)
            assignin('base','Mask_of_sig',mask)
        end
        
    elseif Type == 2
        %--------------------------
        % topoplot
        %--------------------------
        if sum(mask(:)) == 0
            warndlg('no values under threshold','no significant effect');
        else
            EEG.data = M.*mask;
            EEG.setname = 'Lateralization Map';
            EEG.pnts = size(EEG.data,2);
            EEG.xmin = LIMO.data.start/1000; % in sec
            EEG.xmax = LIMO.data.end/1000;   % in sec
            EEG.times =  (LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000); % in sec;
            EEG.trials = 1;
            EEG.chanlocs = LIMO.data.chanlocs;
            EEG.nbchan = size(EEG.data,1);
            pop_topoplot(EEG);
            assignin('base','Plotted_data',EEG.data)
        end
        
    elseif Type == 3
        disp('no ERP Plots to Lateralization'); % we could but nobody asked for it
    end
end % closes if LIMO.level
end % closes the function

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                                   ROUTINES
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% color map
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function color_images_(scale,LIMO)

scale(scale==0) = NaN;
cc              = limo_color_images(scale); % get a color map commensurate to that
colormap(cc);

set(gca,'XMinorTick','on','LineWidth',2)
if isfield(LIMO.data,'expected_chanlocs')
    set(gca,'YTick',1:length(LIMO.data.expected_chanlocs));
else
    set(gca,'YTick',1:length(LIMO.data.chanlocs));
end

ylabel('Electrodes','FontSize',14);
if strcmpi(LIMO.Analysis,'Time')
    xlabel('Time in ms','FontSize',16)
elseif strcmpi(LIMO.Analysis,'Frequency')
    xlabel('Frequency in Hz','FontSize',16)
end

if LIMO.Level == 1
    if strcmpi(LIMO.data.chanlocs,'Components')
        label_electrodes = [];
    else
        for i = length(LIMO.data.chanlocs):-1:1
            if isfield(LIMO.data,'expected_chanlocs')
                label_electrodes{i} = LIMO.data.expected_chanlocs(i).labels;
            else
                label_electrodes{i} = LIMO.data.chanlocs(i).labels;
            end
        end
    end
else
    if isempty(LIMO.design.electrode)
        for i = length(LIMO.data.chanlocs):-1:1
            label_electrodes{i} = LIMO.data.chanlocs(i).labels;
        end
    else
        if length(LIMO.design.electrode) == 1
            label_electrodes = LIMO.design.electrode;
        else
            label_electrodes = ' ';
            ylabel('optimized channel','FontSize',14);
        end
    end
end
set(gca,'YTickLabel', label_electrodes);
end

function call_topolot(EEG,FileName,Domain)

EEG.pnts     = size(EEG.data,2);
EEG.nbchan   = size(EEG.data,1);
EEG.trials   = 1;

if strcmpi(FileName,'R2.mat') || strcmpi(FileName,'R2')
    newname = 'R^2';
else
    newname = [FileName(1:min(strfind(FileName,'_'))-1),FileName(max(strfind(FileName,'_'))+1:end-4)];
end

if strcmpi(Domain,'Time')
    if ~isfield(EEG,'setname')
        if contains(FileName,'con','IgnoreCase',true)
            EEG.setname = sprintf('%s - T values',newname);
        else
            EEG.setname = sprintf('%s - F values',newname);
        end
    end
    pop_topoplot(EEG);
    % set(gca,'Colormap',limo_color_images(EEG.data),'CLim',[min(EEG.data(:)),max(EEG.data(:))])
else % freq
    N = size(EEG.freq,2);
    figure;
    for f=1:N
        if N<=6
            subplot(1,N,f)
        else
            subplot(ceil(N/6),6,f);
        end
        [~,ind] = min(abs(EEG.freq-EEG.freq(f)));
        opt = {'electrodes','on','maplimits','maxmin','verbose','off','colormap', limo_color_images(EEG.data(:,ind))};
        topoplot(EEG.data(:,ind),EEG.chanlocs,opt{:});
        if isfield(EEG,'setname')
            title(sprintf('Frequency %g Hz from \n%s',round(EEG.freq(ind)),EEG.setname));
        else
            title(sprintf('Frequency %g Hz from %s - F values',round(EEG.freq(ind)),newname));
        end
    end
end

function create_interactive_3d_plot(xvect, toplot, channel_labels, mytitle, xlabel_text)
    % Create interactive 3D plot with rotation capability
    
    fig = figure('Name', [mytitle ' - Interactive 3D'], 'Color', 'w');
    
    % Create channel vector
    num_channels = size(toplot,1);
    channel_vector = 1:num_channels;
    
    % Create meshgrid
    [X, Z] = meshgrid(xvect, channel_vector);
    
    % Create surface plot
    h = surf(X, toplot, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    
    % Labels and title
    xlabel(xlabel_text, 'FontSize', 14);
    ylabel('Statistical values', 'FontSize', 14);
    zlabel('Channels', 'FontSize', 14);
    title(mytitle, 'FontSize', 16);
    
    % Colormap and colorbar
    colormap(limo_color_images(toplot));
    c = colorbar;
    ylabel(c, 'F values', 'FontSize', 12);
    
    % Set initial view
    view([-45 30]);
    
    % Channel labels on z-axis
    if num_channels <= 20 && exist('channel_labels','var') && ~isempty(channel_labels)
        set(gca, 'ZTick', channel_vector);
        set(gca, 'ZTickLabel', channel_labels);
    end
    
    % Lighting and material
    light('Position', [-1 -1 2], 'Style', 'local');
    lighting gouraud;
    material dull;
    grid on;
    
    % Add rotation toolbar
    rotate3d on;
    
    % Add custom data cursor
    dcm = datacursormode(fig);
    set(dcm, 'UpdateFcn', @(obj,event_obj) myDataCursorText(obj, event_obj, xvect, channel_vector, toplot, channel_labels, xlabel_text));
end

function txt = myDataCursorText(~, event_obj, xvect, channel_vector, data, channel_labels, xlabel_text)
    % Custom data cursor text
    pos = get(event_obj, 'Position');
    
    % Find nearest indices
    [~, x_idx] = min(abs(xvect - pos(1)));
    [~, z_idx] = min(abs(channel_vector - pos(3)));
    
    % Create text
    if contains(xlabel_text, 'Time')
        x_label = sprintf('Time: %.1f ms', xvect(x_idx));
    else
        x_label = sprintf('Frequency: %.1f Hz', xvect(x_idx));
    end
    
    if exist('channel_labels','var') && ~isempty(channel_labels) && z_idx <= length(channel_labels)
        chan_label = sprintf('Channel: %s', channel_labels{z_idx});
    else
        chan_label = sprintf('Channel: %d', z_idx);
    end
    
    txt = {x_label, chan_label, sprintf('Value: %.3f', data(z_idx, x_idx))};
end

end

