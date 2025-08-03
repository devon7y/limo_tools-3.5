function limo_eeglab_colors()
% LIMO_EEGLAB_COLORS - Apply EEGLAB color scheme to LIMO GUIs
%
% This function reads EEGLAB's color definitions from icadefs.m and applies
% them to all open LIMO GUI windows. Call this function after opening any
% LIMO GUI to inherit EEGLAB's color scheme.
%
% Usage:
%   limo_eeglab_colors()  % Apply colors to all open LIMO figures
%
% The function will:
% - Set figure backgrounds to EEGLAB's GUIBACKCOLOR
% - Set text colors to EEGLAB's GUITEXTCOLOR  
% - Set button colors to EEGLAB's GUIPOPBUTTONCOLOR
%
% EEGLAB Color Variables Used:
%   BACKCOLOR           - Background for analysis/plot figures
%   BACKEEGLABCOLOR     - Main EEGLAB window background
%   GUIBACKCOLOR        - GUI dialog backgrounds
%   GUITEXTCOLOR        - Text color in GUI elements
%   GUIPOPBUTTONCOLOR   - Button colors in GUI windows
%
% Example:
%   limo_gui;           % Open LIMO main GUI
%   limo_eeglab_colors; % Apply EEGLAB colors
%
% Author: Generated to integrate LIMO with EEGLAB color scheme
% Copyright (C) LIMO Team

try
    % Load EEGLAB color definitions
    icadefs;
    
    fprintf('Applying EEGLAB color scheme to LIMO GUIs...\n');
    
    % Get all figure handles
    figs = findobj('Type', 'figure');
    
    % Filter for LIMO-related figures
    limo_figs = [];
    for i = 1:length(figs)
        fig_name = get(figs(i), 'Name');
        if contains(lower(fig_name), 'limo') || ...
           contains(lower(get(figs(i), 'Tag')), 'limo')
            limo_figs(end+1) = figs(i);
        end
    end
    
    if isempty(limo_figs)
        fprintf('No LIMO figures found. Colors will be applied when LIMO GUIs are opened.\n');
        return;
    end
    
    fprintf('Found %d LIMO figure(s) to update.\n', length(limo_figs));
    
    % Apply colors to each LIMO figure
    for i = 1:length(limo_figs)
        fig = limo_figs(i);
        fig_name = get(fig, 'Name');
        
        fprintf('  Updating colors for: %s\n', fig_name);
        
        % Set figure background color
        set(fig, 'Color', GUIBACKCOLOR);
        
        % Update all UI elements in this figure
        apply_colors_to_figure(fig, GUIBACKCOLOR, GUITEXTCOLOR, GUIPOPBUTTONCOLOR);
    end
    
    fprintf('EEGLAB color scheme applied successfully!\n');
    fprintf('Colors used:\n');
    fprintf('  Background: [%.2f %.2f %.2f]\n', GUIBACKCOLOR);
    fprintf('  Text: [%.2f %.2f %.2f]\n', GUITEXTCOLOR);
    fprintf('  Buttons: [%.2f %.2f %.2f]\n', GUIPOPBUTTONCOLOR);
    
catch ME
    fprintf('Error applying EEGLAB colors: %s\n', ME.message);
    fprintf('Make sure EEGLAB is in your path and icadefs.m is accessible.\n');
end
end

function apply_colors_to_figure(fig, bg_color, text_color, button_color)
% Apply color scheme to all elements in a figure

% Get all children of the figure
all_objects = findobj(fig);

for i = 1:length(all_objects)
    obj = all_objects(i);
    obj_type = get(obj, 'Type');
    
    try
        switch obj_type
            case 'uipanel'
                % Panels and containers
                set(obj, 'BackgroundColor', bg_color);
                if isprop(obj, 'ForegroundColor')
                    set(obj, 'ForegroundColor', text_color);
                end
                
            case 'uicontrol'
                % Get the style of the uicontrol
                style = get(obj, 'Style');
                switch style
                    case {'text', 'edit'}
                        % Text and edit boxes
                        if strcmp(style, 'text')
                            set(obj, 'BackgroundColor', bg_color);
                        end
                        set(obj, 'ForegroundColor', text_color);
                        
                    case {'pushbutton', 'togglebutton', 'radiobutton', 'checkbox'}
                        % Buttons
                        set(obj, 'BackgroundColor', button_color);
                        set(obj, 'ForegroundColor', text_color);
                        
                    case {'listbox', 'popupmenu'}
                        % List controls
                        set(obj, 'BackgroundColor', [1 1 1]); % Keep white for readability
                        set(obj, 'ForegroundColor', text_color);
                        
                    case 'slider'
                        % Sliders
                        set(obj, 'BackgroundColor', bg_color);
                end
                
            case 'axes'
                % Axes - only change if they don't contain plots
                children = get(obj, 'Children');
                if isempty(children)
                    set(obj, 'Color', bg_color);
                end
        end
    catch
        % Skip objects that don't support color properties
        continue;
    end
end
end