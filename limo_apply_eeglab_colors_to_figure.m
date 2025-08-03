function limo_apply_eeglab_colors_to_figure(fig, bg_color, text_color, button_color)
% LIMO_APPLY_EEGLAB_COLORS_TO_FIGURE - Apply EEGLAB color scheme to a specific figure
%
% Usage:
%   limo_apply_eeglab_colors_to_figure(fig, bg_color, text_color, button_color)
%
% Inputs:
%   fig         - Figure handle to apply colors to
%   bg_color    - Background color RGB triplet (e.g., GUIBACKCOLOR)
%   text_color  - Text color RGB triplet (e.g., GUITEXTCOLOR)  
%   button_color- Button color RGB triplet (e.g., GUIPOPBUTTONCOLOR)
%
% This function systematically applies EEGLAB's color scheme to all UI
% elements in a LIMO figure window.

if nargin < 4
    % Default to standard MATLAB colors if not provided
    button_color = [0.94 0.94 0.94];
end

% Get all children of the figure recursively
all_objects = findobj(fig);

for i = 1:length(all_objects)
    obj = all_objects(i);
    
    % Skip the figure itself (already handled)
    if obj == fig
        continue;
    end
    
    try
        obj_type = get(obj, 'Type');
        
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
                    case 'text'
                        % Static text
                        set(obj, 'BackgroundColor', bg_color);
                        set(obj, 'ForegroundColor', text_color);
                        
                    case 'edit'
                        % Edit boxes - keep white background for readability
                        set(obj, 'BackgroundColor', button_color);
                        set(obj, 'ForegroundColor', text_color);
                        
                    case {'pushbutton', 'togglebutton'}
                        % Buttons
                        set(obj, 'BackgroundColor', button_color);
                        set(obj, 'ForegroundColor', text_color);
                        
                    case {'radiobutton', 'checkbox'}
                        % Radio buttons and checkboxes
                        set(obj, 'BackgroundColor', bg_color);
                        set(obj, 'ForegroundColor', text_color);
                        
                    case {'listbox', 'popupmenu'}
                        % List controls - keep white background for readability
                        set(obj, 'BackgroundColor', button_color);
                        set(obj, 'ForegroundColor', text_color);
                        
                    case 'slider'
                        % Sliders
                        set(obj, 'BackgroundColor', bg_color);
                        
                    case 'frame'
                        % Frames
                        set(obj, 'BackgroundColor', bg_color);
                end
                
            case 'uibuttongroup'
                % Button groups
                set(obj, 'BackgroundColor', bg_color);
                if isprop(obj, 'ForegroundColor')
                    set(obj, 'ForegroundColor', text_color);
                end
                
            case 'uitable'
                % Tables - keep readable
                set(obj, 'BackgroundColor', [1 1 1]);
                if isprop(obj, 'ForegroundColor')
                    set(obj, 'ForegroundColor', text_color);
                end
                
            case 'axes'
                % Only change axes background if they don't contain data plots
                children = get(obj, 'Children');
                if isempty(children)
                    set(obj, 'Color', bg_color);
                end
        end
        
    catch ME
        % Skip objects that don't support color properties or have other issues
        % This is normal for some object types
        continue;
    end
end

% Force a refresh to show the changes
drawnow;
end