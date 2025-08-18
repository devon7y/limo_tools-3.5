function limo_setup_eeglab_colors()
% LIMO_SETUP_EEGLAB_COLORS - Configure LIMO to use EEGLAB's color scheme
%
% This function applies EEGLAB's color scheme to all currently open LIMO
% windows and sets up automatic color inheritance for future LIMO GUIs.
%
% Usage:
%   limo_setup_eeglab_colors()
%
% The function will:
% 1. Load EEGLAB color definitions from icadefs.m
% 2. Apply colors to all currently open LIMO figures
% 3. Display the color values being used
%
% EEGLAB Colors Applied:
%   - Background: GUIBACKCOLOR (light blue)
%   - Text: GUITEXTCOLOR (dark blue)
%   - Buttons: GUIPOPBUTTONCOLOR (light blue/white)
%
% Requirements:
%   - EEGLAB must be in MATLAB path
%   - icadefs.m must be accessible
%
% Example:
%   % After starting EEGLAB and opening LIMO GUIs:
%   limo_gui;                    % Open main LIMO GUI
%   limo_results;                % Open results GUI
%   limo_setup_eeglab_colors();  % Apply EEGLAB colors to both
%
% Note: Future LIMO GUIs opened after running this function will
% automatically inherit EEGLAB colors if the modified LIMO files
% are used.

fprintf('=== LIMO EEGLAB Color Integration ===\n');

% Check if EEGLAB is available
if ~exist('icadefs', 'file')
    error('EEGLAB not found in MATLAB path. Please start EEGLAB first.');
end

try
    % Load EEGLAB color definitions
    fprintf('Loading EEGLAB color definitions...\n');
    icadefs;
    
    if ~exist('GUIBACKCOLOR', 'var')
        error('EEGLAB color variables not found. Please restart EEGLAB.');
    end
    
    % Display the colors being used
    fprintf('\nEEGLAB Color Scheme:\n');
    fprintf('  Background (GUIBACKCOLOR): [%.3f %.3f %.3f] ', GUIBACKCOLOR);
    fprintf('(RGB: %d %d %d)\n', round(GUIBACKCOLOR*255));
    fprintf('  Text (GUITEXTCOLOR):       [%.3f %.3f %.3f] ', GUITEXTCOLOR);
    fprintf('(RGB: %d %d %d)\n', round(GUITEXTCOLOR*255));
    fprintf('  Buttons (GUIPOPBUTTONCOLOR): [%.3f %.3f %.3f] ', GUIPOPBUTTONCOLOR);
    fprintf('(RGB: %d %d %d)\n', round(GUIPOPBUTTONCOLOR*255));
    
    % Find all open figures
    all_figs = findobj('Type', 'figure');
    
    % Identify LIMO figures
    limo_figs = [];
    limo_names = {};
    
    for i = 1:length(all_figs)
        fig = all_figs(i);
        fig_name = get(fig, 'Name');
        fig_tag = get(fig, 'Tag');
        
        % Check if this is a LIMO figure
        if ~isempty(fig_name) && (contains(lower(fig_name), 'limo') || ...
           contains(lower(fig_tag), 'limo') || ...
           contains(lower(fig_name), 'contrast') || ...
           contains(lower(fig_name), 'batch') || ...
           contains(lower(fig_name), 'random'))
            limo_figs(end+1) = fig;
            if isempty(fig_name)
                limo_names{end+1} = sprintf('Figure %d', fig);
            else
                limo_names{end+1} = fig_name;
            end
        end
    end
    
    if isempty(limo_figs)
        fprintf('\nNo LIMO figures currently open.\n');
        fprintf('Colors will be applied automatically when LIMO GUIs are opened.\n');
        return;
    end
    
    fprintf('\nFound %d LIMO figure(s) to update:\n', length(limo_figs));
    
    % Apply colors to each LIMO figure
    for i = 1:length(limo_figs)
        fig = limo_figs(i);
        name = limo_names{i};
        
        fprintf('  %d. %s\n', i, name);
        
        % Set figure background
        set(fig, 'Color', GUIBACKCOLOR);
        
        % Apply colors to all UI elements
        limo_apply_eeglab_colors_to_figure(fig, GUIBACKCOLOR, GUITEXTCOLOR, GUIPOPBUTTONCOLOR);
    end
    
    fprintf('\n✓ EEGLAB color scheme applied successfully!\n');
    fprintf('\nTip: Future LIMO GUIs will automatically inherit these colors\n');
    fprintf('if you''re using the modified LIMO files.\n');
    
catch ME
    fprintf('\n❌ Error applying EEGLAB colors:\n');
    fprintf('   %s\n', ME.message);
    
    if contains(ME.message, 'icadefs')
        fprintf('\nTroubleshooting:\n');
        fprintf('1. Make sure EEGLAB is started: >> eeglab\n');
        fprintf('2. Check that EEGLAB is in your MATLAB path\n');
        fprintf('3. Verify icadefs.m exists in: functions/sigprocfunc/\n');
    end
end

fprintf('\n=== Color Setup Complete ===\n');
end