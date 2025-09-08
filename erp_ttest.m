function erp_ttest(var1, var2, electrode, time_start, time_end)
% ERP_TTEST Performs paired t-test between two ERP conditions
%
% Inputs:
%   var1        - First ERP variable name (string, e.g., 'SMEhits_ERP')
%   var2        - Second ERP variable name (string, e.g., 'SMEmisses_ERP')
%   electrode   - Electrode number (1-256)
%   time_start  - Start time in milliseconds (e.g., 400)
%   time_end    - End time in milliseconds (e.g., 700)
%
% Output: Results printed to command window
%
% Data format expected: 256 electrodes x 1300 samples x 55 subjects
% Time range: -100 to 2600 ms (samples 1-1300)
    
    % Convert milliseconds to samples
    % ERP data spans -100 to 2600 ms, so we need to add 100ms offset
    % Then divide by sampling factor of 2: sample = (ms + 100) / 2
    sample_start = round((time_start + 100) / 2);
    sample_end = round((time_end + 100) / 2);
    
    % Get the data from workspace variables
    try
        erp1 = evalin('caller', var1);
        erp2 = evalin('caller', var2);
    catch ME
        error('Could not find variables %s or %s in workspace: %s', var1, var2, ME.message);
    end
    
    % Extract data for specified electrode and time range (using converted samples)
    data1_electrode = squeeze(erp1(electrode, sample_start:sample_end, :));    % time x subjects
    data2_electrode = squeeze(erp2(electrode, sample_start:sample_end, :));     % time x subjects
    
    % Calculate mean across time samples for each subject
    data1_mean = mean(data1_electrode, 1);   % 1 x 55 subjects
    data2_mean = mean(data2_electrode, 1);   % 1 x 55 subjects
    
    % Perform paired t-test using MATLAB's built-in function
    [~, p_value, ~, stats] = ttest(data1_mean, data2_mean);
    
    % Extract statistics from the ttest output
    t_stat = stats.tstat;
    
    % Print results
    fprintf('\n=== ERP T-Test Results ===\n');
    fprintf('Electrode: %d\n', electrode);
    fprintf('Time range: %d - %d ms (samples %d - %d)\n', time_start, time_end, sample_start, sample_end);
    fprintf('\n');
    fprintf('t-statistic: %.4f\n', t_stat);
    fprintf('p-value (two-tailed): %.6f\n', p_value);
    fprintf('\n');
    
end