function results_table = erp_ttest_batch(var1, var2, electrode_list, time_windows, sampling_rate_hz, significance_threshold, correction_method)
% ERP_TTEST_BATCH_CORRECTIONS Performs multiple paired t-tests with correction methods
%
% Inputs:
%   var1                   - First ERP variable name (string, e.g., 'SMEhits_ERP')
%   var2                   - Second ERP variable name (string, e.g., 'SMEmisses_ERP')
%   electrode_list         - Vector of electrode numbers (e.g., [21, 36, 101, 224])
%   time_windows           - Matrix of time windows [start_ms, end_ms] (e.g., [400 700; 700 900])
%   sampling_rate_hz       - Sampling rate in Hz (e.g., 500 for 500Hz, 250 for 250Hz)
%   significance_threshold - P-value threshold for significance (e.g., 0.05)
%   correction_method      - String: 'none', 'bonferroni', 'holm', 'fdr_bh', 'fdr_by', 'sidak'
%
% Output:
%   results_table   - Struct array with fields: Electrode, Time_Window, T_Stat, P_Value, 
%                     Corrected_P_Value, Significant, Corrected_Significant
%
% Example:
%   electrodes = [21, 36, 101, 224];
%   windows = [400 700; 700 900; 900 1200];
%   results = erp_ttest_batch_corrections('SMEhits_ERP', 'SMEmisses_ERP', electrodes, windows, 500, 0.05, 'fdr_bh');

    % Set defaults if not provided
    if nargin < 5
        error('sampling_rate_hz is required. Specify the sampling rate in Hz (e.g., 500, 250)');
    end
    if nargin < 6
        significance_threshold = 0.05;
    end
    if nargin < 7
        correction_method = 'none';
    end
    
    % Calculate conversion factor from sampling rate
    % ms_per_sample = 1000 / sampling_rate_hz
    ms_per_sample = 1000 / sampling_rate_hz;
    
    % Initialize results storage
    results = [];
    
    % Get the ERP data from caller's workspace once
    try
        erp1 = evalin('caller', var1);
        erp2 = evalin('caller', var2);
    catch ME
        error('Could not find variables %s or %s in workspace: %s', var1, var2, ME.message);
    end
    
    % Get total number of tests
    n_electrodes = length(electrode_list);
    n_windows = size(time_windows, 1);
    
    % Loop through each electrode and time window combination
    test_count = 0;
    for i = 1:n_electrodes
        electrode = electrode_list(i);
        
        for j = 1:n_windows
            time_start = time_windows(j, 1);
            time_end = time_windows(j, 2);
            
            test_count = test_count + 1;
            
            % Get statistics by calling internal function with the loaded data
            [t_stat, p_value] = erp_ttest_internal(erp1, erp2, electrode, time_start, time_end, ms_per_sample);
            
            % Store results
            result_row.Electrode = electrode;
            result_row.Time_Window = sprintf('%d-%d ms', time_start, time_end);
            result_row.T_Stat = t_stat;
            result_row.P_Value = p_value;
            result_row.Significant = p_value < significance_threshold;
            
            results = [results; result_row];
        end
    end
    
    % Apply multiple comparison correction
    p_values = [results.P_Value]';
    
    switch lower(correction_method)
        case 'none'
            corrected_p = p_values;
        case 'bonferroni'
            corrected_p = bonferroni_correction(p_values);
        case 'holm'
            corrected_p = holm_bonferroni_correction(p_values);
        case 'fdr_bh'
            corrected_p = benjamini_hochberg_correction(p_values);
        case 'fdr_by'
            corrected_p = benjamini_yekutieli_correction(p_values);
        case 'sidak'
            corrected_p = sidak_correction(p_values);
        otherwise
            error('Unknown correction method: %s. Use: none, bonferroni, holm, fdr_bh, fdr_by, sidak', correction_method);
    end
    
    % Add corrected p-values and significance to results
    for i = 1:length(results)
        results(i).Corrected_P_Value = corrected_p(i);
        results(i).Corrected_Significant = corrected_p(i) < significance_threshold;
    end
    
    % Convert to struct array for easier access
    results_table = results;
    
    % Display summary table
    fprintf('\nPaired T-Tests between %s - %s\n', var1, var2);
    fprintf('%-10s %-15s %-10s %-12s %-12s %-12s %-15s\n', 'Electrode', 'Time Window', 'T-Stat', 'P-Value', 'Corrected P', 'Significant', 'Corr Significant');
    
    for i = 1:length(results_table)
        sig_text = 'No';
        if results_table(i).Significant
            sig_text = 'Yes';
        end
        
        corr_sig_text = 'No';
        if results_table(i).Corrected_Significant
            corr_sig_text = 'Yes';
        end
        
        fprintf('%-10d %-15s %-10.4f %-12.6f %-12.6f %-12s %-15s', ...
                results_table(i).Electrode, ...
                results_table(i).Time_Window, ...
                results_table(i).T_Stat, ...
                results_table(i).P_Value, ...
                results_table(i).Corrected_P_Value, ...
                sig_text, ...
                corr_sig_text);
        
        % Add significance indicator based on corrected p-values
        if results_table(i).Corrected_P_Value < 0.001
            fprintf(' ***');
        elseif results_table(i).Corrected_P_Value < 0.01
            fprintf(' **');
        elseif results_table(i).Corrected_P_Value < 0.05
            fprintf(' *');
        end
        fprintf('\n');
    end
    fprintf('\n');
end

function [t_stat, p_value] = erp_ttest_internal(erp1, erp2, electrode, time_start, time_end, ms_per_sample)
    % Internal function that performs t-test on already loaded ERP data
    
    % Convert milliseconds to samples using the provided conversion factor
    % For data with baseline, add 100ms offset before converting
    sample_start = round((time_start + 100) / ms_per_sample);
    sample_end = round((time_end + 100) / ms_per_sample);
    
    % Extract data for specified electrode and time range
    data1_electrode = squeeze(erp1(electrode, sample_start:sample_end, :));
    data2_electrode = squeeze(erp2(electrode, sample_start:sample_end, :));
    
    % Calculate mean across time samples for each subject
    data1_mean = mean(data1_electrode, 1);
    data2_mean = mean(data2_electrode, 1);
    
    % Perform paired t-test using MATLAB's built-in function
    [~, p_value, ~, stats] = ttest(data1_mean, data2_mean);
    
    % Extract t-statistic
    t_stat = stats.tstat;
end

function corrected_p = bonferroni_correction(p_values)
    % Bonferroni correction: multiply each p-value by number of tests
    n = length(p_values);
    corrected_p = min(p_values * n, 1);  % Cap at 1.0
end

function corrected_p = holm_bonferroni_correction(p_values)
    % Holm-Bonferroni step-down correction
    n = length(p_values);
    [sorted_p, sort_idx] = sort(p_values);
    
    % Apply step-down correction
    corrected_p_sorted = zeros(n, 1);
    for i = 1:n
        corrected_p_sorted(i) = min(sorted_p(i) * (n - i + 1), 1);
    end
    
    % Ensure monotonicity (corrected p-values should be non-decreasing)
    for i = 2:n
        corrected_p_sorted(i) = max(corrected_p_sorted(i), corrected_p_sorted(i-1));
    end
    
    % Map back to original order
    corrected_p = zeros(n, 1);
    corrected_p(sort_idx) = corrected_p_sorted;
end

function corrected_p = benjamini_hochberg_correction(p_values)
    % Benjamini-Hochberg False Discovery Rate correction
    n = length(p_values);
    [sorted_p, sort_idx] = sort(p_values);
    
    % Apply step-up correction
    corrected_p_sorted = zeros(n, 1);
    corrected_p_sorted(n) = sorted_p(n);  % Start from the largest
    
    for i = (n-1):-1:1
        corrected_p_sorted(i) = min(corrected_p_sorted(i+1), ...
                                   sorted_p(i) * n / i);
    end
    
    % Ensure corrected p-values don't exceed 1
    corrected_p_sorted = min(corrected_p_sorted, 1);
    
    % Map back to original order
    corrected_p = zeros(n, 1);
    corrected_p(sort_idx) = corrected_p_sorted;
end

function corrected_p = benjamini_yekutieli_correction(p_values)
    % Benjamini-Yekutieli False Discovery Rate correction for dependent tests
    n = length(p_values);
    [sorted_p, sort_idx] = sort(p_values);
    
    % Calculate harmonic sum c(n) for dependent tests
    c_n = sum(1./(1:n));  % Harmonic series: 1 + 1/2 + 1/3 + ... + 1/n
    
    % Apply step-up correction with harmonic sum
    corrected_p_sorted = zeros(n, 1);
    corrected_p_sorted(n) = sorted_p(n) * c_n;  % Start from the largest
    
    for i = (n-1):-1:1
        corrected_p_sorted(i) = min(corrected_p_sorted(i+1), ...
                                   sorted_p(i) * n / i * c_n);
    end
    
    % Ensure corrected p-values don't exceed 1
    corrected_p_sorted = min(corrected_p_sorted, 1);
    
    % Map back to original order
    corrected_p = zeros(n, 1);
    corrected_p(sort_idx) = corrected_p_sorted;
end

function corrected_p = sidak_correction(p_values)
    % Šidák correction for independent tests
    n = length(p_values);
    % Corrected alpha: 1 - (1-alpha)^(1/n)
    % For p-values: 1 - (1-p)^n
    corrected_p = 1 - (1 - p_values).^n;
    corrected_p = min(corrected_p, 1);  % Cap at 1.0
end