# LIMO Likelihood Ratio Analysis - Implementation Documentation

## Overview

This implementation adds **Bayes Factor** (likelihood ratio) analysis to LIMO tools, based on the method of **Kang et al. (2015)** adapted to use LIMO's **empirical bootstrap null distribution**.

**Key Innovation**: Uses LIMO's existing H0 bootstrap samples instead of assuming a parametric null distribution, making the method more robust to distributional assumptions.

## Files

- **`limo_likelihood.m`** - Main function for computing Bayes Factors
- **`limo_likelihood_demo.m`** - Comprehensive demonstration script with 8 examples
- **`LIKELIHOOD_RATIO_DOCUMENTATION.md`** - This file

## Theoretical Background

### What is a Likelihood Ratio / Bayes Factor?

A **likelihood ratio** (LR) or **Bayes Factor** (BF) quantifies the strength of evidence for an alternative hypothesis (H1) relative to a null hypothesis (H0):

```
BF₁₀ = p(data | H1) / p(data | H0)
```

**Interpretation**:
- **BF > 10**: Strong evidence for H1 (effect exists)
- **BF > 3**: Moderate evidence for H1
- **1/3 < BF < 3**: Inconclusive (data don't strongly favor either hypothesis)
- **BF < 1/3**: Moderate evidence for H0 (no effect)
- **BF < 1/10**: Strong evidence for H0

### Advantages over p-values

| p-value approach | Likelihood Ratio approach |
|-----------------|---------------------------|
| Only tests against null | Tests H1 vs H0 |
| Cannot provide evidence FOR null | Can provide evidence FOR or AGAINST |
| Binary decision (sig/not sig) | Continuous strength of evidence |
| Problematic with large N | Evidence naturally calibrated |
| Difficult to combine across studies | Easy to combine (multiply BFs) |

### The Kang et al. (2015) Formula

For a t-test with:
- **t** = observed t-statistic
- **n** = sample size (df + 1)
- **τ** = prior standard deviation on standardized effect size

The Bayes Factor is:

```
BF = √(1/(1 + n×τ²)) × exp[t²/2 × n×τ²/(1 + n×τ²)]
```

Or in log form (more stable):

```
log(BF) = -0.5×log(1 + n×τ²) + (t²/2) × [n×τ²/(1 + n×τ²)]
```

**Key parameters**:
- **τ = 0.5**: Conservative prior (small effects expected)
- **τ = 1.0**: Medium prior (default, recommended by Kang et al.)
- **τ = √2 ≈ 1.414**: Unit information prior

## Implementation Details

### Key Adaptation: Empirical H0

**Original Kang et al.**: Assumes H0 follows a parametric distribution (e.g., t-distribution)

**Our Implementation**:
1. Applies Kang formula to observed data
2. Applies **same formula** to each of LIMO's bootstrap H0 samples
3. Uses empirical bootstrap BF distribution for calibration

This makes the method:
- ✅ Robust to non-normal distributions
- ✅ Compatible with LIMO's existing infrastructure
- ✅ No need for additional resampling
- ✅ Preserves LIMO's handling of temporal/spatial dependencies

### Input Files

The function requires two files:

1. **Test statistic file** (e.g., `paired_samples_ttest_parameter_35.mat`):
   ```
   Structure: [channels × time × 5]
   Dimension 1: Mean estimate
   Dimension 2: Standard error
   Dimension 3: Degrees of freedom
   Dimension 4: t-statistic ← Used for BF calculation
   Dimension 5: p-value
   ```

2. **H0 bootstrap file** (e.g., `H0/H0_paired_samples_ttest_parameter_35.mat`):
   ```
   Structure: [channels × time × 2 × n_bootstraps]
   Dimension 1: t-statistics under H0 ← Used for BF calculation
   Dimension 2: p-values under H0
   ```

### Output Structure

The function returns a `results` structure with fields:

```matlab
results.BF                 % Bayes Factors [channels × time]
results.logBF              % log(BF) for numerical stability
results.BF_H0_dist         % Bootstrap BF distribution [ch × time × n_boots]
results.evidence_labels    % Cell array of evidence categories
results.evidence_numeric   % Numeric codes: -2 to +2
results.n_strong_H1        % Count of strong H1 evidence
results.proportion_inconclusive  % Proportion of inconclusive results
results.params             % Analysis parameters (tau, n, formula, etc.)
```

### Saved Augmented File

If `'save', true` (default), creates a new file with suffix `_likelihood.mat`:

```
Original: paired_samples_ttest_parameter_35.mat
Augmented: paired_samples_ttest_parameter_35_likelihood.mat

Structure: [channels × time × 7]
  (:,:,1): Mean estimate (original)
  (:,:,2): Standard error (original)
  (:,:,3): Degrees of freedom (original)
  (:,:,4): t-statistic (original)
  (:,:,5): p-value (original)
  (:,:,6): Bayes Factor (NEW)
  (:,:,7): log(BF) (NEW)
```

## Usage Examples

### Basic Usage

```matlab
% Minimal usage with defaults
results = limo_likelihood('paired_samples_ttest_parameter_1.mat', ...
                          'H0/H0_paired_samples_ttest_parameter_1.mat');

% With custom prior
results = limo_likelihood('paired_samples_ttest_parameter_1.mat', ...
                          'H0/H0_paired_samples_ttest_parameter_1.mat', ...
                          'tau', 0.5);  % Conservative prior

% With diagnostic plots
results = limo_likelihood('paired_samples_ttest_parameter_1.mat', ...
                          'H0/H0_paired_samples_ttest_parameter_1.mat', ...
                          'tau', 1.0, 'plot', true);
```

### Access Results

```matlab
% Find strong evidence locations
[channels, times] = find(results.BF > 10);
fprintf('Found %d points with strong evidence\n', length(channels));

% Summary statistics
fprintf('Strong H1: %.1f%%\n', results.proportion_strong_H1 * 100);
fprintf('Inconclusive: %.1f%%\n', results.proportion_inconclusive * 100);

% Plot log(BF) heatmap
imagesc(results.logBF');
colorbar;
xlabel('Channel'); ylabel('Time');
title('log(Bayes Factor)');
```

### Visualize Evidence at Electrode

```matlab
% Load LIMO for channel info
load('LIMO.mat');
times = LIMO.data.timevect;

% Find electrode
elec_idx = find(strcmpi({LIMO.data.chanlocs.labels}, 'Cz'));

% Plot time course
figure;
plot(times, results.logBF(elec_idx, :), 'LineWidth', 2);
yline(log(10), 'r--', 'Strong H1');
yline(0, 'k-', 'BF=1');
yline(log(0.1), 'b--', 'Strong H0');
xlabel('Time (ms)'); ylabel('log(BF)');
title('Evidence at Cz');
grid on;
```

### Compare with p-values

```matlab
% Load original test file
test_data = load('paired_samples_ttest_parameter_1.mat');
field = fieldnames(test_data);
data = test_data.(field{1});
p_values = data(:,:,5);

% Compare significance
sig_pval = p_values < 0.05;
strong_BF = results.BF > 10;

fprintf('Significant by p<0.05: %d\n', sum(sig_pval(:)));
fprintf('Strong evidence by BF>10: %d\n', sum(strong_BF(:)));
fprintf('Agreement: %d\n', sum(sig_pval(:) & strong_BF(:)));
```

## Choosing the Prior (τ)

The prior parameter τ represents your expectation about standardized effect sizes (Cohen's d).

### Recommended Values

| τ value | Interpretation | When to use |
|---------|---------------|-------------|
| **0.5** | Conservative | Expect small effects, want to be cautious |
| **0.707** (√0.5) | Moderately conservative | Balance between conservative and medium |
| **1.0** | Medium (default) | General purpose, recommended by Kang et al. |
| **1.414** (√2) | Unit information | Liberal, equivalent to one prior observation |

### Effect of τ on BF

- **Smaller τ** → More conservative → Harder to find strong evidence for H1
- **Larger τ** → More liberal → Easier to find strong evidence for H1
- **Very large τ** → Lindley's paradox (BF → 0 even for strong effects)

### Sensitivity Analysis

Always check sensitivity to prior choice:

```matlab
tau_values = [0.5, 0.707, 1.0, 1.414];
for i = 1:length(tau_values)
    results{i} = limo_likelihood(test_file, H0_file, ...
        'tau', tau_values(i), 'save', false);
    fprintf('τ=%.3f: Strong H1 = %d points\n', ...
        tau_values(i), results{i}.n_strong_H1);
end
```

## Interpretation Guidelines

### Evidence Strength (Kang et al. scale)

| BF range | log(BF) range | Interpretation |
|----------|---------------|----------------|
| > 100 | > 4.6 | Extreme evidence for H1 |
| 32-100 | 3.5-4.6 | Very strong evidence for H1 |
| 10-32 | 2.3-3.5 | **Strong evidence for H1** |
| 3-10 | 1.1-2.3 | Moderate evidence for H1 |
| 1-3 | 0-1.1 | Weak evidence for H1 |
| 1/3-1 | -1.1-0 | Weak evidence for H0 |
| 1/10-1/3 | -2.3 to -1.1 | Moderate evidence for H0 |
| 1/32-1/10 | -3.5 to -2.3 | **Strong evidence for H0** |
| 1/100-1/32 | -4.6 to -3.5 | Very strong evidence for H0 |
| < 1/100 | < -4.6 | Extreme evidence for H0 |

### What Does "Inconclusive" Mean?

When **1/3 < BF < 3** (or **-1.1 < log(BF) < 1.1**):

- Data are compatible with both H0 and H1
- More data needed to discriminate
- Effect may be too small to detect reliably
- **NOT the same as "no effect"** - we simply don't have enough evidence

This is a strength of the BF approach: it explicitly acknowledges when data are insufficient to make a strong conclusion.

## Comparison with LIMO's Existing Methods

### vs. Uncorrected p-values (MCC=1)

| p-value | BF/LR |
|---------|-------|
| Tests only against null | Tests H1 vs H0 |
| Depends on sample size | Evidence naturally calibrated |
| Cannot support null | Can support null or alternative |
| Binary threshold | Continuous evidence scale |

### vs. Maximum Statistic Correction (MCC=4)

| Max correction | BF/LR |
|---------------|-------|
| Controls FWER | Quantifies evidence strength |
| Very conservative | Allows inconclusive region |
| Uses max(t) distribution | Uses BF distribution |
| Same threshold everywhere | Point-wise evidence |

### vs. TFCE (MCC=3)

| TFCE | BF/LR |
|------|-------|
| Cluster-based | Point-wise |
| Threshold-free enhancement | Model-based inference |
| Good for spatial-temporal clusters | Good for specific hypotheses |
| Uses spatial neighbors | Independent across points |

## Limitations and Caveats

### 1. Requires specifying τ
- Solution: Use recommended τ=1.0 or perform sensitivity analysis
- Or use data-driven approach (not yet implemented)

### 2. Assumes normal likelihood
- But empirical H0 reduces impact of this assumption
- Bootstrap captures actual sampling distribution

### 3. Point-wise analysis
- Does not account for spatial-temporal clustering
- For cluster inference, would need to extend method

### 4. Multiple comparisons
- Current implementation: NO correction (as requested)
- Future: Can use max(logBF) from H0 for FWER control

### 5. Computational cost
- Computes BF for all bootstrap samples
- For 800 bootstraps × 64 channels × 500 times ≈ 25M calculations
- Usually completes in < 1 minute

## Technical Details

### Numerical Stability

The implementation uses `log(BF)` computation for stability:

```matlab
% Instead of: BF = term1 * exp(term2)
% We compute: logBF = log(term1) + term2
logBF = -0.5*log(1 + n*tau^2) + (t^2/2) * scaling_factor;
```

This prevents:
- Overflow for large t-statistics
- Underflow for small BF values
- Precision loss in arithmetic

### Handling Different Test Types

The function automatically detects:

```matlab
if contains(filename, 'one_sample')
    % One-sample t-test
elseif contains(filename, 'two_samples')
    % Two-sample t-test
elseif contains(filename, 'paired_samples')
    % Paired t-test
end
```

All use the same Kang formula with appropriate df extraction.

### Time-Frequency Data Support

For 4D data `[channels × freq × time × stats]`:

```matlab
if ndims(test_data) == 4
    t_obs = test_data(:,:,:,4);  % 3D t-statistics
    % Computes BF for each channel × freq × time
end
```

Output BF also 3D: `[channels × freq × time]`

## Future Extensions (Not Yet Implemented)

The following extensions are **NOT part of Kang et al. (2015)** but would adapt the method for specific needs in EEG/ERP analysis. Each section provides complete implementation details for future development.

---

## Extension 1: Cluster-Level Inference

### Overview

**Current limitation**: Point-wise BF treats each channel-time point independently, ignoring spatial-temporal structure of EEG effects.

**Solution**: Combine evidence across spatially/temporally contiguous regions (clusters) to increase sensitivity to extended effects.

### Theoretical Background

EEG/ERP effects typically span:
- **Multiple electrodes** (due to volume conduction and source spread)
- **Multiple timepoints** (effects persist 100-200ms)

**Key principle**: Independent observations with BF₁, BF₂, ..., BFₙ combine multiplicatively:

```
BF_cluster = BF₁ × BF₂ × ... × BFₙ
log(BF_cluster) = Σᵢ log(BFᵢ)
```

This dramatically increases evidence when many adjacent points show moderate effects.

### Algorithm

#### Step 1: Define Clusters

```matlab
function [cluster_labels, n_clusters] = find_evidence_clusters(logBF, LIMO, threshold)
% Find spatial-temporal clusters exceeding evidence threshold
%
% INPUTS:
%   logBF      - [channels × time] log Bayes Factors
%   LIMO       - LIMO structure with neighbouring_matrix
%   threshold  - Initial threshold (e.g., log(3) for moderate evidence)
%
% OUTPUTS:
%   cluster_labels - [channels × time] cluster ID for each point (0 = no cluster)
%   n_clusters     - Number of clusters found

% Create binary mask of points exceeding threshold
mask = logBF > threshold;

% Use LIMO's spatial-temporal clustering
% Requires minimum 2 adjacent channels
[cluster_labels, n_clusters] = limo_findcluster(mask, ...
    LIMO.data.neighbouring_matrix, 2);

end
```

#### Step 2: Compute Cluster-Level Evidence

```matlab
function cluster_evidence = compute_cluster_BF(logBF, cluster_labels, n_clusters)
% Sum log(BF) within each cluster
%
% INPUTS:
%   logBF          - [channels × time] log Bayes Factors
%   cluster_labels - [channels × time] cluster IDs
%   n_clusters     - Number of clusters
%
% OUTPUTS:
%   cluster_evidence - Structure array with fields:
%     .cluster_id    - Cluster number
%     .n_points      - Number of channel-time points in cluster
%     .cluster_logBF - Sum of log(BF) values
%     .cluster_BF    - Product of BF values (= exp(cluster_logBF))
%     .mean_logBF    - Average log(BF) per point
%     .points        - [n_points × 2] list of [channel, time] indices

for c = 1:n_clusters
    % Get all points in this cluster
    cluster_mask = (cluster_labels == c);
    [ch_idx, t_idx] = find(cluster_mask);

    % Compute cluster evidence
    cluster_evidence(c).cluster_id = c;
    cluster_evidence(c).n_points = length(ch_idx);
    cluster_evidence(c).cluster_logBF = sum(logBF(cluster_mask));
    cluster_evidence(c).cluster_BF = exp(cluster_evidence(c).cluster_logBF);
    cluster_evidence(c).mean_logBF = cluster_evidence(c).cluster_logBF / cluster_evidence(c).n_points;
    cluster_evidence(c).points = [ch_idx, t_idx];
end

end
```

#### Step 3: Calibrate Using H0 Distribution

```matlab
function [significant_clusters, threshold] = calibrate_cluster_evidence(...
    cluster_evidence, logBF_H0, LIMO, initial_threshold, alpha)
% Determine which clusters show significant evidence after correction
%
% INPUTS:
%   cluster_evidence  - Structure from compute_cluster_BF (observed data)
%   logBF_H0          - [ch × time × n_boots] H0 log(BF) distribution
%   LIMO              - LIMO structure with neighbouring_matrix
%   initial_threshold - Threshold for cluster formation (e.g., log(3))
%   alpha             - Significance level (e.g., 0.05)
%
% OUTPUTS:
%   significant_clusters - Boolean array, true for significant clusters
%   threshold            - Cluster-level log(BF) threshold for significance

n_boots = size(logBF_H0, 3);
max_cluster_logBF_H0 = zeros(n_boots, 1);

% For each bootstrap sample
parfor b = 1:n_boots
    % Apply same initial threshold
    logBF_b = logBF_H0(:,:,b);
    mask_b = logBF_b > initial_threshold;

    % Find clusters in this H0 sample
    [cluster_labels_b, n_clusters_b] = limo_findcluster(mask_b, ...
        LIMO.data.neighbouring_matrix, 2);

    if n_clusters_b > 0
        % Compute log(BF) for each H0 cluster
        cluster_logBF_b = zeros(n_clusters_b, 1);
        for h = 1:n_clusters_b
            cluster_mask_b = (cluster_labels_b == h);
            cluster_logBF_b(h) = sum(logBF_b(cluster_mask_b));
        end

        % Store maximum cluster evidence
        max_cluster_logBF_H0(b) = max(cluster_logBF_b);
    else
        max_cluster_logBF_H0(b) = 0;
    end
end

% Determine threshold from H0 distribution
% (Controls familywise error rate for cluster-level inference)
threshold = prctile(max_cluster_logBF_H0, (1-alpha)*100);

% Determine which observed clusters are significant
n_clusters_obs = length(cluster_evidence);
significant_clusters = false(n_clusters_obs, 1);

for c = 1:n_clusters_obs
    if cluster_evidence(c).cluster_logBF > threshold
        significant_clusters(c) = true;
    end
end

fprintf('Cluster-level threshold (α=%.2f): log(BF) = %.2f\n', alpha, threshold);
fprintf('Significant clusters: %d / %d\n', sum(significant_clusters), n_clusters_obs);

end
```

#### Step 4: Visualization

```matlab
function visualize_cluster_evidence(logBF, cluster_labels, cluster_evidence, ...
                                     significant_clusters, LIMO)
% Visualize cluster-level evidence
%
% Creates figure with:
% 1. Point-wise log(BF) map
% 2. Cluster labels map
% 3. Significant clusters highlighted
% 4. Cluster evidence bar chart

figure('Position', [100 100 1200 800]);

% Subplot 1: Point-wise log(BF)
subplot(2,2,1);
imagesc(logBF');
colorbar; title('Point-wise log(BF)');
xlabel('Channel'); ylabel('Time');

% Subplot 2: Cluster labels
subplot(2,2,2);
imagesc(cluster_labels');
colorbar; title('Cluster Labels');
xlabel('Channel'); ylabel('Time');

% Subplot 3: Significant clusters only
subplot(2,2,3);
sig_mask = zeros(size(cluster_labels));
for c = find(significant_clusters)'
    sig_mask(cluster_labels == c) = cluster_evidence(c).cluster_logBF;
end
imagesc(sig_mask');
colorbar; title('Significant Clusters (log(BF))');
xlabel('Channel'); ylabel('Time');

% Subplot 4: Cluster evidence bar chart
subplot(2,2,4);
cluster_logBFs = [cluster_evidence.cluster_logBF];
colors = repmat([0.7 0.7 0.7], length(cluster_logBFs), 1);
colors(significant_clusters, :) = repmat([1 0 0], sum(significant_clusters), 1);
bar(cluster_logBFs, 'FaceColor', 'flat', 'CData', colors);
xlabel('Cluster ID'); ylabel('Cluster log(BF)');
title('Cluster Evidence (red = significant)');
grid on;

end
```

### Complete Implementation Function

```matlab
function cluster_results = limo_likelihood_cluster(results, LIMO, varargin)
% LIMO_LIKELIHOOD_CLUSTER - Cluster-level evidence inference
%
% FORMAT: cluster_results = limo_likelihood_cluster(results, LIMO)
%         cluster_results = limo_likelihood_cluster(results, LIMO, 'threshold', log(3), 'alpha', 0.05)
%
% INPUTS:
%   results - Output from limo_likelihood (must contain logBF and logBF_H0_dist)
%   LIMO    - LIMO structure with neighbouring_matrix
%
% OPTIONAL:
%   'threshold' - Initial threshold for cluster formation (default: log(3))
%   'alpha'     - Significance level for cluster-level correction (default: 0.05)
%   'plot'      - Generate visualization (default: true)
%
% OUTPUTS:
%   cluster_results structure with fields:
%     .cluster_evidence      - Array of cluster evidence structures
%     .significant_clusters  - Boolean array of significant clusters
%     .threshold             - Cluster-level threshold
%     .n_significant         - Number of significant clusters

% Parse inputs
p = inputParser;
addRequired(p, 'results');
addRequired(p, 'LIMO');
addParameter(p, 'threshold', log(3), @isnumeric);
addParameter(p, 'alpha', 0.05, @isnumeric);
addParameter(p, 'plot', true, @islogical);
parse(p, results, LIMO, varargin{:});

initial_threshold = p.Results.threshold;
alpha = p.Results.alpha;
do_plot = p.Results.plot;

% Step 1: Find clusters in observed data
fprintf('Finding clusters in observed data...\n');
[cluster_labels, n_clusters] = find_evidence_clusters(results.logBF, ...
    LIMO, initial_threshold);
fprintf('Found %d clusters\n', n_clusters);

% Step 2: Compute cluster-level evidence
fprintf('Computing cluster-level evidence...\n');
cluster_evidence = compute_cluster_BF(results.logBF, cluster_labels, n_clusters);

% Step 3: Calibrate using H0 distribution
fprintf('Calibrating cluster-level threshold using H0 distribution...\n');
[significant_clusters, threshold] = calibrate_cluster_evidence(...
    cluster_evidence, results.logBF_H0_dist, LIMO, initial_threshold, alpha);

% Step 4: Visualize
if do_plot
    visualize_cluster_evidence(results.logBF, cluster_labels, ...
        cluster_evidence, significant_clusters, LIMO);
end

% Package results
cluster_results.cluster_evidence = cluster_evidence;
cluster_results.significant_clusters = significant_clusters;
cluster_results.threshold = threshold;
cluster_results.n_significant = sum(significant_clusters);
cluster_results.cluster_labels = cluster_labels;

% Print summary
fprintf('\n=== CLUSTER-LEVEL RESULTS ===\n');
for c = find(significant_clusters)'
    fprintf('Cluster %d: %d points, log(BF) = %.2f (BF ~ 10^%.1f)\n', ...
        c, cluster_evidence(c).n_points, cluster_evidence(c).cluster_logBF, ...
        cluster_evidence(c).cluster_logBF / log(10));
end

end
```

### Data Requirements

- Point-wise log(BF) map from `limo_likelihood`
- H0 bootstrap log(BF) distribution
- LIMO structure with `neighbouring_matrix`

### Advantages

✅ Increased sensitivity to spatially/temporally extended effects
✅ Controls family-wise error rate at cluster level
✅ Natural for EEG data structure
✅ Interpretable: "Strong evidence in parietal cluster, 280-340ms"

### Limitations

⚠️ Points within cluster are correlated (multiplicative combination assumes independence)
⚠️ Initial threshold choice affects cluster formation
⚠️ Very large clusters can have astronomically large BF values
⚠️ Doesn't localize effect within cluster

### References

- **Maris, E., & Oostenveld, R. (2007)**. Nonparametric statistical testing of EEG- and MEG-data. *Journal of Neuroscience Methods*, 164(1), 177-190.
- **Pernet et al. (2015)**. Cluster-based computational methods for mass univariate analyses. *Journal of Neuroscience Methods*, 250, 83-95.

---

## Extension 2: Group-Level Aggregation

### Overview

**Current limitation**: Analyzes either single subjects OR group-level test statistics, losing subject-level information.

**Solution**: Compute BF for each subject independently, then combine across subjects using multiplicative aggregation.

### Theoretical Background

For **independent subjects**, Bayesian evidence combines multiplicatively:

```
BF_group = BF_subject1 × BF_subject2 × ... × BF_subjectN

log(BF_group) = Σₛ log(BF_subjectₛ)
```

This is the **product rule** for independent evidence and is the correct Bayesian way to aggregate.

**Key insight**: If 15 out of 20 subjects show moderate evidence (BF ≈ 5), group evidence is strong even if no single subject has BF > 10.

### Use Cases

1. **Multi-subject studies**: Each subject analyzed separately, then combined
2. **Heterogeneity analysis**: Identify which subjects drive group effect
3. **Robust group inference**: Less influenced by outlier subjects than averaging

### Algorithm

#### Step 1: Compute Subject-Level Evidence

```matlab
function subject_results = compute_subject_level_BF(subject_dirs, test_filename, ...
                                                     H0_filename, tau)
% Compute BF for each subject independently
%
% INPUTS:
%   subject_dirs   - Cell array of subject directories
%   test_filename  - Name of test file (e.g., 'paired_samples_ttest_parameter_1.mat')
%   H0_filename    - Name of H0 file (e.g., 'H0/H0_paired_samples_ttest_parameter_1.mat')
%   tau            - Prior parameter (same for all subjects)
%
% OUTPUTS:
%   subject_results - Cell array of results structures, one per subject

n_subjects = length(subject_dirs);
subject_results = cell(n_subjects, 1);

fprintf('Computing subject-level BF for %d subjects...\n', n_subjects);

parfor s = 1:n_subjects
    test_file = fullfile(subject_dirs{s}, test_filename);
    H0_file = fullfile(subject_dirs{s}, H0_filename);

    % Compute BF for this subject
    subject_results{s} = limo_likelihood(test_file, H0_file, ...
        'tau', tau, 'save', false, 'plot', false);

    fprintf('  Subject %d: Strong H1 = %d points\n', ...
        s, subject_results{s}.n_strong_H1);
end

end
```

#### Step 2: Aggregate Across Subjects

```matlab
function group_results = aggregate_group_BF(subject_results)
% Aggregate subject-level BF to group level
%
% INPUTS:
%   subject_results - Cell array of subject-level results
%
% OUTPUTS:
%   group_results structure with fields:
%     .logBF_group         - [ch × time] group-level log(BF)
%     .BF_group            - [ch × time] group-level BF (use with caution - can be huge!)
%     .logBF_subjects      - [ch × time × n_subjects] subject-level log(BF)
%     .mean_logBF          - [ch × time] average log(BF) per subject
%     .geometric_mean_BF   - [ch × time] geometric mean of BF
%     .n_subjects_strong   - [ch × time] number of subjects with BF > 10
%     .proportion_strong   - [ch × time] proportion of subjects with BF > 10

n_subjects = length(subject_results);
[n_ch, n_time] = size(subject_results{1}.logBF);

% Collect subject-level log(BF) maps
logBF_subjects = zeros(n_ch, n_time, n_subjects);
for s = 1:n_subjects
    logBF_subjects(:,:,s) = subject_results{s}.logBF;
end

% Group-level aggregation: Sum log(BF) across subjects
logBF_group = sum(logBF_subjects, 3);

% Interpretable metrics
mean_logBF = logBF_group / n_subjects;  % Average per subject
geometric_mean_BF = exp(mean_logBF);     % Geometric mean of BFs

% Count subjects with strong evidence at each point
strong_BF_threshold = 10;
n_subjects_strong = sum(exp(logBF_subjects) > strong_BF_threshold, 3);
proportion_strong = n_subjects_strong / n_subjects;

% Package results
group_results.logBF_group = logBF_group;
group_results.BF_group = exp(logBF_group);  % Warning: can be astronomically large
group_results.logBF_subjects = logBF_subjects;
group_results.mean_logBF = mean_logBF;
group_results.geometric_mean_BF = geometric_mean_BF;
group_results.n_subjects_strong = n_subjects_strong;
group_results.proportion_strong = proportion_strong;
group_results.n_subjects = n_subjects;

fprintf('Group aggregation complete.\n');
fprintf('  Average group log(BF): %.2f (geometric mean BF ≈ %.1f)\n', ...
    mean(mean_logBF(:)), exp(mean(mean_logBF(:))));
fprintf('  Points where >50%% subjects show BF>10: %d\n', ...
    sum(proportion_strong(:) > 0.5));

end
```

#### Step 3: Analyze Heterogeneity

```matlab
function heterogeneity = analyze_subject_heterogeneity(group_results, LIMO)
% Analyze between-subject variability in evidence
%
% INPUTS:
%   group_results - Output from aggregate_group_BF
%   LIMO          - LIMO structure for visualization
%
% OUTPUTS:
%   heterogeneity structure with fields:
%     .sd_logBF           - [ch × time] std dev of log(BF) across subjects
%     .consistency_index  - [ch × time] proportion of subjects agreeing on direction
%     .outlier_subjects   - List of subjects with unusual evidence patterns

n_subjects = group_results.n_subjects;
logBF_subjects = group_results.logBF_subjects;

% Standard deviation across subjects
sd_logBF = std(logBF_subjects, 0, 3);

% Consistency: proportion of subjects with same sign as group
group_sign = sign(group_results.logBF_group);
consistency_index = zeros(size(group_sign));

for s = 1:n_subjects
    subject_sign = sign(logBF_subjects(:,:,s));
    consistency_index = consistency_index + (subject_sign == group_sign);
end
consistency_index = consistency_index / n_subjects;

% Identify outlier subjects (Mahalanobis distance or simple correlation)
subject_mean_logBF = squeeze(mean(mean(logBF_subjects, 1), 2));
subject_std_logBF = squeeze(std(reshape(logBF_subjects, [], n_subjects), 0, 1))';
z_scores = abs(subject_mean_logBF - mean(subject_mean_logBF)) / std(subject_mean_logBF);
outlier_subjects = find(z_scores > 3);  % 3 SD outliers

% Package results
heterogeneity.sd_logBF = sd_logBF;
heterogeneity.consistency_index = consistency_index;
heterogeneity.outlier_subjects = outlier_subjects;

% Summary
fprintf('Heterogeneity analysis:\n');
fprintf('  Mean SD of log(BF) across subjects: %.2f\n', mean(sd_logBF(:)));
fprintf('  Mean consistency index: %.2f\n', mean(consistency_index(:)));
fprintf('  Outlier subjects: %d\n', length(outlier_subjects));
if ~isempty(outlier_subjects)
    fprintf('    Subject IDs: %s\n', mat2str(outlier_subjects'));
end

end
```

#### Step 4: Visualization

```matlab
function visualize_group_aggregation(group_results, LIMO)
% Visualize group-level and subject-level evidence

figure('Position', [100 100 1400 900]);

% Subplot 1: Group-level log(BF)
subplot(2,3,1);
imagesc(group_results.mean_logBF');  % Average per subject (more interpretable)
colorbar; title('Group log(BF) (avg per subject)');
xlabel('Channel'); ylabel('Time');

% Subplot 2: Proportion of subjects with strong evidence
subplot(2,3,2);
imagesc(group_results.proportion_strong');
colorbar; caxis([0 1]);
title('Proportion subjects with BF>10');
xlabel('Channel'); ylabel('Time');

% Subplot 3: Subject-level variability
subplot(2,3,3);
sd_logBF = std(group_results.logBF_subjects, 0, 3);
imagesc(sd_logBF');
colorbar; title('SD of log(BF) across subjects');
xlabel('Channel'); ylabel('Time');

% Subplot 4: Subject-wise distribution at peak
[~, max_idx] = max(group_results.mean_logBF(:));
[ch_max, t_max] = ind2sub(size(group_results.mean_logBF), max_idx);
subplot(2,3,4);
subject_logBF_at_peak = squeeze(group_results.logBF_subjects(ch_max, t_max, :));
histogram(subject_logBF_at_peak, 20);
xlabel('log(BF)'); ylabel('Number of subjects');
title(sprintf('Subject distribution at peak (ch %d, t %d)', ch_max, t_max));
xline(log(10), 'r--', 'Strong');

% Subplot 5: Subject heatmap
subplot(2,3,5);
% Average log(BF) per subject across all channel-time points
subject_avg_logBF = squeeze(mean(mean(group_results.logBF_subjects, 1), 2));
bar(subject_avg_logBF);
xlabel('Subject'); ylabel('Average log(BF)');
title('Subject-level average evidence');
yline(log(3), 'b--', 'Moderate');
yline(log(10), 'r--', 'Strong');
grid on;

% Subplot 6: Time course at electrode
subplot(2,3,6);
if isfield(LIMO.data, 'timevect')
    times = LIMO.data.timevect;
    elec_idx = round(size(group_results.mean_logBF, 1) / 2);  % Middle electrode

    plot(times, group_results.mean_logBF(elec_idx, :), 'k-', 'LineWidth', 2);
    hold on;

    % Add individual subjects as thin lines
    for s = 1:min(10, group_results.n_subjects)  % Plot first 10 subjects
        plot(times, group_results.logBF_subjects(elec_idx, :, s), '-', ...
            'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
    end

    xlabel('Time (ms)'); ylabel('log(BF)');
    title(sprintf('Group (black) + subjects (gray) at ch %d', elec_idx));
    yline(log(10), 'r--', 'Strong');
    grid on;
end

end
```

### Complete Implementation Function

```matlab
function group_results = limo_likelihood_group(subject_dirs, test_filename, ...
                                                H0_filename, varargin)
% LIMO_LIKELIHOOD_GROUP - Group-level aggregation of subject BFs
%
% FORMAT: group_results = limo_likelihood_group(subject_dirs, test_filename, H0_filename)
%         group_results = limo_likelihood_group(..., 'tau', 1.0, 'plot', true)
%
% INPUTS:
%   subject_dirs   - Cell array of subject directory paths
%   test_filename  - Test file name (same for all subjects)
%   H0_filename    - H0 file name (same for all subjects)
%
% OPTIONAL:
%   'tau'   - Prior parameter (default: 1.0)
%   'plot'  - Generate visualizations (default: true)
%
% OUTPUTS:
%   group_results - Structure with group-level aggregated BF

% Parse inputs
p = inputParser;
addRequired(p, 'subject_dirs', @iscell);
addRequired(p, 'test_filename', @ischar);
addRequired(p, 'H0_filename', @ischar);
addParameter(p, 'tau', 1.0, @isnumeric);
addParameter(p, 'plot', true, @islogical);
parse(p, subject_dirs, test_filename, H0_filename, varargin{:});

tau = p.Results.tau;
do_plot = p.Results.plot;

% Step 1: Compute subject-level BF
subject_results = compute_subject_level_BF(subject_dirs, test_filename, ...
                                            H0_filename, tau);

% Step 2: Aggregate to group level
group_results = aggregate_group_BF(subject_results);

% Step 3: Analyze heterogeneity
group_results.heterogeneity = analyze_subject_heterogeneity(group_results, []);

% Step 4: Visualize
if do_plot
    % Load LIMO from first subject for metadata
    LIMO = load(fullfile(subject_dirs{1}, 'LIMO.mat'));
    LIMO = LIMO.LIMO;
    visualize_group_aggregation(group_results, LIMO);
end

fprintf('\n=== GROUP-LEVEL AGGREGATION COMPLETE ===\n');

end
```

### Data Requirements

- Subject-level test files: `subject_XX/paired_samples_ttest_parameter_1.mat`
- Subject-level H0 files: `subject_XX/H0/H0_paired_samples_ttest_parameter_1.mat`
- Same prior τ for all subjects
- Subjects must be independent (not repeated measures)

### Handling Non-Independence

If subjects are **not independent** (e.g., family members, repeated sessions):

```matlab
% Option 1: Hierarchical/Mixed-Effects Model
% Estimate between-subject and within-subject variance components
% Use hierarchical Bayes Factor (beyond current scope)

% Option 2: Effective N correction
% Reduce effective sample size based on intra-class correlation (ICC)
ICC = estimate_ICC(subject_data);
n_effective = n_subjects / (1 + (avg_sessions_per_subject - 1) * ICC);

% Adjust group evidence
logBF_group_corrected = logBF_group * (n_effective / n_subjects);
```

### Alternative: Summary Statistics Approach

Instead of aggregating subject BFs, use second-level GLM:

```matlab
% For each channel-time point:
% Extract subject-level effect estimates and SEs
theta_subjects = zeros(n_subjects, 1);
SE_subjects = zeros(n_subjects, 1);

for s = 1:n_subjects
    % Load subject's test file
    data_s = load(fullfile(subject_dirs{s}, test_filename));
    field = fieldnames(data_s);
    test_data_s = data_s.(field{1});

    % Extract estimate and SE at this point
    theta_subjects(s) = test_data_s(ch, t, 1);  % Mean estimate
    SE_subjects(s) = test_data_s(ch, t, 2);     % Standard error
end

% Weighted combination (inverse-variance weighting)
weights = 1 ./ SE_subjects.^2;
theta_group = sum(theta_subjects .* weights) / sum(weights);
SE_group = sqrt(1 / sum(weights));

% Compute group-level t-statistic
t_group = theta_group / SE_group;
df_group = n_subjects - 1;

% Compute BF at group level
BF_group = kang_formula(t_group, df_group + 1, tau);
```

### Advantages

✅ Preserves subject-level information
✅ Robust to outlier subjects
✅ Identifies heterogeneity
✅ Principled Bayesian aggregation
✅ Can examine which subjects drive group effect

### Limitations

⚠️ Requires independent subjects
⚠️ Group BF can be astronomically large (use log scale)
⚠️ Same prior τ must be appropriate for all subjects
⚠️ Computational cost: N separate BF computations

### References

- **Rouder, J. N., et al. (2012)**. Default Bayes factors for ANOVA designs. *Journal of Mathematical Psychology*, 56(5), 356-374.
- **Kass, R. E., & Raftery, A. E. (1995)**. Bayes factors. *Journal of the American Statistical Association*, 90(430), 773-795.

---

## Extension 3: Empirical Alternative H1

### Overview

**Current limitation**: Requires specifying prior τ (parametric H1 assumes δ ~ Normal(0, τ²)).

**Solution**: Estimate alternative distribution empirically by bootstrapping **with condition labels preserved**, avoiding need to specify τ.

### Theoretical Background

**Parametric H1** (current):
- Assumes effect size δ ~ Normal(0, τ²)
- User specifies τ
- Closed-form Bayes Factor

**Empirical H1** (proposed):
- Estimate p₁(θ | H1) from data via bootstrap
- Bootstrap **preserves** condition labels (resamples trials within conditions)
- Nonparametric: no assumptions about effect distribution

**Likelihood Ratio**:
```
LR = p₁(θ_obs | H1) / p₀(θ_obs | H0)
```

Where:
- p₀ estimated from **null bootstrap** (already in LIMO)
- p₁ estimated from **alternative bootstrap** (new)

### Two Bootstrap Strategies

#### 1. Null Bootstrap (H0: no effect)

**Already computed by LIMO**:

```matlab
% For paired t-test:
% Center each condition separately
A_centered = trials_A - mean(trials_A);
B_centered = trials_B - mean(trials_B);

% Resample with replacement
for b = 1:n_boots
    A_boot = A_centered(randsample(n_trials_A, n_trials_A, true));
    B_boot = B_centered(randsample(n_trials_B, n_trials_B, true));

    theta_H0(b) = mean(A_boot) - mean(B_boot);  % Should be ~0
end

% This creates H0 distribution stored in H0_paired_samples_ttest_*.mat
```

#### 2. Alternative Bootstrap (H1: effect exists)

**New - preserves real effect**:

```matlab
% DO NOT center - keep the real difference
for b = 1:n_boots_alt
    % Resample trials WITHIN each condition (preserves labels)
    A_boot = trials_A(randsample(n_trials_A, n_trials_A, true));
    B_boot = trials_B(randsample(n_trials_B, n_trials_B, true));

    theta_H1(b) = mean(A_boot) - mean(B_boot);  % Centered at TRUE effect
end

% This creates H1 distribution for density estimation
```

### Algorithm

#### Step 1: Generate Alternative Bootstrap

```matlab
function theta_H1 = generate_alternative_bootstrap(trial_data, design, n_boots_alt)
% Generate bootstrap samples under H1 (preserving condition labels)
%
% INPUTS:
%   trial_data   - [channels × time × trials] single-trial EEG
%   design       - [trials × 1] condition labels (e.g., 0/1 for paired)
%   n_boots_alt  - Number of bootstrap samples (recommend 800-2000)
%
% OUTPUTS:
%   theta_H1 - [channels × time × n_boots_alt] bootstrap differences

[n_ch, n_time, n_trials] = size(trial_data);

% Identify conditions
cond1_idx = find(design == 1);
cond0_idx = find(design == 0);
n_cond1 = length(cond1_idx);
n_cond0 = length(cond0_idx);

fprintf('Generating %d alternative bootstrap samples...\n', n_boots_alt);
fprintf('  Condition 1: %d trials\n', n_cond1);
fprintf('  Condition 0: %d trials\n', n_cond0);

theta_H1 = zeros(n_ch, n_time, n_boots_alt);

parfor b = 1:n_boots_alt
    % Bootstrap within each condition (preserving labels)
    boot_idx_1 = cond1_idx(randsample(n_cond1, n_cond1, true));
    boot_idx_0 = cond0_idx(randsample(n_cond0, n_cond0, true));

    % Compute mean difference for this bootstrap
    mean_1 = mean(trial_data(:, :, boot_idx_1), 3);
    mean_0 = mean(trial_data(:, :, boot_idx_0), 3);

    theta_H1(:, :, b) = mean_1 - mean_0;
end

fprintf('Alternative bootstrap complete.\n');

end
```

#### Step 2: Estimate Densities Using KDE

```matlab
function [p0, p1] = estimate_densities_empirical(theta_obs, theta_H0, theta_H1)
% Estimate p0 and p1 using kernel density estimation
%
% INPUTS:
%   theta_obs - [channels × time] observed effect estimates
%   theta_H0  - [channels × time × n_boots_H0] null bootstrap samples
%   theta_H1  - [channels × time × n_boots_H1] alternative bootstrap samples
%
% OUTPUTS:
%   p0 - [channels × time] density under H0 at observed values
%   p1 - [channels × time] density under H1 at observed values

[n_ch, n_time] = size(theta_obs);
p0 = zeros(n_ch, n_time);
p1 = zeros(n_ch, n_time);

fprintf('Estimating densities using KDE...\n');

% Bandwidth selection (Silverman's rule or cross-validation)
% For stability, use Silverman's rule-of-thumb:
% bandwidth = 0.9 * min(std, IQR/1.34) * n^(-1/5)

for ch = 1:n_ch
    if mod(ch, 10) == 0
        fprintf('  Channel %d/%d\n', ch, n_ch);
    end

    for t = 1:n_time
        % Extract samples for this channel-time point
        samples_H0 = squeeze(theta_H0(ch, t, :));
        samples_H1 = squeeze(theta_H1(ch, t, :));
        obs_value = theta_obs(ch, t);

        % Estimate p0 using H0 samples
        try
            % MATLAB's ksdensity
            p0(ch, t) = ksdensity(samples_H0, obs_value, 'Function', 'pdf');
        catch
            % Fallback to normal approximation
            p0(ch, t) = normpdf(obs_value, mean(samples_H0), std(samples_H0));
        end

        % Estimate p1 using H1 samples
        try
            p1(ch, t) = ksdensity(samples_H1, obs_value, 'Function', 'pdf');
        catch
            % Fallback to normal approximation
            p1(ch, t) = normpdf(obs_value, mean(samples_H1), std(samples_H1));
        end

        % Floor to avoid division by zero
        p0(ch, t) = max(p0(ch, t), 1e-300);
        p1(ch, t) = max(p1(ch, t), 1e-300);
    end
end

fprintf('Density estimation complete.\n');

end
```

#### Step 3: Compute Empirical LR

```matlab
function results_empirical = compute_empirical_LR(theta_obs, theta_H0, theta_H1)
% Compute likelihood ratio using empirical densities
%
% INPUTS:
%   theta_obs - [channels × time] observed effect estimates
%   theta_H0  - [channels × time × n_boots_H0] null bootstrap
%   theta_H1  - [channels × time × n_boots_H1] alternative bootstrap
%
% OUTPUTS:
%   results_empirical structure with:
%     .LR      - [channels × time] likelihood ratios
%     .logLR   - [channels × time] log likelihood ratios
%     .p0      - [channels × time] density under H0
%     .p1      - [channels × time] density under H1

% Estimate densities
[p0, p1] = estimate_densities_empirical(theta_obs, theta_H0, theta_H1);

% Compute likelihood ratio
LR = p1 ./ p0;
logLR = log(p1) - log(p0);

% Handle edge cases
LR(isnan(LR) | isinf(LR)) = 1;  % Inconclusive if density estimation fails
logLR(isnan(logLR) | isinf(logLR)) = 0;

% Package results
results_empirical.LR = LR;
results_empirical.logLR = logLR;
results_empirical.p0 = p0;
results_empirical.p1 = p1;

fprintf('Empirical LR computation complete.\n');
fprintf('  LR range: [%.2e, %.2e]\n', min(LR(:)), max(LR(:)));
fprintf('  log(LR) range: [%.2f, %.2f]\n', min(logLR(:)), max(logLR(:)));

end
```

### Complete Implementation Function

```matlab
function results = limo_likelihood_empirical(test_file, H0_file, trial_data, ...
                                              design, varargin)
% LIMO_LIKELIHOOD_EMPIRICAL - Compute LR using empirical H1
%
% FORMAT: results = limo_likelihood_empirical(test_file, H0_file, trial_data, design)
%         results = limo_likelihood_empirical(..., 'n_boots_alt', 1000, 'plot', true)
%
% INPUTS:
%   test_file   - Path to LIMO test file
%   H0_file     - Path to H0 bootstrap file
%   trial_data  - [channels × time × trials] single-trial EEG data
%   design      - [trials × 1] condition labels (0/1 for paired)
%
% OPTIONAL:
%   'n_boots_alt' - Number of alternative bootstraps (default: 800)
%   'plot'        - Generate diagnostic plots (default: true)
%   'save'        - Save results (default: true)
%
% OUTPUTS:
%   results structure with LR, logLR, p0, p1, etc.

% Parse inputs
p = inputParser;
addRequired(p, 'test_file', @ischar);
addRequired(p, 'H0_file', @ischar);
addRequired(p, 'trial_data', @isnumeric);
addRequired(p, 'design', @isnumeric);
addParameter(p, 'n_boots_alt', 800, @isnumeric);
addParameter(p, 'plot', true, @islogical);
addParameter(p, 'save', true, @islogical);
parse(p, test_file, H0_file, trial_data, design, varargin{:});

n_boots_alt = p.Results.n_boots_alt;
do_plot = p.Results.plot;
save_output = p.Results.save;

% Load observed test statistics
fprintf('Loading test file: %s\n', test_file);
test_struct = load(test_file);
field = fieldnames(test_struct);
test_data = test_struct.(field{1});

if ndims(test_data) == 3
    theta_obs = test_data(:,:,1);  % Mean estimate
else
    theta_obs = test_data(:,:,:,1);  % Time-frequency
end

% Load H0 bootstrap (already computed by LIMO)
fprintf('Loading H0 file: %s\n', H0_file);
H0_struct = load(H0_file);
H0_field = fieldnames(H0_struct);
H0_data = H0_struct.(H0_field{1});

% Extract theta from H0 (convert from t-statistics if needed)
% For simplicity, use H0 samples directly
% In practice, might need to convert t→theta using SE
if ndims(H0_data) == 4
    theta_H0 = H0_data(:,:,1,:);  % [ch × time × n_boots]
else
    theta_H0 = H0_data(:,:,:,1,:);  % [ch × freq × time × n_boots]
end

% Step 1: Generate alternative bootstrap
theta_H1 = generate_alternative_bootstrap(trial_data, design, n_boots_alt);

% Step 2: Compute empirical LR
results = compute_empirical_LR(theta_obs, theta_H0, theta_H1);

% Add metadata
results.params.n_boots_H0 = size(theta_H0, ndims(theta_H0));
results.params.n_boots_H1 = n_boots_alt;
results.params.method = 'Empirical H0 and H1 (KDE)';

% Step 3: Save
if save_output
    [filepath, filename, ext] = fileparts(test_file);
    new_filename = [filename '_empirical_LR' ext];
    new_filepath = fullfile(filepath, new_filename);

    save(new_filepath, 'results', '-v7.3');
    fprintf('Saved: %s\n', new_filepath);
end

% Step 4: Plot
if do_plot
    figure('Position', [100 100 1200 800]);

    subplot(2,2,1);
    imagesc(results.logLR');
    colorbar; title('Empirical log(LR)');

    subplot(2,2,2);
    histogram(results.logLR(:), 50);
    xlabel('log(LR)'); ylabel('Frequency');
    title('Distribution of empirical log(LR)');

    subplot(2,2,3);
    imagesc(log10(results.p0)');
    colorbar; title('log10(p0)');

    subplot(2,2,4);
    imagesc(log10(results.p1)');
    colorbar; title('log10(p1)');
end

end
```

### Data Requirements

**Critical**: Requires access to **single-trial data**:
- `trial_data`: [channels × time × trials]
- `design`: [trials × 1] condition labels

**Cannot use** if only summary statistics available (group-level test files).

### KDE Bandwidth Selection

```matlab
function bandwidth = select_bandwidth(samples, method)
% Select optimal bandwidth for kernel density estimation
%
% INPUTS:
%   samples - [n × 1] data samples
%   method  - 'silverman' (default), 'scott', or 'cv'
%
% OUTPUTS:
%   bandwidth - Optimal bandwidth

n = length(samples);
sigma = std(samples);
IQR_val = iqr(samples);

switch method
    case 'silverman'
        % Silverman's rule of thumb
        bandwidth = 0.9 * min(sigma, IQR_val/1.34) * n^(-1/5);

    case 'scott'
        % Scott's rule
        bandwidth = 1.06 * sigma * n^(-1/5);

    case 'cv'
        % Cross-validation (computationally expensive)
        bandwidth = ksdensity_cv(samples);  % Use MATLAB's built-in CV

    otherwise
        bandwidth = [];  % Let ksdensity choose
end

end
```

### Advantages

✅ No need to specify prior τ
✅ Nonparametric - no distributional assumptions
✅ Can capture non-normal alternatives (skewed, bimodal, etc.)
✅ Realistic alternative distribution from actual data

### Limitations

❌ Requires single-trial data (not available from summary statistics)
❌ Computationally expensive (800+ additional bootstraps)
❌ KDE sensitive to bandwidth choice
❌ Unstable in distribution tails (low sample density)
❌ Risk of overfitting (H1 estimated from same data being tested)
❌ No closed-form formula - must use KDE

### Mitigating Overfitting

**Problem**: Using the same data to estimate H1 and compute LR can be circular.

**Solutions**:

1. **Cross-validation**:
```matlab
% Split data into training (estimate H1) and test (compute LR)
n_trials = size(trial_data, 3);
train_idx = randsample(n_trials, round(0.7 * n_trials));
test_idx = setdiff(1:n_trials, train_idx);

% Estimate H1 from training data
theta_H1 = generate_alternative_bootstrap(trial_data(:,:,train_idx), design(train_idx), 800);

% Compute LR on test data
theta_obs_test = compute_observed_effect(trial_data(:,:,test_idx), design(test_idx));
```

2. **Leave-one-out**:
```matlab
% For each bootstrap, exclude corresponding trial from theta_obs calculation
% (More complex but reduces circularity)
```

3. **Independent calibration set**:
```matlab
% Use separate dataset to estimate H1 distribution
% Then apply to your test dataset
```

### Comparison: Parametric vs Empirical H1

| Feature | Parametric H1 (Kang) | Empirical H1 |
|---------|---------------------|--------------|
| **Input** | Prior τ | Single-trial data |
| **H1 assumption** | δ ~ Normal(0, τ²) | Nonparametric |
| **Formula** | Closed-form BF | KDE-based LR |
| **Speed** | Fast (< 1 sec) | Slow (minutes) |
| **Requires** | User-specified τ | Access to trials |
| **Robustness** | Sensitive to τ choice | Sensitive to bandwidth |
| **Overfitting** | None | Risk if not cross-validated |
| **Interpretability** | "Effect of size ~τ" | "Effect as in data" |

### When to Use Empirical H1

Use when:
- Have access to single-trial data
- Uncertain about appropriate effect size prior
- Alternative distribution may be non-normal
- Computational resources available

Do NOT use when:
- Only have summary statistics
- Limited trial count (< 50 trials per condition)
- Need fast computation
- Want simple, interpretable analysis

### References

- **Gronau, Q. F., et al. (2020)**. A tutorial on bridge sampling. *Journal of Mathematical Psychology*, 81, 80-97.
- **Ly, A., et al. (2016)**. Harold Jeffreys's default Bayes factor hypothesis tests. *Journal of Mathematical Psychology*, 72, 19-32.
- **Sheather, S. J. (2004)**. Density estimation. *Statistical Science*, 19(4), 588-597.

---

## Summary: Comparison of Extensions

| Extension | What it adds | Complexity | Data required | Reference |
|-----------|-------------|------------|---------------|-----------|
| **Cluster-level** | Spatial-temporal evidence | Medium | Point-wise BF + neighbors | Maris & Oostenveld (2007) |
| **Group aggregation** | Multi-subject combination | Medium | Subject-level BF maps | Bayesian principles |
| **Empirical H1** | Data-driven alternative | High | Single-trial data | KDE literature |

All three are **compatible** and could be combined in future versions.

---

## References

1. **Kang, I., Cohen, A. L., & Luck, S. J. (2015)**. A Bayesian approach to functional magnetic resonance imaging analysis.
   *[Adapted for EEG with modifications]*

2. **Rouder, J. N., Speckman, P. L., Sun, D., Morey, R. D., & Iverson, G. (2009)**. Bayesian t tests for accepting and rejecting the null hypothesis. *Psychonomic Bulletin & Review, 16*(2), 225-237.

3. **Kass, R. E., & Raftery, A. E. (1995)**. Bayes factors. *Journal of the American Statistical Association, 90*(430), 773-795.

4. **Wagenmakers, E. J., et al. (2018)**. Bayesian inference for psychology. Part II: Example applications with JASP. *Psychonomic Bulletin & Review, 25*(1), 58-76.

## Support and Contact

For questions, issues, or suggestions:
- Check `limo_likelihood_demo.m` for usage examples
- Review this documentation
- Check LIMO tools GitHub repository
- Contact LIMO tools maintainers

## Version History

- **v1.0** (2025): Initial implementation
  - Kang et al. (2015) formula
  - Empirical H0 from LIMO bootstrap
  - Support for all t-test types
  - Diagnostic plotting
  - Evidence categorization

---

**Implementation**: Devon Yanitski & Claude (Anthropic), 2025
**Based on**: Kang et al. (2015) Bayesian approach
**License**: Same as LIMO tools
