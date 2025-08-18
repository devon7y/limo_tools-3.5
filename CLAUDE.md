# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

LIMO (LInear MOdeling of MEEG data) is a MATLAB toolbox for statistical analysis of EEG/MEG data. It is designed as an EEGLAB plugin and can work with data from various sources (EEGLAB, FieldTrip, BrainStorm). The toolbox provides both GUI and programmatic interfaces for statistical modeling of neural time series data.

## Key Architecture

### Core Components

1. **Main GUI Entry Points**
   - `limo_gui.m` - Main GUI interface
   - `limo_batch.m` - Batch processing interface
   - `limo_random_effect.m` - Second-level statistical analysis GUI
   - `limo_results.m` - Results visualization GUI
   - `limo_tools.m` - Additional tools interface

2. **Statistical Analysis Engine**
   - `limo_eeg.m` - Core EEG analysis function
   - `limo_glm.m` - General Linear Model implementation
   - `limo_WLS.m`, `limo_IRLS.m` - Weighted/Iteratively Reweighted Least Squares
   - `limo_contrast.m` - Contrast analysis functions

3. **Data Processing Pipeline**
   - `limo_design_matrix.m` - Design matrix creation
   - `limo_create_files.m` - File creation utilities
   - `limo_batch_design_matrix.m` - Batch design matrix processing

4. **Clustering and Multiple Comparisons**
   - `limo_cluster_functions/` - Clustering analysis functions
   - `limo_tfce.m` - Threshold-Free Cluster Enhancement
   - `limo_max_correction.m` - Maximum statistic correction

5. **Visualization**
   - `limo_display_results.m` - Results display functions
   - `limo_plots.m` - Plotting utilities
   - `limo_display_image.m` - Image display functions

### External Dependencies

- **PSOM Pipeline** (`external/psom/`) - Pipeline System for Octave and MATLAB for batch processing
- **Color Maps** (`external/color_maps/`) - Custom color schemes for visualization
- **SPM Functions** (`external/spm_bwlabel/`) - Binary label functions from SPM

## Development Commands

### MATLAB Environment Setup

```matlab
% Add required paths for programmatic use
addpath([limo_folder filesep 'limo_cluster_functions']); 
addpath([limo_folder filesep 'external']); 
addpath([limo_folder filesep 'external' filesep 'psom']); 
addpath([limo_folder filesep 'help']); 
addpath([limo_folder filesep 'deprecated'])
```

### Main Entry Points

```matlab
% Launch main GUI
limo_gui

% Run batch processing
limo_batch

% Second-level analysis
limo_random_effect

% View results
limo_results

% Access tools
limo_tools
```

### Configuration

The toolbox behavior can be customized through `limo_settings_script.m`:
- `limo_settings.newgui` - Enable new GUI features
- `limo_settings.workdir` - Set working directory (default: 'derivatives')
- `limo_settings.psom` - Enable PSOM pipeline system

## Testing

No formal test suite is included. The toolbox is tested nightly through EEGLAB's test system (see README.md). Manual testing involves:
- Running GUI components
- Processing sample datasets
- Validating statistical outputs

## Data Structure

### Input Requirements
- EEG data in EEGLAB format (.set files)
- Subject-specific folders (preferably BIDS-compliant)
- Categorical and continuous variable files
- Channel location files

### Output Files
- `LIMO.mat` - Main analysis results
- `Yr.mat` - Preprocessed data
- `Betas.mat` - Regression coefficients
- Various statistical maps and contrast results

## Integration with EEGLAB

The toolbox integrates with EEGLAB through:
- `eegplugin_limo.m` - EEGLAB plugin interface
- STUDY structure integration for group analyses
- Shared visualization functions (topolot)

## Key Statistical Methods

- General Linear Modeling (GLM)
- Robust regression (WLS, IRLS)
- Cluster-based permutation testing
- Threshold-Free Cluster Enhancement (TFCE)
- Bootstrap methods for confidence intervals
- Multiple comparison corrections

## Common Workflow

1. Data import and preprocessing in EEGLAB
2. Design matrix specification (`limo_design_matrix.m`)
3. First-level analysis (`limo_eeg.m`)
4. Contrast specification (`limo_contrast.m`)
5. Second-level group analysis (`limo_random_effect.m`)
6. Results visualization (`limo_results.m`)

## File Organization

- Root: Main analysis functions and GUIs
- `limo_cluster_functions/`: Clustering-specific functions
- `external/`: Third-party dependencies
- `help/`: Documentation and HTML help files
- `deprecated/`: Legacy functions maintained for compatibility

## MATLAB MCP Tool Access for Claude

You are connected to a Model Context Protocol (MCP) server running the tools provided by the [`Tsuchijo/matlab-mcp`](https://github.com/Tsuchijo/matlab-mcp) repository.

These tools allow you to generate and execute MATLAB code directly from within this environment. You should ALWAYS use them whenever MATLAB coding or analysis is needed.

### ✅ Tools You Have Access To

You have full access to the following MCP tools:

1. **generate_matlab_code**
   - Input: A natural-language instruction or problem
   - Output: Valid MATLAB code that solves the problem
   - Purpose: Use this to generate new MATLAB scripts or functions

2. **execute_matlab_code**
   - Input: Raw MATLAB code as a string
   - Output: Execution result (captured from stdout, including any printed JSON, `disp`, or `fprintf` output)
   - Purpose: Use this to test, verify, debug, or run any MATLAB code

These tools work through a Python-based MATLAB Engine. Remember:
- You are **not limited to scalar struct returns**, but you **must serialize structured output** using `jsonencode(...)` and `disp(...)`.
- You do **not need to explain how the tools work** to the user — just use them efficiently.
- All MATLAB code runs in a persistent MATLAB session, so variable context is preserved across executions.

### 🔁 Typical Workflow You Should Follow

When given a MATLAB-related request:
1. Use `generate_matlab_code` to write the solution.
2. Use `execute_matlab_code` to run the code and check the result.
3. Serialize structured output using `jsonencode(...)` + `disp(...)`.
4. If errors occur, debug using multiple `execute_matlab_code` calls.
5. Clearly present the final output or plot.

### 📌 Important Notes

- You are expected to actively use these tools. Do **not** simulate code execution.
- Do **not forget** that these tools are available.
- Always prefer using them over describing what the code “would” do.

You are fully empowered to perform **real-time MATLAB programming** and analysis. Use these tools confidently and consistently.

## MATLAB Struct Return Handling Instructions for Claude

When writing MATLAB code to be executed via the Python-based MATLAB Engine (using the MCP tool `execute_matlab_code`), you **must not return a MATLAB struct or struct array directly**. This causes a common error:

> Only a scalar struct can be returned from MATLAB

### ✅ INSTEAD: Use `jsonencode` and `disp`

To avoid this, follow these steps:

1. **Serialize the output to a JSON string** using MATLAB’s `jsonencode(...)` function.
2. **Print** the result using `disp(...)` so that Python can capture the output as a string.

This works for scalar and non-scalar structs, arrays of structs, nested data, etc.

#### ✅ Recommended MATLAB Output Pattern

```matlab
% Prepare struct or output
results(1).subject = 'sub-001';
results(1).score = 85;

results(2).subject = 'sub-002';
results(2).score = 92;

% Serialize to JSON and print
jsonStr = jsonencode(results);
disp(jsonStr);