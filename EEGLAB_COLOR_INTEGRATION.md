# LIMO-EEGLAB Color Integration

This document describes the integration of EEGLAB's color scheme with LIMO GUIs.

## Overview

LIMO has been modified to automatically inherit EEGLAB's color scheme, providing a consistent visual experience across the entire EEGLAB ecosystem.

## Color Mapping

The integration uses the following EEGLAB color variables:

| LIMO Element | EEGLAB Variable | Default RGB | Hex Color |
|--------------|-----------------|-------------|-----------|
| GUI Background | `GUIBACKCOLOR` | [0.66 0.76 1.00] | #A5C2FF |
| Text Color | `GUITEXTCOLOR` | [0.00 0.00 0.40] | #000066 |
| Button Color | `GUIPOPBUTTONCOLOR` | [0.93 0.96 1.00] | #EDF5FF |

## Modified Files

### New Files Added:
- `limo_setup_eeglab_colors.m` - Main setup function
- `limo_apply_eeglab_colors_to_figure.m` - Color application utility
- `limo_eeglab_colors.m` - Standalone color application function

### Modified Files:
- `limo_gui.m` - Added automatic color application in opening function
- `limo_results.m` - Added automatic color application in opening function

## Usage

### Automatic Color Application (Recommended)
The modified LIMO GUIs will automatically apply EEGLAB colors when opened:

```matlab
eeglab;           % Start EEGLAB first
limo_gui;         % Opens with EEGLAB colors automatically
limo_results;     % Opens with EEGLAB colors automatically
```

### Manual Color Application
Apply colors to currently open LIMO windows:

```matlab
% Open LIMO GUIs first
limo_gui;
limo_results;

% Then apply colors manually
limo_setup_eeglab_colors();
```

### Apply Colors to Specific Figure
Apply colors to a single figure:

```matlab
% Get figure handle
fig = gcf;

% Load EEGLAB colors
icadefs;

% Apply colors
limo_apply_eeglab_colors_to_figure(fig, GUIBACKCOLOR, GUITEXTCOLOR, GUIPOPBUTTONCOLOR);
```

## Requirements

1. **EEGLAB must be running**: Start EEGLAB before opening LIMO GUIs
2. **MATLAB path**: EEGLAB must be in MATLAB's search path
3. **icadefs.m access**: The EEGLAB color definition file must be accessible

## Troubleshooting

### Colors Not Applied
- Ensure EEGLAB is started: `eeglab`
- Check MATLAB path includes EEGLAB
- Verify `icadefs.m` exists in `eeglab/functions/sigprocfunc/`

### Inconsistent Colors
- Close all LIMO windows and reopen them
- Run `limo_setup_eeglab_colors()` manually
- Restart MATLAB and EEGLAB if problems persist

### Reverting to Original Colors
To disable EEGLAB color integration:

1. Remove the color application code from the GUI opening functions
2. Or comment out the `try/catch` blocks in:
   - `limo_gui.m` (lines ~44-60)
   - `limo_results.m` (lines ~42-54)

## Technical Details

### How It Works
1. When a LIMO GUI opens, the opening function calls `icadefs`
2. EEGLAB color variables are loaded into workspace
3. `limo_apply_eeglab_colors_to_figure()` systematically applies colors to all UI elements
4. Different element types (buttons, text, panels) receive appropriate colors

### Color Application Strategy
- **Backgrounds**: Use `GUIBACKCOLOR` for panels and static elements
- **Text**: Use `GUITEXTCOLOR` for labels and static text
- **Buttons**: Use `GUIPOPBUTTONCOLOR` for interactive elements
- **Edit boxes/Lists**: Keep white backgrounds for readability
- **Plots/Axes**: Preserve original colors to maintain data visualization integrity

## Integration Benefits

1. **Visual Consistency**: LIMO GUIs match EEGLAB's appearance
2. **User Familiarity**: Consistent color scheme across all tools
3. **Professional Look**: Cohesive visual identity
4. **Accessibility**: Uses EEGLAB's tested color combinations
5. **Automatic Updates**: Changes to EEGLAB colors automatically propagate to LIMO

## Future Enhancements

Potential future improvements:
- Theme switching support
- User-customizable color schemes  
- Dark mode compatibility
- High contrast accessibility options
- Color scheme preferences saving

## Support

For issues with color integration:
1. Check that EEGLAB is properly installed and running
2. Verify all files are in the correct locations
3. Try manual color application using `limo_setup_eeglab_colors()`
4. Report persistent issues to the LIMO development team

## Version History

- **v1.0**: Initial EEGLAB color integration
  - Automatic color inheritance
  - Manual setup function
  - Core GUI modifications