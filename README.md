# LifespanTrajectoryChimpanzees

Analysis code for "Lifespan trajectory of chimpanzee brains characterized by magnetic resonance imaging histology"

## Overview

This repository contains analysis scripts for processing and comparing human and chimpanzee brain MRI data, with a focus on lifespan trajectories and cortical measurements.

## Repository Structure

### In Vivo Human Data Analysis
- `human_maskingSPM_Script.m` - MPM (Magnetization Parameter Mapping) generation
- `human_PreProcScriptMPM_subjList_hMRI.csh` - Segmentation preparation
- `freesurfer_processing_steps_human.txt` - Segmentation steps
- `human_post_freesurfer_pipeline_human_reg_dat` - Post-segmentation processing

### Main Analysis Scripts
- `chimp_to_human_harmonisation.m` - Harmonizes cortical measurements between chimpanzee and human brain data
- `group_freesurfer_analysis_chimps` - FreeSurfer group analysis for chimpanzee data
- `group_freesurfer_analysis_humans` - FreeSurfer group analysis for human data
- `exclude_projection_outliers_for_group_freesurfer` - Outlier detection and correction in cortical surface projections
- `group_freesurfer_analysis_humans_vs_chimps` - Vertex-wise comparisons between human and converted chimp data
- `group_freesurfer_result_plots_including_humans.m` - Visualization of harmonized cortical measurements
- `microstructural_application_pipeline_attempt` - ROI analysis, intracortical profiles, and developmental modeling

### Functions
The `functions` directory contains various utility functions used by the main analysis scripts.

## Related Resources
Postmortem chimpanzee data analysis code is available in a separate repository:
[postmortembrain-mpm](https://github.com/IlonaLipp/postmortembrain-mpm)

## License
This work is licensed under a [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/).

This means you are free to:
- Share — copy and redistribute the material in any medium or format
- Adapt — remix, transform, and build upon the material for any purpose, even commercially

Under the following terms:
- Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made.

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png

## Citation
to come soon

## Contact
Ilona Lipp, ilonalipp@posteo.net