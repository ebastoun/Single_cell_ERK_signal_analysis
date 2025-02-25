# Description of exemplary data and output files for using single cell ERK signal analysis codes

## Experimental data

This example includes an 8-hour recording of MDCK cell monolayer, where each frame was taken in regular intervals of 10 minutes (48 frames in total). The recording was performed with a 40x objective in a Nikon Eclipse Ti2-E Inverted Fluorescence Microscope, with calibration factor of 0.28 μm per pixel, with image size of 1608x1608 pixels. However, due to file size limitation in Github, the exemplary data was resized to an image size of 402x402 pixel with calibration factor of 1.12 μm per pixel.



- Images of nuclei for tracking (CFP channel): Example_CFP.tif
- Images of FRET-based ERK activity biosensor to calculate ERK activity maps (FRET to CFP ratio): Example_FRET.tif



## The files generated during the analysis (e.g., following the STAR protocols) are also provided for comparison.

### Calculation of ERK activity maps from microscopy images in Fiji:

- Images of ERK activity map: Example_ERK.tif

### Nuclei tracking and binary image generation in Fiji:
	
- TrackMate file with calculated tracks: Example_CFP.xml
- Segmented nuclei (magenta contours): Original_segmented_nuclei.tif
- Binary mask of nuclei: MASK_Original_segmented_nuclei.tif

### Extraction of track data from xml file and matching of ERK signal with nuclei position in Matlab (code: Single_cell_ERK_activity_analysis_pipeline.m):

- Data extracted from tracks and single cell ERK activity per nucleus and time point: Analysis.mat

### Results of single cell ERK signal analysis in Matlab (code: Single_cell_ERK_activity_analysis_pipeline.m):

- Single cell ERK signal kymograph: Heatmap.eps
- Single cell ERK signal kymograph statistics: Heatmap_statistics.pdf
- Rose plot of alignment angles: Radial_alignment.pdf
- Boxplot of cell speeds: Cell_speeds.pdf
- Boxplot of linear persistence: Cell_directionality.pdf
- Results of cell migration analysis (alignment angles, cell speeds and linear persistence): Migration_analysis.cvs

### Determination of ERK cell state (active/inactive) per nucleus in Matlab (required for ARCOS analysis) (code: Single_cell_ERK_activity_analysis_pipeline.m):

- Matrix with required information for ARCOS analysis: ARCOS_matrix.cvs 
- Video for comparison of cell state: ERK_waves_median.mp4

### Quantification of collective ERK activation events (activated cell clusters) with ARCOS in R (code: ARCOS_analysis_pipeline.R):

- Folder with plots of cell positions color-coded by activation state and overlaid contours for detected cluster per time point: “Frames” folder.
- File with detected collective ERK activation events: arcos_analysis.csv
