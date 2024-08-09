# Single_cell_ERK_signal_analysis
This repository contains the codes used to extract and analyse single cell ERK activity of MDCK cells. 

## Single cell ERK activity analysis

The matlab file <ins> *Single_cell_ERK_activity_analysis_pipeline.m* </ins> requires the TrackMate file (.xml) with the single cell tracks for the nuclei as well as a .tif file with the ERK signal and a binary mask of the nuclei. The instructions to get this files can be found in the "methods" section of the manuscript. 

This code is divided in four main sections:

1. Extraction and formating of the individual cell tracks.
2. ERK intensity and individual cell tracks matching.
3. Single cell ERK activation calculation and suitable formating for ARCOS algorithm.
4. Single cell ERK signal kymograph and radial alignment figure generation.

## Collective ERK signal analysis with ARCOS algorithm

The R file <ins> *ARCOS_analysis_pipeline.R* </ins> requires the ARCOS_matrix file (.csv) generated in the previous step. 

This code is divided in two main sections:

1. Data reading, detection of collective events by the ARCOS algortihm and results storage.
2. Results reading, filtering and calculation of median events per hour and unit area and median area of events. 


