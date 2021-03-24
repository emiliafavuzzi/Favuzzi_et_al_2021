# MERFISH
this repository contains scripts that were used in GABA receptive micrglia selectively sculpt developing inhibitory circuits, Emilia Favuzzy et al.

First run groupdsmFISH.m to sum smfish images
then run Main_analysis.m to generate count and intensity matrix for repectvely MERFISH encoded genes and smFISH.
to count the smFISH transcripts for 4 genes as done in the paper, run countsmFISHHIGH.m
to recounts Gabbr1 and Gabbr2 in only the Z that are part of the microglia, run recountwithZ.m
After running clustering, edit ExtractCells2plot.m to choose which cluster to generate maps for, then run makecontrolmaps.m to generate maps of all microglia in this cluster for all control slices
