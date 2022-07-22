# FusionSeeker code
This repository contains the source code to reproduce the simulation and benchmark in FusionSeeker paper of Chen et al. 
To install/use FusionSeeker, please visit https://github.com/Maggi-Chen/FusionSeeker/ 



## simulation.py
This script was used for simulating cancer transcriptome with gene fusions. Normal transcripts are randomly assigned to three groups with low/medium/high expression levels. 300 genes are randomly selected and paired to 150 gene fusions. 100 gene fusions are assigned with both breakpoints in exons, and 50 gene fusions have one breakpoint in exon and 1 breakpoint in intron. 150 fused transcripts are also randomly assigned to low/medium/high expression groups. The output file of this script includes three FASTA files, which can then be used for sequencing read simulation. 


## extract_genefusion_jaffal.py
This script was used for extracting gene fusion calls from jaffal's CSV output. Only 'HighConfidence' and 'LowConfidence' gene fusions are kept for comparison. Gene fusions with "PotentialTransSplicing" tag are filtered. 


## extract_genefusion_longgf.py
This script was used for extracting gene fusion calls from LongGF's output in its log file. 


## test.py
This script was used for benchmarking gene fusion discovery accuracy in both simulated and real datasets. Gene fusion calls from each tool are compared to the ground-truth gene fusion list. Two events are considered as 'match' if the two genes are identical. Number of true-positive calls is used for calculation of recall, precision, and F1 score. 
