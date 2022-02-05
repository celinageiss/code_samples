# Code samples

This repository contains code I wrote during some of my former projects at the German Cancer Research Centre (DKFZ) and the University Hospital Heidelberg, structured by programming language:


**R**

- *Data_analysis_sample* (rendered version: [md](https://github.com/celinageiss/code_samples/blob/main/R/Data_anaylsis_sample.md) | source code: [Rmd](https://github.com/celinageiss/code_samples/blob/main/R/Data_anaylsis_sample.Rmd))  
  This analysis pipeline is a part of my current master thesis project, in which I study the effect of the *organellar Ca2+ regulator protein 2 (OCaR2)* during beta-adrenergic induced cardiac arrhythmia in mice. To this end, I apply a linear interaction model on gene expression data and analyse pathways and transcription factor activity based on the model coefficients.

- *plasmid_assembly_app.R* ([here](https://github.com/celinageiss/code_samples/blob/main/R/plasmid_assembly_app.R))  
  This R Shiny app was written to allow an interactive assembly of plasmids for high T cell receptor cloning. The user can choose between multiple building blocks for each plasmid module, filter the library for available pieces and eventually download a fasta file containing the permutations of all selected modules.


**Python**

- *OT-2_miniprep_protocol.py* ([here](https://github.com/celinageiss/code_samples/blob/main/Python/OT2_miniprep_protocol.py))  
  This programme was written for an automation of T cell receptor cloning with the liquid handling robot [OT-2](https://opentrons.com/ot-2/) from Opentrons using their [API](https://docs.opentrons.com/v2/). It takes a 96-well plate and performs a bead-based plasmid purification of the samples.

- *CIViC_annotation.py* ([here](https://github.com/celinageiss/code_samples/blob/main/Python/CIViC_annotation.py))  
  This script takes a tsv table with mutated genes of a cancer patient as an input and annotates hits for the identified variants in the database [CIViC](https://civicdb.org/home) using their API. It was implemented into the  DKFZ's cancer patient analysis pipeline to identify possible tailored therapy options.


**bash**

- *sambamba_subsampling_run.sh* ([here](https://github.com/celinageiss/code_samples/blob/main/bash/sambamba_subsampling_run.sh))  
  This script was used to perform a coverage downsampling of paired tumor-control whole-genome sequencing (WGS) data. The aim of the overall project was to compare the performance of deep whole-exome sequencing and different depths of WGS to find the best cost-effectiveness tradeoff for cancer diagnostics. The run script takes a table with patient IDs and parameters as an input to build a uniform folder structure for each sample and initiates a session on the DKFZ computation cluster by calling a wrapper script.
  
---

Please feel free to also take a look at my [LinkedIn](https://www.linkedin.com/in/celina-geiss/) page and [contact me](mailto:celina.geiss@gmx.de) if you have any questions.
