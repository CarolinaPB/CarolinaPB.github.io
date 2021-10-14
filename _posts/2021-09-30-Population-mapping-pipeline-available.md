---
layout: post
title: Population mapping pipeline is now available 
tags: [snakemake, pipeline, bioinformatics, mapping, bwa, population]
---

The pipeline to map the reads from several samples to a (reference) genome is now available. As input it takes a `.fasta` file of a reference genome, and one or more paths to directories that contain reads. 
The pipeline will go through every subdirectory of the path you gave as input and search for files with `.fq.gz` extension. This means that you can give a path to a directory that has several subdirectories with reads. See example:  

```
/lustre/nobackup/WUR/ABGC/shared/Chicken/Africa/X201SC20031230-Z01-F006_multipath  
├── X201SC20031230-Z01-F006_1  
│   └── raw_data  
│       ├── a109_26_15_1_H  
│       │   ├── a109_26_15_1_H_FDSW202597655-1r_HWFFFDSXY_L3_1.fq.gz  
│       │   ├── a109_26_15_1_H_FDSW202597655-1r_HWFFFDSXY_L3_2.fq.gz  
│       │   └── MD5.txt  
│       └── a20_10_16_1_H  
│           ├── a20_10_16_1_H_FDSW202597566-1r_HWFFFDSXY_L3_1.fq.gz  
│           ├── a20_10_16_1_H_FDSW202597566-1r_HWFFFDSXY_L3_2.fq.gz  
│           └── MD5.txt  
└── X201SC20031230-Z01-F006_2  
    └── raw_data  
        ├── a349_Be_17_1_C  
        │   ├── a349_Be_17_1_C_FDSW202597895-1r_HWFFFDSXY_L3_1.fq.gz  
        │   ├── a349_Be_17_1_C_FDSW202597895-1r_HWFFFDSXY_L3_2.fq.gz  
        │   └── MD5.txt  
        └── a360_Be_05_1_H  
            ├── a360_Be_05_1_H_FDSW202597906-1r_HWFFFDSXY_L3_1.fq.gz  
            ├── a360_Be_05_1_H_FDSW202597906-1r_HWFFFDSXY_L3_2.fq.gz  
            └── MD5.txt  
```
In this case, since the path given was `path1: /lustre/nobackup/WUR/ABGC/shared/Chicken/Africa/X201SC20031230-Z01-F006_multipath`, then all these read files will be used in the pipeline.

The output of this pipeline is the mapped reads per sample - in sorted and indexed `.bam` files and a qualimap report for each sample, as well as a table with a summary of `Mean coverage`, `Mapping rate`, and `Mean mapping quality` for all the samples

You can find the description [here](https://carolinapb.github.io/population-mapping/). There you can also find the tools that the pipeline uses and the output.
The pipeline can be found [here](https://github.com/CarolinaPB/population-mapping).  


To use it, first follow the instructions on [running my snakemake pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/).
