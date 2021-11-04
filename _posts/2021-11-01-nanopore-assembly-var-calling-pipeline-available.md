---
layout: post
title: Nanopore assembly and variant calling pipeline is now available
tags: [snakemake, pipeline, bioinformatics, mapping, short reads, variant calling, freebayes, long reads, nanopore, flye, longshot]
---

A pipeline for creating an assembly from nanopore reads. It creates an assembly, does scaffolding with long reads and polishing with short reads. It also does variant calling with short reads and with long reads separately. 
In addition, this pipeline computes busco scores before and after polishing. From this pipeline you'll get statistics files with assembly statistics, as well as statistics from the variant calling VCFs.
There's an optional step to compute a whole genome alignment between the new assembly and one or more assemblies of your choice.   


The pipeline can be found [here](https://github.com/CarolinaPB/nanopore-assembly).  

To use it, first follow the instructions on [running my snakemake pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/).
