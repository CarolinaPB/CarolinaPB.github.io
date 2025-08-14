---
layout: post
title: Visualise CNAs from EPI2ME wf-human-variation
tags: [bioinformatics, epi2me, CNV, CNA, nanopore, ONT, visualisation, Rshiny]
---

I didn't find an easy way to visualise CNAs from the EPI2ME wf-human-variation, so I wrote a simple Rshiny app to do it. It reads the output files from the workflow and displays the CNAs in an interactive plot.

The code is available on [GitHub](https://github.com/CarolinaPB/visualise_cnv_wf-human-variation)

Two plots are included:

- Multi-chromosome plot with a small overview of all chromosomes - coverage and CNAs
- Single-chromosome plot with coverage and CNAs for the selected chromosome
These plots show the mean coverage per genomic region (blue line) and the CNA predictions from wf-human-variation in green (gain) or red (loss). The small triangles show the start position of the CNA. In addition, the chromosome cytoband is included below the single chromosome plot. The plots are interactive.

The app uses the `<SAMPLE>.wf_cnv.vcf.gz` and `<SAMPLE>.regions.bed.gz` from the wf-human-variation output.

![App overview](images/app.png)
