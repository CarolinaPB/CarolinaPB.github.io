---
layout: post
title: How to fix smoove errors
tags: [bioinformatics,smoove, svtyper, structural variance, SR=0, bwa mem]
---

In this post I will show how to solve some Smoove issues I came across.
I'm using smoove v0.2.8 installed with conda.


## 1. **Duphold error**  
```
fatal.nim(49)            sysFatal
Error: unhandled exception: index -1 not in 0 .. 14130 [IndexDefect]
```

**Solution**
The duphold version in the conda smoove environment is v0.2.1, need to update it to version 0.2.3, which is not available on conda.  
First download duphold from github and check the installation.

```
wget https://github.com/brentp/duphold/releases/download/v0.2.3/duphold
chmod +x ./duphold
./duphold -h # to check if it's installed correctly
```

Add duphold to your path so this version is the version to be used by smoove.
```
PATH=<directory/containing/duphold>:$PATH
```

Since I'm using smoove in a Snakemake pipeline, I add the previous line before my smoove command. 


## 2. **Svtyper error**

```
unrecognized arguments: --max_ci_dist 0
```

**Solution**
The svtyper in the smoove conda environment is quite old, you need to update it.
```
conda install svtyper=0.7.0
```

## 3. **SR=0 for all variants**
When mapping reads with `bwa mem`, if you use the `-M` flag, the split reads are maked as secondary and they will not be used by smoove to get split-read support. This will cause all the variants in your VCF file to have SR=0.
If you don't use split-read support in your analysis, you can run smoove as usual, if you want split-read support values in your structural variant file you can use the following solution.  
This solution was created by Martijn Derks.
> It requires python 2, and the pysam and argparse packages as well as samtools

**Solution**

```
module load samtools

sname=`samtools view -H <sample>.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq`

python bamgroupreads.py -f -M -i <sample>.bam | samblaster --ignoreUnmated -M -a -e -d <sample>.disc.sam -s <smoove call output dir>/<sample>.split.sam -o /dev/null

grep -v "SAMBLASTER" <smoove call output dir>/<sample>.split.sam > $sname.tmp.sam
mv $sname.tmp.sam <smoove call output dir>/<sample>.split.sam
grep -v "SAMBLASTER" <smoove call output dir>/<sample>.disc.sam > $sname.tmp.sam
mv $sname.tmp.sam <smoove call output dir>/<sample>.disc.sam

samtools sort -@ 12 -O bam <smoove call output dir>/<sample>.split.sam > <smoove call output dir>/$sname.split.bam
samtools sort -@ 12 -O bam <smoove call output dir>/<sample>.disc.sam > <smoove call output dir>/$sname.disc.bam

rm <smoove call output dir>/<sample>.split.sam <smoove call output dir>/<sample>.disc.sam
```

The final `split` and `disc` bam files should be in the same directory as the `smoove call` outdir. Smoove will then use these `split` and `disc` bam files for the `smoove call` step.  

For an example on how to use this fix with smoove, see my [Snakemake population level structural variant calling pipeline](https://github.com/CarolinaPB/population-structural-var-calling-smoove)
