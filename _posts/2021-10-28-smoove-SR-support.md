---
layout: post
title: Smoove split read support error
tags: [bioinformatics, smoove, SR, split-reads, bwa, SV, structural variant calling]
comments: true
---

When mapping reads with `bwa mem`, if you use the `-M` flag, the split reads are maked as secondary and they will not be used by smoove to get split-read support. This will cause all the variants in your VCF file to have SR=0.
If you don't use split-read support in your analysis, you can run smoove as usual, if you want split-read support values in your structural variant file you can use the following solution.  
This solution was created by Martijn Derks.
> It requires python 2, and the pysam and argparse packages as well as samtools


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

### Step by step
Get sample name from the sample name read group tag in the bam file. The final files will be named according to the sample name, as expected by smoove.
```
sname=`samtools view -H <sample>.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq`
```

Then small script starts by creating two bam files, one with discordant reads and one with split reads using the python2 script `bamgroupreads.py` and `samblaster`. 
```
python bamgroupreads.py -f -M -i <sample>.bam | samblaster --ignoreUnmated -M -a -e -d <sample>.disc.sam -s <smoove call output dir>/<sample>.split.sam -o /dev/null
```
I tried create sorted bam files from these `split` and `disc` bam files and use them with smoove, but I would always get this error:
```
panic: sam: duplicate program name: line 618: "@PG\tID:SAMBLASTER\tVN:0.1.26\tCL:samblaster -i stdin -o /dev/null -M --acceptDupMarks --excludeDups --ignoreUnmated -d <smoove output dir>/<sample>.disc.sam -s <smoove output dir>/<sample>.split.sam --maxSplitCount 2 --maxUnmappedBases 50 --minIndelSize 50 --minNonOverlap 20"
```

To avoid the previous error, the following steps are needed:
```
grep -v "SAMBLASTER" <smoove call output dir>/<sample>.split.sam > $sname.tmp.sam
mv $sname.tmp.sam <smoove call output dir>/<sample>.split.sam
grep -v "SAMBLASTER" <smoove call output dir>/<sample>.disc.sam > $sname.tmp.sam
mv $sname.tmp.sam <smoove call output dir>/<sample>.disc.sam
```

Then you can create the sorted bam files to be used by smoove
```
samtools sort -@ 12 -O bam <smoove call output dir>/<sample>.split.sam > <smoove call output dir>/$sname.split.bam
samtools sort -@ 12 -O bam <smoove call output dir>/<sample>.disc.sam > <smoove call output dir>/$sname.disc.bam
```

Note that the final `split` and `disc` bams are named according to the sample name present in the read group tag from the bam file.

Remove intermediate files
```
rm <smoove call output dir>/<sample>.split.sam <smoove call output dir>/<sample>.disc.sam
```



For an example on how to use this fix with smoove, see my [Snakemake population level structural variant calling pipeline](https://github.com/CarolinaPB/population-structural-var-calling-smoove)
