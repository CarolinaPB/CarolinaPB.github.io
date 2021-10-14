---
layout: post
title: How to fix smoove errors
tags: [bioinformatics,smoove, svtyper, structural variance]
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

