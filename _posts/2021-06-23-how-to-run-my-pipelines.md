---
layout: post
title: How to run my pipelines 
tags: [snakemake, pipeline, bioinformatics]
---

Most of my pipelines are created with Snakemake and use modules loaded from the HPC and tools installed with conda.

> Click [here](https://github.com/CarolinaPB/snakemake-template/blob/master/Short%20introduction%20to%20Snakemake.pdf) for an introduction to Snakemake

#### Attention!
_Update 17/08/2021:_
Most of my pipelines depend on modules available in WUR's HPC Anunna. These are loaded when necessary using `module load <module>`. In case you're not using Anunna, but another HPC, these programs might already be available on the HPC, but you might need to change their name or how to load them. In any case, most of the programs should be available to install with conda.   
If you want to install them with conda there are two options:
1. Install in your snakemake environment:
2. Install in a separate environment and specify that environment in the Snakefile with:
```
rule example:
(...)
conda:
  "envs/<name of env>.yaml"
shell:
(...)
```

## Clone the repository
#### From github
Go to the repository's page, click the green "Code" button and copy the path
In your terminal go to where you want to download it to and run
```
git clone <path you copied from github>
```

#### From the the WUR HPC (Anunna)
Go to `/lustre/nobackup/WUR/ABGC/shared/PIPELINES/` and choose which pipeline you want to use. 

```
cp -r <pipeline directory> <directory where you want to save it to>
```

First you'll need to do some set up. Go to the pipeline's directory.

## Installation 

Install `conda` if you don't have it  
_Update 05/01/2022:_  
Here I show how to install miniconda in a linux system  
[Download installer](https://docs.conda.io/en/latest/miniconda.html)  
[Installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)  
1. Download the installer to your home directory. Choose the version according to your operating system. You can right click the link, copy and download with
```
wget <link>
```
At the time of writing this update, for me it would be:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

To install miniconda, run:
```
bash <installer name>
```
installer name could be `Miniconda3-latest-Linux-x86_64.sh`

Set up the conda channels in this order:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
### Create conda environment

```
conda create --name <name-of-pipeline> --file requirements.txt
```
> I recommend giving it the same name as the pipeline


This environment contains snakemake and the other packages that are needed to run the pipeline.

### Activate environment
```
conda activate <name-of-pipeline>
```

### To deactivate the environment (if you want to leave the conda environment)
```
conda deactivate
```

## File configuration
### Create HPC config file

Necessary for snakemake to prepare and send jobs.   

#### Start with creating the directory
```
mkdir -p ~/.config/snakemake/<name-of-pipeline>
cd ~/.config/snakemake/<name-of-pipeline>
```

#### Create config.yaml and include the following:
> My pipelines are configured to work with SLURM

```
jobs: 10
cluster: "sbatch -t 1:0:0 --mem=16000 -c 16 --job-name={rule} --exclude=fat001,fat002,fat101,fat100 --output=logs_slurm/{rule}_%j.out --error=logs_slurm/{rule}_%j.err"

use-conda: true
```

> Here you should configure the resources you want to use.


### Go to the pipeline directory and open config.yaml
Configure your paths, but keep the variable names that are already in the config file.

```
OUTDIR: /path/to/output
READS_DIR: /path/to/reads/ 
ASSEMBLY: /path/to/assembly
PREFIX: <output name>
```
<strike>If you want the results to be written to this directory (not to a new directory), open the Snakefile and comment out  `workdir: config["OUTDIR"]` and ignore or comment out the `OUTDIR: /path/to/output` in the config file.</strike>

_Update 10/08/2021:_ If you want the results to be written to this directory (not to a new directory), comment out `OUTDIR: /path/to/output` in the config.yaml file.

**Now the setup is complete**

## How to run the pipeline

Since the pipelines can take a while to run, it's best if you use a [screen session](https://linuxize.com/post/how-to-use-linux-screen/). By using a screen session, Snakemake stays "active" in the shell while it's running, there's no risk of the connection going down and Snakemake stopping.

Start by creating a screen session:

```
screen -S <name of session>
```
You'll need to activate the conda environment again:
```
conda activate <name-of-pipeline>
```

Then run

```
snakemake -np
```

This will show you the steps and commands that will be executed. Check the commands and file names to see if there's any mistake.

If all looks ok, you can now run your pipeline

```
snakemake --profile <name-of-pipeline>
```

If everything was set up correctly, the jobs should be submitted and you should be able to see the progress of the pipeline in your terminal.
