---
layout: post
title: Automatically query Ensembl rapid release homologue page
tags: [bioinformatics, ensembl, homologue, webscraping]
---

Here I share a script to automatically query the ensembl rapid release page for homologues (from a couple of species) of a given list of genes (from one species). 

Usage:
```sh
python ensembl_rapidrelease_homologue_page.py --ids <file_with_ids.txt> --url <rapid release homologue URL> --species1 <species1> --species2 <species2>
```

- ids - text file with one ensembl ID per line. These IDs must correspond to gene IDs from your rapid release annotation.
- url - rapid release homologue page URL, for example, https://rapid.ensembl.org/Meleagris_gallopavo_GCA_905368555.1/Gene/Compara_Homolog?g=
- species1 and species2- name of the species for which you want to find homologues of your genes. This name must be written as in the Ensembl rapid release homologue page.   
  - For example:  `Ruff` or `Green anole` or `Chicken`
  - to add or remove species edit the script

Usage example:
```sh
python ensembl_rapidrelease_homologue_page.py --ids query_genes.txt --url https://rapid.ensembl.org/Meleagris_gallopavo_GCA_905368555.1/Gene/Compara_Homolog?g= --species1 Chicken --species2 Turkey
```

The homologues found correspond to genes in the regular Ensembl release, not in the rapid release.

## Set up
1. Download script `wget https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/ensembl_rapidrelease_homologue_query.py`
2. Download conda environment file `wget https://github.com/CarolinaPB/Bioinfo_scripts/blob/main/envs/ensembl_scraping.yaml`
3. Create conda environment: `conda env create -f ensembl_scraping.yaml`
4. Create Github access token: See instructions [here](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
5. Prepare gene id file, one gene per line
6. Get Ensembl rapid release homologue url: 
    1. Go to https://rapid.ensembl.org/index.html
    2. Find your species and open their page
    3. Click `Example gene`
    4. On the left, click `Homologues`
    5. Copy the URL of this page until `Compara_Homolog?g=`, like this `https://rapid.ensembl.org/Meleagris_gallopavo_GCA_905368555.1/Gene/Compara_Homolog?g=`. 
    6. This will be your input to the `url` parameter.
7. Get name of your two query species  
    1. In the Ensembl `Example gene` page from step 6, look at the homologue table and see if your species are present. If so, copy the names. If not, check for another gene. These will be your input for `species1` and `species2` 

## Run 
```sh
python ensembl_rapidrelease_homologue_page.py --ids <file_with_ids.txt> --url <rapid release homologue URL> --species1 <species1> --species2 <species2>
```

## Output
The output is a file `homologue_table.tsv` with one column per species. The first column is the gene ID from your organism.  
The script also creates a file `processed_genes.txt`, where the IDs of the genes that have been searched are saved. In case you need to restart the script, you can use the same file with genes to be searched. If the gene has already been processed, it will be skipped during the second run. This saves time and resources by avoiding running the same gene several times.
