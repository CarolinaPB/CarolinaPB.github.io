---
layout: post
title: Nanopore assembly and variant calling pipeline is now available
tags: [bioinformatics, freebayes, parallel, snakemake]
---

While trying to use this [snakemake Freebayes wrapper](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/freebayes.html) I came across the following error:

```
File "/path/to/output/.snakemake/scripts/tmp22a3kwl0.wrapper.py", line 71, in <module>
  "({freebayes} {params} -f {snakemake.input.ref}"
File "/path/to/conda/env/lib/python3.9/site-packages/snakemake/shell.py", line 265, in __new__
  raise sp.CalledProcessError(retcode, cmd)
subprocess.CalledProcessError: Command 'set -euo pipefail;  (freebayes-parallel <(fasta_generate_regions.py mother_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa.PolcaCorrected.fa.fai 100000) 2 --use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 -f mother_oneline.k32.w100.ntLink-arks.longstitch-scaffolds.fa.PolcaCorrected.fa mapped/mother_shortreads.mapped.sorted.bam | bcftools view --output-type v - > tmp_var_calling/mother_shortreads.vcf)  2> logs_slurm/freebayes_mother.log' returned non-zero exit status 1
```

```
[E::bcf_hdr_parse_line] Could not parse the header line: "##contig=<ID=scaffold111,364454,f42Z364454,length=364852>"
[E::bcf_hdr_parse_line] Could not parse the header line: "##contig=<ID=scaffold291,3881,f254Z3881,length=3877>"
[W::vcf_parse] Contig 'scaffold284,4471,f158Z4471' is not defined in the header. (Quick workaround: index the file with tabix.)
[E::bcf_hdr_parse_line] Could not parse the header line: "##contig=<ID=scaffold284,4471,f158Z4471>"
[E::vcf_parse] Could not add dummy header for contig 'scaffold284,4471,f158Z4471'
Error: VCF parse error
Traceback (most recent call last):
  File "/path/to/output/.snakemake/conda/b8612e004025a78079693e175b92b161/bin/vcffirstheader", line 16, in <module>
    print(line.strip())
BrokenPipeError: [Errno 32] Broken pipe
Traceback (most recent call last):
  File "/path/to/output/.snakemake/conda/b8612e004025a78079693e175b92b161/bin/fasta_generate_regions.py", line 32, in <module>
    print(chrom_name + ":" + str(region_start) + "-" + str(end))
BrokenPipeError: [Errno 32] Broken pipe
```

The error seems to come from `vcffirstheader`, which is part of `vcflib`. I could not solve this error, so I decided to implement a simplified version of this wrapper in a normal rule using the Freebayes provided scripts [freebayes-parallel](https://github.com/freebayes/freebayes/blob/master/scripts/freebayes-parallel) and [fasta_generate_regions.py](https://github.com/freebayes/freebayes/blob/master/scripts/fasta_generate_regions.py).

This is the resulting rule, but the same principle can be used in a bash script. It uses the Freebayes (v1.3.1), vcflib (v0.00.2019.07.10), samtools (v1.9) and python 2.7.15 modules available in the WUR's HPC Anunna.
```
rule var_calling_freebayes:
    input:
        ref=rules.polish_polca.output.assembly,
        bam='mapped/{prefix}_shortreads.mapped.sorted.bam',
        indexes='mapped/{prefix}_shortreads.mapped.sorted.bam.bai'
    output:
        "results/variant_calling/{prefix}_shortreads.vcf.gz"
    params:
        chunksize=100000, # reference genome chunk size for parallelization (default: 100000)
        scripts_dir = os.path.join(workflow.basedir, "scripts")
    shell:
        """
module load freebayes bcftools vcflib python/2.7.15 samtools

{params.scripts_dir}/freebayes-parallel.sh <({params.scripts_dir}/fasta_generate_regions.py {input.ref}.fai {params.chunksize) 2 \
-f {input.ref} \
--use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2 \
{input.bam} | vcffilter -f 'QUAL > 20' {input} | bgzip -c > {output}
tabix -p vcf {output}
        """
```
Make sure that your freebayes-parallel and fasta_generate_regions.py scripts are executable, you can do it with `chmod +x <filename>`. These scripts should be in a `scripts` directory in your pipeline directory. You can see how it is implemented in my [nanopore assembly and variant calling pipeline](https://github.com/CarolinaPB/nanopore-assembly).  
In the first line of the command ("{params.scripts_dir}/freebayes-parallel.sh <({params.scripts_dir}/fasta_generate_regions.py {input.ref}.fai {params.chunksize) **2**") the number **2** represents the number of threads to be used. Freebayes is also run with a few options `--use-best-n-alleles 4 --min-base-quality 10 --min-alternate-fraction 0.2 --haplotype-length 0 --ploidy 2 --min-alternate-count 2`. You can leave them, remove or adapt to your needs. The resulting VCF is also filtered for `QUAL > 20`.  
Input and output can be chosen freely.

With this rule you avoid the `vcflib` problem in the [Freebayes snakemake wrapper](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/freebayes.html) and use a simplified version of this wrapper. Using Freebayes in parallel significantly decreases computational time.
