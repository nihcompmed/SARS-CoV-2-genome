COVID19 Genome Analysis
=======================
feel free to contact <evancresswell@gmail.com> or <vipulp@niddk.nih.gov > regarding questions on implementation and execution of repo
#### Package Requirements
- Singularity - to run the proper enviornment for code through the singularity container *LADER.simg* (https://hub.docker.com/layers/123851957/evancresswell/erdca/LADER/images/sha256-bf25a591c1682fa1c88eb7cd1b7b9dcbb13848562a0450b5aebfd0b2994e6241?context=explore)
- MAFFT - to align genome sequences

## Table of Contents
- [Genome sequence Alignment](#Align-sequence-alignment)
	- Align full genome
	- extract clades from full alignment
- [Infer interactions using Expectation Reflection](#Infer-interactions-using-Expectation-Reflection)
- [Post-processing and data visualizaion](#Post-processing-and-data-visualization)
	- plotting of interaction maps
	- clade incidence plotting
	- amino acid and nucleotide analysis
- [Results](#Results)
	- Result text files used for paper figures 

# Genome sequence alignment
[Back to Top](#Table-of-Contents)

### Required Files
* **09062020-Cov19gisaid--SeqIdDesc200-aligned-FastaCikti.fa** -- previously aligned fasta file (for testing)
	* downloaded from http://www.dilekbalik.com/SARS-CoV-2_ODOTool/
* **wuhan_ref_fasta** -- covid reference genome: NC_045512.2
* **09062020-Cov19gisaid--SeqIdDesc200-aligned-FastaCikti.fa** -- Full un-aligned fasta file
	* downloaded Oct 20 20202 from https://www.epicov.org/epi3/entities/tmp/tmp_sd_2020_10_20_02_12_qigibl_4j510bafde3/sequences_2020-10-19_07-21.fasta.gz

## Alignment Process
*we want to generate cov_gen_aligned.fasta from our un-aligned fasta file using the Wuhan reference genome*

### Test run

```console
foo@bar:~$ mafft --auto --keeplength --addfragments 09062020-Cov19gisaid--SeqIdDesc200-aligned-FastaCikti.fa wuhan_ref.fasta > 
```

### Full Genome Alignment
- Break up un-aligned fasta files
	- Assumes you have the following files in **/path/to/er_covid19/** 
		- subject genome file: **wuhan_ref.fasta** 
		- aligned genome file: **cov_gen_aligned.fasta** 
		- directory for output: **/path/to/er_covid19/cov_fasta_files/**
```console
foo@bar:~$ singularity exec -B /path/to/er_covid19/biowulf,/path/to/er_covid19/covid_proteins /path/to/er_covid19/LADER.simg python break_up_fasta.py /path/to/er_covid19/ 
```			
	- This creates 
		- broken up fasta files for alignment in *path/to/er_covid19/cov_fasta_files/* 
		- script to align:  **cov_align.swarm** (for cluster computation)

- Run alignment pieces which aligns with mafft (this can be done on a cluster or sequentially).

```console
foo@bar:~$ ./submit_align_swarm.script 
```

	- Individual imulation requirements
		- batches of 15 (lines in **cov_align.swarm**)
		- 10 GB
		- 2 hours
- Finish by concatenating resulting mini-alignments in cov_fasta_files/ (directory created to house mini-alignments) to create full alignment: **covid_genome_full_aligned.fasta**
```console
foo@bar:~$ singularity exec -B /path/to/er_covid19/biowulf,/path/to/er_covid19/covid_proteins /path/to/er_covid19/LADER.simg python concat_fasta.py /path/to/er_covid19/ 
```

### Get Clade Alignments 
- once you have the full genome aligned file you can get clades using get_clades.py
- get_clades.py has subject file hard coded:

``` console
foo@bar:~$  singularity exec -B /path/to/er_covid19/biowulf/,/path/to/er_covid19/covid_proteins /path/to/er_covid19/LADER.simg python get_clades.py /path/to/er_covid19/
```


# Infer interactions using Expectation Reflection
[Back to Top](#Table-of-Contents)
	- Assumes you have the following files in **/path/to/er_covid19/** 
		- subject genome file: **wuhan_ref.fasta** 
		- aligned genome file: **covid_genome_full_aligned.fasta** (created in [Full Genome Alignment][#Full-Genome-Alignment])
		- directory for output: **/path/to/er_covid19/cov_fasta_files/**

- once you have an aligned file you can get DI using ER using run_covGENOME_ER.py
- file generates .pickle files with DI
- the different clade and full sequence run are already defined in run_cov_GENOME_ER file.
	- for all clades see **submit_DI_clade_swarm.script**:
```console
foo@bar:~$ ./submit_DI_clade_swarm.script
```
	- for full genome see **run_covGENOME_ER.py**
		- hardcoded existing full aligned file
```console
foo@bar:~$  singularity exec -B /path/to/er_covid19/biowulf/,/path/to/er_covid19/covid_proteins /path/to/er_covid19/LADER.simg python run_covGENOME_ER.py 
```


# Post-Processing and Data Visualization
[Back to Top](#Table-of-Contents)

## Files Used
* incidence_plotting.py
* gen_incidence_data.py
* make_codon_swarm.py
* make_pair_aa_swarm.py
* codon_mapping.py
* print_pair_counts.py
* covid_aa_bp_variance.py
* plot_covid_genome_CT.py
* ER_covid_PP.py
* expectation_reflection.py

# Results
[Back to Top](#Table-of-Contents)

## Result txt File list:
* s_pairs.txt
* orf7a_pairs.txt
* single_site_bp_aa_freq.txt
* orf7b_pairs.txt
* orf3a_pairs.txt
* orf1ab_pairs.txt
* cov_DI_fullout.txt
* pair_aa_counts.txt
* ORF1ab_bp_vs_aa.txt
* codon_out.txt* 
* cov_DI_G.txt* 
* poster_DI.txt
* cov_DI_S.txt
* cov_DI_GR.txt
* cov_DI_GH.txt
* cov_DI_V.txt
* cov_DI_full1.txt
* cov_DI_full2.txt
* cov_DI_full.txt
* new_er_dis.txt
* covid_sim.txt
