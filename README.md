COVID19 Genome Analysis
=======================


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

	aligned file: cov_gen_aligned.fasta
		- acquired by running: 
			- mafft --auto --keeplength --addfragments 09062020-Cov19gisaid--SeqIdDesc200-aligned-FastaCikti.fa wuhan_ref.fasta > cov_gen_aligned.fasta
		-Files used:
			- 09062020-Cov19gisaid--SeqIdDesc200-aligned-FastaCikti.fa: FILE FROM---> from http://www.dilekbalik.com/SARS-CoV-2_ODOTool/
			- wuhand_ref_fasta: FILE FROM----> covid reference genome: NC_045512.2

			- Downloaded Oct 20 2020
			- https://www.epicov.org/epi3/entities/tmp/tmp_sd_2020_10_20_02_12_qigibl_4j510bafde3/sequences_2020-10-19_07-21.fasta.gz


	Alginment Process	
		- begin by breaking up un-aligned fasta files
			RUN: singularity exec -B /data/cresswellclayec/DCA_ER/biowulf,/data/cresswellclayec/DCA_ER/covid_proteins /data/cresswellclayec/DCA_ER/LADER.simg python break_up_fasta.py #subject_genome_file#.fasta
			- Used  wuhan_ref.fasta as subject 
			 
		- Next run alignment swarm which aligns with mafft.
			RUN: ./submit_align_swarm.script 

		- finally concatenate resulting mini-alignments in cov_fasta_files/
			RUN: singularity exec -B /data/cresswellclayec/DCA_ER/biowulf,/data/cresswellclayec/DCA_ER/covid_proteins /data/cresswellclayec/DCA_ER/LADER.simg python concat_fasta.py 


	Get Clade Alignments 
		- once you have the full genome aligned file you can get clades using get_clades.py
		- get_clades.py has subject file hard coded: 
		RUN: singularity exec -B /data/cresswellclayec/DCA_ER/biowulf/,/data/cresswellclayec/DCA_ER/covid_proteins /data/cresswellclayec/DCA_ER/LADER.simg python get_clades.py





# Infer interactions using Expectation Reflection
[Back to Top](#Table-of-Contents)
- once you have an aligned file you can get DI using ER using run_covGENOME_ER.py
- file generates .pickle files with DI
- the different clade and full sequence run are already defined in run_cov_GENOME_ER file.
	-for all clades:
	RUN: ./submit_DI_clade_swarm.script
	-for full genome run
	- hardcoded existing full aligned file
	RUN: singularity exec -B /data/cresswellclayec/DCA_ER/biowulf/,/data/cresswellclayec/DCA_ER/covid_proteins /data/cresswellclayec/DCA_ER/LADER.simg python run_covGENOME_ER.py 




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
