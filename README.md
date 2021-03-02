Project: Computational pipeline for designing guide RNAs for mismatch-CRISPRi
Authors: Jordi van Gestel, John S. Hawkins, Horia Todor, Carol A. Gross
Date: 1 March 2021

Code description:
generate_sgrna.py generates a list of single-guide RNAs (sgRNAs) for CRISPR interference. For each sgRNA, an additional list of single-nucleotide mismatch sgRNAs is provided, with predicted
knockdown efficacy, which can be used to titrate transcriptional interference. generate_sgrna.py requires installation of pandas and biopython modules. As input, the python script requires a GenBank file 
with genomic features as well as a full genome sequence (.gb,.gbk,.gbff). As output, the computational pipeline generates three CSV-files:
(1) sgRNA.csv: full list of all potential sgRNAs that could be used for a gene of interest
(2) sgRNA_filtered.csv: list of filtered sgRNAs, based on a few filter criteria set by the user (in order to remove potentially weak sgRNAs and avoid off-target hits)
(3) sgRNA_mismatched.csv: list of all single-nucleotide mismatch sgRNAs and their predicted knockdown efficacy relative to their full sgRNAs

Resources:
The python script is provided as part of a STAR Protocols publication, which provides a detail overview of the computational pipeline: van Gestel et al., 2021 STAR Protocols
The STAR Protocols publication follows a previous publication in Cell Systems: Hawkins et al., 2020 Cell Systems

Sample codes:
# Get overview of all parameter settings
python generate_sgRNAs.py -h
# Generate sgRNAs for gene of interest (--locu_tag BSU_00010), using GenBank file GCF_000009045.1_ASM904v1_genomic.gbff as input
python generate_sgrnas.py --genbank GCF_000009045.1_ASM904v1_genomic.gbff --locus_tag BSU_00010
# Generate sgRNAs for all genes in genome
python generate_sgrnas.py --genbank GCF_000009045.1_ASM904v1_genomic.gbff --locus_tag all

Parameters:
--genbank, Genbank file name (default: NA)
--locus_tag, Locus tag or list of locus tags for which to generate sgRNAs. Lists should be separated by commas, not space (default: NA)
--model_param, [Optional] Name of csv file with model parameters (default: model_param.csv)
--step, [Optional] Steps to conduct. When conducting multiple steps list them separated by a comma (default: find,filter,mismatch)
--downstream, [Optional] Number of nucleotides downstream of the PAM sequence to be included in the analysis (default: 15)
--off_target_seed, [Optional] Number of nucleotides of the seed sequence of the sgRNA proximal to the PAM to be included for off-target analysis (default, 9)
--file_find, [Optional] Output file name of step 1 (find) (default: sgRNA.csv)
--file_filter, [Optional] Output file name of step 2 (filter) (default: sgRNA_filter.csv)
--file_mismatch, [Optional] Output file name of step 3 (mismatch) (default: sgRNA_mistmatch.csv)
--sgrna_remove, [Optional] Remove all sgRNAs that match a given sequence, e.g. NNNNNNNNNNNNNNNNNNGG, removes all sgRNAs ending with ‘GG’ (Default: None)
--pam_remove, [Optional] Remove all sgRNAs that are associated with a PAM that matches a given sequence (Default: None)
--downstream_remove, [Optional] Remove all sgRNAs that are associated with an upstream sequence that matches a given sequence, e.g. GNNNNNNNNNNNNNN (Default: None)
--gc_lower, [Optional] Minimal GC content of sgRNA (Default: 0)
--gc_upper, [Optional] Maximum GC content of sgRNA (Default: 1)
--offset_upper, [Optional] Maximum distance (nt) from start codon (i.e. offset) (Default: 1000)
--off_target_upper, [Optional] Maximum number of allowed off-target hits of seed sequence of sgRNA (Default: 10)
