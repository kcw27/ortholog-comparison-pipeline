# Synteny search inputs
Synteny-aware filtering is performed using a wrapper for the [Pynteny](https://robaina.github.io/Pynteny/) tool. Please refer to their documentation if anything is unclear. Below, you will find a description of inputs required by scripts related to synteny search. 

## Inputs for synteny search (synteny_wrapper.sh)
This wrapper encompasses run_pynteny.sh and find_synteny_hits.sh.

### genomeDBsynteny / -g genome_db
Example: 'nextflow/data/synteny_search_inputs/genome_db/'  
This is a database of .gbff.gz files (unzipped is fine). If you created a BLAST database from a genome database, you can repurpose that database for synteny search.

### syntenyInput / -i input_file
Example: 'nextflow/data/synteny_search_inputs/synteny_input.tsv'  
This is a tab-separated table in which column 1 contains [synteny structures](https://robaina.github.io/Pynteny/subcommands/search/) as expected by Pynteny, and column 2 contains output directories. If using this as part of Nextflow, please use relative paths.  
synteny_wrapper.sh automatically cleans the file of Windows-style return characters, but if running scripts individiually, please run the following code to clean your file. Otherwise, return characters may become embedded in your filepaths.
```bash
sed -i 's/\r$//' my_file.txt
```

### hmmsList / -L hmms_list
Example: 'nextflow/data/synteny_search_inputs/hmms_of_interest.txt'  
This is a newline-separated list of all HMMs corresponding to your protein of interest. A protein found in synteny search is considered a match for your protein of interest ONLY IF it contains ALL of the HMMs in this list. Do not list any HMMs that you do not expect to see in your protein of interest.

### hmmsDir / -d hmms_dir
Example: 'nextflow/data/synteny_search_inputs/hmm_db/hmms'  
This is a directory containing .HMM files.

### hmmsMetadata / -m hmms_metadata
Example: 'nextflow/data/synteny_search_inputs/hmm_db/hmms_PA3565.tsv'
This is a tab-separated table. The values in column 1 must match the HMM filenames in hmmsDir.

## Inputs for synteny intersection (intersect_blast_and_synteny_pairwiseBlastApproach.sh)
To explain the intersection step: imagine that your protein of interest is expected to be a LTTR that regulates an adjacent oxidoreductase.
* If you BLAST for your protein, you find genes with sequences similar to your protein, but you don't know whether it actually appears next to an oxidoreductase.
* If you perform a synteny search for an LTTR next to an oxidoreductase, you find a set of proteins that appears in the hits for the synteny search and also contain the HMM(s) corresponding to your protein. However, containing an HMM does not guarantee that this protein actually resembles your protein of interest.
* By intersecting these two datasets, we hope to find the proteins among those in the synteny search hits which actually resemble your protein of interest.
<img width="1440" height="771" alt="image" src="https://github.com/user-attachments/assets/33bb87e7-10b7-40c1-9ad4-e2be763b9b65" />  
  
Intersection is performed by creating a BLAST database from the BLAST hits, blasting the synteny search hits against that database, and filtering out the synteny hit sequences that fall below the specified threshold for pident and/or qcovs, as defined in [BLAST documentation](https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a). As such, the output is a subset of the synteny search hits, not of the original BLAST hits.  
We previously attempted to intersect the two datasets by protein names (intersect_blast_and_synteny_byProteinID.R) to find proteins that unambiguously appeared in both datasets, but this approach lead to very low intersection rates, especially when the BLAST dataset had been obtained by blasting against nr. As such, measures of similarity (pident in conjunction with qcovs) are used to determine whether a protein in the synteny search hits can be considered to have an (almost) identical equivalent among the BLAST hits.  

Defaults (numbers may be anything from 0 to 100):
* intersectPident = "99"
* intersectQcovs = "90"
