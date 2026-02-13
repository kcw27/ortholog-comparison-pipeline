(tutorial for beta testers; public repo will use wspF and mucA for demos)  
Key:  
- $${\color{red}Red\ text:\ .blast\ table}$$  
  * First five columns: (more columns may be added to the right without affecting compatibility with the pipeline)
      1. Genome ID (may be all zeroes, as is typical of blasting against a protein database)
      2. Sequence ID (may be the protein ID or some combination of assembly accession (genome ID), protein ID, and locus tag)
      3. Protein sequence
      4. Evalue
      5. Title(s) of sequence
  * Original BLAST formatting: -outfmt "6 sallgi sallseqid sseq evalue salltitles"
  * Can [convert this to FASTA](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/scripts/blast_processing/convert_blast_to_fasta.sh)    
 
- $${\color{blue}Blue\ text:\ metadata\ table}$$
  * Obtained through either [NCBI esearch](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/scripts/blast_processing/blast2gen.py) or as a byproduct of [synteny search](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/scripts/synteny_wrapper.sh)
  * In the GUI stage of the pipeline, needs to be joined to the contents of a $${\color{green}FASTA}$$  
 
- $${\color{green}Green\ text:\ FASTA}$$

# Pipeline tutorial
TODO
**IMPORTANT:**  
This pipeline is designed to run on Linux. It has not been tested on Mac. A HPC environment is recommended if blasting against a large database, e.g. nr.  
Additionally, this pipeline is intended for use with protein sequences obtained through blastp. [Metadata retrieval through NCBI esearch](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/scripts/blast_processing/blast2gen.py) queries the protein database. It has not been tested on outputs of tblastn or blastn

## Table of contents
- [Obtain ortholog candidates](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#blasting-for-ortholog-candidates)
- [Filter ortholog candidates](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#filter-blast-outputs)
  - [Path 1: basic filtering](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#path-1-using-only-basic-filtering-scripts)
    - [Fetching metadata](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#metadata-retrieval-via-ncbi-esearch)
  - [Path 2: filtering by synteny context](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#path-1-using-only-basic-filtering-scripts)
    - [Intersecting by ID](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#intersection-by-protein-id)
    - [Intersecting by sequence](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#intersection-by-protein-sequence)
- Categorize isolation sources

## BLASTing for ortholog candidates
(need to use run_local_blast.sh, which does the blast and then some postprocessing to have a file in the right format for the rest of the pipeline)
```bash

```

## Filter BLAST outputs
(difference between the paths is the format of the metadata file, but once you reach the categorization stage it doesn't really matter)


### Path 1: using only basic filtering scripts

#### Metadata retrieval via NCBI esearch

### Path 2: filtering by synteny context
(may or may not be preceded by basic filtering)
(should mention that you _can_ use the basic filtering scripts on the data before putting it through the synteny search)
(Metadata is included in the synteny summary file, so it doesn't need to be obtained separately)
(final output of synteny filtering isn't the synteny summary; it's the intersection. The intersection by protein ID script produces a .blast file formatted to convert to FASTA, while the intersection by sequence script directly produces a FASTA.)

#### Intersection by protein ID
(faster, but may fail to overlap a large proportion of sequences that should overlap)
(this script has a utility to filter to the top hit per genome. Come to think of it, the intersection by protein sequence script doesn't have that (not that it would really make sense in the context of that script), but you could extract the protein IDs from that output file, use awk to filter .)

#### Intersection by protein sequence
Use this if the naming conventions aren't consistent between your BLAST file and the synteny search output file. This might happen if you BLAST against a custom database where the protein names aren't the same as the names in the GenBank files queried in the synteny search.

## Process metadata to prepare for isolation source analysis
(wrapper encompasses rescue and categorization. Can use the same wrapper regardless of whether metadata was obtained through path 1 or path 2)

## Multiple sequence alignment
tbd; I use a MUSCLE alignment through the AliView application because it works quickly and doesn't seem to insert as many gaps as 

## GUI
### Utilities
#### Common utilities

#### Utilities exclusive to unaligned data

#### Utilities exclusive to aligned data

### Unaligned data example run

### Aligned data example run
