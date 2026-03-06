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
Please clone this repo to your device and install dependencies as required. Refer to help strings in bash scripts for information on options. If bash scripts do not run, try giving them execute permissions and/or prefix them with "bash".   
It is recommended to use a [tmux](https://www.howtogeek.com/671422/how-to-use-tmux-on-linux-and-why-its-better-than-screen/) window for large jobs so that they run uninterrupted. 

**IMPORTANT:**  
This pipeline is designed to run on Linux. It has not been tested on Mac. A HPC environment is recommended if blasting against a large database, e.g. nr.  
Additionally, this pipeline is intended for use with protein sequences obtained through blastp. [Metadata retrieval through NCBI esearch](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/scripts/blast_processing/blast2gen.py) queries the protein database. It has not been tested on outputs of tblastn or blastn.

## Table of contents
- [Obtain ortholog candidates](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#blasting-for-ortholog-candidates)
- [Filter ortholog candidates](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#filter-blast-outputs)
  - [Path 1: basic filtering](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#path-1-using-only-basic-filtering-scripts)
    - [Fetching metadata](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#metadata-retrieval-via-ncbi-esearch)
  - [Path 2: filtering by synteny context](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#path-1-using-only-basic-filtering-scripts)
    - [Intersecting by ID](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#intersection-by-protein-id)
    - [Intersecting by sequence](https://github.com/kcw27/ortholog-comparison-pipeline/edit/main/example_run.md#intersection-by-protein-sequence)
- Categorize isolation sources
- TODO: rest of table of contents

Prior to running the tutorial, please set these paths in the terminal:
```
rundir="/home/kcw2/data/example_run" # replace with your own path; don't include "/" character at the end
dbfasta="${rundir}/toydb.fasta"
scriptdir="/home/kcw2/ortholog-comparison-pipeline/scripts" # also replace with your own path
blastdir="${scriptdir}/blast_processing"
syntenydir="${rundir}/path2/synteny_search_inputs"
```

## BLASTing for ortholog candidates
In this step, you run a local BLAST with your protein of interest to obtain a table (file with .blast extension) of protein sequences. This is a collection of potential orthologs which you may want to filter (TODO: link to the filtering step).

Dependencies:  
[BLAST](https://anaconda.org/bioconda/blast), installed in conda environment **blast_env**.
```bash
conda install bioconda::blast
conda config --append channels bioconda
conda create -n blast_env blast

conda activate blast_env # please activate blast_env prior to running the rest of this step
```

Inputs:
* A local BLAST database, toydb.
* A FASTA containing a single query sequence, PAO1_PA3565.fasta.

Output:
* A .blast file, PA3565_hits.blast. (TODO: explain columns)

1: Prepare BLAST inputs.  
Prepare BLAST database:
I've provided [toydb](https://github.com/kcw27/ortholog-comparison-pipeline/tree/main/example_run/toydb) for this example run. For your own usage of this pipeline, you may want to [install the nr protein database](https://www.ncbi.nlm.nih.gov/books/NBK569850/) or build your own database from a FASTA. Below, I've included the code used to build toydb from a FASTA. 

```bash
mkdir -p "$rundir/toydb"
cd $"$rundir/toydb" # move to directory that will hold the blast db
conda activate blast_env # if not already active
makeblastdb -in "$dbfasta" -input_type fasta -dbtype prot -title toydb -out toydb
```
Prepare query FASTA: Please refer to [this FASTA](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/example_run/PAO1_PA3565.fasta) as an example. The FASTA should contain a single protein sequence. To analyze multiple proteins, it is best to put each sequence in a separate FASTAs and perform separate runs of the pipeline, starting from separate BLASTs.

2: Run the local BLAST.  
To ensure that the output of the BLAST is compatible with the rest of the pipeline, please use [run_local_blast.sh](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/scripts/blast_processing/run_local_blast.sh).
```bash
$blastdir/run_local_blast.sh -p "$rundir/toydb" -d "toydb" -i "$rundir/PAO1_PA3565.fasta" -o "$rundir/PA3565_hits.blast"
# To limit the max number of hits per query sequence, use the -m flag (default=500,000)
```

An R script is provided to assist with summarizing the BLAST outputs. (Dependencies: R version 4.3.3, tidyverse package, glue package.) This script can be run on any BLAST file in which the first 6 columns are consistent with the expected format of "genome_id", "protein_id", "sequence", "evalue", "protein_title", "organism", so it can be used for raw or filtered BLAST files.
```bash
Rscript $blastdir/summarize_blast.R "$rundir/PA3565_hits.blast"
```

## Filter BLAST outputs
There are two paths of filtering: one that includes filtering by synteny block context, and one that doesn't. The primary difference between these two paths is the metadata: if a synteny search is performed, the summary file doubles as a metadata file; otherwise, metadata must be obtained separately via NCBI esearch.  
The metadata files from the two paths are formatted differently. However, the categorization script (TODO: link to categorization section) and GUI (TODO: also link this) are able to handle both metadata formats.

### Path 1: using only basic filtering scripts
Assumption: within a given genome/organism, the sequence with the lowest evalue is the best match for the query sequence.
Scripts available:
* [filter_blast_by_evalue.sh](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/scripts/blast_processing/filter_blast_by_evalue.sh)
* [get_blast_top_hits_by_organism.sh](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/scripts/blast_processing/get_blast_top_hits_by_organism.sh)
* [get_blast_top_hits_by_genomeID.sh](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/scripts/blast_processing/get_blast_top_hits_by_genomeID.sh)

Dependencies: none.

Input:
* Local BLAST output PA3565_hits.blast

Output: 
* PA3565_hits_filtered.blast. Relative to the input, has the same set of columns but a subset of rows. 

You may notice that PA3565_hits.blast has only 0's in column 1, the genome ID column. This means that we need to retrieve genome IDs, among other metadata, via NCBI esearch. To limit the burden placed on NCBI's server, so we reduce the size of the input file as much as possible prior to retrieving metadata. Aside from using the two scripts shown below, you may also want to filter according to other criteria, e.g. exclude records with "partial" in their titles or limit it to records from a particular species.  
Filtering steps that can be done without genome ID:
```bash
mkdir $rundir/path1
cp $rundir/PA3565_hits.blast $rundir/path1 # copy it to path1 dir so the outputs also save in path1

# sequential filtering- first by evalue threshold, then to top per organism (for ties, only take the first per organism)
# note that you can enter multiple files for the -f option, as described in the help strings
$blastdir/filter_blast_by_evalue.sh -t "1e-100" -f $rundir/path1/PA3565_hits.blast # output: $rundir/path1/PA3565_hits_evalueThreshold_1e-100.blast
$blastdir/get_blast_top_hits_by_organism.sh -f $rundir/path1/PA3565_hits_evalueThreshold_1e-100.blast # output: $rundir/path1/PA3565_hits_evalueThreshold_1e-100_topPerOrganism.blast

# Clean up:
rm $rundir/path1/PA3565_hits_evalueThreshold_1e-100.blast # remove intermediate file
mv $rundir/path1/PA3565_hits_evalueThreshold_1e-100_topPerOrganism.blast  $rundir/path1/PA3565_hits_filtered.blast # give filtered file a shorter name
```

#### Metadata retrieval via NCBI esearch
Recommendations: 
* Please follow [NCBI usage guidelines](https://www.nlm.nih.gov/dataguide/eutilities/utilities.html#:~:text=Timing%3A%20Please%20try%20to%20limit%20large%20jobs%20to,PM%20and%205%3A00%20AM%20Eastern%20time%20during%20weekdays.). In particular, avoid submitting large jobs during business hours. Failure to comply may result in NCBI blocking your IP address.  
* The script will fail if the connection to the server is lost. To mitigate this risk, break large input files should be broken into smaller parts. Do not run metadata retrieval jobs in parallel, as this will overburden the server and cause disconnections.
```bash
# Code example for splitting large file into smaller parts:
bigfile="/home/kcw2/data/blast_outputs/mucA/mucA_wt/mucA_paComplete_genomeIDs.blast"
split -d -a 1 -l 8000 $bigfile /home/kcw2/data/blast_outputs/mucA/mucA_wt/mucA_part_ --additional-suffix=.blast
# -d for numeric suffixes, -a 1 to specify that each suffix should be 1 digit long
# after the filename, put the prefix of the file
# --additional-suffix to give it a file extension
```
* Use a tmux window, even for input files of moderate size. The script features built-in delays to comply with NCBI deadlines, but as a result, it takes a long time to run.

Dependencies: Python 3 (I am using 3.13) and the following Python packages:
```python
import sys
import pandas as pd
from Bio import Entrez, SeqIO
import time
import os
```

Input:
* Filtered BLAST file, PA3565_hits_filtered.blast

Outputs:
* Metadata-containing BLAST file, PA3565_hits_filtered_metadata.blast. The columns from the input are preserved, but metadata columns are added to the right.

```bash
conda deactivate # Python wasn't installed in my blast_env environment
python $blastdir/blast2gen.py $rundir/path1/PA3565_hits_filtered.blast $rundir/path1/PA3565_hits_filtered_metadata.blast & # use a & to submit as a job
```


Now that we have genome ID metadata, we can filter to the top hit (lowest evalue) per genome ID. In the metadata retrieved, the assembly_accession column is treated as the genome ID. To set a different column as genome ID, use the -c flag.  
Input:
* PA3565_hits_filtered_metadata_forGenomeIDfilteringdemo.blast (modified to include multiple protein sequences for some genome IDs)

Output:
* PA3565_hits_filtered_metadata_forGenomeIDfilteringdemo_topPerGenome.blast

```bash
$blastdir/get_blast_top_hits_by_genomeID.sh -r -f $rundir/path1/PA3565_hits_filtered_metadata_forGenomeIDfilteringdemo.blast 
```

### Path 2: filtering by synteny context
TODO
(explain why you might want to use this- e.g. same functional groups occur in multiple genes in the same organism, but only want to find hits that correspond to a protein that occurs in a particular operon)  
You _can_ use the basic filtering scripts on the data before putting it through the synteny search.  



(note: the synteny summary example file has gone through categorization so it has a few extra columns at the end)



(Pynteny, kblin genome download (where did I install this?), blast_env seqtk if using intersection by protein sequence- I installed seqtk in blast_env because the pairwise blast uses blast_env) (I should probably split the dependencies by script- for the synteny wrapper, for intersection by ID which uses R and some R packages, for intersection by sequence)
```bash
TODO
```


(download genome dbs)


(prepare inputs for synteny search)
(important note for -i, input_file: if you prepare the text file in Windows, it may have Windows-style return characters, which show up as ^M. This causes various issues- for example, the synteny wrapper will create directories with names that look normal but actually end with return characters.
Prior to running the synteny search, use code like this to remove Windows-style return characters from your synteny input file.
```bash
sed -i 's/\r$//' "/my/file/path/synteny_input.tsv" # replace with your own filepath
```

If you find that there are Windows-style return characters in your filepath after you've already run the synteny search, use code like this to remove the return characters:
```bash
mv "/my/file/path$(printf '\r')/" "/my/file/path/"
```
)


(run synteny wrapper)


Metadata is included in the synteny summary file, so it doesn't need to be obtained separately like in path 1.

The final output of synteny filtering isn't the synteny summary; it's the intersection of a .blast file (a list of sequences that resemble the protein of interest) and the synteny summary (a list of sequences that occur in the synteny block context of interest, but do not necessarily resemble the protein of interest even though they contain HMMs from that protein). There are two methods for intersection as described in the following two sections. The intersection by protein ID script produces a .blast file formatted to convert to FASTA, while the intersection by sequence script directly produces a FASTA.

#### Intersection by protein ID
(faster, but may fail to overlap a large proportion of sequences that should overlap)
(this script has a utility to filter to the top hit per genome. Come to think of it, the intersection by protein sequence script doesn't have that (not that it would really make sense in the context of that script), but you could extract the protein IDs from that output file, use awk to filter .)

#### Intersection by protein sequence
Use this if the naming conventions aren't consistent between your BLAST file and the synteny search output file. This might happen if you BLAST against a custom database where the protein names aren't the same as the names in the GenBank files queried in the synteny search.

Dependencies:
(seqtk is installed in my base conda environment?? in any case, this will fail to make the fastas at the end if you're in pynteny_env)

Inputs:


Outputs:

```bash
```

## Process metadata to prepare for isolation source analysis
(wrapper encompasses rescue and categorization. Can use the same wrapper regardless of whether metadata was obtained through path 1 or path 2)

## Multiple sequence alignment
TODO; I use a MUSCLE alignment through the AliView application because it works quickly and doesn't seem to insert as many gaps as other methods I've tried, such as MAFFT.
I've tried MUSCLE in the command line (muscle -align "in.fasta" -output out.fasta) and it performed significantly worse than MUSCLE in AliView, so I recommend just using AliView.  
AliView > Align > Realign Everything.

## GUI
TODO. Include screenshots.
### Utilities
#### Common utilities

#### Utilities exclusive to unaligned data

#### Utilities exclusive to aligned data

### Unaligned data example run

### Aligned data example run
