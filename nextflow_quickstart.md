# Pipeline tutorial
The main Nextflow script for this pipeline is main.nf.

## Example run
```bash
# first, cd to the nextflow subdirectory of this repo
nextflow run main.nf -params-file test-params.yaml -output-dir test_results -with-report test_results/report.html
```

## Your own run
### Prepare inputs 
Answer prompts posed by the input wizard to generate a parameter input file (named my-params.yaml below). It may be helpful to refer to the flow chart in the README for this pipeline.
```bash
python input_wizard.py # in the nexflow subdirectory
```
For help with preparing inputs for synteny search, refer to synteny_search_inputs.md.

If you would like to first run the BLAST, then manually examine the BLAST outputs, and _then_ configure the filtering settings, you should respond "y" when the input wizard presents this question, and rerun the input wizard to produce a second parameters file with the filtering settings. You may then run the pipeline with the new parameters file and a **-resume** flag in order to use cached results from the BLAST.

### Run Nextflow
1. [Install Nextflow.](https://www.nextflow.io/docs/latest/install.html)
2. Clone this repository to your device, and CD to the nextflow subdirectory.
3. Run the following (referring to **nextflow run -h** for additional options as desired):
   ```bash
   nextflow run main.nf -params-file my-params.yaml -output-dir my_results -with-report my_results/report.html
   ```
Execution may take several hours. The most time-consuming steps are the initial BLAST, the metadata retrieval (especially if retrieving via NCBI esearch, as rate-limiting is necessary), and the synteny search.

### Additional parameters not configured by the input wizard
* There is a chance that metadata retrieval via NCBI esearch will fail, especially during times of heavy server traffic. Decrease the **--splitSize** parameter (default: split metadata retrieval tasks into files of at most 8000 lines each, then recombine into a single meetadata file) if NCBI esearch fails.  
* If you have already collected metadata for your dataset via this pipeline, you may use the **--existingMetadata** option to load that metadata file and skip metadata retrieval (unless you attempt to filter to top hit per genome with no genome IDs in the BLAST database, which forces NCBI esearch to retrieve genome accessions among other metadata).  
* If you would like to keep the FASTA for the full dataset of hits from the synteny search (as opposed to the fasta for the intersection of BLAST and synteny search, which is a subset of the synteny search dataset), or you would like to keep the BLAST db created from the BLAST hits (which in this pipeline is only used for the intersection step), include the **--keepSyntenyFasta** flag.
* For metadata categorization, the nextflow/data/category_keywords.txt and nextflow/data/subcategory_keywords.txt are used. You may add to these files if you wish. In each section (delimited by ***), the first line is the name of the category or subcategory, and all following lines are the keywords corresponding to that category or subcategory.

## Pipeline outputs
my_results/  
├── benchmarkingFinal/  
│   └── Text file named after whichever filtering step was the last to run  
├── blastFiles/  
│   ├── blastSummaryFigures/  
│   │   ├── evalue_hist.pdf  
│   │   ├── evalue_seqlen_plot.pdf  
│   │   ├── genomeid_hist.pdf  
│   │   └── seqlen_hist.pdf  
│   ├── intermediates/  
│   │   ├── Shortcuts to intermediate .blast files, including the original output of the BLAST of the query sequence against the BLAST database  
│   └── Shortcut to the final filtered .blast file  
├── fasta/  
│   ├── aligned/  
│   │   ├── **One or more aligned FASTA files.** (Multiple FASTA files expected only if performing multiple synteny searches.) If using the pipeline to align, this was aligned by MUSCLE/super5 depending on dataset size.  
│   └── not_aligned/  
│       ├── **One or more FASTA files.** (Multiple FASTA files expected only if performing multiple synteny searches.)  
├── metadata/  
│   ├── **One or more metadata files** with a .blast extension. (They're essentially .tsv files). If synteny search was performed, these are the synteny summary files. Otherwise, only one file is expected.  
└── report.html if the -with-report option is used  

For the GUI, you will need a **FASTA file** (aligned if you'd like to make sequence logos, or not aligned if you'd like to make sequence length histograms) and a **metadata file**.
