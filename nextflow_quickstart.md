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
