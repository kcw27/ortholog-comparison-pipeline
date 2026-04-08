# Metadata Magnet pipeline
<img width="1645" height="337" alt="Metadata Magnet long" src="https://github.com/user-attachments/assets/39370886-ab25-4f8e-a956-8a542dbd2484" />

## Table of Contents
- [Introduction](#introduction)
- [Pipeline diagram](#pipeline-diagram)
- [Nextflow](#stage-2-nextflow)
- [GUI](#stage-4-gui)
- [Dependencies](#dependencies)
- [Acknowledgements](#acknowledgements)

## Introduction
This pipeline facilitates analysis of the ways in which bacterial protein alleles differ across environmental sources.

As of the time of writing, the [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) tool does not include metadata for bacterial proteins. Additionally, [NCBI calculated ortholog datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-ortholog-data-package/) are currently available only for vertebrates and insects. However, by using this pipeline, you can find a set of orthologs for your bacterial protein of interest and collect metadata for these orthologs.

A variety of utilities are provided: you may filter your ortholog dataset in various ways, automatically categorize isolation source metadata into groups such as "host" or "natural environment", generate figures to visualize differences between categories, and run statistical tests on allele properties.

## Pipeline overview
Stages of the pipeline:
1. Obtain a BLAST database of bacterial protein sequences. You may download a prebuilt database such as nr, or use make_blastdb_from_genomes.sh to assemble one from a genome database downloaded using download_databases.sh or directly using [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download).
2. In Nextflow: BLAST for potential orthologs of your bacterial protein sequence of interest, filter the dataset, retrieve the corresponding metadata, and optionally align the sequences. The main outputs of this stage are a **FASTA** and tab-separated **metadata table** for the ortholog dataset.
3. Optionally, further filter the dataset by filtering the metadata. This step may be necessary if you have specific filtering requirements not accounted for by the pipeline. We provide some helper functions in filter_metadata.py, but usually, filtering can be done in the command line. It is not necessary to filter the corresponding FASTA if using the GUI, as the GUI is able to drop sequences lacking metadata. However, if you need the filtered FASTA, use convert_blast_to_fasta.sh to convert the filtered metadata file to FASTA. (If the metadata is a synteny search summary, you will need to reorder columns to match the script's expected format.)
4. In R Shiny GUI: Plug in the **FASTA** and **metadata table** to produce figures and run statistical tests.
<img width="750" height="983" alt="image" src="https://github.com/user-attachments/assets/e0f6993d-b71e-4b17-8bbb-d53fcc2c9ecb" />

## Stage 2: Nextflow
Please refer to [the Nextflow Quickstart page](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/nextflow_quickstart.md) for a tutorial on setting up the input parameter file and running the pipeline.
Currently, Conda is supported for dependency management. We are working on adding Singularity support.

## Stage 4: GUI
### GUI layout
Data input tab:
<img width="2439" height="1351" alt="image" src="https://github.com/user-attachments/assets/f048d26e-cbe4-4e0c-9d0f-10ac6803e79c" />


Figure generation tab:
<img width="3000" height="1747" alt="image" src="https://github.com/user-attachments/assets/6f426a2b-4eb9-4d33-8355-b1c697bb6df1" />

Statistical test tab:
<img width="3000" height="1747" alt="image" src="https://github.com/user-attachments/assets/f0fc154c-7af6-4d6c-8be1-ebb796b3c058" />

### Opening the GUI
If running on a Unix/Linux server, take the following steps to view the GUI:
1. Get your IP address.
   ```bash
   hostname -I
   ```
2. Launch the GUI.
   ```bash
   Rscript "ortholog-comparison-pipeline/scripts/downstream_analysis/R_scripts/pipeline_gui.R"
   ```
3. Replace \<server-ip\> with your IP address, and open the link in a web browser.
   ```text
   http://<server-ip>:3838
   ```

### Example GUI inputs:
(TODO: replace with larger datasets such as mucA and wspF?)  
Use absolute paths. It is okay if there are quotes around the paths entered. If entering the reference sequence ID, please enter it exactly as it appears in the FASTA. 
#### Unaligned FASTA:
```text
"/home/kcw2/ortholog-comparison-pipeline/nextflow/test_results/fasta/not_aligned/results_65_67/results_65_67_filteredSynteny_pident99_qcovs90.fasta"
```
Metadata:
```text
"/home/kcw2/ortholog-comparison-pipeline/nextflow/test_results/metadata/results_65_67_synteny_summary_processedMetadata.blast"
```
Outdir:
```text
/home/kcw2/ortholog-comparison-pipeline/nextflow/test_results/gui_unaligned
```
Reference sequence:
```text
GCA_019434195.1-QYE95824.1-KZ797.12875
```

#### Aligned FASTA:
```text
"/home/kcw2/ortholog-comparison-pipeline/nextflow/test_results/fasta/aligned/results_65_67_filteredSynteny_pident99_qcovs90_aligned.fasta"
```
Metadata:
```text
"/home/kcw2/ortholog-comparison-pipeline/nextflow/test_results/metadata/results_65_67_synteny_summary_processedMetadata.blast"
```
Outdir:
```text
/home/kcw2/ortholog-comparison-pipeline/nextflow/test_results/gui_aligned
```

### GUI tips
If filtering the metadata file after conclusion of the Nextflow stage, you don't need to subset the corresponding FASTA file, but you do need to click the "Apply Subset" button in order to remove the rows with no metadata. (Rows with _NA_ metadata do not show up in the group variable selection menu, but they will be removed if you click "Apply Subset", even if all categories are still selected.)

## Dependencies
It is assumed that the user is running this pipeline on a Unix/Linux server. I have been using Ubuntu 24.04.2 LTS (GNU/Linux 6.8.0-44-generic x86_64).
* Stage 1: If using the scripts in scripts/setup/ to produce a BLAST db, please create the metadata-magnet-setup-env Conda environment from scripts/setup/makeBlastDB.yaml.
* Stage 2: Nextflow handles environment creation automatically, but if you would like to run any of the scripts as standalones, the nextflow/envs/metadata-magnet-env.yaml file will account for most scripts. However, if running synteny search or the subsequent intersection of BLAST with synteny search hits, use nextflow/envs/pynteny-env.yaml instead.
* Stage 4: You will need to install R and the [R Shiny](https://shiny.posit.co/r/getstarted/shiny-basics/lesson1/) package.

## Acknowledgements
* Advisor: Dr. Catherine Armbruster
* Thesis committee:
   * Dr. Irene Kaplow
   * Dr. Phillip Compeau
* Assistance with bioinformatic analysis: Dr. Arkadiy Garber
* Open-source projects used for this pipeline:
   * [bit](https://github.com/AstrobioMike/bit)
   * [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download)
   * [pynteny](https://github.com/Robaina/Pynteny)
   * [seqtk](https://github.com/lh3/seqtk)
# metadata-magnet
