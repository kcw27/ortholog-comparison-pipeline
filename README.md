# Metadata Magnet pipeline
<img width="1645" height="337" alt="Metadata Magnet long" src="https://github.com/user-attachments/assets/39370886-ab25-4f8e-a956-8a542dbd2484" />

## Table of Contents
- [Introduction](#introduction)
- [Pipeline diagram](#pipeline-diagram)
- [GUI](#gui)
- [Dependencies](#dependencies)
- [Example run](#example-run)
- [Acknowledgements](#acknowledgements)

## Introduction
This pipeline facilitates analysis of the ways in which bacterial protein alleles differ across environmental sources.

As of the time of writing, the [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) tool does not include metadata for bacterial proteins. Additionally, [NCBI calculated ortholog datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genes/download-ortholog-data-package/) are currently available only for vertebrates and insects. However, by using this pipeline, you can find a set of orthologs for your bacterial protein of interest and collect metadata for these orthologs.

A variety of utilities are provided: you may filter your ortholog dataset in various ways, automatically categorize isolation source metadata into groups such as "host" or "natural environment", generate figures to visualize differences between categories, and run statistical tests on allele properties.

## Pipeline diagram
<img width="750" height="863" alt="image" src="https://github.com/user-attachments/assets/3a4290eb-8557-43bd-bae3-5224f98338cc" />

Refer to [this page](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/scripts_info.md) for additional information.

## GUI
Data input tab:
<img width="3000" height="1750" alt="image" src="https://github.com/user-attachments/assets/f6f8d807-245c-4ae6-a834-ab7189ff8bba" />

Figure generation tab:
<img width="3000" height="1747" alt="image" src="https://github.com/user-attachments/assets/6f426a2b-4eb9-4d33-8355-b1c697bb6df1" />

Statistical test tab:
<img width="3000" height="1747" alt="image" src="https://github.com/user-attachments/assets/f0fc154c-7af6-4d6c-8be1-ebb796b3c058" />

If running on a Unix server, take the following steps to view the GUI:
1. Get your IP address.
   ```bash
   hostname -I
   ```
2. Launch the GUI.
   ```bash
   Rscript "/home/kcw2/ortholog-comparison-pipeline/scripts/downstream_analysis/R_scripts/pipeline_gui.R"
   ```
3. Replace \<server-ip\> with your IP address, and open the link in a web browser.
   ```text
   http://<server-ip>:3838
   ```

Example inputs to GUI (for debugging):  
### Unaligned FASTA:
```text
"/home/kcw2/data/blast_outputs/pseudomonas_aeruginosa_PA3565_67_synteny_PairwiseBlastIntersected_pident99.fasta"
```
Metadata:
```text
"/home/kcw2/data/results_65_67/synteny_summary.tsv"
```
Outdir:
```text
/home/kcw2/data/testing/foo
```
Reference sequence:
```text
GCF_000006765.1-NP_252255.1
```

### Aligned FASTA:
```text
"/home/kcw2/data/PA3565_align_and_tree_new/alignment/aligned.fasta"
```
Metadata:
```text
"/home/kcw2/data/results_65_67/synteny_summary.tsv"
```
Outdir:
```text
/home/kcw2/data/testing/foo
```
Name map:
```text
"/home/kcw2/data/PA3565_align_and_tree_new/alignment/name_map.tsv"
```

## Dependencies
I am working on getting everything into a single conda environment and setting up a YAML file to facilitate installation. For now, you may simply clone the repository. Note the following main dependencies:
* Command line:
  * [BLAST](https://anaconda.org/bioconda/blast); pipeline assumes the conda environment is called "blast_env"
  * [Pynteny](https://github.com/Robaina/Pynteny); pipeline assumes the conda environment is called "pynteny_env"
  * [seqtk](https://github.com/lh3/seqtk)
  * The pipeline uses ClustalW for alignment and RAxML to produce a phylogenetic tree in [one wrapper script](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/scripts/alignment_and_tree_wrapper.sh). You can use whichever alignment algorithm you prefer, as long as you run fasta_to_phylip.sh on the aligned FASTA before passing it to RAxML. I am considering replacing ClustalW with MAFFT --auto.
* R:
  * reticulate
  * [seqinr](https://www.rdocumentation.org/packages/seqinr/versions/4.2-36)
  * [ggseqlogo](https://github.com/omarwagih/ggseqlogo)
  * [rstatix](https://www.rdocumentation.org/packages/rstatix/versions/0.7.2)
  * shiny
  * pdftools
* Python
* [iTOL](https://itol.embl.de/) (website): you can drag + drop the phylogenetic tree and annotation files to produce a phylogeny annotated by metadata.

## Example run
Refer to [this page](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/example_run.md).

## Acknowledgements
* Advisor: Dr. Catherine Armbruster
* Thesis committee:
   * Dr. Irene Kaplow
   * Dr. Phillip Compeau
* Assistance with bioinformatic analysis: Dr. Arkadiy Garber
