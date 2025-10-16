# AlleleSeeker (working title)

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
(add screenshot)

## Dependencies
I am working on getting everything into a single conda environment and setting up a YAML file to facilitate installation. For now, you can simply clone the repository. Note the following main dependencies:
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
* Python

## Example run
Refer to [this page](https://github.com/kcw27/ortholog-comparison-pipeline/blob/main/example_run.md).

## Acknowledgements
TODO
