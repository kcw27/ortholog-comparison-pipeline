# Pipeline tutorial
Table of contents: TODO
- [Setup](#setup)
- [Running the Nextflow pipeline](#running-the-nextflow-pipeline)
- [Running the R Shiny GUI on Nextflow pipeline outputs](#running-the-r-shiny-gui-on-nextflow-pipeline-outputs)

# Setup
It is assumed that the user is running this pipeline on a Unix/Linux server. I have been using Ubuntu 24.04.2 LTS (GNU/Linux 6.8.0-44-generic x86_64).  
Clone the repository.
```bash
git clone https://github.com/armbrusterlab/metadata-magnet.git
```
[Install Nextflow](https://www.nextflow.io/docs/latest/install.html) to run the pipeline.  
If you would like to use the GUI to produce figures and statistical analyses from the pipeline outputs, install [R](https://www.anaconda.com/docs/getting-started/working-with-conda/packages/using-r-language) and [R Shiny](https://shiny.posit.co/r/getstarted/shiny-basics/lesson1/). Alternatively, activate the metadata-magnet-env conda environment (explained below).  

Nextflow automatically handles conda environment management. However, to use pipeline scripts as standalone functions, you will need to create and activate conda environments from the corresponding YAML files.
```bash
# to run scripts in the metadata-magnet/scripts/setup/ dir:
conda env create -f metadata-magnet/scripts/setup/metadata-magnet-setup-env.yaml

# to run the following scripts: synteny_wrapper.sh, scripts in metadata-magnet/scripts/synteny_search/ dir
conda env create -f metadata-magnet/nextflow/envs/pynteny-env.yaml

# to run all other scripts:
conda env create -f metadata-magnet/nextflow/envs/metadata-magnet-env.yaml

# example of conda activation
conda activate metadata-magnet-setup-env
```

# Running the Nextflow pipeline

The main Nextflow script for this pipeline is main.nf.

## Example run: quickstart
If you would simply like to confirm that the pipeline works, you may run the following code and compare the outputs against those in the metadata-magnet/nextflow/test_results dir.
```bash
cd metadata-magnet/nextflow
nextflow run main.nf -params-file test-params.yaml -output-dir my_test_results -with-report my_test_results/report.html
```

## Example run: setup, Nextflow, R Shiny GUI
### Prepare BLAST database of bacterial proteins
If you are downloading a pre-made BLAST protein database, such as nr, you can skip this step. However, if you would like to curate a set of NCBI genomes to BLAST, you may use the provided scripts to download the genomes and produce a BLAST database.
#### Download genome database
Downloaded genomes for this example run may be found at metadata-magnet/nextflow/example_data/genome_db, in case you are having trouble with the genome download step.  
We are currently working on implementing a script to extract a list of genome IDs from NCBI assembly summaries. For now, we assume that the starting point is a list of NCBI assembly accessions, like [genomeids.txt](https://github.com/armbrusterlab/metadata-magnet/blob/main/nextflow/example_data/genomeids.txt).
```bash
conda activate metadata-magnet-setup-env

datadir="metadata-magnet/nextflow/example_data"
genomedir="${datadir}/genome_db/"

# separate accessions into genbank and refseq
grep "GCA" "${datadir}/genomeids.txt" > "${datadir}/genomeids_genbank.txt"
grep "GCF" "${datadir}/genomeids.txt" > "${datadir}/genomeids_refseq.txt"

# download genomes- this script hasn't been set up with flags yet, so put the genbank IDs first, refseq IDs second, and output directory third
bash metadata-magnet/scripts/setup/download_databases.sh "${datadir}/genomeids_genbank.txt" "${datadir}/genomeids_refseq.txt" $genomedir
```

The output will look something like this:
```
tree "/home/kcw2/metadata-magnet/nextflow/example_data/genome_db/" # to display the directory contents

/home/kcw2/metadata-magnet/nextflow/example_data/genome_db/
├── GCA_003935435.1
│   ├── GCA_003935435.1_ASM393543v1_genomic.gbff.gz
│   └── MD5SUMS
├── GCA_019909815.1
│   ├── GCA_019909815.1_ASM1990981v1_genomic.gbff.gz
│   └── MD5SUMS
├── GCA_019910225.1
│   ├── GCA_019910225.1_ASM1991022v1_genomic.gbff.gz
│   └── MD5SUMS
```
##### Dealing with failed downloads
On occasion, NCBI's rate-limiting will cause some genome downloads to fail. This may be the case if download_databases.sh prints messages like this: 
```
ERROR: No entry for file ending in '_genomic.gbff.gz'
ERROR: Checksum mismatch for 'metadata-magnet/nextflow/example_data/genome_db/refseq/bacteria/GCF_025798945.1/GCF_025798945.1_ASM2579894v1_genomic.gbff.gz'. Expected 'c042891f194b4d756751ffdc42cb491d', got '116e6fc8d3da898cb71c6e082ab15b38'
```
You may use the following code to reattempt downloads, re-iterating until there are 0 genomes still missing.
```bash
# in metadata-magnet-setup-env

# count missing genomes and write their IDs to a file
missing="${datadir}/still_missing_ids.txt"
> $missing

count=0
for dir in "${genomedir}/"*; do
  genome=$(basename "$dir")

  gb=$(ls $dir | grep ".gbff")
  if [[ $gb == "" ]]; then # missing: the gbff file is absent
    count=$(( $count + 1 ))
    echo $genome >> $missing
  elif [[ $(zcat ${dir%/}/$gb | head -n 10 | grep "Service unavailable") ]]; then # missing: the gbff file states that you've been rate-limited
    count=$(( $count + 1 ))
    echo $genome >> $missing
  fi
done
echo Still missing $count .gbff files

# separate the IDs by genbank and refseq
grep "GCA" "$missing" > "${datadir}/missing_genbank.txt"
grep "GCF" "$missing" > "${datadir}/missing_refseq.txt"

# need to remove the dirs corresponding to the missing data, otherwise download will fail
for id in $(cat $missing); do rm -rf $genomedir$id; done # not a typo; genomedir ends with a "/" so there's no slash between the two variables

# reattempt download for the missing data
bash metadata-magnet/scripts/setup/download_databases.sh "${datadir}/missing_genbank.txt" "${datadir}/missing_refseq.txt" $genomedir

# re-count databases with missing data
count=0
for dir in "${genomedir}/"*; do
  genome=$(basename "$dir")

  gb=$(ls $dir | grep ".gbff")
  if [[ $gb == "" ]]; then
    count=$(( $count + 1 ))
  elif [[ $(zcat ${dir%/}/$gb | head -n 10 | grep "Service unavailable") ]]; then
    count=$(( $count + 1 ))
  fi
done

echo Still missing $count .gbff files  # hopefully this is now 0; if not, repeat from the top
```

#### Build a BLAST database from the genome database
```bash
# in metadata-magnet-setup-env
# this is a slow step, so I recommend using a tmux window
bash metadata-magnet/scripts/setup/make_blastdb_from_genomes.sh -g "$genomedir" -d "${datadir}/toydb" -t "toydb"
```

Examine the headers in the FASTA created by this script- the sequence names are in the form of "genomeID-proteinID". **If creating a BLAST db from this script, you need to use the --separateIDs flag when running the Nextflow pipeline, otherwise it will fail to recognize protein names.** 
```bash
head "/home/kcw2/metadata-magnet/nextflow/example_data/toydb/comprehensive.fasta"
>GCA_003935435.1-RRV22772.1 DksA/TraR family C4-type zinc finger protein_EGJ29_10075_RRV22772.1 [Pseudomonas sp. s199]
MATGWANEGAVQEQIDSTVEDAVQRARSRLGQGESLIHCEECGVRIPEARRKALQGVRLCVSCQAELDKQEASFSGYNRRGSKDSQLR
>GCA_003935435.1-RRV22773.1 hypothetical protein_EGJ29_10080_RRV22773.1 [Pseudomonas sp. s199]
MIKRILLPLSIFLLAAGCSTQTKVLEREVGDFDLKLGTAPTRSMAHGLVEPTTAGAFHGGLDLTHASGWYAGQWSPSAGITNGTSLQVNSYAGFLQQPMDDSLGYELGLIHYDFPELENRDRDGYYAGLNFAGSRLGMALNAAPGRTDSTLFLDLGSVTPFGVGVKVKYGSYALENPHYLPGNRSIEMFNDWSLNISRPWLGIQLDLSYTGSSLTGAECEAYSGQNAQCDALVMFRAERQLY
>GCA_003935435.1-RRV22774.1 oxaloacetate decarboxylase_EGJ29_10085_RRV22774.1 [Pseudomonas sp. s199]
MQRRSHHLLRKDFRQLLASDACYHTASVFDPMSARIAADLGFEVGILGGSVASLQVLAAPDFALITLSEFAEQATRIGRVSRLPIIADADHGYGNALNVMRTVIELERAGISALTIEDTSLPAKFGRKSTDLIGMFEAVGKIRAALEARIDPELCIIARTNAGVIGIEEVIARAQAYQKAGADGICLVGIEDFEQLEQVAAGIEIPLMLVTYGNPRLRDNARLAALGVRIVVNGHAAYFAAIKATYDCLREQRQIDASDLNASQLSVKYSTAGEYMVWAEEFMQVKE
>GCA_003935435.1-RRV22775.1 DMT family transporter_EGJ29_10090_RRV22775.1 [Pseudomonas sp. s199]
MDNRHQSLYGVLLILLSGILLASHDGLSKYLTQLYPVFLVVWARYLAQVVLMLGMFAPRMGRRVFHTLRPWPQLLRGLSLVSVSIMFISGLRYIPLAEATAVIFLTPLMVTVASALLGERVSHSQWLAVGVGLLGVMIIVRPGGALFTPAVLLPFGAAISFTVYQLLTRRLSGTDHPVTSNFLSSLVGFLVMSVLVTFNWRTPSVHDAVLMASLGLMAMSGHLVLTQAFRYASAASLAPFTYGQIVFAGIVGFIAFGHIPDVEAIAGMTVIIASGLCMAYVQSRQASRSA
>GCA_003935435.1-RRV22776.1 class I SAM-dependent methyltransferase_EGJ29_10095_RRV22776.1 [Pseudomonas sp. s199]
MPTPSLTIERIYPQQLDAQDHDDRETLRIHMERYDFAAARLIGTRVLDMACGCGYGSDRLAELNPDKVIVGVDIDPAAIAFAQAHYQRPNLSFICADAEAFSAPGSFDTIVSLETIEHLPRPRALIDNCASLLAKGGQIIASVPITPTLDGNPHHLHDFSKRSFFALFEPHGLLPQQRFEQIQWWQFKGLFRRNATKRHRSEGVGNAVLGYYLKHPGYLFKRLGSMLRHGFSNRYLTCQFRGS
```

The FASTA is no longer needed once the BLAST db has been built, so you can remove it to save space.
```bash
rm ${datadir}/toydb/comprehensive.fasta
```

Because the files for this database are large, I was unable to include them in the repository. I have, however, included them in a zipped form in case you are unable to run make_blastdb_from_genomes.sh.
```bash
# two stages of unzipping
gunzip "metadata-magnet/nextflow/example_data/toydb.zip.gz"
unzip "metadata-magnet/nextflow/example_data/toydb2.zip"
```

### Prepare Nextflow params file
For your convenience, we have included an input wizard to gather pipeline inputs and produce a YAML file which may then be used with Nextflow.  
If you are interested in performing synteny-aware filtering of your ortholog dataset, please refer to [synteny_search_inputs.md](https://github.com/armbrusterlab/metadata-magnet/blob/main/synteny_search_inputs.md).
```bash
python metadata-magnet/nextflow/input_wizard.py
```
The params file used for this example run is [metadata-magnet/nextflow/example_data/input.yaml](https://github.com/armbrusterlab/metadata-magnet/blob/main/nextflow/example_data/input.yaml).
Here are the responses I entered to produce that file:
<img width="2058" height="1278" alt="image" src="https://github.com/user-attachments/assets/68615a90-882d-4f26-b5b9-940fed7671fa" />  
If you would like to first run the BLAST, then manually examine the BLAST outputs, and _then_ configure the filtering settings, you should respond "y" when the input wizard presents this question, and rerun the input wizard to produce a second parameters file with the filtering settings. You may then run the pipeline with the new parameters file and a **-resume** flag in order to use cached results from the BLAST.  
Please provide a real email address so that NCBI can contact you about issues instead of IP-banning your institution.

#### Additional parameters not configured by the input wizard
* There is a chance that metadata retrieval via NCBI esearch will fail, especially during times of heavy server traffic. Decrease the **--splitSize** parameter (default: split metadata retrieval tasks into files of at most 8000 lines each, then recombine into a single meetadata file) if NCBI esearch fails.  
* If you have already collected metadata for your dataset via this pipeline, you may use the **--existingMetadata** option to load that metadata file and skip metadata retrieval (unless you attempt to filter to top hit per genome with no genome IDs in the BLAST database, which forces NCBI esearch to retrieve genome accessions among other metadata).  
* If you would like to keep the FASTA for the full dataset of hits from the synteny search (as opposed to the fasta for the intersection of BLAST and synteny search, which is a subset of the synteny search dataset), or you would like to keep the BLAST db created from the BLAST hits (which in this pipeline is only used for the intersection step), include the **--keepSyntenyFasta** flag.
* For metadata categorization, the nextflow/data/category_keywords.txt and nextflow/data/subcategory_keywords.txt are used. You may add to these files if you wish. In each section (delimited by ***), the first line is the name of the category or subcategory, and all following lines are the keywords corresponding to that category or subcategory.

### Run the Nextflow pipeline
Make sure that Nextflow is installed. You don't need to activate any Conda environments; Nextflow will handle it. Currently, only Conda is supported for dependency management, but we are working on adding Singularity support.
Execution may take several hours. The most time-consuming steps are the initial BLAST, the metadata retrieval (especially if retrieving via NCBI esearch, as rate-limiting is necessary), and the synteny search.
```bash
tmux new -s nextflow # optional but recommended for longer runs
cd metadata-magnet/nextflow/ # need to be in the same directory as main.nf and nextflow.config
mv example_results/ example_results_old/ # outdir files persist even when the pipeline is rerun with the same outdir, so it's recommended to move/delete the old outputs or use a new outdir name

nextflow run main.nf -params-file example_data/input.yaml -output-dir example_results -with-report example_results/report.html
```
Outputs appear in metadata-magnet/nextflow/example_results. Key outputs:
* metadata/ dir
* fasta/ dir, with subdirs for aligned/ (if you choose to align using the pipeline) and not_aligned/ fastas
If synteny search is conducted, there may be multiple files in these dirs depending on how many synteny search outdirs you have specified. Otherwise, there's a single file in each dir.  
The full output tree may be found [later in this document](#pipeline-outputs).

### Optional: additional filtering/processing of metadata
The GUI is capable of filtering by isolation source category, isolation source subcategory, genus, or species. It also adds a column for sequence length if loading unaligned sequences. If you would like to perform additional filtering of the metadata, you may use helper functions such as those in [filter_metadata.py](https://github.com/armbrusterlab/metadata-magnet/blob/main/scripts/downstream_analysis/filter_metadata.py). You may also add columns to the metadata, e.g. using [find_pq_repeats.py](https://github.com/armbrusterlab/metadata-magnet/blob/main/scripts/downstream_analysis/find_pq_repeats.py) (TODO: refactor the script to take CLIs).  
If filtering the metadata file after conclusion of the Nextflow stage, you don't need to subset the corresponding FASTA file, as the GUI is able to drop sequences lacking metadata. However, when you load the data in the GUI, you do need to click the "Apply Subset" button in order to remove the rows with no metadata. (Rows with _NA_ metadata do not show up in the group variable selection menu, but they will be removed if you click "Apply Subset", even if all categories are still selected.) If you need the filtered FASTA, use convert_blast_to_fasta.sh to convert the filtered metadata file to FASTA. (If the metadata is a synteny search summary, you will need to reorder columns to match the script's expected format.)

### Analyze Nextflow outputs in the GUI
#### Open the GUI
If running on a Unix/Linux server, take the following steps to view the GUI:
1. Get your IP address.
   ```bash
   hostname -I
   ```
2. Launch the GUI in the terminal.
   ```bash
   conda activate metadata-magnet-env
   Rscript metadata-magnet/scripts/downstream_analysis/R_scripts/pipeline_gui.R
   ```
3. Replace \<server-ip\> with your IP address, and open the link in a web browser.
   ```text
   http://<server-ip>:3838
   ```  

#### Load data
<img width="582" height="709" alt="image" src="https://github.com/user-attachments/assets/d745cf29-eb65-475f-9d61-0163e42a5d8c" />  

Use absolute paths. Replace with the paths to the respective files on your own device. It is okay for the paths and output directory to be surrounded by quotes.   
The reference sequence ID, if you provide it, must match exactly with the sequence's name in the FASTA. Note that if you filter out the reference sequence from your dataset, it won't appear in figures.

##### Unaligned data
FASTA:
```
"/home/kcw2/metadata-magnet/nextflow/example_results/fasta/not_aligned/negBackward/negBackward_filteredSynteny_pident99_qcovs90.fasta"
```
Metadata:
```
"/home/kcw2/metadata-magnet/nextflow/example_results/metadata/negBackward_synteny_summary_processedMetadata.blast"
```
Outdir:
```
/home/kcw2/metadata-magnet/nextflow/example_results/gui/unaligned
```
Reference sequence ID (arbitrarily picking the first one in the FASTA):
```
GCF_021060745.1-WP_003092534.1-K4W69.RS10915
```

##### Aligned data
FASTA:
```
"/home/kcw2/metadata-magnet/nextflow/example_results/fasta/aligned/negBackward_filteredSynteny_pident99_qcovs90_aligned.fasta"
```
Metadata:
```
"/home/kcw2/metadata-magnet/nextflow/example_results/metadata/negBackward_synteny_summary_processedMetadata.blast"
```
Outdir:
```
/home/kcw2/metadata-magnet/nextflow/example_results/gui/aligned
```
(ignore the name map file field; it was for a feature that hasn't yet been implemented)  

Reference sequence ID (arbitrarily picking the first one in the FASTA):
```
GCF_021060745.1-WP_003092534.1-K4W69.RS10915
```

It may take a few seconds to load the data, especially if it's a large dataset.  
After loading data: you can subset the data according to the following types of metadata: isolation source category, isolation source subcategory, genus, species. Click the "Apply Subset" button after checking or unchecking buttons. If you want to export the metadata table at the current level of filtering, click "Save current data to TSV". To undo all subsetting and return to the original file, click "Clear Subset".
<img width="2404" height="1247" alt="image" src="https://github.com/user-attachments/assets/9e5ba740-deb1-44f7-a592-267f3af1f397" />  
I will filter out the "no category" records, and filter to only "aeruginosa" species. Filtering can be stacked in this manner. Make sure to hit "Apply subset" before navigating away from the current group variable. 
<img width="1015" height="58" alt="image" src="https://github.com/user-attachments/assets/69ecfa91-4410-4da8-aca8-80dba267c8c8" />  

#### Visualize the data
##### Histograms
You can make histograms on any numeric variable, such as the sequence_length column which is created if loading unaligned sequences. Histograms are produced on several levels: a histogram for the entire dataset, a set of histograms for each group variable category, and a set of histograms for each subcategory within the category. If the selected group variable is (isolation source) category, then the subcategory will be (isolation source) subcategory. If the selected group variable is genus, then the subcategory will be species.  
If the reference sequence appears in the current subset of the data, its length will also be marked on the histograms.  
You can select histogram bin widths. Try adjusting this if no bins appear in your histogram.  
Consistent Y-axis across plots: you may want to uncheck this if one category has a much larger number of observations than the others, so that the axes will be scaled appropriately for individual histograms.  
Histograms for each isolation source category:  
<img width="2439" height="1461" alt="image" src="https://github.com/user-attachments/assets/458213a1-3fb5-4194-accb-eff1ee8235f9" />  
Histograms for subcategories within a category (in this case, "natural environment" category, but such plots are made for all categories):  
<img width="1384" height="699" alt="image" src="https://github.com/user-attachments/assets/623837d1-e1d5-46a6-9775-6374b1e810e9" />  
Histogram for the full dataset:  
<img width="1037" height="723" alt="image" src="https://github.com/user-attachments/assets/de051e0c-8145-4a85-aeaa-3357d6d59519" />  

##### Violin plots
Violin plots are available as a more compact alternative to the histograms.  
Violin plots for the overall dataset:  
<img width="729" height="726" alt="image" src="https://github.com/user-attachments/assets/d262b162-a531-4b68-acfa-701320df6358" />  

Violin plots within subcategories (of "natural environment" category, but such plots are made for all categories):  
<img width="686" height="684" alt="image" src="https://github.com/user-attachments/assets/1bfe8dd6-c30b-4396-8e0a-bee7e4955734" />  

##### Sequence logos
Sequence logos are exclusive to aligned FASTA inputs. If creating a sequence logo, you must provide a reference sequence ID. Loci of interest are specified relative to your original reference sequence, so you don't need to worry about manually calibrating the position of an amino acid of interest to its position in the alignment.  
Sequence logo for the overall dataset:  
<img width="2439" height="1460" alt="image" src="https://github.com/user-attachments/assets/538367ca-8b58-4805-ab6a-b8258fb0abed" />  
Sequence logos for each isolation source category:  
<img width="1421" height="715" alt="image" src="https://github.com/user-attachments/assets/b21a9b70-b8cb-4799-9eb4-b785b3f9010d" />  
Sequence logos for each subcategory in the "natural environment" category:  
<img width="1427" height="721" alt="image" src="https://github.com/user-attachments/assets/6336441e-d465-4ba6-af40-94e86557d35a" />  

#### Run statistical analyses
TODO

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
