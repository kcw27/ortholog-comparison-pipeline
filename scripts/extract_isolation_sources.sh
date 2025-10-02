#!/bin/bash

### Inputs:
## $1: a genome database directory $genome_db with the following structure:
# genome_db="~/data/genome_db/" # its direct children are subdirectories named after accessions

# - $genome_db
# -- GCF_001646865.1
# --- GCF_001646865.1_rk21_genomic.gbff.gz 
# -- GCF_001646945.1
# --- GCF_001646945.1_AVMART05_1.0_genomic.gbff.gz 
# -- GCF_001647025.1
# --- GCF_001647025.1_EnAtr_2.0_genomic.gbff.gz 

## $2: name of output file (tsv) with the following columns: assembly_accession, organism, isolation_source, title
# title is included so that isolation source can be rescued using rescue_source() from metadata_processing.py

### Example run:
# bash "/home/kcw2/ortholog-comparison-pipeline/scripts/extract_isolation_sources.sh" "/home/share/pa/pa_genomes/refseq/bacteria/" "/home/kcw2/data/testing/pa_genomes_isoSources.tsv"

genome_db=$1
fname=$2
> "$fname"  # Create a blank summary file
export fname

# Function to extract metadata
write_iso_sources() {
  genome=${1%/}
  genome_id=$(basename "$genome")
  database_genome="${genome_db}/${genome_id}"
  
  echo "Processing: $genome"

  # In $database_genome, search one level down for first matching .gbff* file
  genbank_file=$(find "$genome" -mindepth 1 -maxdepth 1 -name "*.gbff*" -print -quit)

  # First, get the information that is the same across all loci in locus_tags.
  organism=$(zcat $genbank_file | grep "\bORGANISM\b" | head -n 1 | awk '{$1=$1;print}' | cut -d " " -f 2-)
  
  # further-updated code that is supposed to filter out extraneous lines and remove duplicates:
  # after finding matches and collapsing into a single line, reintroduce newlines (only between separate entries, not within entries)
  # then grep for the desired header (e.g. "^/isolation_source" to only keep lines that start with /isolation_source;
  # for some reason, sed include lines that came directly after the intended isolation source match), and use sort and uniq.
  isolation_source=$(zcat $genbank_file | sed -n '/isolation_source="/,/"/p' | tr -d '\n' | tr -s " " \
    | sed 's|" /|"\n/|g' | grep "^/isolation_source" | sort | uniq | tr '\n' ' ')
  # remove all newlines, then squeeze spaces, then reintroduce newlines only between separate entries to use sort and uniq,
  # then remove newlines again by replacing with spaces
  titles=$(zcat $genbank_file | sed -n '/\bTITLE\b/,/\bJOURNAL\b/p' | grep -v "JOURNAL" | tr -d '\n' | tr -s " " \
    | sed 's|TITLE|\nTITLE|g' | grep "^TITLE" | sort | uniq | tr '\n' ' ')
  # finds matches between lines featuring TITLE and JOURNAL
  # (with word boundaries), then uses inverse grep to remove the lines featuring JOURNAL, then removes newlines and squeezes spaces
  # then reintroduce newlines only between separate entries to use sort and uniq,
  # then remove newlines again by replacing with spaces

  echo "$genome_id	$organism	$isolation_source	$titles" >> "$fname"

}


wrapperdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # get location of current script
log="${wrapperdir}/write_iso_sources.log"
echo "CLIs: ${genome_db}, ${input_file}, ${hmms}, ${hmms_dir}, ${hmms_metadata}" > $log

start_time="$(date -u +%s)"


# Parallelize metadata extraction

export -f write_iso_sources
  
# Go through directories two levels down from the Pynteny output directory to access outputs for each genome,
# because the directory one level down is just the "genomes" directory created to organize them.
# But also pass $outdir itself as an argument so that the corresponding synteny_summary.tsv file can be found.
find "$genome_db" -mindepth 1 -maxdepth 1 -type d | parallel -j $(nproc)  --eta --progress --joblog $log 'write_iso_sources {}'

end_time="$(date -u +%s)"


elapsed="$(($end_time-$start_time))"
echo "Finished! Elapsed: ${elapsed} seconds." >> $log
echo "Start time: ${start_time}" >> $log
echo "End time: ${end_time}" >> $log