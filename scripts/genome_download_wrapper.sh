#!/bin/bash

# Wrapper script to take a BLAST hits file with -outfmt "6 sgi sseqid sseq evalue stitle"
# and download relevant genome databases (.gbff files) that will be used for synteny search.

### CLIs:
# $1: name of BLAST hits file (in long format)
# #2: directory to which the genome databases should be downloaded

### Outputs:
# Gets the list of accessions corresponding to organisms from the BLAST hits
# Sorts the accessions into GenBank and RefSeq categories
# Calls ncbi-genome-download to download the databases corresponding to these accessions
# (download_databases.sh assumes that the genomes are in the bacterial section; modify download_databases.sh if that's not the case)

# Example run:
# bash genome_download_wrapper.sh "${HOME}/data/Fha1_orthologs.txt" "${HOME}/data/fha1_dbs" &

start_time="$(date -u +%s)"

blast=$1
dbdir=$2

wrapperdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # get location of current script
scriptsdir="${wrapperdir}/blast_processing" # this is where I put the scripts
log="${wrapperdir}/genome_download_wrapper_log.txt"
echo "CLIs for genome_download_wrapper.sh: ${blast}, ${dbdir}" > $log

name_dir=$(dirname $blast)
name_root=$(basename $blast | cut -d "." -f 1)

# run extract_organisms.sh
orgs="${name_dir}/${name_root}_organisms.txt"
bash "${scriptsdir}/extract_organisms.sh" "$blast" "$orgs" &
pid0=$! # capture process ID of the job
echo "Run extract_organisms.sh with CLIs ${blast} ${orgs}; Job ID ${pid0}, time $(date -u +%s)" >> $log
wait $pid0  # Stall execution until command finishes

# run convert_organisms_to_accessions.sh
accessions="${name_dir}/${name_root}_accessions.txt"

bash "${scriptsdir}/convert_organisms_to_accessions.sh" "$orgs" "$accessions" &
pid1=$! # capture process ID of the job
echo "Run convert_organisms_to_accessions.sh with CLIs ${orgs} ${accessions}; Job ID ${pid1}, time $(date -u +%s)" >> $log
wait $pid1  # Stall execution until command finishes

# run sort_accessions.sh
bash "${scriptsdir}/sort_accessions.sh" "$accessions" "$name_dir" &
pid2=$! # capture process ID of the job
echo "Run sort_accessions.sh with CLIs ${accessions} ${name_dir}; Job ID ${pid2}, time $(date -u +%s)" >> $log
wait $pid2  # Stall execution until command finishes

gb="${name_dir}/${name_root}_accessions_gb.txt"
ref="${name_dir}/${name_root}_accessions_ref.txt"
#gb="${name_dir}/${name_root}_accessions_genbank.txt"
#ref="${name_dir}/${name_root}_accessions_refseq.txt"

# Commented out for now, but I should set up a flag so that the user can require confirmation before downloading genomes.
#read -p "Please review ${gb} and ${ref}, then decide whether to proceed with downloading genomes (y/n). " choice
#case "$choice" in 
#  y|Y ) echo "Proceeding...";;
#  n|N ) echo "Terminating."; exit 1;;
#  * ) echo "Invalid choice."; exit 1;;
#esac


# run download_databases.sh
bash "${wrapperdir}/setup/download_databases.sh" "$gb" "$ref" "$dbdir" &
pid3=$! # capture process ID of the job
echo "Run download_databases.sh with CLIs ${gb} ${ref} ${dbdir}; Job ID ${pid3}, time $(date -u +%s)" >> $log
wait $pid3  # Stall execution until command finishes

end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"
echo "Finished! Elapsed: ${elapsed} seconds." >> $log
echo "Start time: ${start_time}" >> $log
echo "End time: ${end_time}" >> $log