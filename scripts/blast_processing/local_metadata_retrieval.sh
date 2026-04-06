#!/bin/bash

# integrating this into Nextflow; not implementing flags unless it becomes a user-facing script
### Inputs:
# $1: genome db to search, e.g.
# my_genome_db
# ├── GCA_003691345.1
# │   ├── GCA_003691345.1_ASM369134v1_genomic.gbff.gz
# ├── GCF_040931635.1
# │   ├── GCF_040931635.1_ASM4093163v1_genomic.gbff.gz
# ├── GCF_044759015.1
# │   ├── GCF_044759015.1_6536_genomic.gbff.gz

# $2: BLAST file with legitimate genome IDs in column 1 and protein IDs in col 2
# It is assumed that there is no header row, and the 6 existing columns are sallgi sallseqid sseq evalue salltitles organism

### Output: adds metadata columns

# need to add a header at some point
# and make sure that the columns are in the right order, matching the outputs of blast2gen.py in the first 6 columns
# (can't find the other two kinds of accessions)
# To save time, sort by genome ID, and within the genome find all the proteins pertaining to that genome... or not necessarily; I think this can be parallelized like find_synteny_hits.sh
# if there are multiple matches for the protein ID (i.e. it occurs at multiple loci) then just take the first... the bulk of the metadata's the same, but the sequence might differ between loci

# Write to temp files in current working directory? Or rather, subdirs within it so that they don't clobber each other...? Wait, do we even need to have multiple temp files? I don't think so

export genome_db=${1%/}
export blast=$2

export metadata_file="local_metadata.blast" # this is for Nextflow, so just save it locally with a name matching the expected output format
echo "genome_id	protein_id	sequence	evalue	title_old	locus_tag	organism	isolation_source	titles	seq_tech	assembly_method	genome_coverage	isolation_site	sequence_old	organism_old" > $metadata_file

get_metadata() {
    # Input: a single row of the BLAST file
    # From that, gets the genome accession (col 1) and protein accession (col 2) and writes metadata for these to the output file
    input="$1"

    # existing metadata
    genome_id=$(printf '%s\n' "$input" | cut -f1)
    protein_id=$(printf '%s\n' "$input" | cut -f2)
    sequence_old=$(printf '%s\n' "$input" | cut -f3)
    evalue=$(printf '%s\n' "$input" | cut -f4)
    title_old=$(printf '%s\n' "$input" | cut -f5)
    organism_old=$(printf '%s\n' "$input" | cut -f6)

    # get new metadata
    database_genome="${genome_db}/${genome_id}"
    if [[ ! -d $database_genome ]]; then
        echo "Error: ${genome_id} not found in database"
        return 0
    fi

    echo "Processing ${genome_id}..."
    # metadata that applies to the entire genome
    genbank_file=$(find "$database_genome" -mindepth 1 -maxdepth 1 -name "*.gbff*" -print -quit)
    
    organism=$(zcat $genbank_file | grep "\bORGANISM\b" | head -n 1 | awk '{$1=$1;print}' | cut -d " " -f 2-)

    isolation_source=$(zcat $genbank_file | sed -n '/isolation_source="/,/"/p' | tr -d '\n' | tr -s " " \
        | sed 's|" /|"\n/|g' | grep "^/isolation_source" | sort | uniq | tr '\n' ' ')

    titles=$(zcat $genbank_file | sed -n '/\bTITLE\b/,/\bJOURNAL\b/p' | grep -v "JOURNAL" | tr -d '\n' | tr -s " " \
        | sed 's|TITLE|\nTITLE|g' | grep "^TITLE" | sort | uniq | tr '\n' ' ')

    isolation_site=$(zcat "$genbank_file" |
        awk '
            /Isolation Site/ {
            if (done) next   # ignore additional matches
            inblock = 1

            line = $0
            sub(/^[[:space:]]+/, "", line)
            sub(/.*::[[:space:]]*/, "", line)
            print line
            next
            }

            inblock {
            if ($0 ~ /^ {12}[^[:space:]]/) {
                inblock = 0
                done = 1
                next
            }

            line = $0
            sub(/^[[:space:]]+/, "", line)
            sub(/.*::[[:space:]]*/, "", line)
            print line
            }
        ' |
        sed ':a;N;$!ba;s/\n/ /g' |
        tr '\t' ' ' |      # <-- new: convert tabs to spaces
        tr -s ' '          # squeeze spaces
        )

    seq_tech=$(zcat $genbank_file | awk '$0 ~/Sequencing Technology/ { sub(/^[^:]*[: ]+/, ""); print }' | head -n 1 )
    asm_method=$(zcat $genbank_file | awk '$0 ~/Assembly Method/ { sub(/^[^:]*[: ]+/, ""); print }' | head -n 1 )
    genome_coverage=$(zcat $genbank_file | awk '$0 ~/Genome Coverage/ { sub(/^[^:]*[: ]+/, ""); print }' | head -n 1)

    # metadata that pertains to this particular sequence
    # if there are multiple matches to $protein_id (which might be the case- the same protein may appear on multiple loci), this only returns the first.
    # while it is possible to get metadata for all loci with this protein, I am wary of overrepresenting the isolation source metadata for genomes that have multiple copies of a protein.
    sequence=$(zcat "$genbank_file" | awk -v tag="$protein_id" '$0 ~ tag, /  gene  /' | \
        awk '/translation="/ {flag=1; sub(/.*translation="/, "")} flag && /"/ {flag=0; sub(/".*/, ""); print} flag' | \
        tr -d '[:space:]') # in the .gbff, the sequence was split over multiple lines; this combines it into a single line

    locus_tag=$(zcat "$genbank_file" |
        awk -v pid="$protein_id" '
            /locus_tag="/ {
                # extract the value inside quotes
                match($0, /locus_tag="([^"]+)"/, m)
                if (m[1] != "") last = m[1]
            }
            $0 ~ pid {
                print last
                exit
            }
        '
    )

    echo "$genome_id	$protein_id	$sequence	$evalue	$title_old	$locus_tag	$organism	$isolation_source	$titles	$seq_tech	$asm_method	$genome_coverage	$isolation_site	$sequence_old	$organism_old" >> $metadata_file
}

export -f get_metadata
parallel -j $(nproc) --eta --progress get_metadata :::: "$blast"