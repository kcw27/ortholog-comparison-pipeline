/*
 * BLAST the query FASTA against the provided BLAST database
 */
process runBlast {

    input:
    path queryFasta
    path blastPath
    val blastName

    output:
    path "*_hits.blast", emit: blast
    path "benchmarking.txt", emit: benchmark
    // path "has_genome_ids.tmp", emit: hasGenomeIDsFile

    script:
    """
    
    """

    stub:
    """
    # to find the script, need to search relative to the project directory
    projDir="${workflow.projectDir}"
    head "\$projDir/../scripts/blast_processing/run_local_blast.sh"

    echo "${queryFasta}"
    echo "${blastPath}"
    echo "${blastName}"

    inputf="${queryFasta}"
    outname="\${inputf%.*}_hits.blast"

    # write a genome ID for testing how the filter-to-top-per-genome process responds
    #genomeID="0"
    genomeID="1"

    echo "\$genomeID	runBlast	stub	output	file" > "\$outname" # to use $ to denote variables made in bash (as opposed to being passed in), you need to escape the dollar sign
    echo "\$genomeID	here's	yet	another	row" >> "\$outname" # to use $ to denote variables made in bash (as opposed to being passed in), you need to escape the dollar sign

    cat "/home/kcw2/test_scripts/nextflow_tests/data/myInput.txt" >> "\$outname"

    lineCount=\$(wc -l "\$outname" | cut -f 1 -d " ")
    echo "Number of hits from original BLAST: \$lineCount" > "benchmarking.txt"

    ### decided to move the GID check to the filter-to-top-per-genome process unless I can make it work here
    # assume that if the first genome ID is 0, all genome IDs are 0 and therefore you need to get genome IDs from NCBI esearch
    # and that if the first genome ID isn't 0, then you do have genome IDs
    #firstGID=\$(head -n 1 queryHits.blast | cut -f 1)

    #if [[ "\$firstGID" != "0" ]]; then
    #    echo "true" > has_genome_ids.tmp
    #else
    #    echo "false" > has_genome_ids.tmp
    #fi
    """
}