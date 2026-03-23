/*
 * BLAST the query FASTA against the provided BLAST database
 */
process runBlast {

    input:
    path queryFasta
    path blastPath
    val blastName

    output:
    path "queryHits.blast", emit: blast
    path "has_genome_ids.tmp", emit: hasGenomeIDsFile

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
    outname="queryHits.blast"
    #genomeID="0"
    genomeID="1"
    echo "\$genomeID	runBlast	stub	output	file" > "\$outname" # to use $ to denote variables made in bash (as opposed to being passed in), you need to escape the dollar sign
    echo "\$genomeID	here's	yet	another	row" >> "\$outname" # to use $ to denote variables made in bash (as opposed to being passed in), you need to escape the dollar sign

    cat "/home/kcw2/test_scripts/nextflow_tests/data/myInput.txt" >> "\$outname"

    # assume that if the first genome ID is 0, all genome IDs are 0 and therefore you need to get genome IDs from NCBI esearch
    # and that if the first genome ID isn't 0, 
    firstGID=\$(head -n 1 queryHits.blast | cut -f 1)

    if [[ "\$firstGID" != "0" ]]; then
        echo "true" > has_genome_ids.tmp
    else
        echo "false" > has_genome_ids.tmp
    fi
    """
}

// process checkForGenomeIDs {
//     input:
//     path blastFile
    
//     output:
//     val hasGenomeIDs
    
//     script:
//     """
//     first_gid=\$(head -n 1 ${blastFile} | cut -f 1)
    
//     if [[ "\$first_gid" != "0" ]]; then
//         hasGenomeIDs="true"
//     else
//         hasGenomeIDs="false"
//     fi
    
//     echo "\$hasGenomeIDs"
//     """
    
//     stub:
//     """
//     first_gid=\$(head -n 1 ${blastFile} | cut -f 1)
    
//     if [[ "\$first_gid" != "0" ]]; then
//         hasGenomeIDs="true"
//     else
//         hasGenomeIDs="false"
//     fi
    
//     echo "\$hasGenomeIDs"
//     """
// }