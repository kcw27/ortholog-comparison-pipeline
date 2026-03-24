process foo {
    output:
    val "foo"

    script:
    """
    echo "foo"
    """
}

process bar {
    output:
    val "bar"

    script:
    """
    echo "bar"
    """
}

process filterToEvalue {
    input:
    path inputFile
    val evalue

    output:
    path "*_evalueThreshold_${evalue}.blast"
    
    script:
    """
    """

    stub:
    """
    projDir="${workflow.projectDir}"
    head "\$projDir/../scripts/blast_processing/filter_blast_by_evalue.sh"

    sleep 3 # to confirm that later processes WILL wait for this to finish

    # assign the variables to string variables in bash to be safe
    inputf="${inputFile}"
    thr="${evalue}"
    output="\${inputf%.*}_evalueThreshold_\${thr}.blast"
    cat "\${inputf}" > "\$output"
    echo "1: file saved to \$output" >> "\$output"

    sleep 3
    """
}

process filterToOrganism {
    input:
    path inputFile

    output:
    path "*_topPerOrganism.blast"
    
    script:
    """
    """

    stub:
    """
    projDir="${workflow.projectDir}"
    head "\$projDir/../scripts/blast_processing/get_blast_top_hits_by_organism.sh"

    sleep 3 # to confirm that later processes WILL wait for this to finish

    # assign the variables to string variables in bash to be safe
    inputf="${inputFile}"
    output="\${inputf%.*}_topPerOrganism.blast"
    cat "\${inputf}" > "\$output"
    echo "2: file saved to \$output" >> "\$output"
    """
}

process filterToGenome {
    input:
    path inputFile
    val splitSize
    val email

    output:
    path "*_topPerGenome.blast", emit: blast
    path "file_containing_metadata_path.txt", emit: fileOfMetadataPath 
    // determine the value within the script- either a real path, or still ''
    // you'll need to read the value from this file when making a channel
    
    script:
    """
    """

    stub:
    """
    set -euo pipefail
    
    retrieveMetadata() {
        MAX_RETRIES=3
        BACKOFF=30   # seconds
        retrievalLog='retrieveMetadata.log'
        > \$retrievalLog

        local blastSegment="\$1"
        local metadataSegment="\${blastSegment%.*}_metadata.blast" # the output

        echo "Starting \$blastSegment" >> \$retrievalLog
        echo "=== Log for \$blastSegment started at \$(date +"%Y-%m-%d %H:%M:%S") ===" >> \$retrievalLog

        local attempt=1
        while (( attempt <= MAX_RETRIES )); do
            echo "Attempt \$attempt for \$blastSegment" >> \$retrievalLog

            awk '{print \$0 "_METADATA-TAG AND \${metadataSegment} AND ${email}"}' \$blastSegment > \$metadataSegment & # in real script, replace this with the actual blast2gen.py call, with email param

            pid=\$!
            echo "PID: \$pid" >> \$retrievalLog

            # Wait for job to finish
            if wait "\$pid"; then # indicates that the job has exited with an error code of 0
                # once job is done, check if file exists; treat as a failure otherwise
                # theoretically, if the job exits successfully then 
                if [[ -f "\$metadataSegment" ]]; then
                    echo "Success on attempt \$attempt" >> \$retrievalLog
                    echo "Finished \$blastSegment successfully" >> \$retrievalLog
                    return
                else 
                    echo "Failure on attempt \$attempt; file was not written to \$metadataSegment" >> \$retrievalLog
                    (( attempt++ ))
                    if (( attempt <= MAX_RETRIES )); then
                        echo "Retrying after \$BACKOFF seconds..." >> \$retrievalLog
                        sleep "\$BACKOFF"
                    fi
                fi
            else # reach this block if the job has failed, e.g. from an unhandled exception in the script
                echo "Failure on attempt \$attempt; job has crashed" >> \$retrievalLog
                (( attempt++ ))
                if (( attempt <= MAX_RETRIES )); then
                    echo "Retrying after \$BACKOFF seconds..." >> \$retrievalLog
                    sleep "\$BACKOFF"
                fi
            fi
        done

        echo "FAILED after \$MAX_RETRIES attempts: \$blastSegment" >> \$retrievalLog
        echo "=== FAILED after \$MAX_RETRIES attempts ===" >> \$retrievalLog
    }

    projDir="${workflow.projectDir}"
    head "\$projDir/../scripts/blast_processing/get_blast_top_hits_by_genomeID.sh"

    sleep 3 # to confirm that later processes WILL wait for this to finish

    # assign the variables to string variables in bash to be safe
    inputf="${inputFile}"
    output="\${inputf%.*}_topPerGenome.blast"

    # assume that if the first genome ID is 0, all genome IDs are 0 and therefore you need to get genome IDs from NCBI esearch
    # and that if the first genome ID isn't 0, then you do have genome IDs
    firstGID=\$(head -n 1 queryHits.blast | cut -f 1)

    cat "\${inputf}" > "\$output"

    if [[ "\$firstGID" == "0" ]]; then # do need to get metadata
        head "\$projDir/../scripts/blast_processing/blast2gen.py" 
        
        # file splitting
        split -d -a 2 -l \$splitSize \$inputFile splitFile_part --additional-suffix=.blast

        # sequential retrieval
        # it's guaranteed to run on only one file at a time because the wait pid in retrieveMetadata blocks, and it won't return until completed
        for file in splitFile_part*; do
            echo "Processing: \$file"
            retrieveMetadata "\$file"
        done

        # re-joining all the metadata parts, saving to metadataFilePath
        # first write the header of the first file, splitFile_part00.blast
        head -n 1 "splitFile_part00.blast" > "temp.tmp"

        # then write everything after the header in all the files matching the pattern
        tail -n +2 -q splitFile_part*.blast >> "temp.tmp" # q option so it doesn't write filenames to the file, and no quotes around input filenames so the wildcard works

        echo "\$output" > "file_containing_metadata_path.txt" # is it okay to overwrite the input BLAST file with blast2gen.py? probably 
        mv "temp.tmp" \$metadataFilePath"
        echo "3.1: metadata retrieval WAS performed; metadataFilePath is \$metadataFilePath" >> "\$output"
    else # in the real script, still need to keep the else block to assign metadataFilePath
        > "file_containing_metadata_path.txt" # writes nothing to the file
        echo "3.1: metadata retrieval NOT performed; metadataFilePath is \$metadataFilePath" >> "\$output"
    fi

    echo "3.2: file saved to \$output" >> "\$output"
    """
}

process filterSynteny {
    input:
    tuple path(genomeDBsynteny), path(syntenyInput), path(hmmsList), path(hmmsDir), path(hmmsMetadata)

    // for the intersection
    val keepSyntenyFasta
    val intersectPident
    val intersectQcovs
    val initialBlast 

    output:
    path "*/synteny_input.tsv", emit: metadataFiles // synteny output dirs will be saved as subdirs of this
    path "*/*.fasta", emit: fastaFiles
    script:
    """
    """

    stub:
    """
    projDir="${workflow.projectDir}"
    head "\$projDir/../scripts/synteny_wrapper.sh"

    # man, but what do I do about this taking multiple outdirs? Then there will be multiple metadata files???
    # actually, I think it's fine to have multiple values in the metadata channel?? all we do with the metadata files afterward is publish them
    # though it is concerning that we have to check whether the metadata channel contains "" vs something else
    # and come to think of it, intersecting blast with synteny is going to produce multiple fastas
    # so do we align multiple fastas??
    # let me think about whether I want to restrict this to one synteny search or not

    # also run synteny intersection

    # for now, do a very simple stub
    # later comment out the out2 and out3 lines to test with a single synteny output
    out1="outdir1"
    echo "4: synteny file saved to \${out1}/synteny_input.tsv" >> "\${out1}/synteny_input.tsv"
    echo "4: fasta file saved to \${out1}/\${out1}.fasta" >> "\${out1}/\${out1}.fasta"

    out2="outdir2"
    echo "4: synteny file saved to \${out2}/synteny_input.tsv" >> "\${out2}/synteny_input.tsv"
    echo "4: fasta file saved to \${out2}/\${out2}.fasta" >> "\${out2}/\${out2}.fasta"

    out3="outdir3"
    echo "4: synteny file saved to \${out3}/synteny_input.tsv" >> "\${out3}/synteny_input.tsv"
    echo "4: fasta file saved to \${out3}/\${out3}.fasta" >> "\${out3}/\${out3}.fasta"
    """
}