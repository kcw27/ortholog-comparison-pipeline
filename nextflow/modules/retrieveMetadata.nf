process retrieveMetadata {
    conda "${workflow.projectDir}/envs/metadata-magnet-env.yaml"

    input:
    val signal
    path inputFile // a blast file
    val splitSize // segmenting for metadata retrieval
    val email // for NCBI esearch
    path genomeDBmetadata // for local metadata retrieval

    output:
    path "*_metadata.blast" 

    when:
    signal == "true"

    script:
    """    
    retrieveMetadata() {
        MAX_RETRIES=3
        BACKOFF=30   # seconds

        local blastSegment="\$1"
        local metadataSegment="\$PWD/\${blastSegment%.*}_metadata.temp" # the output; need the dir (PWD) on it so that fetch_info has an outdir to write to
        echo "blastSegment \$blastSegment, metadataSegment \$metadataSegment" > \$retrievalLog

        echo "Starting \$blastSegment" >> \$retrievalLog
        echo "=== Log for \$blastSegment started at \$(date +"%Y-%m-%d %H:%M:%S") ===" >> \$retrievalLog

        local attempt=1
        while (( attempt <= MAX_RETRIES )); do
            echo "Attempt \$attempt for \$blastSegment" >> \$retrievalLog

            projDir="${workflow.projectDir}"
            python "\$projDir/../scripts/blast_processing/blast2gen.py" "\$blastSegment" "\$metadataSegment" "${email}" &

            pid=\$!
            echo "PID: \$pid" >> \$retrievalLog

            # Wait for job to finish
            if wait "\$pid"; then # indicates that the job has exited with an error code of 0, i.e. successfully
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

    inputf="${inputFile}"
    output="\${inputf%.*}_metadata.blast"
    export retrievalLog='retrieveMetadata.log'
    
    firstGID=\$(head -n 1 ${inputFile} | cut -f 1) 
    if [[ "\$firstGID" != "0" && -d ${genomeDBmetadata} ]]; then # assume we have genome IDs, so we can perform local database search
        echo "Performing local database search..."
        projDir="${workflow.projectDir}"
        bash "\$projDir/../scripts/blast_processing/local_metadata_retrieval.sh" ${genomeDBmetadata} ${inputFile}
    else # we do need to get metadata from blast2gen.py
        echo "Performing NCBI esearch..."
        # file splitting
        split -d -a 2 -l ${splitSize} ${inputFile} "splitFile_part" --additional-suffix=.blast

        # sequential retrieval
        # it's guaranteed to run on only one file at a time because the wait pid in retrieveMetadata blocks, and it won't return until completed
        for file in splitFile_part*; do
            echo "Processing: \$file"
            retrieveMetadata "\$file"
        done

        # re-joining all the metadata parts
        # first write the header of the first metadata file. There's guaranteed to be at least one file.
        head -n 1 "splitFile_part00_metadata.temp" > "temp.tmp"

        # then write everything after the header in all the files matching the pattern
        tail -n +2 -q splitFile_part*_metadata.temp >> "temp.tmp" # q flag so it doesn't write filenames to the file, and no quotes around input filenames so the wildcard works

        mv "temp.tmp" "\$output" # rename to the expected output name
    fi
    """

    stub:
    """
    projDir="${workflow.projectDir}"
    head "\$projDir/../scripts/blast_processing/blast2gen.py"

    inputf="${inputFile}"
    output="\${inputf%.*}_metadata.blast"
    cat "${inputFile}" > "\$output"
    echo "Called blast2gen.py with arguments: inputFile ${inputFile}, output \$output, email ${email}" >> "\$output"
    """
}

process foo_metadata {
    input:
    val signal
    path input_file

    output:
    path "foo_metadata.txt"

    when:
    signal == "true"

    script:
    """
    cat "/home/kcw2/ortholog-comparison-pipeline/example_run/path1/PA3565_hits_filtered_metadata.blast" > "foo_metadata.txt"
    cat "${input_file}" >> "foo_metadata.txt"
    echo "foo metadata" >> "foo_metadata.txt"
    wc -l "foo_metadata.txt" >> "foo_metadata.txt"
    """
}

process bar_metadata {
    input:
    path blastFile
    path hasGenomeIDsFile

    output:
    path "report.txt"

    script:
    """
    cat "${blastFile}" > report.txt
    echo "" >> report.txt
    echo "***" >> report.txt
    cat "${hasGenomeIDsFile}" >> report.txt
    """
}