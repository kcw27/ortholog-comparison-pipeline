process produceFasta {
    input:
    val signal // run this process only if the fasta channel is empty
    path input_file

    output:  
    path "${inputBasename}.fasta"
    when:
    signal == "true"

    script:
    inputBasename = input_file.simpleName
    """
    # this process is only supposed to run on metadata files not obtained through synteny search; signal is only true if no synteny search was done
    # when you pass values to this, make sure it gets one file as input
    projDir="${workflow.projectDir}"
    bash "\$projDir/../scripts/blast_processing/convert_blast_to_fasta.sh" -b ${input_file} -r -g -c 'locus' 
    # in this case, removing gap characters is fine because  
    """

    stub:
    inputBasename = input_file.simpleName
    """
    cat "${input_file}" > "${inputBasename}.fasta"
    echo "made into a fasta" >> "${inputBasename}.fasta"
    wc -l "${inputBasename}.fasta" >> "${inputBasename}.fasta"
    """
}

process alignFasta {
    conda "${workflow.projectDir}/envs/metadata-magnet-env.yaml"

    input:
    path input_file // fasta

    output:  
    path "${inputBasename}_aligned.fasta"

    script:
    inputBasename = input_file.simpleName
    """
    seqCount=\$(grep '^>' ${input_file} | wc -l) # check for lines starting with >
    echo "There are \$seqCount sequences in the input FASTA."

    if [ "\$seqCount" -le 1 ]; then # can't align fasta with 0 or 1 sequences
        echo "Fewer than 2 sequences; not running MUSCLE."
        cp ${input_file} "${inputBasename}_aligned.fasta"
    elif [ "\$seqCount" -le 1000 ]; then # <= 1000 seqs: use full MUSCLE algorithm, since the size is reasonably small
        muscle -align ${input_file} -output "${inputBasename}_aligned.fasta"
    else # use super5 for larger datasets
        muscle -super5 ${input_file} -output "${inputBasename}_aligned.fasta"
    fi
    """

    stub:
    inputBasename = input_file.simpleName
    """
    cat "${input_file}" > "${inputBasename}_aligned.fasta"
    echo "fasta is now aligned" >> "${inputBasename}_aligned.fasta"
    wc -l "${inputBasename}_aligned.fasta" >> "${inputBasename}_aligned.fasta"
    """
}