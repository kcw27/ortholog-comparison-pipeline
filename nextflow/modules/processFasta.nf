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
    # this process is only supposed to run on metadata files not obtained through synteny search
    # when you pass values to this, make sure it gets one file as input
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
    input:
    path input_file

    output:  
    path "${inputBasename}_aligned.fasta"

    script:
    inputBasename = input_file.simpleName
    """
    """

    stub:
    inputBasename = input_file.simpleName
    """
    cat "${input_file}" > "${inputBasename}_aligned.fasta"
    echo "fasta is now aligned" >> "${inputBasename}_aligned.fasta"
    wc -l "${inputBasename}_aligned.fasta" >> "${inputBasename}_aligned.fasta"
    """
}