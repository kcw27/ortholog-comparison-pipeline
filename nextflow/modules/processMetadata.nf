process processMetadata {
    input:
    path input_file

    output:  
    path "${inputBasename}_processedMetadata.txt"

    script:
    inputBasename = input_file.simpleName
    """
    # you can now pass an outname to metadata_processing_wrapper.sh
    # use -o "processedMetadata.txt" or whatever the output filename is
    """

    stub:
    inputBasename = input_file.simpleName
    """
    cat "${input_file}" > "${inputBasename}_processedMetadata.txt"
    echo "metadata processing process" >> "${inputBasename}_processedMetadata.txt"
    wc -l "${inputBasename}_processedMetadata.txt" >> "${inputBasename}_processedMetadata.txt" 
    """
    
    // output:
    // path "foo.txt"

    // script:
    // """
    // """

    // stub:
    // """
    // echo "foo" > "foo.txt"
    // """
}