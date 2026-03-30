process processMetadata {
    conda "${workflow.projectDir}/envs/metadata-magnet-env.yaml"

    input:
    path input_file
    path categoryFile
    path subcategoryFile

    output:  
    path "${inputBasename}_processedMetadata.blast"

    script:
    inputBasename = input_file.simpleName
    """
    projDir="${workflow.projectDir}"
    outname="${inputBasename}_processedMetadata.blast"

    bash "\$projDir/../scripts/metadata_processing_wrapper.sh" -f ${input_file} -o "\$outname" -c ${categoryFile} -s ${subcategoryFile}
    """

    stub:
    inputBasename = input_file.simpleName
    """
    cat "${input_file}" > "${inputBasename}_processedMetadata.blast"
    echo "metadata processing process" >> "${inputBasename}_processedMetadata.blast"
    wc -l "${inputBasename}_processedMetadata.blast" >> "${inputBasename}_processedMetadata.blast" 
    """
}