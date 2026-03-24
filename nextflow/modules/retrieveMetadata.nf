// eventually I'll set it up so that there are two processes: one to run blast2gen.py and one to run the 

process foo_metadata {
    input:
    val signal // doesn't do anything, just blocks until the upstream process executes
    path input_file

    output:
    path "foo_metadata.txt"

    script:
    """
    cat "${input_file}" > "foo_metadata.txt"
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