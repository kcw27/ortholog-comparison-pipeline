process foo {
    output:
    val "foo"

    script:
    """
    echo "foo"
    """
}

process bar {
    conda "${workflow.projectDir}/envs/metadata-magnet-env.yaml"

    output:
    val "bar"

    script:
    """
    #!/usr/bin/env Rscript

    installed.packages()
    """
}

process filterToEvalue {
    input:
    path inputFile
    val evalue
    path benchmark_file

    output:
    path "*_evalueThreshold_${evalue}.blast", emit: blast
    path "benchmarking_evalue.txt", emit: benchmark
    
    script:
    """
    # use evalue for -t 
    """

    stub:
    """
    projDir="${workflow.projectDir}"
    head "\$projDir/../scripts/blast_processing/filter_blast_by_evalue.sh"

    # sleep 3 # to confirm that later processes WILL wait for this to finish

    # assign the variables to string variables in bash to be safe
    inputf="${inputFile}"
    thr="${evalue}"
    output="\${inputf%.*}_evalueThreshold_\${thr}.blast"
    cat "\${inputf}" > "\$output"
    echo "1: file saved to \$output" >> "\$output"

    lineCount=\$(wc -l "\$output" | cut -f 1 -d " ")
    newBenchmarkFile="benchmarking_evalue.txt"
    cat ${benchmark_file} > \$newBenchmarkFile
    echo "Number of hits remaining after filtering to evalue <= ${evalue}: \$lineCount" >> \$newBenchmarkFile

    # sleep 3
    """
}

process filterToOrganism {
    input:
    path inputFile
    path benchmark_file

    output:
    path "*_topPerOrganism.blast", emit: blast
    path "benchmarking_organism.txt", emit: benchmark
    
    script:
    """
    """

    stub:
    """
    projDir="${workflow.projectDir}"
    head "\$projDir/../scripts/blast_processing/get_blast_top_hits_by_organism.sh"

    # sleep 1 # to confirm that later processes WILL wait for this to finish

    # assign the variables to string variables in bash to be safe
    inputf="${inputFile}"
    output="\${inputf%.*}_topPerOrganism.blast"
    cat "\${inputf}" > "\$output"
    echo "2: file saved to \$output" >> "\$output"

    lineCount=\$(wc -l "\$output" | cut -f 1 -d " ")
    newBenchmarkFile="benchmarking_organism.txt"
    cat ${benchmark_file} > \$newBenchmarkFile
    echo "Number of hits remaining after filtering to top hit per organism: \$lineCount" >> \$newBenchmarkFile
    """
}

process filterToGenome {
    input:
    val signal
    path inputFile
    path benchmark_file
    val hasHeader

    output:
    path "*_topPerGenome.blast", emit: blast 
    path "benchmarking_genome.txt", emit: benchmark

    
    script:
    """
    # genome ID col will be col 1
    # determine whether to use -r based on hasHeader
    # -f is inputFile
    # output is auto-nsmed.
    """

    stub:
    """
    projDir="${workflow.projectDir}"
    head "\$projDir/../scripts/blast_processing/get_blast_top_hits_by_genomeID.sh"

    # sleep 3 # to confirm that later processes WILL wait for this to finish

    # assign the variables to string variables in bash to be safe
    inputf="${inputFile}"
    output="\${inputf%.*}_topPerGenome.blast"
    cat "\${inputf}" > "\$output"
    echo "3: file saved to \$output" >> "\$output"

    lineCount=\$(wc -l "\$output" | cut -f 1 -d " ")
    newBenchmarkFile="benchmarking_genome.txt"
    cat ${benchmark_file} > \$newBenchmarkFile
    echo "Number of hits remaining after filtering to top hit per genome: \$lineCount" >> \$newBenchmarkFile

    hasHeaderStr="${hasHeader}"
    echo "hasHeaderStr: \$hasHeaderStr" >> "\$output"
    if [[ \${hasHeaderStr} == "true" ]]; then
        echo "Input file has a header (which should be the case if metadata is retrieved, otherwise no header)" >> "\$output"
    else
        echo "Input file does NOT have a header (i.e. metadata has not yet been retrieved)" >> "\$output"
    fi
    """
}

process filterSynteny {
    conda "${workflow.projectDir}/envs/pynteny-env.yaml"

    input:
    // for the synteny search
    tuple path(genomeDBsynteny), path(syntenyInput), path(hmmsList), path(hmmsDir), path(hmmsMetadata)

    // for the intersection
    tuple val(keepSyntenyFasta), val(intersectPident), val(intersectQcovs), val(blastFiles) 

    output:
    // path "*/synteny_summary.tsv", emit: metadataFiles // synteny output dirs will be saved as subdirs of this
    path "*/*_synteny_summary.tsv", emit: syntenySummaries // synteny output dirs will be saved as subdirs of this
    path "*.fasta", emit: fastaFiles
    path "*/*_benchmarking.txt", emit: benchmarkFiles

    script:
    """
    # run synteny wrapper

    # then run the intersection script
    tempDir=$(pwd) # to keep temp files in the work directory
    # pass this as -t to the intersection script

    # only after the intersection script is done do you rename the synteny summary files, since the intersection script expects them to have unchanged names
    """

    stub:
    """
    #pynteny --help # check .command.out
    #pynteny build -i "/home/kcw2/ortholog-comparison-pipeline/nextflow/data/GCA_000967305.2_ASM96730v2_genomic.gbff" -o "my_outfile.txt"

    projDir="${workflow.projectDir}"
    head "\$projDir/../scripts/synteny_wrapper.sh"

    # This function features both the synteny search and synteny intersection

    # I've also successfully tested this with a single synteny output
    out1="outdir1"
    mkdir -p \$out1
    echo "4: synteny file saved to \${out1}/synteny_summary.tsv" > "\${out1}/synteny_summary.tsv" 
    echo "4: fasta file saved to \${out1}.fasta" > "\${out1}.fasta"
    echo "4: benchmarking file saved to \${out1}/\${out1}_benchmarking.txt" > "\${out1}/\${out1}_benchmarking.txt"
    echo "print the number of lines minus 1 in the synteny_summary.tsv file, and the number of lines divided by 2 in the fasta" >> "\${out1}/\${out1}_benchmarking.txt"

    out2="outdir2"
    mkdir -p \$out2
    echo "4: synteny file saved to \${out2}/synteny_summary.tsv" > "\${out2}/synteny_summary.tsv"
    echo "4: fasta file saved to \${out2}.fasta" > "\${out2}.fasta"
    echo "4: benchmarking file saved to \${out2}/\${out2}_benchmarking.txt" > "\${out2}/\${out2}_benchmarking.txt"
    echo "print the number of lines minus 1 in the synteny_summary.tsv file, and the number of lines divided by 2 in the fasta" >> "\${out2}/\${out2}_benchmarking.txt"

    out3="outdir3"
    mkdir -p \$out3
    echo "4: synteny file saved to \${out3}/synteny_summary.tsv" > "\${out3}/synteny_summary.tsv" 
    echo "4: fasta file saved to \${out3}.fasta" > "\${out3}.fasta"
    echo "4: benchmarking file saved to \${out3}/\${out3}_benchmarking.txt" > "\${out3}/\${out3}_benchmarking.txt"
    echo "print the number of lines minus 1 in the synteny_summary.tsv file, and the number of lines divided by 2 in the fasta" >> "\${out3}/\${out3}_benchmarking.txt"

    for outdir in \${PWD}/*/; do # PWD to get the working directory of the current process
        dirName="\$(basename "\$outdir")" # extract the name of the deepest directory, e.g. "outdir1"
        mv "\${outdir}/synteny_summary.tsv" "\${outdir}/\${dirName}_synteny_summary.tsv" # give unique names so they don't clobber each other
    done
    """
}

process collectSyntenyBenchmarking {
    input:
    path benchmark_file
    path syntenyBenchmarking

    output:
    path "benchmarking_synteny.txt"

    script:
    """
    cat "${benchmark_file}" > "benchmarking_synteny.txt"
    echo "No longer filtering the file from the original BLAST. Instead, for each synteny search, a subset is taken from the hits that closely resemble sequences in the original BLAST." >> "benchmarking_synteny.txt"
    cat ${syntenyBenchmarking} >> "benchmarking_synteny.txt" # don't put quotes around syntenyBenchmarking; let it expand if there are multiple, and cat them all
    """
}