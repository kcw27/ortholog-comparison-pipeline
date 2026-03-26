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
    path benchmark_file

    output:
    path "*_evalueThreshold_${evalue}.blast", emit: blast
    path "benchmarking_evalue.txt", emit: benchmark
    
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

    sleep 1 # to confirm that later processes WILL wait for this to finish

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

    output:
    path "*_topPerGenome.blast", emit: blast 
    path "benchmarking_genome.txt", emit: benchmark

    
    script:
    """
    """

    stub:
    """
    projDir="${workflow.projectDir}"
    head "\$projDir/../scripts/blast_processing/get_blast_top_hits_by_genomeID.sh"

    sleep 3 # to confirm that later processes WILL wait for this to finish

    # assign the variables to string variables in bash to be safe
    inputf="${inputFile}"
    output="\${inputf%.*}_topPerGenome.blast"
    cat "\${inputf}" > "\$output"
    echo "3: file saved to \$output" >> "\$output"

    lineCount=\$(wc -l "\$output" | cut -f 1 -d " ")
    newBenchmarkFile="benchmarking_genome.txt"
    cat ${benchmark_file} > \$newBenchmarkFile
    echo "Number of hits remaining after filtering to top hit per genome: \$lineCount" >> \$newBenchmarkFile
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