#!/usr/bin/env nextflow
include { runBlast } from './modules/runBlast.nf'
// include { foo; bar} from './modules/filtering.nf'
include { filterToEvalue; filterToOrganism; filterToGenome as filterToGenome_noMetadata; filterToGenome as filterToGenome_withMetadata; filterSynteny; collectSyntenyBenchmarking } from './modules/filtering.nf' 
// include { foo_metadata as foo_metadata_beforeFiltering; foo_metadata as foo_metadata_afterFiltering; bar_metadata } from './modules/retrieveMetadata.nf'
include { retrieveMetadata as retrieveMetadata_beforeFiltering; retrieveMetadata as retrieveMetadata_afterFiltering } from './modules/retrieveMetadata.nf'
include { processMetadata } from './modules/processMetadata.nf' 
include { produceFasta; alignFasta } from './modules/processFasta.nf' 
// include { FILTER1 } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter1.nf'
// include { FILTER2 } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter2.nf'
// include { filter3NeedSignal as filter3_noMetadata; filter3NeedSignal as filter3_withMetadata } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter3.nf'
// include { filter4WithMultipleOutputs } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter4.nf'
// include { FINAL } from '/home/kcw2/test_scripts/nextflow_tests/modules/finalStep.nf'

params {
    // required inputs for initial blast
    queryFasta: Path
    blastPath: Path
    blastName: String

    // optional input for initial blast
    separateIDs: String // technically Boolean, but use string so that 

    // filtering parameters
    stopBeforeFiltering: Boolean

    evalueThreshold: String // even though it's a number, it should be passed as a string to bash
    filterByOrganism: Boolean
    filterByGenome: Boolean

    // // params related to synteny search
    genomeDBsynteny: Path
    syntenyInput: Path
    hmmsList: Path
    hmmsDir: Path
    hmmsMetadata: Path

    keepSyntenyFasta: String // even though it's a Boolean, it should be passed as a string to bash
    intersectPident: String // even though it's a number, it should be passed as a string to bash
    intersectQcovs: String // even though it's a number, it should be passed as a string to bash

    // retrieving metadata
    genomeDBmetadata: Path
    entrezEmail: String
    splitSize: String

    // alternatively, provide a metadata file from a previous run and skip the retrieval
    existingMetadata: Path

    // metadata categorization
    category: Path
    subcategory: Path

    // align using this pipeline?
    align: Boolean
}

workflow {
    main:    
    def metadata = Channel.empty() // May contain multiple values if performing synteny search
    def fasta = Channel.empty() // AMay contain multiple values if performing synteny search
    def alignedFasta = Channel.empty() // may or may not align FASTA, but either way, initialize it so there are no errors with publishing

    // for publication, initialize empty channels specifically for the blast outputs
    def filteredBlastEvalue_output = Channel.empty()
    def filteredBlastOrganism_output = Channel.empty()
    def filteredBlastGenome_output = Channel.empty()

    // Get the BLAST file
    runBlast(params.queryFasta, params.blastPath, params.blastName, params.separateIDs)
    def currentBlast = runBlast.out.blast
    def benchmarking = runBlast.out.benchmark

    // only proceed if not instructed to stop early
    if (!params.stopBeforeFiltering) {

        // filtering gauntlet
        // Filter by evalue
        if (params.evalueThreshold != null) {
            println("Filtering BLAST output to lines with evalue <= ${params.evalueThreshold}")
            filteredBlastEvalue = filterToEvalue(currentBlast, params.evalueThreshold, benchmarking)
            filteredBlastEvalue_output = filteredBlastEvalue.blast
            currentBlast = filteredBlastEvalue.blast
            benchmarking = filteredBlastEvalue.benchmark
        } 
        
        // Filter by organism
        if (params.filterByOrganism) {
            println("Filtering BLAST output to top hit (lowest evalue) per organism")
            filteredBlastOrganism = filterToOrganism(currentBlast, benchmarking)
            filteredBlastOrganism_output = filteredBlastOrganism.blast
            currentBlast = filteredBlastOrganism.blast
            benchmarking = filteredBlastOrganism.benchmark
        } 
        
        if (params.filterByGenome) {
            println("Filtering BLAST output to top hit (lowest evalue) per genome")
            
            // To determine whether we need to get metadata now or leave it until after filtering is complete, read the first genome ID from the BLAST file
            // if the first genome ID isn't 0, assume it has genome IDs
            def blastWithCheck = currentBlast.map { blastFile ->
                def firstLine = blastFile.text.readLines().first() // Get the first BLAST hit record
                def firstGID = firstLine.split()[0]  // The first value in the row represents the genome ID, which is either 0 (representing no data) or a legitimate genome ID
                def hasGenomeIDs = (firstGID != "0")  
                return tuple(blastFile, hasGenomeIDs)
            }.view { blastFile, hasGenomeIDs -> "Blast: $blastFile, Has genome IDs: $hasGenomeIDs" }
            
            // Branch based on the whether we have enough info to filter to top per genome yet
            def hasGIDs_result = blastWithCheck.branch { blastFile, hasGenomeIDs ->
                yes: hasGenomeIDs == true
                no: hasGenomeIDs == false
            }
            
            // YES branch: we already have genome IDs, so we don't need to get metadata yet; we can get it after we're done filtering
            def yesOutput = filterToGenome_noMetadata(
                    hasGIDs_result.yes,
                    currentBlast,
                    benchmarking,
                    hasGIDs_result.yes.map { "false" }
                )

            
            // NO branch: don't have genome IDs yet, so we have no choice but to get metadata now (even though it's generally better to retrieve metadata for fewer records)
            def metadataOut = retrieveMetadata_beforeFiltering(
                    hasGIDs_result.no.map { "true" } ,
                    currentBlast,
                    params.splitSize,
                    params.entrezEmail,
                    params.genomeDBmetadata
                )
                
                noOutput = filterToGenome_withMetadata(
                    metadataOut.map { "true" },
                    metadataOut,
                    benchmarking,
                    metadataOut.map { "true" }
                )
                metadata = metadata.mix(noOutput.blast)//.view{ content -> "Contents of metadata channel after filter3: $content" }
            
            // Combine outputs
            filteredBlastGenome_output = yesOutput.blast.mix(noOutput.blast)
            currentBlast = filteredBlastGenome_output
            benchmarking = yesOutput.benchmark.mix(noOutput.benchmark)
        }
        
        // determining whether to do synteny search
        def syntenyInputList = [params.genomeDBsynteny, params.syntenyInput, params.hmmsList, params.hmmsDir, params.hmmsMetadata] // mandatory synteny sesarch inputs
        if (syntenyInputList.any { it == null }) {
            syntenySearchFlag = false
        } else {
            syntenySearchFlag = true
        }

        if (syntenySearchFlag) {
            println("Performing synteny search")

            // Though you could technically pass the arguments as lists, this causes issues with abiility to cache because of something to do with serialization
            // It's better to create channel tuples 
            Channel
                .of(params.genomeDBsynteny)
                .combine(Channel.of(params.syntenyInput))
                .combine(Channel.of(params.hmmsList))
                .combine(Channel.of(params.hmmsDir))
                .combine(Channel.of(params.hmmsMetadata))
                .set { syntenyInputChannel }
            
            Channel
                .of(params.keepSyntenyFasta)
                .combine(Channel.of(params.intersectPident))
                .combine(Channel.of(params.intersectQcovs))
                .combine(currentBlast)
                .set { intersectionInputChannel }
            
            current = filterSynteny(syntenyInputChannel, intersectionInputChannel)

            // no longer updating currentBlast
            // synteny search outputs are saved to the metadata and fasta channels
            metadata_individual_values = current.syntenySummaries.flatMap { it }//.view()
            metadata = metadata.mix(metadata_individual_values)//.view{ content -> "Contents of metadata channel after synteny search: $content" }

            fasta = fasta.mix(current.fastaFiles.flatMap { it })//.view{ content -> "Contents of fasta channel after synteny search: $content" }

            benchmarking = collectSyntenyBenchmarking(benchmarking, current.benchmarkFiles.collect())
        }
        
        // after filtering, retrieve metadata (if it hasn't already been obtained in steps 3 and/or 4 of filtering)
        // metadata.view{v -> "metadata channel contents: $v"}

        if (params.existingMetadata != null) {
            println("Loading the provided metadata.")
            def addThisMetadata = channel.fromPath(params.existingMetadata)
            metadata = metadata.mix(addThisMetadata)
        }

        def metadataCheck = metadata
        .ifEmpty { "true" }
        // metadataCheck.view{v-> "Value of metadataCheck: $v"}

        def metadataIsEmpty = metadataCheck.filter { it == "true" }
        // metadataIsEmpty.view{v-> "Value of metadataIsEmpty: $v"}

        retrieveMetadata_afterFiltering(metadataIsEmpty.map { "true" },
            currentBlast, 
            params.splitSize, 
            params.entrezEmail,
            params.genomeDBmetadata
        )


        // Process the final metadata
        def finalMetadata = metadata.mix(retrieveMetadata_afterFiltering.out)//.view{ content -> "Contents of finalMetadata channel: $content" } 
        metadata = processMetadata(finalMetadata, params.category, params.subcategory) // will publish this

        // produce a FASTA from the current blast only if synteny search hasn't been performed (i.e. only if the fasta channel is still empty)
        // the synteny search may produce multiple FASTAs, but if you don't do synteny search, currentBlast is guaranteed to have 1 thing in it, so you just make a FASTA based on that
        def fastaCheck = fasta
        .ifEmpty { "true" }
        // fastaCheck.view{v-> "Value of fastaCheck: $v"}

        def fastaIsEmpty = fastaCheck.filter { it == "true" }
        // fastaIsEmpty.view{v-> "Value of fastaIsEmpty: $v"}
        produceFasta(fastaIsEmpty, metadata) // aligning metadata channel instead of currentBlast channel because seqs from metadata are more complete than seqs from the original BLAST alignment
        fasta = fasta.mix(produceFasta.out)//.view{ content -> "Contents of fasta channel: $content" } 

        // and then, if params.align is true, we align the FASTA(s)- anything in the fasta channel. 
        if (params.align) {
            alignFasta(fasta)
            alignedFasta = alignedFasta.mix(alignFasta.out)//.view{ content -> "Contents of alignedFasta channel: $content" } 
        }
    }

    publish:
    currentBlastFile = currentBlast
    blast_output = runBlast.out.blast
    blastSummaryFigures = runBlast.out.figures
    filteredBlastEvalue_outputFile = filteredBlastEvalue_output
    filteredBlastOrganism_outputFile = filteredBlastOrganism_output
    filteredBlastGenome_outputFile = filteredBlastGenome_output

    fastaFile = fasta
    alignedFastaFile = alignedFasta

    metadataFile = metadata //processMetadata.out

    benchmarkFile = benchmarking
}


output {
    // only the final outputs will be copied; the rest are published as soft links.
    currentBlastFile {
        path { "blastFiles" }
    }

    blastSummaryFigures {
        path { "blastFiles" }
        mode 'copy'
    }

    fastaFile {
        path { "fasta/not_aligned" }
        mode 'copy'
    }

    alignedFastaFile {
        path { "fasta/aligned" }
        mode 'copy'
    }

    metadataFile {
        path { "metadata" }
        mode 'copy'
    }

    benchmarkFile {
        path { "benchmarkingFinal" }
        mode 'copy'
    }

    blast_output {
        path { "blastFiles/intermediates" }
    }

    filteredBlastEvalue_outputFile {
        path { "blastFiles/intermediates" }
    }

    filteredBlastOrganism_outputFile {
        path { "blastFiles/intermediates" }
    }

    filteredBlastGenome_outputFile {
        path { "blastFiles/intermediates" }
    }
}