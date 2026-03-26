#!/usr/bin/env nextflow
include { runBlast } from './modules/runBlast.nf'
// include { foo; bar} from './modules/filtering.nf'
include { filterToEvalue; filterToOrganism; filterToGenome as filterToGenome_noMetadata; filterToGenome as filterToGenome_withMetadata; filterSynteny } from './modules/filtering.nf' 
include { foo_metadata as foo_metadata_beforeFiltering; foo_metadata as foo_metadata_afterFiltering; bar_metadata } from './modules/retrieveMetadata.nf'
include { retrieveMetadata as retrieveMetadata_beforeFiltering; retrieveMetadata as retrieveMetadata_afterFiltering } from './modules/retrieveMetadata.nf'
include { processMetadata } from './modules/processMetadata.nf' 
include { produceFasta; alignFasta } from './modules/processFasta.nf' 
include { FILTER1 } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter1.nf'
include { FILTER2 } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter2.nf'
include { filter3NeedSignal as filter3_noMetadata; filter3NeedSignal as filter3_withMetadata } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter3.nf'
include { filter4WithMultipleOutputs } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter4.nf'
include { FINAL } from '/home/kcw2/test_scripts/nextflow_tests/modules/finalStep.nf'

params {
    // required inputs
    queryFasta: Path
    blastPath: Path
    blastName: String

    // temp; delete later
    gate1: Boolean
    gate2: Boolean
    gate3: Boolean
    gate4: Boolean

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

    keepSyntenyFasta: Boolean
    intersectPident: String // even though it's a number, it should be passed as a string to bash
    intersectQcovs: String // even though it's a number, it should be passed as a string to bash

    // retrieving metadata
    genomeDBmetadata: Path
    entrezEmail: String
    splitSize: String

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
    // // TODO: initialize empty channels for the outputs of the synteny filtering step

    // Get the BLAST file
    runBlast(params.queryFasta, params.blastPath, params.blastName)
    def currentBlast = runBlast.out.blast
    def benchmarking = runBlast.out.benchmark


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
                benchmarking
            )

        
        // NO branch: don't have genome IDs yet, so we have no choice but to get metadata now (even though it's generally better to retrieve metadata for fewer records)
        def metadataOut = retrieveMetadata_beforeFiltering(
                hasGIDs_result.no.map { "true" } ,
                currentBlast,
                params.splitSize,
                params.entrezEmail
            )
            
            noOutput = filterToGenome_withMetadata(
                metadataOut.map { "true" },
                metadataOut,
                benchmarking
            )
            metadata = metadata.mix(noOutput.blast).view{ content -> "Contents of metadata channel after filter3: $content" }
        
        // Combine outputs
        filteredBlastGenome_output = yesOutput.blast.mix(noOutput.blast)
        currentBlast = filteredBlastGenome_output
        benchmarking = yesOutput.benchmark.mix(noOutput.benchmark)
    }
    
    // Filter 4
    // // move this over in a bit: detremining whether to do synteny search
    // def syntenyInputList = [params.genomeDBsynteny, params.syntenyInput, params.hmmsList, params.hmmsDir, params.hmmsMetadata]
    // if (syntenyInputList.any { it == null }) {
    //     println "List contains at least one null value"
    //     syntenySearchFlag = false
    // } else {
    //     println "List contains no null values"
    //     syntenySearchFlag = true
    // }
    // println("syntenySearchFlag: '${syntenySearchFlag}'")

    // if (syntenySearchFlag) {
    //     foo()
    // } else {
    //     bar()
    // }
    if (params.gate4) {
        println("Filter 4 applied")
        current = filter4WithMultipleOutputs(currentBlast)

        currentBlast = current.syntenySummaries
        current_individual_values = current.syntenySummaries.flatMap { it }.view()
        metadata = metadata.mix(current_individual_values).view{ content -> "Contents of metadata channel after filter4: $content" }

        fasta = fasta.mix(current.fastas.flatMap { it }).view{ content -> "Contents of fasta channel after filter4: $content" }
    }
    
    // after filtering, retrieve metadata (if it hasn't already been obtained in steps 3 and/or 4 of filtering)
    metadata.view{v -> "metadata channel contents: $v"}

    def metadataCheck = metadata
    .ifEmpty { "true" }
    metadataCheck.view{v-> "Value of metadataCheck: $v"}

    def metadataIsEmpty = metadataCheck.filter { it == "true" }
    def metadataIsNotEmpty = metadataCheck.filter { it == "false" } // I don't think we need this?
    metadataIsEmpty.view{v-> "Value of metadataIsEmpty: $v"}
    metadataIsNotEmpty.view{v-> "Value of metadataIsNotEmpty: $v"}

    foo_metadata_afterFiltering(metadataIsEmpty, currentBlast)


    // Process the final metadata
    def finalMetadata = metadata.mix(foo_metadata_afterFiltering.out).view{ content -> "Contents of finalMetadata channel: $content" } 
    processMetadata(finalMetadata) // will publish this

    // produce a FASTA from the current blast only if synteny search hasn't been performed (i.e. only if the fasta channel is still empty)
    // the synteny search may produce multiple FASTAs, but if you don't do synteny search, currentBlast is guaranteed to have 1 thing in it, so you just make a FASTA based on that
    def fastaCheck = fasta
    .ifEmpty { "true" }
    fastaCheck.view{v-> "Value of fastaCheck: $v"}

    def fastaIsEmpty = fastaCheck.filter { it == "true" }
    fastaIsEmpty.view{v-> "Value of fastaIsEmpty: $v"}
    produceFasta(fastaIsEmpty, currentBlast)
    fasta = fasta.mix(produceFasta.out).view{ content -> "Contents of fasta channel: $content" } 

    // and then, if params.align is true, we align the FASTA(s)- anything in the fasta channel. 
    if (params.align) {
        alignFasta(fasta)
        alignedFasta = alignedFasta.mix(alignFasta.out).view{ content -> "Contents of alignedFasta channel: $content" } 
    }

    publish:
    currentBlastFile = currentBlast
    blast_output = runBlast.out.blast
    filteredBlastEvalue_outputFile = filteredBlastEvalue_output
    filteredBlastOrganism_outputFile = filteredBlastOrganism_output
    filteredBlastGenome_outputFile = filteredBlastGenome_output
    // TODO: publish the channels for outputs of synteny filtering step

    fastaFile = fasta
    alignedFastaFile = alignedFasta
    metadataFile = processMetadata.out

    benchmarkFile = benchmarking
}

//  // test run workflow
// workflow {

//     main:
//     def metadata = Channel.empty() // May contain multiple values if performing synteny search
//     def fasta = Channel.empty() // AMay contain multiple values if performing synteny search
//     def alignedFasta = Channel.empty() // may or may not align FASTA, but either way, initialize it so there are no errors with publishing

//     // // move this over in a bit: detremining whether to do synteny search
//     // def syntenyInputList = [params.genomeDBsynteny, params.syntenyInput, params.hmmsList, params.hmmsDir, params.hmmsMetadata]
//     // if (syntenyInputList.any { it == null }) {
//     //     println "List contains at least one null value"
//     //     syntenySearchFlag = false
//     // } else {
//     //     println "List contains no null values"
//     //     syntenySearchFlag = true
//     // }
//     // println("syntenySearchFlag: '${syntenySearchFlag}'")

//     // if (syntenySearchFlag) {
//     //     foo()
//     // } else {
//     //     bar()
//     // }

//     // Get the BLAST file
//     runBlast(params.queryFasta, params.blastPath, params.blastName)
//     currentBlast = runBlast.out.blast


//     // test filtering gauntlet
//     // Filter 1
//     if (params.gate1) {
//         println("Filter 1 applied")
//         currentBlast = FILTER1(currentBlast)
//     }
    
//     // Filter 2
//     if (params.gate2) {
//         println("Filter 2 applied")
//         currentBlast = FILTER2(currentBlast)
//     }
    
//     if (params.gate3) {
//         println("Filter 3 applied")
        
//         // Read the first value from each BLAST file to determine the branch
//         def blastWithCheck = currentBlast.map { blastFile ->
//             def firstLine = blastFile.text.readLines().first()
//             def firstGID = firstLine.split()[0]  // Get first column
//             def hasGenomeIDs = (firstGID != "0")  // if the first genome ID isn't 0, assume it has genome IDs
//             return tuple(blastFile, hasGenomeIDs)
//         }.view { blastFile, hasGenomeIDs -> "Blast: $blastFile, Has genome IDs: $hasGenomeIDs" }
        
//         // Branch based on the check
//         def hasGIDs_result = blastWithCheck.branch { blastFile, hasGenomeIDs ->
//             yes: hasGenomeIDs == true
//             no: hasGenomeIDs == false
//         }
        
//         // Process YES branch
//         def yesOutput = filter3_noMetadata(
//                 hasGIDs_result.yes,
//                 currentBlast
//             )
        
//         // Process NO branch  
//         def metadataOut = foo_metadata_beforeFiltering(
//                 hasGIDs_result.no.map { "true" } ,
//                 currentBlast
//             )
            
//             noOutput = filter3_withMetadata(
//                 metadataOut.map { "true" },
//                 metadataOut
//             )
//             metadata = metadata.mix(noOutput).view{ content -> "Contents of metadata channel after filter3: $content" }
        
//         // Combine outputs
//         currentBlast = yesOutput.mix(noOutput)
//     }
    
//     // Filter 4
//     if (params.gate4) {
//         println("Filter 4 applied")
//         current = filter4WithMultipleOutputs(currentBlast)

//         // // for use when all the outputs are in the same channel and you have to split them
//         // currentBlast = current // for the synteny search, you don't update the current blast, you just update the metadata channel. But for testing I'll use this.
//         // current_individual_values = current.flatMap { it }.view()
//         // metadata = metadata.mix(current_individual_values).view{ content -> "Contents of metadata channel after filter4: $content" }

//         currentBlast = current.syntenySummaries
//         current_individual_values = current.syntenySummaries.flatMap { it }.view()
//         metadata = metadata.mix(current_individual_values).view{ content -> "Contents of metadata channel after filter4: $content" }

//         fasta = fasta.mix(current.fastas.flatMap { it }).view{ content -> "Contents of fasta channel after filter4: $content" }
//     }
    
//     // Final step
//     currentBlast = FINAL(currentBlast)

//     metadata.view{v -> "metadata channel contents: $v"}

//     def metadataCheck = metadata
//     .ifEmpty { "true" }
//     metadataCheck.view{v-> "Value of metadataCheck: $v"}

//     def metadataIsEmpty = metadataCheck.filter { it == "true" }
//     def metadataIsNotEmpty = metadataCheck.filter { it == "false" } // I don't think we need this?
//     metadataIsEmpty.view{v-> "Value of metadataIsEmpty: $v"}
//     metadataIsNotEmpty.view{v-> "Value of metadataIsNotEmpty: $v"}

//     foo_metadata_afterFiltering(metadataIsEmpty, currentBlast)


//     // Process the final metadata
//     def finalMetadata = metadata.mix(foo_metadata_afterFiltering.out).view{ content -> "Contents of finalMetadata channel: $content" } 
//     processMetadata(finalMetadata)

//     // produce a FASTA from the current blast only if synteny search hasn't been performed (i.e. only if the fasta channel is still empty)
//     // the synteny search may produce multiple FASTAs, but if you don't do synteny search, currentBlast is guaranteed to have 1 thing in it, so you just make a FASTA based on that
//     def fastaCheck = fasta
//     .ifEmpty { "true" }
//     fastaCheck.view{v-> "Value of fastaCheck: $v"}

//     def fastaIsEmpty = fastaCheck.filter { it == "true" }
//     fastaIsEmpty.view{v-> "Value of fastaIsEmpty: $v"}
//     produceFasta(fastaIsEmpty, currentBlast)
//     fasta = fasta.mix(produceFasta.out).view{ content -> "Contents of fasta channel: $content" } 

//     // and then, if params.align is true, we align the FASTA(s)- anything in the fasta channel. 
//     if (params.align) {
//         alignFasta(fasta)
//         alignedFasta = alignedFasta.mix(alignFasta.out).view{ content -> "Contents of alignedFasta channel: $content" } 
//     }

//     publish:
//     blast_output = runBlast.out.blast
//     currentBlastFile = currentBlast
//     fastaFile = fasta
//     alignedFastaFile = alignedFasta
//     metadataFile = processMetadata.out
// }

output {
    // only the final outputs will be copied; the rest are published as soft links.
    currentBlastFile {
        path { "blastFiles" }
        // mode 'copy' // for debug
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
    //  TODO: put other filteirng outputs here

}