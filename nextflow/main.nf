#!/usr/bin/env nextflow
include { runBlast } from './modules/runBlast.nf'
include { foo; bar; filterToEvalue; filterToOrganism; filterToGenome; filterSynteny } from './modules/filtering.nf'
include { foo_metadata as foo_metadata_beforeFiltering; foo_metadata as foo_metadata_afterFiltering; bar_metadata } from './modules/retrieveMetadata.nf'
include { processMetadata } from './modules/processMetadata.nf' 
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

    // templ delete later
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
    // println("hi")

    main:
    // def currentBlast = Channel.empty() // not strictly necessary because it's guaranteed to populate right after, but this is for the sake of clarity
    def metadata = Channel.empty() // May contain multiple values.
    // currentFasta = Channel.empty() // Actually, I should keep this empty so that if we opt in for alignment, it's restricted until after we have a FASTA at all
    // either make the FASTA from the metadata file (maybe I shouldn't initialize that then...), or if synteny search was performed, use the FASTA from that 
    def alignedFasta = Channel.empty() // may or may not align FASTA, but either way, initialize it so there are no errors with publishing

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

    // Get the BLAST file
    runBlast(params.queryFasta, params.blastPath, params.blastName)
    currentBlast = runBlast.out.blast
    
    // runBlast.out.hasGenomeIDsFile.view { "File path: $it" }
    // runBlast.out.hasGenomeIDsFile
    //     .map { file -> file.text.trim() }
    //     .view { contents -> "Current contents: $contents" }


    // test filtering gauntlet
    // def current = Channel.empty()
    // Filter 1
    if (params.gate1) {
        println("Filter 1 applied")
        current = FILTER1(currentBlast)
        currentBlast = current
    }
    
    // Filter 2
    if (params.gate2) {
        println("Filter 2 applied")
        current = FILTER2(currentBlast)
        currentBlast = current
    }
    
    if (params.gate3) {
        println("Filter 3 applied")
        
        // Read the first value from each BLAST file to determine the branch
        def blastWithCheck = currentBlast.map { blastFile ->
            def firstLine = blastFile.text.readLines().first()
            def firstGID = firstLine.split()[0]  // Get first column
            def hasGenomeIDs = (firstGID != "0")  // if the first genome ID isn't 0, assume it has genome IDs
            return tuple(blastFile, hasGenomeIDs)
        }.view { blastFile, hasGenomeIDs -> "Blast: $blastFile, Has genome IDs: $hasGenomeIDs" }
        
        // Branch based on the check
        def hasGIDs_result = blastWithCheck.branch { blastFile, hasGenomeIDs ->
            yes: hasGenomeIDs == true
            no: hasGenomeIDs == false
        }
        
        // Process YES branch
        def yesOutput = filter3_noMetadata(
                hasGIDs_result.yes,
                currentBlast
            )
        
        // Process NO branch  
        def metadataOut = foo_metadata_beforeFiltering(
                hasGIDs_result.no.map { "true" } ,
                currentBlast
            )
            
            noOutput = filter3_withMetadata(
                metadataOut.map { "true" },
                metadataOut
            )
            metadata = metadata.mix(noOutput).view{ content -> "Contents of metadata channel after filter3: $content" }
        
        // Combine outputs
        current = yesOutput.mix(noOutput)
        currentBlast = current
    }
    
    // Filter 4
    if (params.gate4) {
        println("Filter 4 applied")
        current = filter4WithMultipleOutputs(currentBlast)

        currentBlast = current // for the synteny search, you don't update the current blast, you just update the metadata channel
        current_individual_values = current.flatMap { it }.view()
        metadata = metadata.mix(current_individual_values).view{ content -> "Contents of metadata channel after filter4: $content" } //  filter4WithMultipleOutputs.out

    }
    
    // Final step
    current = FINAL(currentBlast)
    currentBlast = current

    // test out the metadata processes
    // Create the metadata channel by mixing outputs
    // def metadata = noOutput[0].mix(filter4outputs_ch[0]).view{ content -> "Contents of metadata channel: $content" }

    // // Split into two branches based on whether metadata is empty
    // def (metadataWithData, metadataEmpty) = metadata
    //     .branch {
    //         hasData: it != null  // or some condition that indicates it has data
    //         empty: true  // fall-through for empty
    //     }

    // metadataWithData.view{ content -> "Contents of metadataWithData channel: $content" }
    // metadataEmpty.view{ content -> "Contents of metadataEmpty channel: $content" }

    // // Run foo_metadata_afterFiltering only on the empty branch
    // def fooMetadataOut = foo_metadata_afterFiltering(metadataEmpty.map {true}, metadataEmpty)

    // // Combine the non-empty metadata with the processed empty branch
    // def finalMetadata = metadataWithData.mix(fooMetadataOut).view{ content -> "Contents of finalMetadata channel: $content" }

    // Try to get metadata with a default if empty
    // def call_foo_metadata_on_this = Channel.empty()

    // def finalMetadata = metadata
    //     .ifEmpty { 
    //         // If empty, create metadata from foo_metadata_afterFiltering
    //         // foo_metadata_afterFiltering(Channel.of(true), Channel.empty())
    //         call_foo_metadata_on_this = currentBlast
    //     }
    //     .ifEmpty {
    //         // This shouldn't happen, but just in case
    //         error "Failed to generate metadata"
    //     }

    // def fooMetadataOut = foo_metadata_afterFiltering(call_foo_metadata_on_this.map {true}, call_foo_metadata_on_this)


    metadata.view{v -> "metadata channel contents: $v"}

    def metadataCheck = metadata
    .ifEmpty { "true" }
    metadataCheck.view{v-> "Value of metadataCheck: $v"}

    def metadataIsEmpty = metadataCheck.filter { it == "true" }
    def metadataIsNotEmpty = metadataCheck.filter { it == "false" } // I don't think we need this?
    metadataIsEmpty.view{v-> "Value of metadataIsEmpty: $v"}
    metadataIsNotEmpty.view{v-> "Value of metadataIsNotEmpty: $v"}

    foo_metadata_afterFiltering(metadataIsEmpty, currentBlast)

    // // the code below works, but I'm trying something else
    // metadataIsEmpty = metadata | ifEmpty { "true" } // if not empty, it'll just contain the contents of the metadata channel
    // metadataIsEmpty.view{v -> "Contents of metadataIsEmpty channel: $v"}
    // foo_metadata_afterFiltering(metadataIsEmpty, currentBlast)


    // Process the final metadata
    def finalMetadata = metadata.mix(foo_metadata_afterFiltering.out).view{ content -> "Contents of finalMetadata channel: $content" } // fix this in a bit
    processMetadata(finalMetadata)

    // metadata = noOutput.out[0]
    //     .mix(filter4WithMultipleOutputs.out[0])
    //     .view{ content -> "Contents of metadata channel: $content" }

    // metadata
    //     .ifEmpty {
    //         println("Metadata channel is empty; using the output of FILTER1 for this test")
    //         def metadata_now = FILTER1.out.view()
    //     }
    
    // foo_metadata_afterFiltering(metadata_now)
    // updated_metadata = metadata[0].mix(foo_metadata_afterFiltering.out[0]).view()

    // processMetadata(updated_metadata)


    // check whether metadata has been retrieved yet
    // metadata = // mix the channels from foo_metadata.out and filter4.out
    // and if this channel is empty, run foo_metadata_afterFiltering and then the categorization process; else, skip to the categorization process

    // Run through the filtering gauntlet
    // if (params.evalueThreshold ) {
    //     currentBlast = FILTER1(currentBlast)
    //     filter1_out = current
    // }







    
    // println "Resolved params:"
    // println "  character = ${params.character}"
    // println "  batch     = ${params.batch}"
    // println "  input_tsv = ${params.input_tsv}"

    // if (params.extraLineStr) {
    //     System.out.println("In workflow main: printing an extra line")
    //     System.out.println("Here is the extra line")
    // } else {
    //     System.out.println("In workflow main: NOT printing an extra line")
    // }
    
    // // read from TSV with header
    // // when you have a header, you can pull the column using foo.colName where foo is the variable you're using in map
    // // the column names are keys
    // greeting_ch = channel.fromPath(params.input_tsv)
    //     .splitCsv( header:true, sep:'\t' ) // the rows become arrays
    //     .map { row -> row.greeting } // instead of flatten, which would pass every element of the array, use map which is much more flexible

    // extraLine_ch = channel.fromPath(params.input_tsv)
    //     .splitCsv( header:true, sep:'\t' ) // the rows become arrays
    //     .map { row -> row.extraLine } // instead of flatten, which would pass every element of the array, use map which is much more flexible

    // // emit a greeting
    // sayHello(greeting_ch, extraLine_ch.map { it.toString() })

    // // convert greeting to uppercase
    // convertToUpper(sayHello.out) // takes the output of the previous process; sayHello.out is a channel.

    // // combine the outputs
    // convertToUpper.out.view { contents -> "Before collect: $contents" }
    // convertToUpper.out.collect().view { contents -> "After collect: $contents" }
    // // collectGreetings(convertToUpper.out.collect()) //.view { contents -> "After collect: $contents" } // the output of the prev process came in 3 parts, so the channel operator .collect() is used to combine them into one
    // // collectGreetings(convertToUpper.out) // .view { contents -> "After collect: $contents" } // if you don't collect the inputs, three separate calls are made to collectGreetings 

    // // passing multiple parameters to a function:
    // // VERY important that the order in which you pass the arguments to the process matches the order in which they're specified in the input block of that process; it uses positional matching
    // collectGreetings(
    //     convertToUpper.out.collect(), 
    //     params.batch
    // )

    // // run cowpy
    // cowpy(collectGreetings.out.outfile, params.character)

    // // run the process I designed to test running Python code in Nextflow
    // runPython(params.myString)

    publish:
    blast_output = runBlast.out.blast
    currentBlastFile = currentBlast
    // final_results = filtered_blast
    alignedFastaFile = alignedFasta
    metadataFile = processMetadata.out
}

output {
    // only the final outputs will be copied; the rest are published as soft links.
    blast_output {
        path { "blastFiles/intermediates" }
    }

    // publish the BLASTS from simple filtering in blastFiles/intermediates too

    currentBlastFile {
        path { "blastFiles" }
        mode 'copy' // for debug
    }

    // publish FASTA later, once you've actually initialized the channel in the script
    // with mode 'copy'

    alignedFastaFile {
        mode 'copy'
    }

    metadataFile {
        mode 'copy'
    }

    // final_results {
    //     path { "final" }
    //     mode 'copy'
    // }

}