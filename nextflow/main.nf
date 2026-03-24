#!/usr/bin/env nextflow
include { runBlast } from './modules/runBlast.nf'
include { foo; bar; filterToEvalue; filterToOrganism; filterToGenome; filterSynteny } from './modules/filtering.nf'
include { foo_metadata as foo_metadata_beforeFiltering; foo_metadata as foo_metadata_afterFiltering; bar_metadata } from './modules/retrieveMetadata.nf' 
include { FILTER1 } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter1.nf'
include { FILTER2 } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter2.nf'
include { filter3NeedSignal as filter3_noMetadata; filter3NeedSignal as filter3_withMetadata } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter3.nf'
include { FILTER4 } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter4.nf'
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

    main:
    // def currentBlast = Channel.empty() // not strictly necessary because it's guaranteed to populate right after, but this is for the sake of clarity
    def metadata = Channel.empty() // will later decide whether we still need to get metadata based on whether this is still empty or not. (Can't initialize with empty strings)
    // May contain multiple values.
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
    
    // Filter 3
    // runBlast.out.hasGenomeIDsFile
    //     .branch { file ->
    //         def content = file.text.trim()
    //         yes: content == "true"
    //         no: content == "false"
    //         other: true
    //     }
    //     .set { result }

    // if (params.gate3) {
    //     println("Filter 3 applied")
    //     result_y = result.yes.view { v -> "$v is yes" }
    //     result_n = result.no.view { v -> "$v is no" }
    //     result_o = result.other.view { v -> "$v is other" }

    //     currentBlast.view { blastFile -> 
    //         def contents = blastFile.text  // This reads the file content
    //         println "Current contents of blastFile: $contents"
    //     }

    //     // yes branch: just run FILTER3
    //     f3_noMetadata_out = filter3_noMetadata(
    //         result.yes,
    //         currentBlast
    //     )
        

    //     // no branch: run foo_metadata and then FILTER3
    //     fooMetadata_out = foo_metadata_beforeFiltering(
    //         result.no,
    //         currentBlast
    //     )

    //     f3_withMetadata_out = filter3_withMetadata(
    //         fooMetadata_out.map { true },
    //         fooMetadata_out
    //     )

    //     // either way, assign to out
    //     current = f3_noMetadata_out.out[0]
    //         .mix(
    //             f3_withMetadata_out.out[0]
    //         )
    //         .view()

    //     currentBlast = current
    // }
    if (params.gate3) {
        println("Filter 3 applied")
        
        // Read the first value from each BLAST file to determine the branch
        def blastWithCheck = currentBlast.map { blastFile ->
            def firstLine = blastFile.text.readLines().first()
            def firstValue = firstLine.split()[0]  // Get first column
            def hasGenomeIDs = (firstValue != "0")  // true if not 0
            return tuple(blastFile, hasGenomeIDs)
        }.view { blastFile, hasGenomeIDs -> "Blast: $blastFile, Has genome IDs: $hasGenomeIDs" }
        
        // Branch based on the check
        def result = blastWithCheck.branch { blastFile, hasGenomeIDs ->
            yes: hasGenomeIDs == true
            no: hasGenomeIDs == false
        }
        
        // Process YES branch
        def yesOutput = filter3_noMetadata(
                result.yes,
                currentBlast
            )
        
        
        // Process NO branch  
        def metadataOut = foo_metadata_beforeFiltering(
                result.no,
                currentBlast
            )
            
            noOutput = filter3_withMetadata(
                metadataOut.map { true },
                metadataOut
            )
        
        // Combine outputs
        current = yesOutput.mix(noOutput)
        currentBlast = current
    }
    
    // Filter 4
    if (params.gate4) {
        println("Filter 4 applied")
        current = FILTER4(currentBlast)
        currentBlast = current
    }
    
    // Final step
    current = FINAL(currentBlast)
    currentBlast = current



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
    metadataFile = metadata
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