#!/usr/bin/env nextflow
include { runBlast } from './modules/runBlast.nf'
include { foo; bar; handleGenomeFlag } from './modules/filtering.nf'
include { foo_metadata; bar_metadata } from './modules/retrieveMetadata.nf' 
include { FILTER1 } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter1.nf'
include { FILTER2 } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter2.nf'
include { FILTER3 } from '/home/kcw2/test_scripts/nextflow_tests/modules/filter3.nf'
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

    // metadata categorization
    category: Path
    subcategory: Path

    // align using this pipeline?
    align: Boolean
}

// This one almost worked, but hasGenomeIDs wasn't actually a Boolean, but rather a channel
workflow applyFiltersWithMetadata {
    take:
    blastWithFlag   // each element is [hasGenomeIDs (bool), blastFile]

    main:
    // Split the flag out so conditional logic can use it
    hasGenomeIDs = blastWithFlag.map { flag, _blast -> flag }.first()
    hasGenomeIDs.view{ contents -> "hasGenomeIDs contents: $contents" }
    def current  = blastWithFlag.map { _flag, blast -> blast }
    // This will wait for the value and assign it as a string
    def myString = hasGenomeIDs.value as String
    
    
    // Filter 1
    if (params.gate1) {
        println("Filter 1 applied")
        current = FILTER1(current)
    }
    
    // Filter 2
    if (params.gate2) {
        println("Filter 2 applied")
        current = FILTER2(current)
    }
    
    // Filter 3
    if (params.gate3) {
        println("Filter 3 applied")
        println "The value is: ${myString}"
        if (myString == "false") {
            println("No genome IDs; running foo_metadata")
            current = foo_metadata(current)
            hasMetadata = true
        }
        current = FILTER3(current)
    }
    
    // Filter 4
    if (params.gate4) {
        println("Filter 4 applied")
        current = FILTER4(current)
    }
    
    // Final step
    current = FINAL(current)
    
    emit:
    filtered_blast = current
}

// // this one seems to work, but doesn't guarantee that foo_metadata runs in between FILTER2 and FILTER3
// // Also, as I've discovered, this runs foo_metadata even if filter3 isn't applied...
// workflow applyFiltersWithMetadata {
//     take:
//     blastWithFlag   // each element is [hasGenomeIDs (bool), blastFile]

//     main:
//     // Split into two branches based on the flag
//     def (branchWithGenomes, branchWithoutGenomes) = blastWithFlag.branch {
//         withGenomes: it[0] == true
//         withoutGenomes: it[0] == false
//     }
    
//     // Process branch without genomes through foo_metadata
//     // Extract just the blast files for each branch
//     def blastsWithGenomes = branchWithGenomes.map { flag, blast -> blast }
//     def blastsWithoutGenomes = branchWithoutGenomes.map { flag, blast -> blast }
    
//     // Apply foo_metadata process to the branch without genome IDs
//     def processedWithoutGenomes = foo_metadata(blastsWithoutGenomes)
    
//     // Combine the branches back
//     def current = blastsWithGenomes.mix(processedWithoutGenomes)
    
//     // Now apply filters that don't depend on the flag
//     if (params.gate1) {
//         println("Filter 1 applied")
//         current = FILTER1(current)
//     }
    
//     if (params.gate2) {
//         println("Filter 2 applied")
//         current = FILTER2(current)
//     }
    
//     if (params.gate3) {
//         println("Filter 3 applied")
//         current = FILTER3(current)
//     }
    
//     if (params.gate4) {
//         println("Filter 4 applied")
//         current = FILTER4(current)
//     }
    
//     current = FINAL(current)
    
//     emit:
//     filtered_blast = current
// }




workflow {

    main:
    currentBlast = Channel.empty() // not strictly necessary because it's guaranteed to populate right after, but this is for the sake of clarity
    // currentMetadata = Channel.empty() // make metadataFile only after the filtering gauntlet is complete, and only if synteny search wasn't run.
    // currentFasta = Channel.empty() // Actually, I should keep this empty so that if we opt in for alignment, it's restricted until after we have a FASTA at all
    // either make the FASTA from the metadata file (maybe I shouldn't initialize that then...), or if synteny search was performed, use the FASTA from that 

    hasMetadata = Channel.of("false") // keeping it as a string because processes are buggy with Booleans
    // hasFasta = Channel.of("false") // actually, the better way to ensure this is just to rely on the presence/absence of the currentFasta channel

    syntenyInputList = [params.genomeDBsynteny, params.syntenyInput, params.hmmsList, params.hmmsDir, params.hmmsMetadata]
    if (syntenyInputList.any { it == null }) {
        println "List contains at least one null value"
        syntenySearchFlag = false
    } else {
        println "List contains no null values"
        syntenySearchFlag = true
    }
    println("syntenySearchFlag: '${syntenySearchFlag}'")

    if (syntenySearchFlag) {
        foo()
    } else {
        bar()
    }

    // Get the BLAST file
    runBlast(params.queryFasta, params.blastPath, params.blastName)
    currentBlast = runBlast.out.blast

    // test modification to the current BLAST
    bar_metadata(runBlast.out.blast, runBlast.out.hasGenomeIDsFile)
    currentBlast = bar_metadata.out

    hasGenomeIDs = runBlast.out.hasGenomeIDsFile
        .map { file -> file.text.trim() == "true" }

    runBlast.out.hasGenomeIDsFile.view { hasGenomeIDsFile -> 
        def contents = hasGenomeIDsFile.text  // This reads the file content
        println "Current contents: $contents"
    }

    // First, check what's actually in the channel
    hasGenomeIDs.view { "Channel contains: $it" }

    // Also check if the file exists and has content
    runBlast.out.hasGenomeIDsFile.view { "File path: $it" }

    // reading this file: need strip to remove the newline character baked into the file 
    runBlast.out.hasGenomeIDsFile.text.strip()
        .branch { v ->
            yes: v == "true"
            no: v == "false"
            other: true //fall-through
        }
        .set { result }

    result_y = result.yes.view { v -> "$v is yes" }
    result_n = result.no.view { v -> "$v is no" }
    result_o = result.other.view { v -> "$v is other" }
    // hasGenomeIDs = 

        
    // Print the actual file content after the process runs
    runBlast.out.hasGenomeIDsFile
        .map { file ->
            println "=== FILE CONTENTS ==="
            println "Content: '${file.text}'"
            println "Trimmed: '${file.text.trim()}'"
            println "Length: ${file.text.trim().length()}"
            println "====================="
            return file.text.trim()
        }
        .view()

    // // Create a value channel with the string
    // def genomeFlagChannel = runBlast.out.hasGenomeIDsFile
    //     .map { file -> file.text.trim() }
    
    // // Now capture it as a string using .value inside the workflow
    // def genomeFlagString = genomeFlagChannel.value
    
    // // Now you have it as a regular string
    // println "The flag is: '$genomeFlagString'"  // Should print: The flag is: 'false'
    
    
    // // Use it in conditional logic
    // if (genomeFlagString == "true") {
    //     println "Genome IDs present"
    //     // Do something with genome IDs
    // } else {
    //     println "No genome IDs found"
    //     // Handle the false case
    // }
 
    handleGenomeFlag(runBlast.out.hasGenomeIDsFile.map { file -> file.text.trim() })

    // Combine blast records with the flag as a tuple: [hasGenomeIDs, blastFile]
    blastWithFlag = hasGenomeIDs.combine(currentBlast)
    // blastWithFlag = runBlast.out.hasGenomeIDsFile.combine(currentBlast)

    filtered_blast = applyFiltersWithMetadata(blastWithFlag)

    println("Has metadata:")
    println(hasMetadata)

    // currentBlast.view { blastFile -> 
    //     def contents = blastFile.text  // This reads the file content
    //     println "Current contents: $contents"
    // }

    // use .branch for branching paths


    
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
    final_results = filtered_blast

}

output {
    // only the final outputs will be copied; the rest are published as soft links.
    blast_output {
        path { "blastFiles/intermediates" }
    }

    // publish the BLASTS from simple filtering in blastFiles/intermediates too

    currentBlastFile {
        path { "blastFiles" }
    }

    final_results {
        path { "final" }
        mode 'copy'
    }

}