#!/usr/bin/env python3

### Define functions
yes = ["y", "yes"]
def format_input(param_name, param_value):
    # returns the string to write to write to the file.
    return f"{param_name}: '{param_value.strip().strip('"')}'\n"

def preamble():
    print("Welcome to the Metadata Magnet input file setup wizard. " \
        "\nCompletion of this form will generate a params file (YAML format) for your pipeline run in Nextflow. " \
        "(nextflow run option: -params-file my_params_file.yaml)" \
        "\nThe inputs you enter here will be used by the pipeline to produce a FASTA that represents orthologs of " \
        "your protein of interest (with various kinds of optional filtering) as well as a metadata file. " \
        "\nThe FASTA and metadata outputs may then be used as inputs to the GUI for statistical tests and visualization.")

    outfile = input("First, please enter the file path to which you would like to save the params YAML file: ").strip()
    with open(outfile, "w") as file:
        pass # create an empty file, overwriting the file at outfile if it already exists 
    print(f"Params will be saved to {outfile}." \
          "For the rest of these prompts, when entering filepaths, if entering relative paths, please enter them relative to the nextflow directory where main.nf is found. " \
          "(Alternatively, use absolute paths.)\n")
    return outfile

def blast(outfile):
    print("Collecting inputs for BLAST.")
    with open(outfile, "a") as file:
        query = input("Please enter the path to the FASTA containing your query sequence. Ensure that the FASTA contains only one sequence. ")
        file.write(format_input("queryFasta", query))

        db = input("Please enter the path to the directory containing the BLAST db files. ")
        file.write(format_input("blastPath", db))

        title = input("Please enter the name/title of your BLAST database. This is the part of the filename that is shared across the BLAST db files, not including extensions. ")
        file.write(format_input("blastName", title))

        separateIDsString = input("In the BLAST file, do you anticipate protein IDs to be in the format genomeID-proteinID? " \
            "This would be the case if you use the provided script to make a blast DB from a genome DB. (y/n): ")
        if separateIDsString.lower() in yes: # if not, just use the default, which is false
            file.write(format_input("separateIDs", "true"))
            separateIDs = True
        else:
            separateIDs = False

        file.write(f"\n")
    print()
    return separateIDs

def filter(outfile):
    print("Setting filtering settings now.")
    with open(outfile, "a") as file:
        evalueBool = input("Would you like to filter by evalue? (y/n): ")
        if evalueBool.lower() in yes:
            evalueThreshold = input("Please enter an evalue threshold. Sequences with evalue <= this threshold will be kept. For scientific notation, please use e or E, e.g. 1e-30: ")
            file.write(format_input("evalueThreshold", evalueThreshold))

        organismBool = input("Would you like to filter to the top hit (sequence with lowest evalue, keeping only one in the event of a tie) per ORGANISM? (Different strains of the same bacteria are considered separate organisms.) (y/n): ")
        if organismBool.lower() in yes:
            file.write(format_input("filterByOrganism", "true"))

        genomeBool = input("Would you like to filter to the top hit (sequence with lowest evalue, keeping only one in the event of a tie) per GENOME? (y/n): ")
        if genomeBool.lower() in yes:
            file.write(format_input("filterByGenome", "true"))
            filterByGenome = True
        else:
            filterByGenome = False

        syntenyBool = input("Would you like to filter by synteny context? (y/n): ")
        if syntenyBool.lower() in yes:
            print("Please prepare inputs after referring to Pynteny documentation and this pipeline's provided example synteny search inputs.")
            filterBySynteny = True
            syntenyGenomeDB = input("Please enter the path to the genome database directory in which to perform the synteny searches. " \
                "This is a directory in which all subdirectories are named after assembly accessions, "
                "and each subdirectory contains a .gbff.gz file (unzipped is fine too).: ")
            file.write(format_input("genomeDBsynteny", syntenyGenomeDB))

            syntenyInput = input("Please enter the path to a tab-separated table in which column 1 contains synteny structures (refer to documentation for 'pynteny search'), "
                "and column 2 contains the corresponding output directory names as RELATIVE paths. ")
            file.write(format_input("syntenyInput", syntenyInput))

            hmms_list = input("Please enter the path to a file that lists HMMs from your query sequence (do NOT include HMMs from the surrounding genes): ")
            file.write(format_input("hmmsList", hmms_list))

            hmms_dir = input("Please enter the path to the HMM database directory: ")
            file.write(format_input("hmmsDir", hmms_dir))

            hmms_meta = input("Please enter the path to the HMM metadata table: ")
            file.write(format_input("hmmsMetadata", hmms_meta))

            print("After the synteny context search is performed, the outputs of that search will be intersected with the set of sequences from the BLAST "
                "(with any filtering specified earlier), and the intersection will return the subset of sequences from the synteny context search that " \
                "closely resemble sequences among the BLAST hits. More specifically, the synteny hits will be blasted against a database made from the BLAST hits. " \
                "You may now specify filtering thresholds for pident (percent identity) and qcovs (query coverage per subject) for this intersection blast. " \
                "Based on these thresholds, a subset of the hits from the synteny search (i.e. orthologs that appear in the specified synteny context) will be kept." \
                "Refer to BLAST documentation: https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a")
            
            pidentBool = input("Would you like to specify a pident threshold different from the default of 99? (y/n): ")
            if pidentBool.lower() in yes:
                pident = input("Please specify a pident threshold as a number from 0 to 100. " \
                "The filtering will keep only synteny hit sequences with pident >= this threshold in the intersection BLAST. ")
                file.write(format_input("intersectPident", pident))

            qcovsBool = input("Would you like to specify a qcovs threshold different from the default of 90? (y/n): ")
            if qcovsBool.lower() in yes:
                qcovs = input("Please specify a qcovs threshold as a number from 0 to 100. " \
                "The filtering will keep only synteny hit sequences with qcovs >= this threshold in the intersection BLAST. ")
                file.write(format_input("intersectQcovs", qcovs))

            # keepBool = input("Would you like to keep a FASTA of the full dataset of hits that occur in the synteny context? " \
            # "and/or" \
            # "Would you like to keep a copy of the temporary BLAST db made from the input BLAST hits? (Recommended: no) (y/n): ")
            # if keepBool.lower() in yes:
            #     print("Look for these output files in the work dir corresponding to the filterSynteny process. They won't be published to the final output directory.")
            #     file.write(format_input("keepSyntenyFasta", "true"))
        else:
            filterBySynteny = False

        file.write(f"\n")

    print()

    return filterByGenome, filterBySynteny # need this as part of the decision on whether to get metadata
    
def retrieveMetadata(outfile, hasGenomeIDs):
    print("Collecting inputs for metadata retrieval.")
    with open(outfile, "a") as file:
        if hasGenomeIDs:
            dbBool = input("You previously indicated that genome IDs are encoded in your BLAST db, " \
            "which means that it may be possible to perform a local database search. " \
            "This is generally preferable to using NCBI esearch, which is slower due to rate-limiting and has the chance to disconnect and fail." \
            "\nDo you have a genome database that you can use for this search? This might be the database you built the blast DB from using our provided script," \
            "or a collection of genomes downloaded using https://github.com/kblin/ncbi-genome-download encompassing all genomes assembly accessions in your blast DB."
            "\nExpected format: a directory in which all subdirectories are named after assembly accessions, " \
            "and each subdirectory contains a .gbff.gz file (unzipped is fine too)." \
            "\n(y/n): ")
            if dbBool.lower() in yes:
                metadataDB = input("Please enter the path to the genome DB: ")
                file.write(format_input("genomeDBmetadata", metadataDB))
            else:
                print("In that case, you will need to retrieve metadata through NCBI esearch.")
                email = input("Please enter your email so that NCBI can contact you if there's a problem with your requests: ")
                file.write(format_input("entrezEmail", email))
        else:
            print("You will need to retrieve metadata through NCBI esearch.")
            email = input("Please enter your email so that NCBI can contact you if there's a problem with your requests: ")
            file.write(format_input("entrezEmail", email))
        file.write(f"\n")
        
def categorize(outfile):
    print("Collecting inputs for categorization of isolation source metadata.")
    catBool = input("By default, categorization is performed using the keyword maps found in the nextflow/data dir, at category_keywords.txt and subcategory_keywords.txt." \
    "Would you like to provide your own keyword maps instead of using the defaults? (y/n): ")
    if catBool.lower() in yes:
        with open(outfile, "a") as file:
            cat = input("Please enter the path to your category keywords file: ")
            file.write(format_input("category", cat))

            subcat = input("Please enter the path to your subcategory keywords file: ")
            file.write(format_input("subcategory", subcat))

            file.write(f"\n")

def align(outfile):
    print("Collecting inputs for alignment of FASTAs.")
    alignBool = input("Would you like to align your FASTAs automatically using this pipeline? " \
    "MUSCLE/super5 is selected based on the size of the final FASTA." \
    "If you require greater control over the alignment method, please select 'n' and find the FASTAs in the fasta/not_aligned subdirectory of the output directory. " \
    "(y/n): ")
    if alignBool.lower() in yes:
        with open(outfile, "a") as file:
            file.write(format_input("align", "true"))
            # file.write(f"\n")

### Run wizard
def run_wizard():
    outfile = preamble()
    hasGenomeIDs = blast(outfile) # needing to separate IDs is equivalent to stating that genome IDs are encoded in the IDs

    filterNow = input("You may specify in advance the filtering thresholds for the BLAST hits for your query sequence according to certain criteria. " \
        "But would you like to inspect the BLAST hits manually before filtering? (y/n): ")
    if filterNow.lower() in yes:
        with open(outfile, "a") as file:
            file.write(format_input("stopBeforeFiltering", "true"))
        print(f"YAML has been exported to {outfile}. " \
            "First run the pipeline with this inputs yaml file, and inspect the output BLAST file. " \
            "When you are ready to continue, rerun this wizard with the exact same set of arguments as defined previously, "
            "and then rerun the pipeline with the new inputs YAML file and the -resume flag to use cached outputs.")
        return
    
    print("Continuing to retrieve inputs.")
    filteredByGenome, filteredBySynteny = filter(outfile)

    # decide now: do we need to retrieve metadata?
    if (filteredByGenome and not hasGenomeIDs):
        print("Earlier, you indicated that you would like to filter by genome ID, and that there are no genome IDs encoded in your BLAST database. " \
        "As a result, you will need to retrieve metadata for this dataset. (If you have previously retrieved metadata for your dataset using this pipeline, " \
        "you can provide a path to it using the --existingMetadata option when calling nextflow run, but please provide responses to the following questions.)")
        retrieveMetadata(outfile, hasGenomeIDs)
    elif filteredBySynteny:
        print("Synteny search provides metadata, so it is not necessary to retrieve metadata separately under most circumstances. Proceeding to the next set of inputs.\n")
    else:
        retrieveMetadata(outfile, hasGenomeIDs)
        
    # categorize(outfile) # the average user probably won't make their own keyword maps, so use defaults unless using the --category and --subcategory flags

    align(outfile)

    print(f"The input parameters file is complete! {outfile}")
    print("Example run: CD to the nextflow dir (where main.nf is located) and run the following:")
    print(f"nextflow run main.nf -params-file {outfile}")
    print("(Adding the -resume flag is recommended on subsequent runs, and there is no disadvantage to using it even when it does not apply.)")

run_wizard()