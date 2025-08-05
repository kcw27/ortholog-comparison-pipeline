### Katherine Wang, 6/20/25
# For a given set of locus tags, produces metadata tables which (in the R segment of the pipeline) are
# joined with sequences from the corresponding FASTA, and used to generate annotations for iTOL trees.
# For datasets that have gone through synteny filtering, we recommend using the raw metadata file
# produced by find_synteny_hits.sh (synteny_summary.tsv in the specified outdir)
# rather than querying NCBI, as the raw metadata file is more informative than the NCBI query results.

import pandas as pd # working with dataframes
import os # for working with directories
import sys # raise error, then exit

# for NCBI queries
from Bio import Entrez, SeqIO
from urllib.error import HTTPError
import time

# for isolation source rescue and categorization
import re


def metadata_wrapper(locus_tags_file, outname, category_file, subcategory_file, raw_metadata_file="", db_to_search="protein", locus_tag_col="protein_id", email=""):
	'''
	Inputs: file containing newline-separated list of locus tags (protein or nucleotide IDs), 
	output filename, keyword file for category terms, keyword file for subcategory terms,
	(optional; no default value) raw metadata file (produced by find_synteny_hits.sh,
	which is run as part of synteny_wrapper.sh) to process into a clean metadata file, 
	(optional; default value "protein") the NCBI database to search if fetching metadata from NCBI,
	(optional; default value "protein_id") the locus tag col by which to filter the metadata to unique records,
	(optional) an email address to use for Entrez.
	Keyword file format: a newline-separated file in which sections are delimited by "***",
	and within each section, the first line is the category/subcategory
	while the subsequent lines are keywords associated with that category/subcategory.
	If find_synteny_hits.sh was run upstream of producing locus_tags_file, we recommend using the 
	raw metadata file rather than querying NCBI, as the raw metadata file is more informative.

	Output: a metadata file with the following columns: protein_id, nucleotide_id, genome_id,
	organism, isolation_source, titles, rescued_source (added by rescue_source()),
	category (added by categorize()), subcategory (added by categorize())
	If raw metadata file is provided, processes that file.
	If no raw metadata file is provided, queries NCBI to produce a metadata file.
	Performs isolation source rescue, then categorizes the isolation sources.
	Saves the metadata file to the output filename.
	'''
	# get outdir from outname
	outdir = os.path.dirname(outname)

	# make outdir if it doesn't already exist
	os.makedirs(outdir, exist_ok=True)

	# get a metadata file to build upon
	if raw_metadata_file == "": # if no raw metadata was provided, get metadata from NCBI
		metadata_origin = "fetched"
		fetch_metadata(locus_tags_file, outname, db_to_search)
		# outname will be overwritten to add columns for rescued_source, category, subcategory
	else: # raw metadata was provided; reformat it for the sake of consistency
		metadata_origin = "synteny_summary"
		process_raw_metadata(locus_tags_file, outname, raw_metadata_file)
		# outname will be overwritten to add columns for rescued_source, category, subcategory

	# rescue isolation source from the titles column where applicable
	rescue_source(outname, metadata_origin)

	# categorize isolation sources, using rescued_source as backup when isolation_source not available
	categorize(outname, category_file, subcategory_file)
	print(f"Wrote processed metadata file to {outname}")

	# New addition: make a copy of the metadata file with only one row per record
	print(f"Filtering {outname} to unique records by the {locus_tag_col} column...")
	prefix, extension = os.path.splitext(outname) # prefix is everything before the last "."
	unique_name = f"{prefix}_unique{extension}"
	get_unique_records(outname, locus_tag_col, unique_name)
	print(f"Wrote processed metadata file to {unique_name}")

def fetch_metadata(locus_tags_file, outname, db_to_search="protein", email=""):
	'''
	Inputs: file containing newline-separated list of locus tags, output filename, (optional) the NCBI database to search (either protein or nucleotide)
	(e.g. if the initial BLAST was on a protein database, then the locus tag list is of protein IDs, so search the protein database)
	Output: metadata with the following columns: protein_id, nucleotide_id, genome_id, organism, isolation_source, titles
	Note that there's no rescued_source column; this is handled separately by rescue_source()
	'''
	# 7/11/25: added sequencing technology metadata; still need to test this function

	start_time = time.time()

	if db_to_search not in ["protein", "nucleotide"]:
		sys.exit("db_to_search must be either 'protein' or 'nucleotide'. Exiting fetch_metadata().")

	Entrez.email = email # replace with your own email address

	locus = []
	organism = []
	isolation_source = []
	titles = []
	source_genome = [] # in protein records, db_source indicates the genome the protein came from
	genome_id = [] # in some nucleotide records, the GCF can be found
	sequencing_technology = []

	with open(locus_tags_file, 'r') as f:
		txt = f.read()
	locus_tags = txt.split() # split by whitespace, including newlines; txt is newline-delimited

	print(f"Will write any errors to fetch_metadata_log.txt in the current directory: {os.getcwd()}")
	log = open("fetch_metadata_log.txt", "w")
	log.write(f"Inputs to fetch_metadata() from metadata_processing.py: {locus_tags_file}, {outname}, {db_to_search}\n")
	log.write(f"Current time: {time.ctime(time.time())}\n") # time.ctime() for human-readable time

	for locus_tag in locus_tags:
		print(f"Locus tag: {locus_tag}")
		try:
			handle = Entrez.esearch(db=db_to_search, term=locus_tag)
			record = Entrez.read(handle)
			handle.close()

		except HTTPError as e:
			if e.code == 400:
				log.write(f"Locus: {locus_tag} Bad Request Error: {e}\n")
			else:
				log.write(f"Locus: {locus_tag} HTTP Error: {e}\n")

		except Exception as e:
			log.write(f"Locus: {locus_tag} An error occurred: {e}\n")

		finally:
			time.sleep(1) # Add delay to respect NCBI rate limits

		if record["IdList"]:
			for accession in record["IdList"]:
				try:
					handle = Entrez.efetch(db=db_to_search, id=accession, rettype="gb", retmode="text")
					gb_record = SeqIO.read(handle, "genbank")
					handle.close()

					# now add a row (i.e. append to each column that will go in the dataframe)
					locus.append(locus_tag)
					
					anno = gb_record.annotations # some records lack certain keys like 'references', so every time you try to access a value in anno, use get()
					info = gb_record.features[0].qualifiers

					organism.append(info["organism"][0]) # it's a list containing a single string, so access that string

					# if isolation source(s) found, it's formatted as a list.
					# Parse it to turn it into a string.
					iso_sources = info.get("isolation_source", "") # it's a string if not found
					if isinstance(iso_sources, list): # i.e. isolation source was found
						iso_sources = ", ".join(iso_sources)
					isolation_source.append(iso_sources) # isolation source might not exist as a key in the dict; if not, just write an empty string

					titles.append(anno.get('references', ''))
					source_genome.append(anno.get("db_source", "")) # if protein record, this could be a way to get the source genome

					# looking for information that is two levels deep in a nested dict if it exists,
					# hence the two levels of the get method
					# can't assume that 'structured_comment' always exists as a key in anno
					structured_comment = anno.get('structured_comment', {}) # default to empty dict to ensure that get() works on structured_comment
					genome_anno = structured_comment.get('Genome-Annotation-Data', {}) # default to empty dict to ensure that get() works on genome_anno
					genome_id.append(genome_anno.get('Annotation Name', ''))

					# also get the sequencing technology
					genome_assembly = structured_comment.get('Genome-Assembly-Data', {}) # default to empty dict to ensure that get() works on genome_assembly
					sequencing_technology.append(genome_assembly.get('Sequencing Technology', 'NA'))

				except HTTPError as e:
					if e.code == 400:
						log.write(f"Accession: {accession} Bad Request Error: {e}")
					else:
						log.write(f"Accession: {accession} HTTP Error: {e}")
				except Exception as e:
					log.write(f"Locus: {locus_tag} An error occurred: {e}\n")
		else: # didn't find any records matching the locus tag.
			log.write(f"No {db_to_search} records found for {locus_tag}\n")

	# write data to file
	empty_list = [""] * len(locus) # use an empty list to stand in for unobtainable info depending on the db searched

	if db_to_search=="protein":
		df = pd.DataFrame(list(zip(locus, empty_list, source_genome, organism, isolation_source, titles, sequencing_technology)), 
			columns=['protein_id', 'nucleotide_id', 'genome_id', 'organism', 'isolation_source', 'titles', 'sequencing_technology'])
		list_of_lists=[locus, empty_list, source_genome, organism, isolation_source, titles, sequencing_technology]
	elif db_to_search=="nucleotide":
		df = pd.DataFrame(list(zip(empty_list, locus, genome_id, organism, isolation_source, titles, sequencing_technology)), 
			columns=['protein_id', 'nucleotide_id', 'genome_id', 'organism', 'isolation_source', 'titles', 'sequencing_technology'])
		list_of_lists=[empty_list, locus, genome_id, organism, isolation_source, titles, sequencing_technology]

	# At the end, check that all lists are the same length; print an error message otherwise
	if len(set(len(x) for x in list_of_lists)) > 1:
		print("Lists are not all the same length. Please double-check your data.")
		log.write("Lists are not all the same length. Please double-check your data.\n")
		for t in list_of_lists:
			print(len(t), t[:5])
			log.write(f"{len(t)} {t[:5]} \n")

	df.to_csv(outname, sep='\t', index=False)

	print(f"Finished writing metadata to {outname}!")
	end_time = time.time()
	print(f"Time elapsed in seconds: {end_time - start_time}")
	log.write(f"Time elapsed in seconds: {end_time - start_time}")
	log.close()


def process_raw_metadata(locus_tags_file, outname, raw_metadata_file):
	'''
	Inputs: newline-separated list of locus tags, output filename, tab-separated raw metadata file.
	Raw metadata is assumed to be produced by find_synteny_hits.sh (please refer to synteny_summary.tsv),
	or otherwise has columns ordered like this:
	genome_id, contig_id, organism, isolation_source, titles, locus_tag, protein_id, protein_sequence, sequencing_technology

	Output: metadata with the following columns: protein_id, nucleotide_id, genome_id, organism, isolation_source, titles, sequencing_technology
	Note that there's no rescued_source column; this is handled separately by rescue_source()
	'''
	# Read input files
	with open(locus_tags_file, 'r') as f:
		txt = f.read()
	locus_tags = txt.split() # split by whitespace, including newlines; txt is newline-delimited
	print(f'There are {len(locus_tags)} locus tags in the locus tags file. Here are the first few: {locus_tags[:10]}')

	raw = pd.read_csv(raw_metadata_file, sep='\t', 
		names=["genome_id", "contig_id", "organism", "isolation_source", 
		"titles", "nucleotide_id", "protein_id", "protein_sequence", "sequencing_technology"])
		# The input file is assumed to follow the format specified by the names list

	# Keep rows if the protein_id is in the locus_tags list.
	# Clarification: locus_tags list contains protein IDs,
	# hence the usage of the protein_id column from the raw metadata.
	filtered = raw[raw["protein_id"].isin(locus_tags)]
	print(f'There were {len(raw)} lines in the original metadata file, and {len(filtered)} remaining after filtering to rows for which the protein_id column contained values in the locus tags list.')

	# count each protein id in the raw metadata
	protein_id_freqs = {}
	for id in raw["protein_id"]:
		protein_id_freqs[id]=protein_id_freqs.get(id, 0)+1
	non_unique = {id for id in protein_id_freqs.keys() if protein_id_freqs[id]>1}
	d = {id:protein_id_freqs[id] for id in non_unique}

	# count each protein id in the filtered metadata
	protein_id_freqs = {}
	for id in filtered["protein_id"]:
		protein_id_freqs[id]=protein_id_freqs.get(id, 0)+1
	non_unique = {id for id in protein_id_freqs.keys() if protein_id_freqs[id]>1}
	d = {id:protein_id_freqs[id] for id in non_unique}

	# next step: rearrange the columns in the filtered df
	# columns from raw metadata that are left out: contig_id, protein_sequence
	new_order = ["protein_id", "nucleotide_id", "genome_id", "organism", "isolation_source", "titles", "sequencing_technology"]
	filtered = filtered[new_order]

	# Finally, save the filtered df to a file with the name outname.
	filtered.to_csv(outname, sep="\t", index=False) # tab separator in case of commas; no rownames

	print(f"Saved file to {outname}")


def rescue_source(metadata_file, metadata_origin, outname = ""):
	'''
	Isolation source rescue function.
	Inputs: tab-separated processed metadata file (obtained either from process_raw_metadata() 
	for raw metadata obtained as a synteny search summary, or from fetch_metadata() 
	for metadata obtained directly from NCBI), metadata_origin ("synteny_summary" or "fetched"),
	which determines the approach used to pull isolation sources out of the metadata,
	and optionally an output filename (will overwrite the input metadata file if none provided).
	Output: adds a rescued_source column to the metadata, using string following "isolated from"
	in the titles column of the metadata.
	'''
	# check if value of metadata_origin is legitimate
	if metadata_origin not in ["synteny_summary", "fetched"]:
		sys.exit("metadata_origin must be either 'synteny_summary' or 'fetched'. Exiting rescue_source().")

	# if no outname provided, overwrite the input metadata file
	if outname == "":
		outname = metadata_file

	# initialize rescued_source column
	df = pd.read_csv(metadata_file, sep="\t") #
	rescued_source = [""] * len(df.index) # give it the same number of elements as rows in df

	# The string may be either "Isolated from" or "isolated from"
	# so start the match with a case-insensitive match to these strings
	if metadata_origin == "synteny_summary":
		# titles are separated by the string "TITLE"
		# so the match should go either until TITLE or to the end of the string
		pattern = re.compile(r'(?i)isolated from\s+(.*?)(?=\s+TITLE|$)')
		# (?i): Enables case-insensitive matching for just "isolated from".
		# \s+: Matches one or more whitespace characters after "isolated from".
		# The whitespace before and after the match are not included in the match.
		# (.*?): Non-greedy capture of any character (i.e. the match target).
		# (?=\s+TITLE|$): Positive lookahead; stop matching before " TITLE" or at end of string.

	elif metadata_origin == "fetched":
		# titles in reference column reliably end with ', even the last title
		# because the last title is followed by ', ...)
		pattern = re.compile(r'(?i)isolated from\s+(.*?)(?=\',|$)') # goes until ',

	# iterate thru df rows and apply the pattern to look for isolation sources within titles
	for i in range(len(df)):
		refs = safe_str(df["titles"][i])

		# use regex to find an isolation source, if any
		match = re.search(pattern, refs) # re.search to find the first match, if any
		if match: # if a match was found
			rescued_source[i] = match.group() # get the string from the match

	# add new column to df
	df["rescued_source"] = rescued_source

	# save this new version of df to outname
	df.to_csv(outname, sep="\t", index=False)
	print("Metadata with rescued isolation sources written to file", outname)



def categorize(metadata_file, category_file, subcategory_file, outname = ""):
	'''
	Source categorization function.
	Inputs: processed metadata file which may or may not have a rescued_source column,
	keyword file for category terms, keyword file for subcategory terms,
	and optionally an output filename (will overwrite the input metadata file if none provided).
	Keyword file format: a newline-separated file in which sections are delimited by "***",
	and within each section, the first line is the category/subcategory
	while the subsequent lines are keywords associated with that category/subcategory.

	Output: adds a category column and subcategory column to the metadata,
	using regular expressions to match terms from isolation_source (and rescued_source if
	available and if nothing in isolation_source) to category/subcategory keywords.
	'''
	df = pd.read_csv(metadata_file, sep="\t")

	# if no outname provided, overwrite the input metadata file
	if outname == "":
		outname = metadata_file

	# if rescued_source column exists, assume that it should be sed to categorize
	# in the case that isolation_source is empty for that row
	rescue = hasattr(df, "rescued_source") # Boolean- does df have a rescued_source column?
	if rescue:
		print(f"Will rescue missing categories and subcategories for {metadata_file} using rescued_source column.")

	# get dictionaries for category and subcategory
	cat_dict = parse_keywords(category_file)
	subcat_dict = parse_keywords(subcategory_file)

	# combine keywords (keys of the input dictionaries) into regex expressions
	cat_kw = list(cat_dict.keys())
	cat_pattern = re.compile(r'\b(' + '|'.join(map(re.escape, cat_kw)) + r')') # all strings are converted to lowercase for consistency anyway, so no need to ignore case

	subcat_kw = list(subcat_dict.keys())
	subcat_pattern = re.compile(r'\b(' + '|'.join(map(re.escape, subcat_kw)) + r')') # all strings are converted to lowercase for consistency anyway, so no need to ignore case

	# initialize category columns with empty strings; will be overwritten if regex finds a match
	category = [""] * len(df.index)
	subcategory = [""] * len(df.index)

	# the empty strings might have been causing issues, so use non-empty instead
	category = ["no category"] * len(df.index)
	subcategory = ["no subcategory"] * len(df.index)

	# iterate thru df to determine categories and subcategories
	for i in range(len(df)):
		source = safe_str(df["isolation_source"][i]).lower() # must convert to lowercase since all keys in the dict are lowercase

		# update the category
		# if you want this to be more accurate, you could get the whole list and assign whichever associated
		# category is most abundant among these keywords as the new value at category[i]
		cat_match = cat_pattern.search(source)

		if cat_match: # if a match was found
			category[i] = cat_dict[cat_match.group()] # general category corresponding to the keyword found
		elif rescue: # try to rescue category with rescued_source
			rescued_str = safe_str(df["rescued_source"][i]).lower() # must convert to lowercase since all keys in the dict are lowercase
			cat_rescue = cat_pattern.search(rescued_str)

			if cat_rescue: # successfully rescued category from rescued_str
				category[i] = cat_dict[cat_rescue.group()] # overwrite the empty string
				#print(f"\tRescued category '{category[i]}' from '{rescued_str}'; isolation_source was {source}.")

		# update the specific category
		subcat_match = subcat_pattern.search(source)

		if subcat_match:
			subcategory[i] = subcat_dict[subcat_match.group()] # specific category corresponding to the keyword found
		elif rescue: # try to rescue subcategory with rescued_source
			rescued_str = safe_str(df["rescued_source"][i]).lower() # must convert to lowercase since all keys in the dict are lowercase
			subcat_rescue = subcat_pattern.search(rescued_str)

			if subcat_rescue: # successfully rescued subcategory from rescued_str
				subcategory[i] = subcat_dict[subcat_rescue.group()] # overwrite the empty string
				#print(f"\tRescued subcategory '{subcategory[i]}' from '{rescued_str}'; isolation_source was {source}.")

	# add new columns to df
	df["category"] = category
	df["subcategory"] = subcategory

	# save this new version of df to outname
	df.to_csv(outname, sep="\t", index=False)
	print("Isolation source categories and subcategories written to file", outname)

def parse_keywords(keywords_file):
	'''
	Helper function for categorize().
	Input: a newline-separated file in which sections are delimited by "***",
	and within each section, the first line is the category/subcategory
	while the subsequent lines are keywords associated with that category/subcategory.
	Output: returns a dictionary in which the keys are the keywords 
	and the values are the corresponding descriptor (i.e. category/subcategory).
	'''
	d = {}

	# read in the file
	with open(keywords_file, 'r') as f:
		txt = f.read().lower() # convert all to lowercase for consistency

	# separate into sections
	sections = txt.split("***")

	# separate each section into descriptor and keywords
	for section in sections:
		section = section.strip() # remove flanking whitespace, including newlines
		terms = section.split("\n")
		descriptor = terms[0]
		keywords = terms[1:]

		# add entries to dictionary
		for kw in keywords:
			d[kw] = descriptor

	return d

def safe_str(x):
	'''
	Helper function to ensure that when a value is a null of some kind, 
	it returns an empty string rather than the null value.
	Used for safe regular expression searching.
	'''
	return str(x) if pd.notnull(x) else ""

def get_locus_tags_from_fasta(fasta_file, outname):
	'''
	If you're starting off with a FASTA file and you want to fetch its metadata from NCBI,
	use this function to produce a text file containing the list of locus tags from that FASTA
	This is used as input for fetch_metadata().
	'''
	with open(fasta_file, 'r') as f:
		txt = f.read()
	records = txt.split(">") # each record in the FASTA begins with a caret
	locus_tags = [r.split()[0] for r in records if r != ""]
	# assume the locus tag is everything between the caret and the first space
	# remove empty items from the list to avoid errors

	with open(outname, 'w') as f:
		for tag in locus_tags:
			f.write(f"{tag}\n")
	print(f"Wrote locus tags from {fasta_file} to {outname}")

def get_unique_records(metadata_file, locus_tag_col, outname):
	'''
	Input: metadata file, name of column to treat as locus tags, output file name.
	Output: writes the first record for each locus tag to outname.
	'''
	df = pd.read_csv(metadata_file, sep="\t")
	original_col_order = list(df.columns)

	# unique_locus_tags = set(df[locus_tag_col])
	# print(unique_locus_tags)

	grouped_df = df.groupby(locus_tag_col).first()

	# Reset indices to match format
	first_values = grouped_df.reset_index()
	# print(first_values)

	# The operations above may have changed the order of the columns, so reset the order
	first_values = first_values[original_col_order]

	# Save the filtered df
	first_values.to_csv(outname, sep="\t", index=False)
	print(f"Wrote unique values of {metadata_file} to {outname}")

# Main script
if __name__ == "__main__":
	'''
	# Move from the scripts directory to the data directory to access input data
	os.chdir("../data") # TODO: works on my local device, but change this to work on any device
	
	# run metadata_wrapper(), saving outputs in data/processed_metadata/
	# for the two synteny-filtered datasets:
	metadata_wrapper("locus_tags/locus_tags_65_66_67.txt", "processed_metadata/metadata_65_66_67.tsv",
		"keywords/category_keywords.txt", "keywords/subcategory_keywords.txt",
		raw_metadata_file="raw_metadata/synteny_summary_65_66_67.tsv")

	metadata_wrapper("locus_tags/locus_tags_65_67.txt", "processed_metadata/metadata_65_67.tsv",
		"keywords/category_keywords.txt", "keywords/subcategory_keywords.txt",
		raw_metadata_file="raw_metadata/synteny_summary_65_67.tsv")
	'''
	'''
	# for the synteny-unfiltered datset: I ran fetch_metadata() on the cepacia server,
	# saving the output to "processed_metadata/metadata_unfiltered.tsv",
	# and I'm currently running this part of the analysis on my local device.
	# I just need to run rescue_source() and categorize().
	rescue_source("processed_metadata/metadata_unfiltered.tsv", "fetched")
	categorize("processed_metadata/metadata_unfiltered.tsv",
		"keywords/category_keywords.txt", "keywords/subcategory_keywords.txt")

	# test the fetch version of metadata_wrapper() using a smaller dataset
	metadata_wrapper("locus_tags/locus_tags_test.txt", "processed_metadata/metadata_unfiltered_test.tsv",
		"keywords/category_keywords.txt", "keywords/subcategory_keywords.txt", db_to_search="protein", email="kcw2@andrew.cmu.edu")
	'''

	'''
	# run this pipeline on POG002101.fasta to produce an up-to-date metadata file
	os.chdir("C:/Users/achro/OneDrive/Desktop/CMU/Spring 2025/Armbruster Lab research/ortholog_comparison_refactor/test_data")
	fasta_file = "orthologs/POG002101.fasta"
	locus_tags_file = "locus_tags.txt"
	outname = "processed_metadata/metadata.tsv"
	category_file = "keywords/category_keywords.txt"
	subcategory_file = "keywords/subcategory_keywords.txt" # initially used category_keywords.txt on accident

	# get locus tags from POG002101.fasta, and save to file
	get_locus_tags_from_fasta(fasta_file, locus_tags_file)

	# call metadata_wrapper()
	metadata_wrapper(locus_tags_file, outname, category_file, subcategory_file, db_to_search="nucleotide", email="kcw2@andrew.cmu.edu")
	'''
	'''
	# filter POG002101 metadata to the first record for each locus tag (nucleotide_id)
	os.chdir("C:/Users/achro/OneDrive/Desktop/CMU/Spring 2025/Armbruster Lab research/ortholog_comparison_refactor/test_data")
	category_file = "keywords/category_keywords.txt"
	subcategory_file = "keywords/subcategory_keywords.txt"
	metadata_file = "processed_metadata/metadata.tsv"
	locus_tag_col = "nucleotide_id"
	outname = "processed_metadata/metadata_unique.tsv"

	categorize(metadata_file, category_file, subcategory_file)
	get_unique_records(metadata_file, locus_tag_col, outname)
	'''
	'''
	# get unique records for the synteny unfiltered dataset, since at the time I retrieved the metadata, the pipeline didn't include get_unique_records()
	get_unique_records("processed_metadata/metadata_unfiltered.tsv", "protein_id", "processed_metadata/metadata_unfiltered_unique.tsv")
	# redo categorization because the initial categorization was done with an earlier version of the function
	categorize("processed_metadata/metadata_unfiltered_unique.tsv", "keywords/category_keywords.txt", "keywords/subcategory_keywords.txt")
	'''

	'''
	# Process raw metadata from Fha1 synteny search
	# run metadata_wrapper(), saving outputs in data/processed_metadata/
	metadata_wrapper("../data/locus_tags/locus_tags_Fha1_orthologs_synteny_filtered.txt", "../data/processed_metadata/metadata_Fha1.tsv",
		"../data/keywords/category_keywords.txt", "../data/keywords/subcategory_keywords.txt",
		raw_metadata_file="../data/raw_metadata/synteny_summary_Fha1.tsv")
	'''

	# categorizing the isolation sources in the fha1 Pa-only raw metadata output from my modified blast2gen.py
	# (I forgot to retrieve titles so there's no source rescue here)
	#metadata_file = "../data/raw_metadata/fha1_genome_info_paOnly.tsv"
	metadata_file = "../data/raw_metadata/fha1_genome_info_paOnly_top.tsv"
	category_file = "../data/keywords/category_keywords.txt"
	subcategory_file = "../data/keywords/subcategory_keywords.txt"
	#outname = "../data/processed_metadata/fha1_genome_info_paOnly_categorized.tsv"
	outname = "../data/processed_metadata/fha1_genome_info_paOnly_categorized_top.tsv"
	categorize(metadata_file, category_file, subcategory_file, outname)