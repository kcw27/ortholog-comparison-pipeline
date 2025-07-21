# 7/3/25
# Later, combine this with functions from find_inverted_repeats.py

import re
from Bio import SeqIO
import os

def find_matches(p, string):
	"""
	Input: a compiled regex pattern p, and a string to search.
	Output: a tuple containing the matches (list[str]), starts (list[int]), and ends (list[int]) of matches.
	"""
	#print(string)
	matches = []
	starts = []
	ends = []
	for m in p.finditer(string):
		matches.append(m.group())
		starts.append(m.start())
		ends.append(m.end())

	return (matches, starts, ends)

def summarize_pq_repeats(multifasta, outname, pattern=r"(?:PQ)+"):
	'''
	New version that writes the coordinates of only the longest PQ repeat per sequence. 
	Inputs: a multifasta file in which to find PQ repeats, an output file name,
	(optional) a regular expression pattern in case you want to apply this to a
	motif that isn't a PQ repeat.
	Output: To outname, writes a table with columns locus_tag, start, end.
	Information for each PQ repeat ("PQ" or longer) for each sequence is written to a row.
	If a sequence has no instances of PQ at all, start=0 and end=0 for that locus tag.
	'''
	p = re.compile(pattern) # set up the regex object

	make_outdir(outname) # make outdir if necessary

	with open(outname, "w") as f:
		# open the fasta and get the seqs with locus tags
		with open(multifasta, "r") as fasta_file:
			for record in SeqIO.parse(multifasta, "fasta"):
				# get info from the fasta
				locus_tag = record.id
				seq = str(record.seq)

				# find the pattern in the sequence
				matches, starts, ends = find_matches(p, seq)

				if len(matches) == 0: # no matches found
					# still include the locus tag, but "start" and "end"
					# of the "pattern match" are both 0
					f.write(f"{locus_tag}\t0\t0\n")
				else: # matches were found!
					# get the index of the longest pattern match
					longest = max(range(len(matches)), key=lambda i: len(matches[i]))
						# range(len(matches)) generates a list of indices
						# max() with the key of the length of matches[i] picks index i with the longest matches[i]
					f.write(f"{locus_tag}\t{starts[longest]}\t{ends[longest]}\n") # only write info for the longest match

	print(f"Done writing to {outname}!")


def summarize_pq_repeats_old(multifasta, outname, pattern=r"(?:PQ)+"):
	'''
	Old version of the function that wrote the coordinates of every single PQ repeat to the file,
	rather than the coordinates of only the longest PQ repeat per sequence.
	Inputs: a multifasta file in which to find PQ repeats, an output file name,
	(optional) a regular expression pattern in case you want to apply this to a
	motif that isn't a PQ repeat.
	Output: To outname, writes a table with columns locus_tag, start, end.
	Information for each PQ repeat ("PQ" or longer) for each sequence is written to a row.
	If a sequence has no instances of PQ at all, start=0 and end=0 for that locus tag.
	'''
	'''
	# open the fasta and get the seqs with locus tags
	with open(multifasta, "r") as f:
		for record in SeqIO.parse(multifasta, "fasta"):
			print(record.id, len(record), str(record.seq)[:10])


	make_outdir(outname)
	with open(outname, "w") as f:
		# do all the stuff
		print("hi")
		f.write("hello")
	'''
	p = re.compile(pattern) # set up the regex object

	make_outdir(outname) # make outdir if necessary

	with open(outname, "w") as f:
		# open the fasta and get the seqs with locus tags
		with open(multifasta, "r") as fasta_file:
			for record in SeqIO.parse(multifasta, "fasta"):
				# get info from the fasta
				locus_tag = record.id
				seq = str(record.seq)

				# find the pattern in the sequence
				matches, starts, ends = find_matches(p, seq)

				if len(matches) == 0: # no matches found
					# still include the locus tag, but "start" and "end"
					# of the "pattern match" are both 0
					f.write(f"{locus_tag}\t0\t0\n")
				else: # matches were found!
					for i in range(len(starts)): # starts and ends have same len
						f.write(f"{locus_tag}\t{starts[i]}\t{ends[i]}\n")

	print(f"Done writing to {outname}!")


def make_outdir(outname):
	'''
	General helper function.
	Input: name of output file.
	Output: if the directory in outname doesn't exist, create it.
	'''
	# get outdir path from outname
	directory_path = os.path.dirname(outname)

	# Check if the directory exists, and create it if it doesn't
	if not os.path.exists(directory_path):
		os.makedirs(directory_path)
		print(f"Directory '{directory_path}' created.")
	else:
		print(f"Directory '{directory_path}' already exists.")



if __name__ == "__main__":
	'''
	pattern = r"(?:PQ)+"
	compiled = re.compile(pattern)
	foo=["PQPQPQ", "PQAPQPQ", "PPPQQQ", "PPPAQQQ", "PQPQAPQ", "PQP", "pqpq", "PPPPQQQQ"]

	for s in foo:
    		print(find_matches(compiled, s))
	'''

	# multifasta1="../data/orthologs/Fha1_orthologs_synteny_filtered.fasta"
	# outname1="../data/motifs/Fha1_pq_repeats.tsv"
	# summarize_pq_repeats(multifasta1, outname1)


	# multifasta2="../data/orthologs/Fha1_orthologs.fasta"
	# outname2="../data/motifs/Fha1_and_Fha2_pq_repeats.tsv"
	# summarize_pq_repeats(multifasta2, outname2)

	# multifasta3="../data/orthologs/Fha1_orthologs_Pa_only.fasta"
	# outname3="../data/motifs/Fha1_and_Fha2_Pa_only_pq_repeats.tsv"
	# summarize_pq_repeats(multifasta3, outname3)

	# multifasta4="../data/orthologs/Fha1_orthologs_alt.fasta"
	# outname4="../data/motifs/Fha1_and_Fha2_alt_pq_repeats.tsv"
	# summarize_pq_repeats(multifasta4, outname4)

	# multifasta5="../data/orthologs/Fha1_orthologs_alt_Pa_only.fasta"
	# outname5="../data/motifs/Fha1_and_Fha2_alt_Pa_only_pq_repeats.tsv"
	# summarize_pq_repeats(multifasta5, outname5)

	multifasta6="../data/orthologs/fha1_genome_info_paOnly_top.fasta"
	outname6="../data/motifs/Fha1_Pa_top_pq_repeats.tsv"
	summarize_pq_repeats(multifasta6, outname6)