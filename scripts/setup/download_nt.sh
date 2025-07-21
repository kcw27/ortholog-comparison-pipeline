# to download the nt nucleotide database for BLAST searches

source activate blast_env # so that it can find update_blastdb.pl; it's part of https://anaconda.org/bioconda/blast
cd /home/share # download to the shared directory
# ideally, you'd create /home/share/nt and then cd to it; I didn't realize that it wouldn't make a subdirectory for the files
update_blastdb.pl --decompress nt # make sure the conda blast_env environment is active so you have access to update_blastdb.pl