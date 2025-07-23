#!/bin/bash

# CLIs:
# $1: long-format BLAST file with organisms in column 6
# $2: name of output file, to which the list of unique organisms (genus + species combinations; cut off the strain) is written
# Strain is excluded from the organism name because this script is intended to produce an input for convert_organisms_to_accessions.sh,
# and querying each organism with strain name attached seems to only result in redundant results and an increase in runtime.

# Example run:
# bash "${blastscripts}/extract_organisms.sh" "${datadir}/PA3565_with_orgs_long.blast" "${datadir}/PA3565_orgs.txt"

cut -f 6 $1 | cut -d " " -f 1,2 | sort | uniq > $2