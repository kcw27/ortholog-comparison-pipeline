#!/bin/bash

# --- Parse command-line flags ---
print_help() {
    echo "Usage: $0 -f <metadata_file> -o <out> -c <category_file> -s <subcategory_file>"
    echo ""
    echo "  -f    Path to metadata file"
    echo "  -o    Output filename (optional, default "" to overwrite the input file)"
    echo "  -c    Path to category file"
    echo "  -s    Path to subcategory file"
    echo "  -h    Show help message and exit"
}

# Set default output filename
out=""

while getopts ":f:o:c:s:h" opt; do
  case $opt in
    f) metadata_file="$OPTARG" ;;
    o) out="$OPTARG" ;;
    c) category_file="$OPTARG" ;;
    s) subcategory_file="$OPTARG" ;;
    h) print_help; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; print_help; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; print_help; exit 1 ;;
  esac
done

if [[ -z "$metadata_file" || -z "$category_file" || -z "$subcategory_file" ]]; then
    echo "Error: Flags -f, -c, and -s are required." >&2
    print_help
    exit 1
fi

# --- Run embedded Python ---
wrapperdir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # get location of current script
scriptsdir="${wrapperdir}/downstream_analysis"

python3 - <<EOF
import importlib.util
import sys

script_path = "$scriptsdir/metadata_processing.py"
module_name = "mp"

spec = importlib.util.spec_from_file_location(module_name, script_path)
mp = importlib.util.module_from_spec(spec)
sys.modules[module_name] = mp
spec.loader.exec_module(mp)

# guess the metadata origin
origin = mp.determine_origin("$metadata_file")
print(f"Inferred metadata origin: {origin}")
if origin not in ['synteny_summary', 'fetched']:
  raise Exception(f"Invalid origin {origin}")

mp.rescue_source("$metadata_file", origin, outname="$out")
mp.categorize("$metadata_file", "$category_file", "$subcategory_file", outname="$out")
EOF


# --- Run R script for benchmarking ---
# Rscript "$scriptsdir/R_scripts/benchmarking.R" "$metadata_file"
Rscript -e "source('$scriptsdir/R_scripts/benchmarking.R'); benchmark('$out')"