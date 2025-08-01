#!/bin/bash

# --- Parse command-line flags ---
print_help() {
    echo "Usage: $0 -f <metadata_file> -o <metadata_origin> -c <category_file> -s <subcategory_file>"
    echo ""
    echo "  -f    Path to metadata file"
    echo "  -o    Origin type (must be 'synteny_summary' or 'fetched')"
    echo "  -c    Path to category file"
    echo "  -s    Path to subcategory file"
    echo "  -h    Show help message and exit"
}

while getopts ":f:o:c:s:h" opt; do
  case $opt in
    f) metadata_file="$OPTARG" ;;
    o) metadata_origin="$OPTARG" ;;
    c) category_file="$OPTARG" ;;
    s) subcategory_file="$OPTARG" ;;
    h) print_help; exit 0 ;;
    \?) echo "Invalid option: -$OPTARG" >&2; print_help; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; print_help; exit 1 ;;
  esac
done

# --- Validate required arguments ---
if [[ -z "$metadata_file" || -z "$metadata_origin" || -z "$category_file" || -z "$subcategory_file" ]]; then
    echo "Error: All four flags -f, -o, -c, and -s are required." >&2
    print_help
    exit 1
fi

# --- Validate metadata_origin ---
if [[ "$metadata_origin" != "synteny_summary" && "$metadata_origin" != "fetched" ]]; then
    echo "Error: -o must be either 'synteny_summary' or 'fetched'" >&2
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

mp.rescue_source("$metadata_file", "$metadata_origin")
mp.categorize("$metadata_file", "$category_file", "$subcategory_file")
EOF
