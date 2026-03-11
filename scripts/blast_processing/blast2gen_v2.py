#!/usr/bin/env python3

# Implement caching and retries

import sys
import pandas as pd
from Bio import Entrez, SeqIO # I added the SeqIO import
import time
import os # to get the path
import math
import numpy as np

def entrez_retry(func, log, max_retries=3, delay=2, **kwargs):
    for attempt in range(1, max_retries+1):
        try:
            return func(**kwargs)
        except Exception as e:
            msg = f"[Retry {attempt}/{max_retries}] Error: {e}\n"
            print(msg.strip())
            log.write(msg)
            if attempt == max_retries:
                raise
            time.sleep(delay)

def extract_accession(subject):
    return subject.strip().split()[0]

def wp_to_nucleotide(wp_acc, log):
    """Get nucleotide accession via IPG database for WP proteins."""
    try:
        #handle = Entrez.esearch(db="ipg", term=wp_acc)
        handle = entrez_retry(Entrez.esearch, log, db="ipg", term=wp_acc)
        record = Entrez.read(handle)
        handle.close()
        ipg_ids = record['IdList']
        if not ipg_ids:
            return None
        #handle = Entrez.esummary(db="ipg", id=ipg_ids[0])
        handle = entrez_retry(Entrez.esummary, log, db="ipg", id=ipg_ids[0])
        summary = Entrez.read(handle)
        handle.close()
        nuc_acc = summary['DocumentSummarySet']['DocumentSummary'][0]['NucleotideAccession']
        return nuc_acc
    except:
        return None

def protein_to_nucleotide(protein_acc, log):
    """Direct protein->nucleotide query."""
    try:
        #handle = Entrez.elink(dbfrom="protein", db="nuccore", id=protein_acc, linkname="protein_nuccore")
        handle = entrez_retry(Entrez.elink, log, dbfrom="protein", db="nuccore", id=protein_acc, linkname="protein_nuccore")
        records = Entrez.read(handle)
        handle.close()
        links = records[0]['LinkSetDb'][0]['Link']
        return links[0]['Id'] if links else None
    except:
        return None

def nucleotide_to_bioproject_assembly(nuc_id, log):
    """Nucleotide accession to BioProject and Assembly."""
    bioproject_acc, assembly_acc = 'NA', 'NA'
    try:
        # Fetch nucleotide record
        #handle = Entrez.efetch(db='nuccore', id=nuc_id, rettype='gb', retmode='xml')
        handle = entrez_retry(Entrez.efetch, log, db='nuccore', id=nuc_id, rettype='gb', retmode='xml')
        records = Entrez.read(handle)
        handle.close()

        rec = records[0]
        # Get BioProject
        for xref in rec.get('GBSeq_xrefs', []):
            if xref['GBXref_dbname'] == 'BioProject':
                bioproject_acc = xref['GBXref_id']
                break
        # Get Assembly
        #handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=nuc_id)
        handle = entrez_retry(Entrez.elink, log, dbfrom="nuccore", db="assembly", id=nuc_id)
        asm_records = Entrez.read(handle)
        handle.close()
        if asm_records[0]['LinkSetDb']:
            assembly_id = asm_records[0]['LinkSetDb'][0]['Link'][0]['Id']
            #handle = Entrez.esummary(db="assembly", id=assembly_id)
            handle = entrez_retry(Entrez.esummary, log, db="assembly", id=assembly_id)
            asm_summary = Entrez.read(handle)
            handle.close()
            assembly_acc = asm_summary['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
    except:
        pass
    return bioproject_acc, assembly_acc

def load_cache(cache_file):
    if not os.path.exists(cache_file):
        return {}

    try:
        df = pd.read_csv(cache_file, sep="\t")
    except Exception as e:
        print(f"Cache file unreadable: {e}")
        print("Deleting last row and retrying.")
        return _truncate_last_row(cache_file)

    cache = {}
    for idx, row in df.iterrows():
        try:
            acc = row["accession"]
        except KeyError:
            print("Cache row missing 'accession'. Truncating last row.")
            return _truncate_last_row(cache_file)
        cache[acc] = row.to_dict()

    return cache

def _truncate_last_row(cache_file):
    df = pd.read_csv(cache_file, sep="\t", on_bad_lines="skip")

    # If the file has only one row, delete it entirely
    if len(df) <= 1:
        os.remove(cache_file)
        return {}

    # Drop the last row
    df = df.iloc[:-1]

    # Rewrite the cache cleanly
    df.to_csv(cache_file, sep="\t", index=False)

    # Return the cleaned cache as a dict
    cache = {}
    for _, row in df.iterrows():
        cache[row["accession"]] = row.to_dict()
    return cache


def append_cache_row(cache_file, row_dict):
    df = pd.DataFrame([row_dict])
    header = not os.path.exists(cache_file)
    df.to_csv(cache_file, sep="\t", mode="a", header=header, index=False)

#def cache_row_incomplete(row_dict):
#    """
#    Returns True if the cached row is missing real metadata.
#    We treat 'NA' or empty strings as incomplete.
#    """
##    keys_to_check = [
##        "nucleotide_acc", "bioproject_acc", "assembly_acc",
##        "organism", "protein_seq"
##    ]
##    for k in keys_to_check:
##        if k not in row_dict or row_dict[k] in ("NA", "", None):
##            return True
#    print("Running cache_row_incomplete")
#    print(row_dict)
#    for value in row_dict.values(): # just check all of them
#        if value is pd.NA:
#            print("Value was pd.NA")
#            return True
#        if isinstance(value, float) and math.isnan(value):
#            print("Value was a mathematical nan")
#            return True
#        if value is np.nan:
#            print("Value was np.nan")
#            return True
#        if value in ("NA", "", None):
#            print("Value was a stringlike nan")
#            return True
#    return False

def cache_row_incomplete(row_dict):
    """
    Returns True if ANY value in the cached row is missing or incomplete.
    Treats pd.NA, numpy.nan, float nan, None, and 'NA' as incomplete.
    """
    for value in row_dict.values():

        # Pandas NA
        if value is pd.NA:
            return True

        # numpy.nan or float nan
        if isinstance(value, float) and math.isnan(value):
            return True

        # numpy.nan (identity check doesn't work, so use isnan)
        if isinstance(value, np.generic) and np.isnan(value):
            return True

        # String-like missing values
        if value in ("NA", "", None):
            return True

    return False

def normalize_row(row_dict):
    """
    Replace any NaN-like value with the string 'NA'.
    """
    clean = {}
    for k, v in row_dict.items():

        # Pandas NA
        if v is pd.NA:
            clean[k] = "NA"
            continue

        # numpy.nan or float nan
        try:
            if isinstance(v, float) and math.isnan(v):
                clean[k] = "NA"
                continue
        except:
            pass

        # numpy.nan
        try:
            if isinstance(v, np.generic) and np.isnan(v):
                clean[k] = "NA"
                continue
        except:
            pass

        # None or empty string
        if v in ("", None):
            clean[k] = "NA"
            continue

        clean[k] = v

    return clean


def fetch_info(accessions, outdir):
    start_time = time.time()
    print(f"Will write any errors to fetch_info.log in the output directory: {outdir}")
    log = open(f"{outdir}/fetch_info.log", "w")
    log.write(f"Current time: {time.ctime(time.time())}\n") # time.ctime() for human-readable time
    
    cache_file = os.path.join(outdir, "cache.tsv")
    cache = load_cache(cache_file)

    # If cache exists, check if the last row is incomplete
    if cache:
        # Get the last accession written
        last_acc = list(cache.keys())[-1]
        last_row = cache[last_acc]
    
        if cache_row_incomplete(last_row):
            print(f"Last cached row for {last_acc} is potentially incomplete. Re-fetching...")
            log.write(f"Last cached row for {last_acc} potentially incomplete; retrying.\n")
    
            # Remove the bad row from the cache dict
            del cache[last_acc]
    
            # Rewrite cache.tsv without the bad row
            df = pd.DataFrame(cache.values())
            df.to_csv(cache_file, sep="\t", index=False)
    
            # Now force re-fetching this accession
            accessions.insert(0, last_acc)


    info_dict = {}
    total_accs = len(accessions)
    for i, acc in enumerate(accessions):
        print(f"Processing accession {acc}; {i+1} out of {total_accs}...")
        
        if acc in cache:
          print(f"Using cached result for {acc}")
          info_dict[acc] = cache[acc]
          continue
          
        nucleotide_acc, bioproject_acc, assembly_acc, organism = 'NA', 'NA', 'NA', 'NA'
        isolation_source, sequencing_technology, title = 'NA', 'NA', 'NA' # adding three pieces of metadata I'm interested in

        # Fetch organism name first
        try:
            #handle = Entrez.esummary(db='protein', id=acc)
            #record = Entrez.read(handle)
            
            # trying this out with SeqIO.read on Entrez.efetch instead
            # because that's what I wrote my isolation source and sequencing technology code to work with
            #handle = Entrez.efetch(db='protein', id=acc, rettype="gb", retmode="text")
            handle = entrez_retry(Entrez.efetch, log, db='protein', id=acc, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank") 
            handle.close()
            #organism = record[0].get('Title', 'NA').split('[')[-1].rstrip(']')
            info = record.features[0].qualifiers
            organism = info["organism"][0] # it's a list containing a single string, so access that string
            protein_seq = str(record.seq)

            # protein_id and locus_tag
            #protein_id = "NA"
            locus_tag = "NA"
            
            for feature in record.features:
                if feature.type == "CDS": # CDS is for gene-level features
                    q = feature.qualifiers
                    #protein_id = q.get("protein_id", ["NA"])[0]
                    #protein_id = record.id
                    # or
                    #protein_id = record.annotations.get("accessions", ["NA"])[0]

                    locus_tag = q.get("locus_tag", ["NA"])[0]
                    break

            
            # Fetch isolation source and sequencing technology too, while we have this record
            # (copying and pasting some code from my own function, fetch_metadata() in metadata_processing.py
            # Isolation source(s):

            iso_sources = info.get("isolation_source", "NA") # it's a string if not found
            if isinstance(iso_sources, list): # i.e. isolation source was found
                iso_sources = ", ".join(iso_sources)
            isolation_source = iso_sources
            
            # Sequencing technology:
            anno = record.annotations
            # looking for information that is two levels deep in a nested dict if it exists,
            # hence the two levels of the get method
            # can't assume that 'structured_comment' always exists as a key in anno
            structured_comment = anno.get('structured_comment', {}) # default to empty dict to ensure that get() works on structured_comment
            genome_assembly = structured_comment.get('Genome-Assembly-Data', {}) # default to empty dict to ensure that get() works on genome_assembly
            sequencing_technology = genome_assembly.get('Sequencing Technology', 'NA')
            assembly_method = genome_assembly.get('Assembly Method', 'NA')
            genome_coverage = genome_assembly.get('Genome Coverage', 'NA')
            
            migs = structured_comment.get('MIGS-Data', {}) # once again default to empty dict
            isolation_site = migs.get('Isolation Site', 'NA')
            
            # Title:
            title = anno.get('references', 'NA')
        except Exception as e:
            #organism = 'NA'
            print(f"An exception occurred: {e}")
            log.write(f"An exception occurred for {acc}: {e}\n")

        # WP proteins need IPG lookup
        if acc.startswith('WP_'):
            try:
                nucleotide_acc = wp_to_nucleotide(acc, log)
            except:
                nucleotide_acc = 'NA'
        else:
            try:
                nuc_id = protein_to_nucleotide(acc, log)
                if nuc_id:
                    #handle = Entrez.esummary(db="nuccore", id=nuc_id)
                    handle = entrez_retry(Entrez.esummary, log, db="nuccore", id=nuc_id)
                    nuc_record = Entrez.read(handle)
                    handle.close()

                    nucleotide_acc = nuc_record[0]['AccessionVersion']
            except Exception as e:
                print(f"An exception occurred while retrieving nucleotide_id: {e}")
                log.write(f"An exception occurred while retrieving nucleotide_id for {acc}: {e}\n")
                nucleotide_acc = 'NA'

        # Now fetch BioProject and Assembly from nucleotide accession
        if nucleotide_acc:
            try: # adding a try-except block because I got a runtime error here
                #handle = Entrez.esearch(db="nuccore", term=nucleotide_acc)
                handle = entrez_retry(Entrez.esearch, log, db="nuccore", term=nucleotide_acc)
                search_record = Entrez.read(handle)
                handle.close()
                nuc_ids = search_record['IdList']
                if nuc_ids:
                    bioproject_acc, assembly_acc = nucleotide_to_bioproject_assembly(nuc_ids[0], log)
            except RuntimeError as e:
                print(f"Error with accession {nucleotide_acc}: {e}")
                log.write(f"Error with accession {nucleotide_acc} for {acc}: {e}\n")


        info_dict[acc] = {
            'nucleotide_acc': nucleotide_acc if nucleotide_acc else 'NA',
            'bioproject_acc': bioproject_acc,
            'assembly_acc': assembly_acc,
            'organism': organism,
            'isolation_source': isolation_source,
            'sequencing_technology': sequencing_technology,
            'title': title,
            'assembly_method': assembly_method,
            'genome_coverage': genome_coverage,
            'protein_seq': protein_seq,
            #'protein_id': protein_id,
            'locus_tag': locus_tag,
            'isolation_site': isolation_site
        }
        
        row = {"accession": acc}
        row.update(info_dict[acc])
        
        #row = normalize_row(row) # so that no NaNs end up in the cache
        
        append_cache_row(cache_file, row)


        print(f"Processed {acc}: Nuc={nucleotide_acc}, BioProj={bioproject_acc}, Asm={assembly_acc}, Org={organism}, Source={isolation_source}, Seqtech={sequencing_technology}, Title={title}, Assembly method={assembly_method}, Coverage={genome_coverage}, Isolation site={isolation_site}")
        time.sleep(0.4)  # Respect rate limit. If needed, increase it so NCBI doesn't complain about too many requests
        
    end_time = time.time()
    print(f"Time elapsed in seconds: {end_time - start_time}")
    log.write(f"Time elapsed in seconds: {end_time - start_time}")
    log.close()
    return info_dict

def main(input_blast, output_tsv, email="example@mail.com"):
    Entrez.email = email
    print(f"Entrez email: {Entrez.email}")
    outdir = os.path.dirname(output_tsv)

    cols = ['genome_id_old', 'subject', 'sequence_old', 'evalue', 'title_old', 'organism_old']
    blast_df = pd.read_csv(input_blast, sep='\t', header=None, names=cols, engine='python')
    blast_df['clean_accession'] = blast_df['subject'].apply(extract_accession)
    unique_accs = blast_df['clean_accession'].unique().tolist()

    print(f"Fetching detailed info for {len(unique_accs)} proteins from NCBI...")
    info_dict = fetch_info(unique_accs, outdir)

    # Map cached/fetched metadata back into the dataframe
    for col in [
        'nucleotide_acc', 'bioproject_acc', 'assembly_acc', 'organism',
        'isolation_source', 'sequencing_technology', 'title',
        'assembly_method', 'genome_coverage', 'protein_seq',
        'locus_tag', 'isolation_site'
    ]:
        blast_df[col] = blast_df['clean_accession'].map(lambda x: info_dict[x][col])

    final_df = blast_df[['assembly_acc', 'subject', 'protein_seq', 'evalue', 'title_old', 'locus_tag', 'organism',
                         'nucleotide_acc', 'bioproject_acc', 'genome_id_old', 'organism_old',
                         'isolation_source', 'sequencing_technology', 'title', 'assembly_method',
                         'genome_coverage', 'sequence_old', 'isolation_site']]

    final_df.rename(columns={'assembly_acc': 'genome_id', 'protein_seq': 'sequence'}, inplace=True)
    final_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Annotated output saved to {output_tsv}")

if __name__ == '__main__':
    if len(sys.argv) not in (3, 4):
        print(f"Usage: {sys.argv[0]} input_blast.tsv output_annotated.tsv [user_email@mail.com]")
        sys.exit(1)

    input_blast = sys.argv[1]
    output_tsv  = sys.argv[2]
    email       = sys.argv[3] if len(sys.argv) == 4 else "example@mail.com"

    main(input_blast, output_tsv, email)
