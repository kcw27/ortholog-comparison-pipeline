#!/usr/bin/env python3

# Original script by Dr. Arkadiy Garber at ASU
# Modifying this to work with my long-format nr blastp outputs
# I think that just entails changing the column names around in main

# Example run: only run one at a time, otherwise it might have an issue with too many requests being made to NCBI
# python $blastscripts/blast2gen.py $datadir/PA3565_with_orgs_long.blast $datadir/PA3565_with_orgs_long_annotatedByBlast2gen.blast
# python blast2gen.py /home/kcw2/data/blast_outputs/fha1_nr_orgs_long_paOnly.txt /home/kcw2/data/blast_outputs/fha1_genome_info_paOnly.tsv &

# python blast2gen.py /home/kcw2/data/blast_outputs/fha1_nr_orgs_long.txt /home/kcw2/data/blast_outputs/fha1_genome_info.tsv &
# python blast2gen.py /home/kcw2/data/blast_outputs/PA3565_nr_orgs_long.txt /home/kcw2/data/blast_outputs/PA3565_genome_info.tsv &
# python blast2gen.py /home/kcw2/data/testing/out/PA3565_nr_small_orgs_long.txt /home/kcw2/data/testing/out/PA3565_nr_small_genome_info.tsv

import sys
import pandas as pd
from Bio import Entrez, SeqIO # I added the SeqIO import
import time
import os # to get the path

def extract_accession(subject):
    return subject.strip().split()[0]

def wp_to_nucleotide(wp_acc):
    """Get nucleotide accession via IPG database for WP proteins."""
    try:
        handle = Entrez.esearch(db="ipg", term=wp_acc)
        record = Entrez.read(handle)
        handle.close()
        ipg_ids = record['IdList']
        if not ipg_ids:
            return None
        handle = Entrez.esummary(db="ipg", id=ipg_ids[0])
        summary = Entrez.read(handle)
        handle.close()
        nuc_acc = summary['DocumentSummarySet']['DocumentSummary'][0]['NucleotideAccession']
        return nuc_acc
    except:
        return None

def protein_to_nucleotide(protein_acc):
    """Direct protein->nucleotide query."""
    try:
        handle = Entrez.elink(dbfrom="protein", db="nuccore", id=protein_acc, linkname="protein_nuccore")
        records = Entrez.read(handle)
        handle.close()
        links = records[0]['LinkSetDb'][0]['Link']
        return links[0]['Id'] if links else None
    except:
        return None

def nucleotide_to_bioproject_assembly(nuc_id):
    """Nucleotide accession to BioProject and Assembly."""
    bioproject_acc, assembly_acc = 'NA', 'NA'
    try:
        # Fetch nucleotide record
        handle = Entrez.efetch(db='nuccore', id=nuc_id, rettype='gb', retmode='xml')
        records = Entrez.read(handle)
        handle.close()

        rec = records[0]
        # Get BioProject
        for xref in rec.get('GBSeq_xrefs', []):
            if xref['GBXref_dbname'] == 'BioProject':
                bioproject_acc = xref['GBXref_id']
                break
        # Get Assembly
        handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=nuc_id)
        asm_records = Entrez.read(handle)
        handle.close()
        if asm_records[0]['LinkSetDb']:
            assembly_id = asm_records[0]['LinkSetDb'][0]['Link'][0]['Id']
            handle = Entrez.esummary(db="assembly", id=assembly_id)
            asm_summary = Entrez.read(handle)
            handle.close()
            assembly_acc = asm_summary['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
    except:
        pass
    return bioproject_acc, assembly_acc

def fetch_info(accessions, outdir):
    start_time = time.time()
    print(f"Will write any errors to fetch_info.log in the output directory: {outdir}")
    log = open(f"{outdir}/fetch_info.log", "w")
    log.write(f"Current time: {time.ctime(time.time())}\n") # time.ctime() for human-readable time

    info_dict = {}
    total_accs = len(accessions)
    for i, acc in enumerate(accessions):
        print(f"Processing accession {acc}; {i+1} out of {total_accs}...")
        nucleotide_acc, bioproject_acc, assembly_acc, organism = 'NA', 'NA', 'NA', 'NA'
        isolation_source, sequencing_technology, title = 'NA', 'NA', 'NA' # adding three pieces of metadata I'm interested in

        # Fetch organism name first
        try:
            #handle = Entrez.esummary(db='protein', id=acc)
            #record = Entrez.read(handle)
            
            # trying this out with SeqIO.read on Entrez.efetch instead
            # because that's what I wrote my isolation source and sequencing technology code to work with
            handle = Entrez.efetch(db='protein', id=acc, rettype="gb", retmode="text")
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
                nucleotide_acc = wp_to_nucleotide(acc)
            except:
                nucleotide_acc = 'NA'
        else:
            try:
                nuc_id = protein_to_nucleotide(acc)
                if nuc_id:
                    handle = Entrez.esummary(db="nuccore", id=nuc_id)
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
                handle = Entrez.esearch(db="nuccore", term=nucleotide_acc)
                search_record = Entrez.read(handle)
                handle.close()
                nuc_ids = search_record['IdList']
                if nuc_ids:
                    bioproject_acc, assembly_acc = nucleotide_to_bioproject_assembly(nuc_ids[0])
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

        print(f"Processed {acc}: Nuc={nucleotide_acc}, BioProj={bioproject_acc}, Asm={assembly_acc}, Org={organism}, Source={isolation_source}, Seqtech={sequencing_technology}, Title={title}, Assembly method={assembly_method}, Coverage={genome_coverage}, Isolation site={isolation_site}")
        time.sleep(0.4)  # Respect rate limit. If needed, increase it so NCBI doesn't complain about too many requests
        
    end_time = time.time()
    print(f"Time elapsed in seconds: {end_time - start_time}")
    log.write(f"Time elapsed in seconds: {end_time - start_time}")
    log.close()
    return info_dict

def main(input_blast, output_tsv, email = "example@mail.com"):
    Entrez.email = email
    print(f"Entrez email: {Entrez.email}")
    outdir = os.path.dirname(output_tsv)

    cols = ['genome_id_old', 'subject', 'sequence_old', 'evalue', 'title_old', 'organism_old'] # I believe 'subject' is meant to be the protein accession

    blast_df = pd.read_csv(input_blast, sep='\t', header=None, names=cols, engine='python')
    blast_df['clean_accession'] = blast_df['subject'].apply(extract_accession)
    unique_accs = blast_df['clean_accession'].unique().tolist()
    print(len(unique_accs))

    print(f"Fetching detailed info for {len(unique_accs)} proteins from NCBI...")
    info_dict = fetch_info(unique_accs, outdir)

    blast_df['nucleotide_accession'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['nucleotide_acc'])
    blast_df['bioproject_accession'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['bioproject_acc'])
    blast_df['assembly_accession'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['assembly_acc'])
    blast_df['organism'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['organism'])
    # adding metadata I'm interested in
    blast_df['isolation_source'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['isolation_source'])
    blast_df['sequencing_technology'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['sequencing_technology'])
    blast_df['titles'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['title'])
    blast_df['assembly_method'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['assembly_method'])
    blast_df['genome_coverage'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['genome_coverage'])
    blast_df['protein_seq'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['protein_seq'])
    #blast_df['protein_id'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['protein_id'])
    blast_df['locus_tag'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['locus_tag'])
    blast_df['isolation_site'] = blast_df['clean_accession'].map(lambda x: info_dict[x]['isolation_site'])

    # Columns: same set of columns as the original, with new info written at the end
    # Update: swap some old columns with columns obtained from metadata retrieval
    # and place them such that the output can be plugged directly into convert_blast_to_fasta.sh
    final_df = blast_df[['assembly_accession', 'subject', 'protein_seq', 'evalue', 'title_old', 'locus_tag', 'organism',
                         'nucleotide_accession', 'bioproject_accession', 'genome_id_old', 'organism_old',
                         'isolation_source', 'sequencing_technology', 'titles', 'assembly_method',
                         'genome_coverage', 'sequence_old', 'isolation_site']] # removed 'protein_id' because it was redundant
                         
    # rename columns for compatibility with downstream analysis functions in the GUI
    final_df.rename(columns={'assembly_accession': 'genome_id', 'protein_seq': 'sequence'}, inplace=True)

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
