# Utility functions to filter the output of blast2gen.py.
import pandas as pd

def filter_by_aligned_proportion(metadata_file, outname, threshold=0.9):
  '''
  Filters the input metadata file down to rows in which the length of the sequence_old column
  (representing the BLAST alignment of query against subject), not counting gap characters,
  is at least the specified proportion of the length of the sequence column (representing the
  sequence obtained from the GenBank file during NCBI esearch).
  '''
  df = pd.read_csv(metadata_file, sep='\t')
  print(f"Length of input file: {len(df)} rows")
  
  df["alignment_length"] = [len(s.replace("-", "")) for s in df["sequence_old"]]
  df["subject_length"] = [len(s.replace("-", "")) for s in df["sequence"]]
  
  df = df[df["alignment_length"] / df["subject_length"] >= threshold]
  
  df = df.drop(["alignment_length", "subject_length"], axis=1)
  
  print(f"Length of output file: {len(df)} rows")
  
  df.to_csv(outname, index=False, sep='\t')
  
  
def filter_by_length(metadata_file, outname, lower=0, upper=float('inf')):
  '''
  Filters the input metadata file down to rows in which the length of the sequence column
  is between the specified upper and lower bounds.
  '''
  df = pd.read_csv(metadata_file, sep='\t')
  print(f"Length of input file: {len(df)} rows")
  
  df["subject_length"] = [len(s.replace("-", "")) for s in df["sequence"]]
  
  df = df[(df["subject_length"] >= lower) & (df["subject_length"] <= upper)]
  
  df = df.drop(["subject_length"], axis=1)
  
  print(f"Length of output file: {len(df)} rows")
  
  df.to_csv(outname, index=False, sep='\t')
  

# could add a function to filter by evalue, since there's a column for that.