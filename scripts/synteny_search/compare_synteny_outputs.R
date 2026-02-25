# for a quick comparison of the four wspF synteny outputs

library(dplyr)
library(glue)

compare <- function(df1, df2) {
  df1_size <- nrow(df1) 
  df2_size <- nrow(df2) 
  print(glue("Size of df1: {df1_size}"))
  print(glue("Size of df2: {df2_size}"))
  print(glue("Number of rows shared: {nrow(merge(df1, df2))}"))
  #print(merge(df1, df2) |> select(genome_id, contig, locus_tag, protein_id))
#  smaller_df1 <- df1 |> select(genome_id, contig, locus_tag, protein_id)
#  smaller_df2 <- df2 |> select(genome_id, contig, locus_tag, protein_id)
#  smaller_df1 <- df1 |> select(genome_id, protein_id)
#  smaller_df2 <- df2 |> select(genome_id, protein_id)
#  print(glue("Number of rows shared if only using genome_id and protein_id columns: {nrow(merge(smaller_df1, smaller_df2))}"))
  
  print("Locus tags shared (with genome ID so you can look them up):")
  print(merge(df1, df2) |> select(genome_id, locus_tag))
  
#  print(glue("Number of unique genome_id in df1 and df2: {length(unique(df1$genome_id))}, {length(unique(df2$genome_id))}"))
#  print(glue("Intersection of genome_id: {length(intersect(df1$genome_id, df2$genome_id))}"))
#  
#  print(glue("Number of unique contig in df1 and df2: {length(unique(df1$contig))}, {length(unique(df2$contig))}"))
#  print(glue("Intersection of contig: {length(intersect(df1$contig, df2$contig))}"))
#  
#  print(glue("Number of unique locus_tag in df1 and df2: {length(unique(df1$locus_tag))}, {length(unique(df2$locus_tag))}"))
#  print(glue("Intersection of locus_tag: {length(intersect(df1$locus_tag, df2$locus_tag))}"))
#  
#  print(glue("Number of unique protein_id in df1 and df2: {length(unique(df1$protein_id))}, {length(unique(df2$protein_id))}"))
#  print(glue("Intersection of protein_id: {length(intersect(df1$protein_id, df2$protein_id))}"))
  
  print("***")
}

nb <- read.csv("/home/kcw2/data/blast_outputs/wspF/synteny_search_outputs/negBackward/synteny_summary.tsv", header=TRUE, sep="\t")
nf <- read.csv("/home/kcw2/data/blast_outputs/wspF/synteny_search_outputs/negForward/synteny_summary.tsv", header=TRUE, sep="\t")
pb <- read.csv("/home/kcw2/data/blast_outputs/wspF/synteny_search_outputs/posBackward/synteny_summary.tsv", header=TRUE, sep="\t")
pf <- read.csv("/home/kcw2/data/blast_outputs/wspF/synteny_search_outputs/posForward/synteny_summary.tsv", header=TRUE, sep="\t")

print("nb, nf")
compare(nb, nf)
print("nb, pb")
compare(nb, pb)
print("nb, pf")
compare(nb, pf)
print("nf, pb")
compare(nf, pb)
print("nf, pf")
compare(nf, pf)
print("pb, pf")
compare(pb, pf)