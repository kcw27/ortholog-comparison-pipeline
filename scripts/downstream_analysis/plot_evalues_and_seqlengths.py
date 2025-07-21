from Bio import SeqIO
# import os
import matplotlib.pyplot as plt
import re
import pandas as pd
#import seaborn as sns

def plot_evals_and_seqlens(multifasta, outname):
    '''
    Input: multifasta file, name to save the output plot to
    Assumption: fasta headers are in the form >locus_tag title|evalue (which is the output format of convert_blast_to_fasta.sh)
    Output: scatter plot with sequence length as x axis, evalue as y axis (log-scaled)
    Points are colored based on whether they correspond to a P. aeruginosa record.
    '''
    with open(multifasta, "r") as f:
        # can't assign the SeqIO parse output to a variable and use it multiple times; must create a new one each time
        seq_lens = [int(len(record.seq)) for record in SeqIO.parse(multifasta, "fasta")] # sequence lengths as ints
        evalues = [float(record.description.split("|")[1]) for record in SeqIO.parse(multifasta, "fasta")] # parse out the evalues and convert from str to float

        pattern = r"\[(.*?)\]" # when ? comes after a quantifier like *, it makes it non-greedy
        organisms = [re.search(pattern, record.description).group() for record in SeqIO.parse(multifasta, "fasta")]

        #pattern2 = r"Pseudomonas aeruginosa"
        pattern2 = r"Pseudomonas aeruginosa." # add the . so it allows for matches containing strain names; doesn't seem to make a difference though
        is_pa = [re.search(pattern2, org) != None for org in organisms]
        pa_color_map = ['red' if species_is_pa else 'blue' for species_is_pa in is_pa]

        # print(seq_lens[len(seq_lens)-10:])
        # print(evalues[len(seq_lens)-10:])
        # print(max(evalues))
        # print(len(seq_lens), len(evalues), len(organisms), len(is_pa))
    
    df = pd.DataFrame({'seq_len': seq_lens, 'evalue': evalues, 'organism': organisms, 'is_pa': is_pa})
    #print(df)
    # print(sum(df.is_pa))
    # print(pa_color_map[:10])



    # df.plot.scatter('seq_len', # x-axis
    #             'evalue', # y-axis
    #             hue='is_pa'
    #             #c=pa_color_map # color by whether it's P. aeruginosa
    #            )

    #sns.relplot(data=df, x='seq_len', y='evalue', hue='is_pa')

    # groups = df.groupby('is_pa')
    # for name, group in groups:
    #     plt.plot(group.seq_len, group.evalue, marker='o', linestyle='', markersize=12, label=name)

    
    # for name, group in df.groupby('is_pa'):
    #     plt.scatter(
    #         group['seq_len'],
    #         group['evalue'],
    #         color=colors[name],
    #         label=str(name),
    #         s=10,  # marker size
    #         alpha=0.7
    #     )

    # plt.legend()

    # # plt.xlabel('Sequence Length')
    # # plt.ylabel('E-value')
    # plt.xlabel('Sequence length')
    # plt.ylabel('evalue (log10)')
    # plt.title('BLAST hit length vs evalue (log10)')
    # plt.yscale('log')  # Switch to log scale
    # plt.legend(title='Is P. aeruginosa')
    # plt.tight_layout()
    # #plt.figure(figsize=(10, 6))
    # plt.show()

    # # plt.figure(figsize=(10, 6))
    # # # plt.scatter(seq_lens, evalues, marker='o', color='b')
    # # plt.yscale('log')  # Set y-axis to log scale
    # # plt.xlabel('Sequence length')
    # # plt.ylabel('evalue (log10)')
    # # plt.title('BLAST hit length vs evalue (log10)')
    # # #plt.legend()
    # # #plt.grid(True)
    # #plt.show()
    # plt.savefig(outname) 

    fig, ax = plt.subplots(figsize=(10, 6))  # Create a new figure and set size

    colors = {True: 'red', False: 'blue'} # color map for is_pa

    for name, group in df.groupby('is_pa'):
        ax.scatter(
            group['seq_len'],
            group['evalue'],
            color=colors[name],
            label=str(name),
            s=10,
            alpha=0.7
        )

    ax.set_xlabel('Sequence length')
    ax.set_ylabel('evalue (log10)')
    ax.set_title(f'BLAST hit length vs evalue (log10)\n{sum(df.is_pa)} Pseudomonas aeruginosa records out of {len(df)} total records')
    ax.set_yscale('log')
    ax.legend(title='Is P. aeruginosa')
    fig.tight_layout()

    fig.savefig(outname)
    plt.close(fig)  # Ensures clean slate before next iteration, otherwise it'll save the same figure to multiple runs of this function

def plot_evals_and_seqlens_v2(multifasta, subset_multifasta, outname):
    '''
    Input: multifasta file, multifasta file containing a subset of the first, name to save the output plot to
    Assumption: fasta headers are in the form >locus_tag title|evalue (which is the output format of convert_blast_to_fasta.sh)
    Output: scatter plot based on multifasta with sequence length as x axis, evalue as y axis (log-scaled).
    Points are colored according to whether they're in subset_multifasta.
    '''
    seq_lens = []
    evalues = []
    organisms = []
    in_subset = []

    pattern = r"\[(.*?)\]" # when ? comes after a quantifier like *, it makes it non-greedy

    subset_ids = {record.id for record in SeqIO.parse(subset_multifasta, "fasta")}

    with open(multifasta, "r") as f:
        for record in SeqIO.parse(multifasta, "fasta"):
            seq_lens.append(int(len(record.seq))) # sequence lengths as ints
            evalues.append(float(record.description.split("|")[1])) # parse out the evalues and convert from str to float
            organisms.append(re.search(pattern, record.description).group())
            in_subset.append(record.id in subset_ids)

    df = pd.DataFrame({'seq_len': seq_lens, 'evalue': evalues, 'organism': organisms, 'in_subset': in_subset})
    # print(subset_ids)
    # print(df.head())

    # create the plot
    fig, ax = plt.subplots(figsize=(10, 6))  # Create a new figure and set size

    colors = {True: 'red', False: 'blue'} # color map for in_subset

    for name, group in df.groupby('in_subset'):
        ax.scatter(
            group['seq_len'],
            group['evalue'],
            color=colors[name],
            label=str(name),
            s=10,
            alpha=0.7
        )

    ax.set_xlabel('Sequence length')
    ax.set_ylabel('evalue (log10)')
    ax.set_title(f'BLAST hit length vs evalue (log10)\n{sum(df.in_subset)} records out of {len(df)} total records appeared in the fha1 synteny search')
    ax.set_yscale('log')
    ax.legend(title='Appeared in synteny search for fha1')
    fig.tight_layout()

    fig.savefig(outname)
    plt.close(fig)  # Ensures clean slate before next iteration, otherwise it'll save the same figure to multiple runs of this function

def subset_fasta_to_pa(multifasta, outname):
    '''
    Input: multifasta file, name to save the output plot to
    Assumption: fasta headers are in the form >locus_tag title|evalue (which is the output format of convert_blast_to_fasta.sh)
    Output: fasta file containing only records for which the species is P. aeruginosa 
    '''

    pattern = r"\[(.*?)\]" # when ? comes after a quantifier like *, it makes it non-greedy
    pattern2 = r"Pseudomonas aeruginosa."

    with open(outname, "w") as outfile:
        with open(multifasta, "r") as f:
            for record in SeqIO.parse(multifasta, "fasta"):
                organism = re.search(pattern, record.description).group()
                if re.search(pattern2, organism): # if the organism is P. aeruginosa
                    outfile.write(f">{record.description}\n")
                    outfile.write(f"{str(record.seq)}\n")





if __name__ == "__main__":
    multifasta = "../data/orthologs/Fha1_orthologs.fasta"
    outname = "../Fha1_outputs/Fha1_seqlens_evalues_plot.png"
    #plot_evals_and_seqlens(multifasta, outname)

    multifasta2 = "../data/orthologs/Fha1_orthologs_synteny_filtered.fasta"
    outname2 = "../Fha1_outputs/Fha1_synteny_filtered_seqlens_evalues_plot.png"
    #plot_evals_and_seqlens(multifasta2, outname2)

    outname3 = "../Fha1_outputs/Fha1_seqlens_evalues_plot_colored_by_synteny.png"
    #plot_evals_and_seqlens_v2(multifasta, multifasta2, outname3)

    #subset_fasta_to_pa(multifasta, "../data/orthologs/Fha1_orthologs_Pa_only.fasta")

    # there's no multifasta3; I'm just keeping the number consistent with outname's number
    multifasta4 = "../data/orthologs/Fha1_orthologs_alt.fasta"
    outname4 = "../Fha1_outputs/Fha1_seqlens_evalues_plot_new.png"
    #plot_evals_and_seqlens(multifasta4, outname4)

    #subset_fasta_to_pa(multifasta4, "../data/orthologs/Fha1_orthologs_alt_Pa_only.fasta")

    multifasta5 = "../data/orthologs/fha1_genome_info_paOnly_top.fasta"
    outname5 = "../Fha1_outputs_v3/Fha1_seqlens_evalues_plot_paOnly_top.png"
    plot_evals_and_seqlens(multifasta5, outname5)

    # i=0
    # with open(multifasta, "r") as f:
    #         # for record in SeqIO.parse(multifasta, "fasta"):
    #         #         while i<10:
    #         #             #print(record.id, len(record), str(record.seq)[:10])
    #         #             print(record.description)
    #         #             i+=1
    #         # print(list(SeqIO.parse(multifasta, "fasta"))[:10])
    #     desc = [record.description for record in SeqIO.parse(multifasta, "fasta")]
    #     print(len(desc), desc[:10])
