import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# ======= Functions to Process Input Data ======= #
def readAAsequence(filename):
    aa_seq = ""
    with open(filename) as file:
        next(file)
        for line in file:
            aa_seq += line.replace('\n', '')

    return aa_seq


def AAtypes(aa_seq):
    polar = ['Y', 'W', 'H', 'K', 'R', 'E', 'D', 'N', 'Q', 'S', 'T', 'C']
    small = ['P', 'A', 'G', 'S', 'N', 'D', 'T', 'C', 'V']
    hydro = ['I', 'L', 'V', 'M', 'F', 'Y', 'W', 'H', 'K', 'T', 'C', 'A']

    p_count = 0
    s_count = 0
    h_count = 0

    for aa in aa_seq:
        if aa in polar:
            p_count += 1
        if aa in small:
            s_count += 1
        if aa in hydro:
            h_count += 1

    len_seq = len(aa_seq)
    aa_usage = [p_count/len_seq, s_count/len_seq, h_count/len_seq]

    return aa_usage


def plot_aa_usage(df, file_save_loc, file_ext):
    save_name = file_save_loc.replace(file_ext, '.jpg')
    fig = plt.figure()
    sns.barplot(x="type", y="aa_usage", data=df, palette="Set2")
    plt.xlabel("\nAmino Acid Classification", fontweight='bold')
    plt.ylabel("Amino Acid Usage\n", fontweight='bold')
    fig.savefig(save_name, dpi=300, bbox_inches='tight')


def AAtypetable(file_save_loc, file_ext):
    # Create empty data frame to populate:
    df = pd.DataFrame(columns=['type', 'aa_usage'], dtype=object)
    df['type'] = ['Polar', 'Small', 'Hydro']

    aa_seq = readAAsequence(file_save_loc)
    aa_usage = AAtypes(aa_seq)

    # adding information to created data frame (rounding the amino acid usage to 2 d.p's):
    df['aa_usage'] = aa_usage

    # Writing to output file with tab-separated values:
    #save_name = file_save_loc.replace(file_ext, '.csv')
    #df.to_csv(save_name, index=False, sep='\t')

    # Plotting usage stats:
    plot_aa_usage(df, file_save_loc, file_ext)

