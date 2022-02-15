#### GENETIC DIVERSITY ####
def GeneticDiversity(seg_pos, ac_seg, pop1, start, stop):
    # Calculating the sequence diversity
    pi = allel.sequence_diversity(seg_pos, ac_seg[pop1][:], start, stop)

    return pi


def plot_GenDiv(seg_pos, ac_seg, pop1, plot_name='GenDiv_graph'):
    if len(seg_pos) >= 2001:
        wind_size = 2000
    elif len(seg_pos) >= 101:
        wind_size = 100
    elif len(seg_pos) >= 11:
        wind_size = 10
    else:
        wind_size = 1

    windows = allel.moving_statistic(seg_pos, statistic=lambda v: [v[0], v[-1]], size=wind_size)

    # Compute Genetic Diversity
    y = allel.windowed_diversity(seg_pos, ac_seg[pop1][:], windows=windows)[0]

    # use the block centres as the X coordinate
    x = windows

    # plot
    fig, ax = plt.subplots(figsize=(12, 4))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y, 'k-', lw=.5)
    plt.title('Slidng window showing the Genetic Diversity through chromosome 22')
    ax.set_ylabel('Genetic Diversity')
    ax.set_xlabel('Chromosome %s position (bp)' % chrom)
    ax.set_xlim(seg_pos.min(), seg_pos.max())
    plt.savefig(f'{plot_name}.png')

    return fig



#### TAJIMAS D ####
def plot_tajimas(seg_pos, ac_seg, pop1, pop2, plot_name='tajimasd_graph'):
    if len(seg_pos) >= 2001:
        wind_size = 2000
    elif len(seg_pos) >= 100:
        wind_size = 100
    elif len(seg_pos) >= 10:
        wind_size = 10
    else:
        wind_size = 1

    windows = allel.moving_statistic(pos, statistic=lambda v: [v[0], v[-1]], size=wind_size)
    x = np.asarray(windows).mean(axis=1)

    # compute Tajima's D
    y1, _, _ = allel.windowed_tajima_d(seg_pos, ac_seg[pop1][:], windows=windows)
    y2, _, _ = allel.windowed_tajima_d(seg_pos, ac_seg[pop2][:], windows=windows)

    # plot
    fig, ax = plt.subplots(figsize=(12, 4))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y1, lw=.5, label='AFR')
    ax.plot(x, y2, lw=.5, label='EUR')
    plt.title("Plot of Tajima's $D$ in Window over Chromosome %s" % chrom)
    ax.set_ylabel("Tajima's $D$")
    ax.set_xlabel('Chromosome %s position (bp)' % chrom)
    ax.set_xlim(seg_pos.min(), seg_pos.max())
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1));
    plt.savefig(f'{plot_name}.png')

    return fig



#### FST VALUES AND PLOT ####
def plot_fst(ac1, ac2, pos, plot_name='fst_graph', blen=2000):
    fst, fst_se, vb, _ = allel.blockwise_hudson_fst(ac1, ac2, blen=blen)

    # use the per-block average Fst as the Y coordinate
    y = vb

    # use the block centres as the X coordinate
    x = allel.moving_statistic(pos, statistic=lambda v: (v[0] + v[-1]) / 2, size=blen)

    # plot
    fig, ax = plt.subplots(figsize=(12, 4))
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y, 'k-', lw=.5)
    plt.title("Plot of Fst in Window over Chromosome %s\n For Populations AFR and EUR" % chrom)
    ax.set_ylabel('$F_{ST}$')
    ax.set_xlabel('Chromosome %s position (bp)' % chrom)
    ax.set_xlim(0, pos.max())
    plt.savefig(f'{plot_name}.png')
    # F string, aka function string. Can add any number or list

    return fig



#### MASSIVE USER INPUT FUNCTION ####
def UserInput(query_id, start, stop, pop1, pop2, chrom='22', m_allele=3):
    # Still need to add gene name, rsID

    # hardcoded the zarr path
    zarr_path = 'Documents/GitHub/Software_Dev_Group_Project/website/data/ALL_30X_Chr22_GR38.zarr'

    ### GENOTYPE DATA ###
    # Opening the zarr file
    callset = zarr.open_group(zarr_path, mode='r')

    # Creating variants table
    variants = allel.VariantChunkedTable(callset[chrom]['variants'],
                                         names=['POS', 'REF', 'ALT', 'GENE', 'RS_VAL', ],
                                         index='POS')

    # Extract variant positions and store in array:
    pos = variants['POS'][:]

    print('Positions extracted.')

    ## POPULATION DATA ###
    # need to hardcode the population data
    pop_path = 'Documents/integrated_call_samples_v3.20130502.ALL.panel'
    pop_data = pd.read_csv(pop_path, sep='\t')

    # Drop any Unnamed columns from population data:
    pop_data = pop_data.loc[:, ~pop_data.columns.str.contains('^Unnamed')]

    print('Population data processed')

    # Create Boolean array classifying population membership:
    sample_selection = pop_data.super_pop.isin({pop1, pop2}).values

    # Create a subset of the pop_data table based on this selection:
    samples_subset = pop_data[sample_selection].reset_index(drop=True)

    # Subsetting genotypes
    calldata = callset[chrom]['calldata']
    genotypes = allel.GenotypeChunkedArray(calldata['GT'])
    genotypes_subset = genotypes.subset(sel0=None, sel1=sample_selection)

    # Allele count within the two populations listed, creating a dictionary:
    subpops = {
        'ALL': list(range(len(samples_subset))),
        pop1: samples_subset[samples_subset.super_pop == pop1].index.tolist(),
        pop2: samples_subset[samples_subset.super_pop == pop2].index.tolist(),
    }

    # Counting the alleles within the sub-populations:
    ac_subpops = genotypes_subset.count_alleles_subpops(subpops, max_allele=3)

    print("Allele count done")

    # Boolean array of those SNPs that do segregate:
    is_seg = ac_subpops['ALL'].is_segregating()[:]

    # Subset the two population genotypes subset to keep only the segregating SNPs:
    genotypes_seg = genotypes_subset.compress(is_seg, axis=0)

    # Subset the allele counts as well:
    ac_seg = ac_subpops.compress(is_seg)

    # And the variants data:
    variants_seg = variants.compress(is_seg)

    # Extract variant positions and store in array:
    seg_pos = variants_seg['POS'][:]

    # Calling the Genetic Diversity function:
    pop1_pi = GeneticDiversity(seg_pos, ac_seg, pop1, start, stop)
    pop2_pi = GeneticDiversity(seg_pos, ac_seg, pop2, start, stop)
    # Genetic Diversity for sliding window:
    GenDiv_plot = plot_GenDiv(seg_pos, ac_seg, pop1, plot_name='GenDiv_graph' + query_id)

    # Calculating fst:
    fst, fst_se, _, _ = allel.blockwise_hudson_fst(ac_seg[pop1], ac_seg[pop2], blen=100000)
    # Fst sliding window:
    fst_plot = plot_fst(ac_seg[pop1], ac_seg[pop2], variants_seg['POS'][:], plot_name='fst_graph' + query_id)

    # Calculating Tajima's D:
    # taj_pop1 = allel.windowed_tajima_d(seg_pos, ac_seg[pop1][:], windows=windows)
    # taj_pop2 = allel.windowed_tajima_d(seg_pos, ac_seg[pop2][:], windows=windows)
    # Tajimas slidng window:
    tajimas_plot = plot_tajimas(seg_pos, ac_seg, pop1, pop2, plot_name='tajimasd_graph' + query_id)

    return fst, fst_se, fst_plot, pop1_pi, pop2_pi, tajimas_plot

    # return haplotype

