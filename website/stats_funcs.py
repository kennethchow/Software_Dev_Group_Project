import allel
import numpy as np
import pandas as pd
import itertools


# ======= Dictionaries used to Create Tables of Data to Display ======= #

pop_dict = {'AFR': 'African', 'AMR': 'American', 'EAS': 'East Asian',
            'EUR': 'European', 'SAS': 'South Asian'}

stats_dict = {'seq_div': 'Nucleotide Diversity', 'watt_thet': 'Watterson\'s Theta', 'hap_div': 'Haplotype Diversity',
              'taj_d': 'Tajima\'s D', 'fst': 'Fst', 'daf': 'Derived Allele Frequency'}


# ======= Functions used to process the Chromosome VCF Data and create the summary stats ======= #

def subset_G_array(data, start_idx, stop_idx, mask, subset_type):
    """ Function to retain chunked compression during subsetting """

    # To subset by index:
    if subset_type == 'index':
        data = allel.chunked.core.subset(data, range(start_idx, stop_idx))
    # To subset according to a mask:
    elif subset_type == 'mask':
        data = allel.chunked.core.compress(mask, data)

    # Convert back to a chunked Genotype array:
    data = allel.GenotypeChunkedArray(data)

    return data


def create_index_mask(data, start_idx, stop_idx):
    """ Creates a Boolean mask from start and stop indices """

    sel = np.asanyarray(range(start_idx, stop_idx))
    mask = np.zeros(len(data), dtype=bool)

    # Setting the values between the indices equal to True:
    mask[sel] = True

    return mask


def haplotype_diversity(genotype_chunked_array):
    """ Calculting haplotype diversity from a genotype chunked array """

    # converting genotype chunked array to genotype array
    genotype_array = genotype_chunked_array[:]

    # computing haplotype array from genotype array
    haplotype_array = genotype_array.to_haplotypes()

    # computing haplotype diversity
    hap_div = allel.haplotype_diversity(haplotype_array)

    return hap_div


def SummaryStats(stats, seg_pos, ac_seg, pop, pop_data, phased_genotypes):
    """ Calculate Summary Stats """
    comb_stats = []

    if 'seq_div' in stats:
        seq_div = allel.sequence_diversity(seg_pos, ac_seg[pop][:])
        comb_stats.append(seq_div)

    if 'watt_thet' in stats:
        watt_thet = allel.watterson_theta(seg_pos, ac_seg[pop][:])
        comb_stats.append(watt_thet)

    if 'taj_d' in stats:
        taj_d = allel.tajima_d(ac_seg[pop][:], pos=seg_pos)
        comb_stats.append(taj_d)

    if 'hap_div' in stats:
        # If there is no matching phased genotype then hap diversity cannot be calculated:
        if not isinstance(phased_genotypes, str):
            # Create Boolean array classifying population membership:
            sample_selection = pop_data.super_pop.isin({pop}).values

            # Subset the phased haplotype data based on population membership:
            ph_gen_subset = allel.GenotypeChunkedArray(allel.chunked.core.subset(
                phased_genotypes, sel0=None, sel1=sample_selection))

            # Calculate haplotype diversity from that subset:
            hap_div = haplotype_diversity(ph_gen_subset)
            comb_stats.append(hap_div)

        else:
            comb_stats.append('**')

    return comb_stats


def create_data_table(pops, stats, stats_data):
    """ Function to calculate the output table to display the summary stats """

    full_pop_names = []
    full_stats_names = []

    # Append which pops and stats have been selected by the user:
    for key in pops:
        full_pop_names.append(pop_dict[key])

    for st in stats:
        full_stats_names.append(stats_dict[st])

    # Removing fst and daf as they are not calculated here:
    try:
        full_stats_names.remove('Fst')
        full_stats_names.remove('Derived Allele Frequency')
    except:
        pass

    # Create a pandas dataframe using the summary stats:
    stats_df = pd.DataFrame(stats_data, columns=full_stats_names)
    stats_df.insert(loc=0, column='Population', value=full_pop_names)
    stats_df = stats_df.fillna('*')
    stats_df = stats_df.round(4)

    return stats_df


def PopulationFiltering(pop_data, stats, pops, genotypes, phased_genotypes, variants):
    """ Subsetting the data based on population membership """

    # Create Boolean array classifying population membership:
    sample_selection = pop_data.super_pop.isin(pops).values

    # Create a subset of the pop_data table based on this selection:
    samples_subset = pop_data[sample_selection].reset_index(drop=True)

    # Subsetting genotypes based on the selection mask (i.e. which populations to keep):
    # sel0 indicates variants (i.e. rows) to drop, sel1 indicates samples (i.e. columns) to drop:
    genotypes = allel.GenotypeChunkedArray(allel.chunked.core.subset(
        genotypes, sel0=None, sel1=sample_selection))

    """ Calculate allele count data for these populations"""

    # Create a dictionary to hold information about whether each sample is in a population:
    subpops = {'ALL': list(range(len(samples_subset)))}

    # Add new key and value for each population:
    for p in range(len(pops)):
        subpops[pops[p]] = samples_subset[samples_subset.super_pop == pops[p]].index.tolist()

    # Counting the alleles within the sub-populations (max alternate allele set to 5):
    ac_subpops = genotypes.count_alleles_subpops(subpops, max_allele=5)

    """ Subsetting Data Based on Whether the SNPs Segregate """

    # Boolean array of those SNPs that do segregate among the populations selected:
    is_seg = ac_subpops['ALL'].is_segregating()[:]

    #print("Is Seg Contents: ", is_seg)
    # If there are no segregating variants, then cannot calculate summary stats.
    if not any(is_seg):
        pop_stats = 'No segregating variants at this SNP.'

    else:
        # Subset the population data to keep only the segregating SNPs:
        ac_seg = ac_subpops.compress(is_seg)
        variants_seg = variants.compress(is_seg)
        seg_pos = variants_seg['POS'][:]

        """ Calculating Summary Stats """
        pop_stats = []

        for pop in pops:
            comb_stats = SummaryStats(stats, seg_pos, ac_seg, pop, pop_data, phased_genotypes)
            pop_stats.append(comb_stats)

    # A string message is returned in the case there are no segregating variants:
    if isinstance(pop_stats, str):
        stats_df = pop_stats
        fst_df = ""
        ac_seg = ""
        seg_pos = ""

    else:
        stats_df = create_data_table(pops, stats, pop_stats)

        """ Calc. Fst for Combinations of Populations """
        # Create a list of all the possible population combinations from those selected by the user:
        pop_combs = list(itertools.combinations(pops, 2))

        if 'fst' in stats:
            # Working out fst for each pair of populations:
            fst_comps = []
            fst_col_vals = []
            for c in range(len(pop_combs)):
                num, den = allel.hudson_fst(ac_seg[pop_combs[c][0]],
                                            ac_seg[pop_combs[c][1]])
                try:
                    fst = np.sum(num) / np.sum(den)
                except:
                    fst = 0
                fst_comps.append(fst)
                fst_col_vals.append("%s vs. %s" % (pop_dict[pop_combs[c][0]], pop_dict[pop_combs[c][1]]))

            # Create a pandas dataframe using the fst stats:
            fst_df = pd.DataFrame({'Populations Compared': fst_col_vals,
                                   'Fst Value': fst_comps},
                                  columns=['Populations Compared', 'Fst Value'])
            fst_df = fst_df.fillna('***')

        else:
            fst_df = ""

    return stats_df, fst_df, ac_seg, seg_pos