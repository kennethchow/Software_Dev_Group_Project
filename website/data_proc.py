from website.stats_funcs import *
import zarr
import allel

# Use disk as storage rather than memory:
allel.chunked.storage_registry["default"] = allel.chunked.storage_zarr.ZarrTmpStorage()

# Data Paths:
zarr_path = 'website/data/FINAL_30x_GR38_NoBiIndel.zarr'
pop_data_path = 'website/data/pop_data.panel'


def filter_data(chrom, start_pos, stop_pos, rs_val, gene_name, stats, pops):
    """ Loading the 1000 Genomes Data (Stored in dask array in zarr file) """

    # Opening the top level folder of the stored zarr file:
    callset = zarr.open_group(zarr_path, mode='r')

    # Determine which population-specific data to load:
    pop_data_dict = {'AFR': ['AF_AFR', 'DAF_AFR', 'GF_HET_AFR', 'GF_HOM_REF_AFR', 'GF_HOM_ALT_AFR'],
                     'AMR': ['AF_AMR', 'DAF_AMR', 'GF_HET_AMR', 'GF_HOM_REF_AMR', 'GF_HOM_ALT_AMR'],
                     'EAS': ['AF_EAS', 'DAF_EAS', 'GF_HET_EAS', 'GF_HOM_REF_EAS', 'GF_HOM_ALT_EAS'],
                     'EUR': ['AF_EUR', 'DAF_EUR', 'GF_HET_EUR', 'GF_HOM_REF_EUR', 'GF_HOM_ALT_EUR'],
                     'SAS': ['AF_SAS', 'DAF_SAS', 'GF_HET_SAS', 'GF_HOM_REF_SAS', 'GF_HOM_ALT_SAS'],
                     }
    pop_var_data = []
    for pop in pops:
        [pop_var_data.append(x) for x in pop_data_dict[pop]]

    # Specifying default data to always load:
    def_data = ['CHROM', 'POS', 'REF', 'ALT', 'GENE', 'RS_VAL']
    var_data = def_data + pop_var_data

    # Loading SNP Variants Data:
    variants = allel.VariantChunkedTable(callset[chrom]['variants'],
                                         names=var_data,
                                         index=('CHROM', 'POS'))

    # Extract variant positions and store in array:
    pos = variants['POS'][:]

    # Loading Genotypes Data:
    genotypes = allel.GenotypeChunkedArray(callset[chrom]['calldata']['GT_BI'])

    # Loading Phased Genotypes Data:
    phased_genotypes = allel.GenotypeChunkedArray(callset[chrom]['calldata']['PH_GT_BI'])

    # Loading SNP Position data for Phased Genotypes:
    ph_pos = np.array(callset[chrom]['variants']['PH_POS'])

    # Loading population information:
    pop_data = pd.read_csv(pop_data_path, sep='\t')

    # Drop any Unnamed columns:
    pop_data = pop_data.loc[:, ~pop_data.columns.str.contains('^Unnamed')]

    """ Filtering the data according to the user specifications """

    # ---- POSITIONAL INFORMATION (Chromosome, start, stop) ---- #
    # Check if user has entered positional information:
    pos_info = [chrom, start_pos, stop_pos]

    # If they have, then reduce the data based on this criteria:
    if all(s.strip() for s in pos_info):

        # Returning indices of where to slice the data:
        # Except is used when the locations are out of the range of the data.
        start_array = (pos - int(start_pos))
        try:
            start_idx = list(start_array).index(min(start_array[start_array >= 0]))
        except:
            start_idx = 0
        stop_array = (pos - int(stop_pos))
        try:
            stop_idx = list(stop_array).index(max(stop_array[stop_array <= 0]))
        except:
            stop_idx = 0

        # Similarly for the phased data:
        ph_start_array = (ph_pos - int(start_pos))
        try:
            ph_start_idx = list(ph_start_array).index(min(ph_start_array[ph_start_array >= 0]))
        except:
            ph_start_idx = 0
        ph_stop_array = (ph_pos - int(stop_pos))
        try:
            ph_stop_idx = list(ph_stop_array).index(max(ph_stop_array[ph_stop_array <= 0]))
        except:
            ph_stop_idx = 0

        # Reducing Data According to Positions specified:
        pos_mask = create_index_mask(variants, start_idx, stop_idx + 1)
        variants = variants.compress(pos_mask)
        pos = pos[start_idx:stop_idx + 1]
        genotypes = subset_G_array(genotypes, start_idx, stop_idx + 1, None, 'index')
        phased_genotypes = subset_G_array(phased_genotypes, ph_start_idx, ph_stop_idx + 1, None, 'index')
        ph_pos = ph_pos[ph_start_idx:ph_stop_idx + 1]

    # ----       GENE NAME        ---- #
    # Check if user has entered gene name information:
    if any(gene_name.strip()):

        # Create a Boolean mask array:
        gene_mask = variants['GENE'][:] == gene_name

        # If mask has >0 True values, then a match has been found:
        if np.count_nonzero(gene_mask) > 0:
            variants = variants.compress(gene_mask)
            pos = pos.compress(gene_mask)
            genotypes = subset_G_array(genotypes, None, None, gene_mask, 'mask')

            # Return the positions of the first and last SNPs in that gene:
            gene_start_pos = pos.min()
            gene_stop_pos = pos.max()

            # Subset phased data by matching to the positions in the original data:
            ph_gn_start_arr = (ph_pos - gene_start_pos)

            # Ensuring these positions fall within the range of this data:
            try:
                ph_gene_start = list(ph_gn_start_arr).index(min(ph_gn_start_arr[ph_gn_start_arr >= 0]))
            except:
                ph_gene_start = 0
            ph_gn_stop_arr = (ph_pos - gene_stop_pos)
            try:
                ph_gene_stop = list(ph_gn_stop_arr).index(max(ph_gn_stop_arr[ph_gn_stop_arr <= 0]))
            except:
                ph_gene_stop = -1

            # Filtering the phased data:
            phased_genotypes = subset_G_array(phased_genotypes, ph_gene_start, ph_gene_stop + 1, None, 'index')
            ph_pos = ph_pos[ph_gene_start:ph_gene_stop + 1]

        else:
            variants = ""
            pos = ""
            genotypes = ""
            phased_genotypes = ""
            ph_pos = ""

    # ---- SNP INFORMATION (Rs Value/(s)) ---- #
    # Check if user has entered rs value information:
    if any(rs_val.strip()):

        # Create a Boolean mask array:
        rs_val_mask = variants['RS_VAL'][:] == rs_val

        # If mask has >0 True values, then a match has been found:
        if np.count_nonzero(rs_val_mask) > 0:
            variants = variants.compress(rs_val_mask)
            pos = pos.compress(rs_val_mask)
            genotypes = subset_G_array(genotypes, None, None, rs_val_mask, 'mask')

            # Find the matching phased genotypes index position:
            pos_loc = pos[0]
            # For where there is a matched phased genotype:
            try:
                ph_pos_idx = np.where(ph_pos == pos_loc)[0][0]
                phased_genotypes = subset_G_array(phased_genotypes, ph_pos_idx, ph_pos_idx + 1, None, 'index')
                ph_pos = ph_pos[ph_pos_idx]
            except:
                phased_genotypes = ""
                ph_pos = ""

        else:
            variants = ""
            pos = ""
            genotypes = ""
            phased_genotypes = ""
            ph_pos = ""

    """ Filtering and calculating according to populations and statistics chosen"""

    # A string is returned for these dataframes if no matching SNPs for query:
    if isinstance(genotypes, str) | isinstance(variants, str):
        stats_df = "Query returned no matching SNPs."
        fst_df = ""
        ac_seg = ""
        seg_pos = ""

    else:
        stats_df, fst_df, ac_seg, seg_pos = PopulationFiltering(pop_data, stats,
                                                                pops, genotypes, phased_genotypes, variants)

    return stats_df, fst_df, ac_seg, seg_pos, variants



