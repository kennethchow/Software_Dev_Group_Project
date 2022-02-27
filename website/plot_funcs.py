import allel
import numpy as np
import pandas as pd


def plot_data_seqdiv(seg_pos, ac_seg, pop1, pop2, w_size, step):
    """Compute DataFrame of Nucleotide Diversity for one/two populations"""
    windows = get_windows(seg_pos, w_size, step)
    # use the block centres as the X coordinate
    x = np.asarray(windows).mean(axis=1)

    # Nucleotide Diversity
    y_1 = allel.windowed_diversity(seg_pos, ac_seg[pop1][:], windows=windows)[0]
    # Reduce data to only windows and values where SNPs are present:
    x_1, y_1 = remove_nans(x, y_1)
    # Create dictionary of x and y values:
    d_pop1 = {'chrom_pos': x_1, 'seq_div': y_1, 'population': pop1}

    # If there is a second population selected, return a df containing both values:
    if pop2 is not None:
        y_2 = allel.windowed_diversity(seg_pos, ac_seg[pop2][:], windows=windows)[0]
        x_2, y_2 = remove_nans(x, y_2)
        d_pop2 = {'chrom_pos': x_2, 'seq_div': y_2, 'population': pop2}
        df = pd.concat(
            [pd.DataFrame(data=d_pop1), pd.DataFrame(data=d_pop2)] )

    else:
        df = pd.DataFrame(data=d_pop1)

    return df


def plot_data_watt_thet(seg_pos, ac_seg, pop1, pop2, w_size, step):
    """Compute DataFrame of Watterson's Theta for one/two populations"""
    windows = get_windows(seg_pos, w_size, step)
    # use the block centres as the X coordinate
    x = np.asarray(windows).mean(axis=1)

    # Watterson's Theta
    y_1 = allel.windowed_watterson_theta(seg_pos, ac_seg[pop1][:], windows=windows)[0]
    # Reduce data to only windows and values where SNPs are present:
    x_1, y_1 = remove_nans(x, y_1)
    # Create dictionary of x and y values:
    d_pop1 = {'chrom_pos': x_1, 'watt_thet': y_1, 'population': pop1}

    # If there is a second population selected, return a df containing both values:
    if pop2 is not None:
        y_2 = allel.windowed_watterson_theta(seg_pos, ac_seg[pop2][:], windows=windows)[0]
        x_2, y_2 = remove_nans(x, y_2)
        d_pop2 = {'chrom_pos': x_2, 'watt_thet': y_2, 'population': pop2}
        df = pd.concat(
            [pd.DataFrame(data=d_pop1), pd.DataFrame(data=d_pop2)] )

    else:
        df = pd.DataFrame(data=d_pop1)

    return df


def plot_data_tajimas(seg_pos, ac_seg, pop1, pop2, w_size, step):
    """Compute DataFrame of Tajima's D for one/two populations"""
    windows = get_windows(seg_pos, w_size, step)
    # use the block centres as the X coordinate
    x = np.asarray(windows).mean(axis=1)

    # Compute Tajima's D
    y_1, _, _ = allel.windowed_tajima_d(seg_pos, ac_seg[pop1][:], windows=windows)
    # Reduce data to only windows and values where SNPs are present:
    x_1, y_1 = remove_nans(x, y_1)
    # Create dictionary of x and y values:
    d_pop1 = {'chrom_pos': x_1, 'tajima_d': y_1, 'population': pop1}

    # If there is a second population selected, return a df containing both values:
    if pop2 is not None:
        y_2, _, _ = allel.windowed_tajima_d(seg_pos, ac_seg[pop2][:], windows=windows)
        x_2, y_2 = remove_nans(x, y_2)
        d_pop2 = {'chrom_pos': x_2, 'tajima_d': y_2, 'population': pop2}
        df = pd.concat(
            [pd.DataFrame(data=d_pop1), pd.DataFrame(data=d_pop2)])

    else:
        df = pd.DataFrame(data=d_pop1)

    return df


def plot_data_fst(seg_pos, ac_seg, pop1, pop2, w_size, step):
    """Compute DataFrame of Fst for comparing two populations"""
    # use the block centres as the X coordinate
    windows = get_windows(seg_pos, w_size, step)
    # use the block centres as the X coordinate
    x = np.asarray(windows).mean(axis=1)

    # Compute Fst:
    y, _, _ = allel.windowed_hudson_fst(seg_pos, ac_seg[pop1][:], ac_seg[pop2][:],
                                          windows=windows, fill=np.nan)

    # Reduce data to only windows and values where SNPs are present:
    x_1, y_1 = remove_nans(x, y)

    # make into df for plotting
    d_pops = {'chrom_pos': x_1, 'fst': y_1, 'population': pop1 + ' vs. ' + pop2}
    df = pd.DataFrame(data=d_pops)

    return df


def get_windows(seg_pos, w_size, step):
    # Return the start and end positions of the variants selected:
    min_pos = min(seg_pos)
    max_pos = max(seg_pos)

    # Create an array between these positions:
    pos = np.array(range(min_pos, max_pos + 1))

    # Calculate windows based on size specified by user:
    windows = allel.moving_statistic(pos, statistic=lambda v: [v[0], v[-1]], size=w_size, step=step)

    return windows


def remove_nans(x, y):
    """ Removes nan values when no SNPs are present in a given window """
    # Create a mask of the nan values returned by the summary stat:
    nan_mask = np.isnan(y)

    x_no_nan = x[nan_mask == False]
    y_no_nan = y[nan_mask == False]

    return x_no_nan, y_no_nan


