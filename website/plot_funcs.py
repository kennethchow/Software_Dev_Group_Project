import allel
import numpy as np
import pandas as pd


def plot_data_seqdiv(seg_pos, ac_seg, pop1, pop2, w_size):
    """Compute DataFrame of Genetic Diversity for two populations"""
    windows = allel.moving_statistic(seg_pos, statistic=lambda v: [v[0], v[-1]], size=w_size)
    # use the block centres as the X coordinate
    x = np.asarray(windows).mean(axis=1)
    # x = [x[0] for x in windows]

    # Sequence Diversity
    y_1 = allel.windowed_diversity(seg_pos, ac_seg[pop1][:], windows=windows)[0]
    y_2 = allel.windowed_diversity(seg_pos, ac_seg[pop2][:], windows=windows)[0]

    # make into df for plotting
    d_pop1 = {'chrom_pos': x, 'seq_div': y_1, 'population': pop1}
    d_pop2 = {'chrom_pos': x, 'seq_div': y_2, 'population': pop2}
    df = pd.concat(
        [pd.DataFrame(data=d_pop1), pd.DataFrame(data=d_pop2)]
    )
    return df


def plot_data_tajimas(seg_pos, ac_seg, pop1, pop2, w_size):
    """Compute DataFrame of Tajima's D for two populations"""
    windows = allel.moving_statistic(seg_pos, statistic=lambda v: [v[0], v[-1]], size=w_size)
    # use the block centres as the X coordinate
    x = np.asarray(windows).mean(axis=1)

    # Compute Tajima's D
    y_1, _, _ = allel.windowed_tajima_d(seg_pos, ac_seg[pop1][:], windows=windows)
    y_2, _, _ = allel.windowed_tajima_d(seg_pos, ac_seg[pop2][:], windows=windows)

    # make into df for plotting
    d_pop1 = {'chrom_pos': x, 'tajima_d': y_1, 'population': pop1}
    d_pop2 = {'chrom_pos': x, 'tajima_d': y_2, 'population': pop2}
    df = pd.concat(
        [pd.DataFrame(data=d_pop1), pd.DataFrame(data=d_pop2)]
    )
    return df


def plot_data_fst(seg_pos, ac_seg, pop1, pop2, w_size):
    """Compute DataFrame of Fst for comparing two populations"""
    # use the block centres as the X coordinate
    x = allel.moving_statistic(seg_pos, statistic=lambda v: (v[0] + v[-1]) / 2, size=w_size)

    # Compute Fst:
    fst, fst_se, vb, _ = allel.blockwise_hudson_fst(ac_seg[pop1][:], ac_seg[pop2][:], blen=w_size)

    # use the per-block average Fst as the Y coordinate
    y = vb

    # make into df for plotting
    d_pops = {'chrom_pos': x, 'fst': y, 'population': pop1 + ' vs. ' + pop2}
    df = pd.DataFrame(data=d_pops)

    return df
