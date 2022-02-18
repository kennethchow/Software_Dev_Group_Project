import allel
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def plot_seqdiv(seg_pos, ac_seg, pop1, plot_name='GenDiv_graph', chrom='22'):

    if len(seg_pos) >= 2001:
        wind_size = 2000
    elif len(seg_pos) >= 101:
        wind_size = 100
    elif len(seg_pos) >= 11:
        wind_size = 10
    else:
        wind_size = 1

    windows = allel.moving_statistic(seg_pos, statistic=lambda v: [v[0], v[-1]], size=wind_size)

    # Compute Sequence Diversity
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


def plot_tajimas(seg_pos, ac_seg, pop1, pop2, plot_name='tajimasd_graph', chrom='22'):
    if len(seg_pos) >= 2001:
        wind_size = 2000
    elif len(seg_pos) >= 100:
        wind_size = 100
    elif len(seg_pos) >= 10:
        wind_size = 10
    else:
        wind_size = 1

    windows = allel.moving_statistic(seg_pos, statistic=lambda v: [v[0], v[-1]], size=wind_size)
    x = np.asarray(windows).mean(axis=1)

    # Compute Tajima's D
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


def plot_fst(ac1, ac2, pos, plot_name='fst_graph', blen=2000, chrom='22'):
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

