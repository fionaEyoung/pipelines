# Summarise segmentation comparison metrics in various plots
#
# Fiona Young
# University College London
# February 2023

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from os import path, getlogin
import os, sys, argparse
import json
from itertools import combinations
import pandas as pd
from math import factorial
import re


DEBUG = False
METHODS = ['TG', 'TF', 'AT', 'TSX', 'TSD']
TRACTS = ['cst', 'or', 'af']
HEMS = ['l', 'r']

def set_box_colours(bp, c, index=None):

    for prop in bp.keys():
        if not index:
            plt.setp(bp[prop], color=c)
        else:
            plt.setp(bp[prop][index], color=c)

def main():

    P = argparse.ArgumentParser()
    P.add_argument('dataset', type=str)
    args = P.parse_args(sys.argv[1:])

    if args.dataset == 'tractoinferno':
        METHODS.append('TGR')

    try:
        root_dir = os.environ["DATAROOT"]
    except KeyError:
        root_dir = f"/Users/{getlogin()}/Documents/{'UCL/CDT' if getlogin()=='fiona' else 'Research'}/tractfinder"

    data_dir = path.join(root_dir, 'images_and_data')
    results_dir = path.join(root_dir, 'analyses')

    datafile = f'all_results_{args.dataset}.pkl'
    all_results = pd.read_pickle(path.join(results_dir, datafile))

    if 'TGR' in METHODS:
        lut = {'TF': {'TG': 0, 'TGR' : 1, 'TSD' : 2, 'TSX': 3, 'AT': 4},
               'TG': {'TF': 0, 'TGR' : 5, 'TSD' : 6, 'TSX': 7, 'AT': 8},
               'TGR':{'TF': 1, 'TG' : 5, 'TSD' : 9, 'TSX': 10, 'AT': 11},
               'TSD' : {'TF': 2, 'TG': 6, 'TGR' : 9, 'TSX': 12, 'AT': 13},
               'TSX' : {'TF': 3, 'TG': 7, 'TGR' : 10, 'TSD': 12, 'AT': 14},
               'AT' : {'TF': 4, 'TG': 8, 'TGR' : 11, 'TSD': 13, 'TSX' : 14}}
    else:
        lut = {'TF': {'TG': 0, 'TSD' : 1, 'TSX': 2, 'AT': 3},
               'TG': {'TF': 0, 'TSD' : 4, 'TSX': 5, 'AT': 6},
               'TSD' : {'TF': 1, 'TG': 4, 'TSX': 7, 'AT': 8},
               'TSX' : {'TF': 2, 'TG': 5, 'TSD': 7, 'AT': 9},
               'AT' : {'TF': 3, 'TG': 6, 'TSD': 8, 'TSX' : 9}}
    n = len(METHODS)

    lut_ = [None] * int(factorial(n) / factorial(2) / factorial(n - 2))
    for m1 in lut:
        for m2 in lut[m1]:
            if not lut_[lut[m1][m2]]:
                lut_[lut[m1][m2]] = [m1, m2]


    ## Specify some plot settings
    # Colourmaps
    cmap_generalised = cm.get_cmap('inferno')
    cmap_binary = cm.get_cmap('viridis')
    # Boxplot parameters
    c = 'w'
    boxprops = dict(linestyle='-', color=c, facecolor=c)
    flierprops = dict(marker='o', markerfacecolor=c, markeredgecolor=c,
                       markersize=1.5)
    medianprops = dict(linestyle='-', linewidth=1, color='k')
    whiskerprops = dict(color=c)
    flierprops=dict(marker='.', markerfacecolor='k', markeredgecolor='none', alpha=0.5)
    meanpointprops = dict(marker='o', markeredgecolor='k', markersize=10,
                      markerfacecolor='k')
    meanlineprops = dict(linestyle='--', linewidth=2.5, color='purple')

    ## Which metrics to include
    include_metrics = ['binary_dice', 'weighted_dice',
                       'density_correlation',
                       'ba_schilling_volume', 'ba_schilling_signed_m1']
    n_m = len(include_metrics) # nr of metrics
    n_c = n-1 # nr of comparisons
    n_t = len(TRACTS) # nr of tracts

    ## For one against everyone else comparison: which method to compare
    if args.dataset=='tractoinferno':
        compare_against = 'TGR'
    else:
        compare_against = 'TG'

    #https://davidmathlogic.com/colorblind/#%23000000-%233D840B-%23F38B03-%235D32B1
    tract_colours = {'af' :'#F38B03',
                     'or' :'#3D840B',
                     'cst':'#5D32B1'}
    method_colours = {'TSX': '#FF8100',
                      'TSD': '#C80000',
                      'TG' : '#7797A7',
                      'TGR': '#1DAFFF',
                      'TF' : '#FFC818',
                      'AT' : '#F046A8'}
    method_shapes = {'TSX': 'v',
                      'TSD': '^',
                      'TG' : 's',
                      'TGR': 's',
                      'TF' : 'o',
                      'AT' : 'd'}
    metric_params = {'binary_dice': {'lims': [0,1],
                                     'title': 'Binary DSC'},
                   'hd':            {'lims': [0,40],
                                    'title': 'Hausdorff Distance'},
                   'ba_schilling_boundary': {'lims': [0,15],
                                    'title': 'Bundle Adjacency (Boundary)'},
                   'ba_schilling_volume': {'lims': [0,15],
                                    'title': 'Bundle Distance'},
                   'weighted_dice': {'lims': [0,1],
                                     'title': 'Weidthed DSC'},
                   'density_correlation': {'lims': [0,1],
                                     'title': 'Density Correlation'},
                   'ba_schilling_signed_m1': {'lims': [-15,15],
                                     'title': 'Signed Bundle Distance'},}

    order = METHODS.copy()
    order.remove(compare_against)

    ## All metrics box plots, grouped by method
    all_metrics_fig, all_metrics_axs = plt.subplots(nrows=1,
                                                    ncols=len(include_metrics),
                                                    figsize=(5*len(include_metrics), 4))

    for metric, ax in zip(include_metrics, all_metrics_axs):

        x = np.arange(1, (n_t+1)*n_c, (n_t+1))

        one_of_each = []

        for i, tract in enumerate(TRACTS):

            # Left and right all combined
            mask = ((all_results.methods.apply(lambda x: compare_against in x)) &
                    (all_results['tract'] == tract))

            if not metric == 'ba_schilling_signed_m1':

                box = all_results[mask].assign(
                        group=lambda df: df.methods.apply(lambda x: x.difference({compare_against}).pop()
                                                          ).astype(pd.CategoricalDtype(order, ordered=True))
                        ).boxplot(column=metric, by='group', ax=ax, grid=False,
                                          return_type='dict', patch_artist=True,
                                          flierprops=flierprops,
                                          boxprops=dict(alpha=0.7),
                                          positions=x+i)
                set_box_colours(box[metric], tract_colours[tract])
                one_of_each.append(box[metric]['boxes'][0])

            else:
                # special case for signed
                box = all_results[mask].assign(
                    group=lambda df: df.methods.apply(lambda x: x.difference({compare_against}).pop()
                                                      ).astype(pd.CategoricalDtype(order, ordered=True)),
                    signed=lambda df: np.where(df.method1.apply(lambda x: x['name']) == df.group, df.ba_schilling_signed_m1, df.ba_schilling_signed_m2)
                ).boxplot(column='signed', by='group', ax=ax, grid=False,
                          return_type='dict', patch_artist=True,
                          flierprops=flierprops,
                          boxprops=dict(alpha=0.7),
                          positions=x+i)
                set_box_colours(box['signed'], tract_colours[tract])
                one_of_each.append(box['signed']['boxes'][0])

        if metric == 'ba_schilling_signed_m1':
            ax.axhline(y=0, c='#8a8a8a', ls='--')

        ax.set_xticks(x+1)
        # ax.set_xticklabels(METHODS[1:])
        ax.set_ylabel(metric_params[metric]['title'])
        ax.set_xlabel('')
        ax.tick_params(axis='x', length=0)
        ax.set_ylim(metric_params[metric]['lims'])

        ax.legend(one_of_each, TRACTS)
        ax.set_title("")

    all_metrics_fig.suptitle("")
    all_metrics_fig.savefig(path.join(results_dir, f'all_metrics_{args.dataset}.png'),
                transparent=False, dpi=120, bbox_inches="tight")

    ## All metrics box plots, grouped by tract
    all_metrics_fig, all_metrics_axs = plt.subplots(nrows=1,
                                                    ncols=len(include_metrics),
                                                    figsize=(5*len(include_metrics), 4))

    for metric, ax in zip(include_metrics, all_metrics_axs):

        x = np.arange(1, (n_c+1)*n_t, (n_c+1))

        one_of_each = []

        for i, method in enumerate(order):

            # Left and right all combined
            mask = (all_results['methods'] == {compare_against, method})

            if not metric == 'ba_schilling_signed_m1':
                box = all_results[mask].boxplot(column=metric, by='tract', ax=ax, grid=False,
                                          return_type='dict', patch_artist=True,
                                          flierprops=flierprops,
                                          boxprops=dict(alpha=0.7),
                                          positions=x+i)
                set_box_colours(box[metric], method_colours[method])
                one_of_each.append(box[metric]['boxes'][0])

            else:
                # special case for signed
                box = all_results[mask].assign(
                    signed=lambda df: np.where(df.method1.apply(lambda x: x['name']) == method, df.ba_schilling_signed_m1, df.ba_schilling_signed_m2)
                ).boxplot(column='signed', by='tract', ax=ax, grid=False,
                          return_type='dict', patch_artist=True,
                          flierprops=flierprops,
                          boxprops=dict(alpha=0.7),
                          positions=x+i)
                set_box_colours(box['signed'], method_colours[method])
                one_of_each.append(box['signed']['boxes'][0])

        if metric == 'ba_schilling_signed_m1':
            ax.axhline(y=0, c='#8a8a8a', ls='--')

        ax.set_xticks(x+2)
        ax.set_xticklabels([t.get_text().upper() for t in ax.get_xticklabels()])
        ax.set_ylabel(metric_params[metric]['title'])
        ax.set_xlabel('')
        ax.tick_params(axis='x', length=0)
        ax.set_ylim(metric_params[metric]['lims'])

        ax.legend(one_of_each, order, ncol=len(order),
                columnspacing=1, markerscale=0.5, handlelength=1)
        ax.set_title("")

    all_metrics_fig.suptitle("")
    all_metrics_fig.savefig(path.join(results_dir, f'all_metrics_by_tract_{args.dataset}.png'),
                transparent=False, dpi=120, bbox_inches="tight")

    ## DSC matrices
    if args.dataset == 'clinical':
        all_results['set']=all_results.method1.apply(lambda x: re.findall(r'\w*imri', x['path'])[0])
    split_on = 'hem' # 'hem' or 'set'
    left_right = HEMS # HEMS or ['gosh_imri', 'imri']
    for tract in TRACTS:

        fig = plt.figure(figsize=(8, 8))
        gs = fig.add_gridspec(n, n, hspace=0, wspace=0)
        ax = gs.subplots(sharex=True, sharey=True)

        inds1 = np.triu_indices(n,1)

        for i, (m1, m2) in enumerate(combinations(METHODS,2)):

            mask = ((all_results['methods'] == {m1, m2} ) &
                   (all_results['tract'] == tract ) )#&
                   # (all_results['hem'] == h))

            for s in left_right:
                yi = all_results[mask & (all_results[split_on] == s)]['weighted_dice']
                yj = all_results[mask & (all_results[split_on] == s)]['binary_dice']

                idx = lut[m1][m2]

                for j, y in enumerate((yi, yj)):
                    row, col = inds1[j][idx], inds1[0 if j else 1][idx]
                    ax_ = ax[row, col]

                    ax_.axvspan(-2 if s==left_right[0] else 2, 0,
                                      color = cmap_binary(y.mean()) if j else cmap_generalised(y.mean()))

                    box = ax_.boxplot(y, positions=[-0.5 + (s==left_right[1])], widths=[0.8],
                                      patch_artist=True, showcaps=0,
                                      boxprops=boxprops, flierprops=flierprops,
                                      medianprops=medianprops,
                                      whiskerprops=whiskerprops)

                    # If first column or top row
                    if row == 0:
                            ax_.set_title(lut_[idx][1])
                            #print(row, col, m1, m2)
                    if col == 0:
                            ax_.set_ylabel(lut_[idx][1])
                            #print(row, col, m1, m2)

        ax[0][0].set_ylabel(lut_[0][0])
        ax[0][0].set_title(lut_[0][0])

        plt.setp(ax, xlim=[-2,2], ylim=[0,1], yticks=[0.1,0.5,0.9], xticks=[])

        cbar_ax = fig.add_axes([.95, 0.3, 0.05, .577])
        cb1 = plt.colorbar(cm.ScalarMappable(norm=None, cmap=cmap_generalised),
                     cax=cbar_ax, label='generalised Dice score')
        cb1.set_label('mean generalised Dice score', fontsize=15)
        cbar_ax = fig.add_axes([.125, 0, 0.577, .05])
        cb2 = plt.colorbar(cm.ScalarMappable(norm=None, cmap=cmap_binary),
                     cax=cbar_ax, orientation='horizontal')
        cb2.set_label('mean binary Dice score', fontsize=15)


        fig.suptitle(f'{tract.upper()}', fontsize=20)
        fig.savefig(path.join(results_dir, f'{tract}_pairwise_{args.dataset}.png'),
                    transparent=False, dpi=120, bbox_inches="tight")
        plt.cla()

    ## Scatter plots
    for tract in TRACTS:

        distance_metric = 'ba_schilling_signed_m1'
        volume_metric = 'binary_dice'

        mask_l = ((all_results['tract'] == tract) &
                (all_results['hem'] == 'l'))
        mask_r = ((all_results['tract'] == tract) &
                (all_results['hem'] == 'r'))

        fig, ax = plt.subplots(figsize=(8,8))

        for m in METHODS:

            if distance_metric == 'ba_schilling_signed_m1':
                x = np.where(all_results.method1.apply(lambda x: x['name']) == m,
                             all_results.ba_schilling_signed_m1,
                             all_results.ba_schilling_signed_m2)
                x_l = x[mask_l & (all_results['methods'] == {m, compare_against})]
                x_r = x[mask_r & (all_results['methods'] == {m, compare_against})]
            else:
                x_l = all_results[mask_l & (all_results['methods'] == {m, compare_against})
                                ][distance_metric]
                x_r = all_results[mask_r & (all_results['methods'] == {m, compare_against})
                                ][distance_metric]

            y_l = all_results[mask_l & (all_results['methods'] == {m, compare_against})
                             ][volume_metric]
            y_r = all_results[mask_r & (all_results['methods'] == {m, compare_against})
                             ][volume_metric]


            ax.plot(x_l,y_l,method_shapes[m], color=method_colours[m], fillstyle='full',
                    label=m+' left')
            ax.plot(x_r,y_r,method_shapes[m], color=method_colours[m], fillstyle='none',
                    label=m+' right')

        ax.set_title(f'{tract.upper()}', fontsize=20)
        ax.set_xlabel(distance_metric)
        ax.set_ylabel(volume_metric)
        ax.set_xlim(metric_params[distance_metric]['lims'])
        #ax.set_ylim(metric_params[volume_metric]['lims'])
        ax.legend()

        fig.savefig(path.join(results_dir, f'{tract}_scatter_{args.dataset}.png'))


if __name__ == '__main__':
    main()
