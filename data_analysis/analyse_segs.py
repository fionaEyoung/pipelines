# Compare segmentations for given dataset
#
# Fiona Young
# University College London
# February 2023

import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm

from skimage.metrics import hausdorff_distance
from skimage.segmentation import find_boundaries
from scipy.ndimage import distance_transform_edt as edt

from itertools import combinations
from math import factorial

from os import path, getlogin
from glob import glob
import os, sys, json, argparse

METHODS = ['TF', 'TSX', 'TSD', 'TG', 'AT']
TRACTS = ['cst', 'or', 'af']
HEMS = ['l', 'r']

class Score:
    # Just a structure for holding computed metrics
    def __init__(self, method1, method2, tract, hem):
        if isinstance(method1, str):
            self.method1 = {'name' : method1}
        elif isinstance(method1, dict):
            assert 'name' in method1
            self.method1 = method1

        if isinstance(method2, str):
            self.method2 = {'name' : method1}
        elif isinstance(method2, dict):
            assert 'name' in method2
            self.method2 = method2

        self.tract = tract.lower()
        self.hem = hem.lower()[0]

def dsc(A,B):
    assert A.dtype=='bool' and B.dtype=='bool'
    return 2*np.logical_and(A,B).sum()/(np.count_nonzero(A)+np.count_nonzero(B))

def w_dsc(A,B):
    # As defined in 10.1016/j.nicl.2017.07.020
    assert ((0 <= A) & (A <= 1)).all() & ((0 <= B) & (B <= 1)).all()
    return (A[np.logical_and(A,B)].sum() + B[np.logical_and(A,B)].sum())/(A.sum() + B.sum())

def jon_dsc(A,B):
    assert ((0 <= A) & (A <= 1)).all() & ((0 <= B) & (B <= 1)).all()
    return 2 * np.sqrt((A*B)).sum() / (A.sum() + B.sum())

def g_dsc(A,B):
    # Alternative generalisation
    return 2 * (A*B).sum() / ((A**2).sum() + (B**2).sum())

def ba_schilling(A, B, vox=(1,1), boundary_only=False, signed=False):
    # Bundle adjacency as defined (interpreted) in Schilling et al. 2021 Neuroimage
    # (10.1016/j.neuroimage.2021.118502)
    # It is sort of like Hausdorff distance, except considering all non-overlapping
    # voxels, not just surface voxels, and taking the average, not max
    # Variations from original: using only boundary voxels, and signed version
    assert A.dtype=='bool' and B.dtype=='bool'

    # get distances map of background to foreground
    dist_A = edt(np.logical_not(A), sampling=vox[0])
    dist_B = edt(np.logical_not(B), sampling=vox[1]) * (-1 if signed else 1)

    if boundary_only:
        return np.hstack((dist_A[np.logical_and(find_boundaries(B, mode='inner'),
                                                np.logical_not(A))],
                          dist_B[np.logical_and(find_boundaries(A, mode='inner'),
                                                np.logical_not(B))])
                         ).mean()
    else:
        return np.hstack((dist_A[np.logical_and(B, np.logical_not(A))],
                          dist_B[np.logical_and(A, np.logical_not(B))])).mean()

def get_segmentation_results(dirname, filename_format, tracts=['cst', 'af', 'or'], hems=HEMS):

    assert all(s in filename_format for s in ['{h}', '{t}'])

    try_paths = {t.lower()+h[0].lower() : path.join(dirname, filename_format.format(h=h,t=t))
                for t in tracts for h in hems}

    if all( [path.exists(p) for p in try_paths.values()] ):
        return try_paths
    else:
        raise FileNotFoundError("""Unable to locate all segmentations with file name format """
                                         +filename_format+" in "
                                         +dirname)

def compare_segs(S, weighted_dice_version='Clayden'):

    if not isinstance(S, Score):
        raise ValueError("Input must be score object")

    I = [{'info' : m,
          'img'  : nib.load(m['path'])}
        for m in (S.method1, S.method2)]

    for img in I:
        img['data'] = img['img'].get_fdata()
        img['bin']  = img['data'] > img['info']['t_lower']
        img['vox']  = img['img'].header.get_zooms()[:3]

    ## Hausdorff distance
    S.hd = hausdorff_distance(find_boundaries(I[0]['bin'], mode='inner'),
                              find_boundaries(I[1]['bin'], mode='inner'))
    ## Density correlation
    S.density_correlation = np.corrcoef(I[0]['data'].flatten(),
                                        I[1]['data'].flatten()
                                        )[0,1]

    ## Distance (BA)
    S.ba_schilling_volume    = ba_schilling(I[0]['bin'], I[1]['bin'],
                                    vox=[img['vox'] for img in I])
    S.ba_schilling_boundary  = ba_schilling(I[0]['bin'], I[1]['bin'],
                                    vox=[img['vox'] for img in I], boundary_only=True)
    S.ba_schilling_signed_m1 = ba_schilling(I[0]['bin'], I[1]['bin'],
                                    vox=[img['vox'] for img in I], signed=True)
    S.ba_schilling_signed_m2 = -S.ba_schilling_signed_m1

    ## Binary dice
    S.binary_dice      = dsc(I[0]['bin'], I[1]['bin'])
    ## 'generalised' dice
    S.generalised_dice = g_dsc(I[0]['data'], I[1]['data'])
    ## Weighted dice (involves data rescaling, so must always be last)
    # GROSS DATA WRANGLING ENSUES
    for img in I:
        if img['info']['name'] in ('TF' , 'AT'):
            # rescale aribtrary forage values to [0,1]
            img['data'] /= img['data'].max()
            # numerical errors introduced by e.g. transformation give noisy values -> clip these
            img['data'][img['data']<0] = 0
        if 'TG' in img['info']['name']:
            if weighted_dice_version == 'Cousineau':
                img['data'] = np.clip((img['data'] - img['info']['t_lower'])/(img['info']['t_upper'] - img['info']['t_lower']),0,1)
            else:
                img['data'] /= img['data'].max()
        if img['info']['name'] in ('TSX', 'TSD'):
            if weighted_dice_version == 'Cousineau':
                img['data'][img['data']<0] = 0
            else:
                img['data'][img['data']<1e-3] = 0

    S.wdsc_version = weighted_dice_version
    if weighted_dice_version == 'Cousineau':
        S.weighted_dice = w_dsc(I[0]['data'], I[1]['data'])
    elif weighted_dice_version == 'Clayden':
        S.weighted_dice = jon_dsc(I[0]['data'], I[1]['data'])

def plot_results(ax, results, title, n=5):

    dsc_mat = np.zeros((n,n))
    dsc_mat[np.triu_indices(n,1)] = np.array([d.weighted_dice for d in results])
    dsc_mat[np.triu_indices(n,1)[::-1]] = np.array([d.binary_dice for d in results])
    labels = [results[0].method1['name']] + [d.method2['name'] for d in results if d.method1 == results[0].method1]

    ax.matshow(dsc_mat, cmap=plt.cm.plasma, vmin=0, vmax=1)
    for i in range(dsc_mat.shape[0]):
        for j in range(dsc_mat.shape[1]):
            ax.text(x=j, y=i,s=round(dsc_mat[i, j],5), va='center', ha='center')

    ax.set_xlabel('upper triangle: wDSC, lower triangle: DSC')
    ax.set_title(title)

    plt.setp(ax, xticks=np.arange(max(dsc_mat.shape)), xticklabels=labels,
                 yticks=np.arange(max(dsc_mat.shape)), yticklabels=labels)


def main():

    P = argparse.ArgumentParser()
    P.add_argument('--dataset', '-d', type=str, choices=['HCP', 'clinical', 'hcp', 'test', 'tractoinferno'], required=True)
    P.add_argument('--from_saved', '-s', action='store_true')
    P.add_argument('--weighted_dice_version', '-w', type=str, choices=['Cousineau', 'Clayden', 'Generalised'], default='Clayden')
    args = P.parse_args(sys.argv[1:])

    try:
        root_dir = os.environ["DATAROOT"]
    except KeyError:
        root_dir = f"/Users/{getlogin()}/Documents/{'UCL/CDT' if getlogin()=='fiona' else 'CDT'}/tractfinder"

    data_dir = path.join(root_dir, 'images_and_data')
    results_dir = path.join(root_dir, 'analyses')

    ## Collect data paths for different datasets
    if args.dataset == 'clinical':
        subject_list = ["gosh_imri/1/0", "gosh_imri/1/1",
                        "gosh_imri/3/0", "gosh_imri/3/1",
                        "gosh_imri/4/0", "gosh_imri/4/1",
                        "gosh_imri/ks/1",
                        "imri/1/0", "imri/1/4",
                        "imri/5/0", "imri/5/2",
                        "imri/6/0", "imri/6/3",
                        "imri/7/0", "imri/7/2"]
    elif args.dataset.lower() == 'hcp':
        subject_list =  [path.relpath(x, data_dir) for x in glob(path.join(data_dir, 'hcp/??????', 'dwi3'))]
    elif args.dataset == 'tractoinferno':
        METHODS.append('TGR')
        subject_list = [path.basename(x) for x in glob(path.join(data_dir, 'tractoinferno/sub-*'))]
        # Subjects to exclude because of failed tractography
        remove = ["sub-1245", "sub-1264", "sub-1265", "sub-1271",
                  "sub-1160", "sub-1167", "sub-1169", "sub-1221", "sub-1274"]
        subject_list = [s for s in subject_list if s not in remove]
    elif args.dataset == 'test':
        subject_list = ["hcp/100307/dwi3"]
    else:
        raise ValueError("Invalid dataset name. Options are 'clinical', 'HCP' or 'tractoinferno'")

    ## Threshold values for binarising continuous segmentations
    # Can be modified
    t_lower = {'TF' : 0.05,
               'TSX': 0.5,
               'TSD': 0.5,
               'TG' : 10,
               'TGR': 0,
               'AT' : 0.1}
    t_upper = {'TF' : None,
               'TSX': None,
               'TSD': None,
               'TG' : 100,
               'TGR': 100,
               'AT' : None}

    ## Construct elaborate lookup tables so it doesn't matter which method is 'first'
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

    ## List of comparisons, each containing all metrics
    all_results = []
    for s in subject_list:

        dir_path = path.join(data_dir, s)

        ## Inconsistent naming has delivered me here
        if args.dataset.lower() in [ 'hcp','test' ]:
            segmentation_paths = {
                'TSD': get_segmentation_results(path.join(dir_path,'tractseg_dkfz_ssst', 'bundle_segmentations'),
                            '{t}_{h}.nii.gz', tracts=[t.upper() for t in TRACTS], hems=['left', 'right']),
                'TSX': get_segmentation_results(path.join(dir_path,'tractseg_xtract_ssst', 'bundle_segmentations'),
                            '{t}_{h}.nii.gz'), # default args for hems and tracts
                'TG' : get_segmentation_results(path.join(dir_path,'rois_and_tracks'),
                            '{h}{t}.nii.gz'),
                'TF' : get_segmentation_results(dir_path,
                            '{h}_{t}_tractmap_msmt1.nii.gz'),
                'AT' : get_segmentation_results(path.join(dir_path,'..'),
                            '{h}_{t}_atlas.nii.gz')
                }
        elif args.dataset == 'clinical':
            segmentation_paths = {
                'TSD': get_segmentation_results(path.join(dir_path,'tractseg_dkfz_ssst', 'bundle_segmentations'),
                            '{t}_{h}.nii.gz', tracts=[t.upper() for t in TRACTS], hems=['left', 'right']),
                'TSX': get_segmentation_results(path.join(dir_path,'tractseg_xtract_ssst', 'bundle_segmentations'),
                            '{t}_{h}.nii.gz'),
                'TG' : get_segmentation_results(path.join(dir_path,'rois_and_tracks'),
                            '{h}{t}.nii.gz'),
                'TF' : get_segmentation_results(path.join(dir_path, 'tractfinder'),
                            '{h}_{t}_tractmap.nii.gz'),
                'AT' : get_segmentation_results(path.join(dir_path, 'tractfinder'),
                            '{h}_{t}_atlas.nii.gz')
                }
        elif args.dataset == 'tractoinferno':
            segmentation_paths = {
                'TSD': get_segmentation_results(path.join(dir_path,'tractseg', 'bundle_segmentations'),
                            '{t}_{h}.nii.gz', tracts=[t.upper() for t in TRACTS], hems=['left', 'right']),
                'TSX': get_segmentation_results(path.join(dir_path,'tractseg', 'bundle_segmentations'),
                            '{t}_{h}.nii.gz'),
                'TG' : get_segmentation_results(path.join(dir_path, 'tractography-ich'),
                            '{h}{t}.nii.gz'),
                'TF' : get_segmentation_results(path.join(dir_path, 'forage'),
                            '{h}_{t}_tractmap.nii.gz'),
                'AT' : get_segmentation_results(path.join(dir_path, 'forage'),
                            '{h}_{t}_atlas.nii.gz'),
                'TGR' : get_segmentation_results(path.join(dir_path, 'tractography'),
                            '{h}{t}.nii.gz')
                }
        else:
            raise ValueError("Invalid dataset name.")

        if not args.from_saved:
            fig, ax = plt.subplots()

        for tract in TRACTS:
            for h in HEMS:
                save_name = f'{tract}_{h}'
                outpath = path.join(dir_path, 'analysis', f'{save_name}.json')
                if not path.isdir(path.dirname(outpath)):
                    os.makedirs(path.dirname(outpath))

                if not args.from_saved:

                    results = [None]*len(lut_)

                    ## Compute metrics

                    for i, (m1, m2) in enumerate(combinations(METHODS,2)):

                        results[lut[m1][m2]] = Score(*[{'name':m, 't_lower':t_lower[m],
                                                        't_upper':t_upper[m],
                                                        'path':segmentation_paths[m][tract+h[0]]}
                                                        for m in (m1, m2)],
                                                     tract,
                                                     h)
                        results[lut[m1][m2]].methods = {m1, m2}
                        compare_segs(results[lut[m1][m2]], weighted_dice_version=args.weighted_dice_version)

                    with open(outpath, 'w') as f:
                        json.dump([{i:x.__dict__[i] for i in x.__dict__ if i!='methods'}  for x in results], f)

                    all_results += [x.__dict__ for x in results]

                    ## Plotting
                    plot_results(ax, results, save_name, n=n)
                    fig.savefig(outpath.replace('json', 'png'),
                                transparent=False, dpi=80, bbox_inches="tight")
                    plt.cla()

                else: # Aggregate precomputed results only
                    with open(outpath, 'r')  as f:
                        results = json.load(f)
                        for r in results:
                            r.update({'methods': {r['method1']['name'],
                                                  r['method2']['name']}})
                        all_results += results

        if not args.from_saved:
            plt.close(fig)

    ## Save all results for entire dataset
    with open(path.join(results_dir, f'all_results_{args.dataset.lower()}.json'), 'w') as f:
        json.dump({
            'thresholds':t_lower,
            'paths':segmentation_paths,
            'weighted DSC':args.weighted_dice_version,
            'attributes':list(all_results[0].keys()),
            }, f)

    all_results = pd.DataFrame(all_results)
    all_results.to_pickle(path.join(results_dir, f'all_results_{args.dataset.lower()}.pkl'))

if __name__ == '__main__':
    main()
