#!/bin/bash -l

# Wallclock time (format hours:minutes:seconds).
#$ -l h_rt=1:00:0

# Cores
#$ -pe smp 36

# RAM (must be an integer followed by M, G, or T)
#$ -l mem=4G

# TMPDIR space (default is 10 GB - remove if cluster is diskless)
# avg 5MB per .tck x 3 tracts x 2 hems x (2 + 3x3) x 49 sbj = ca. 16.2GB total, 330MB per subject. 66 tractograms per subject
#$ -l tmpfs=5G

# Set up the job array task number range
#$ -t 1-49

# Set the name of the job.
#$ -N hcp-tractography

# Set the working directory to somewhere in scratch space.
#$ -wd /home/zcqsfey/Scratch/output

#Local2Scratch

# ------------------------ Job Setup -------------------------

# Parse parameter file to get variables.
number=$SGE_TASK_ID
# TODO: check subject number list that 100208 removed
paramfile=/home/zcqsfey/Scratch/input/49.txt
subject=""$(sed -n ${number}p $paramfile )""

SUBJECT_DIR=/home/zcqsfey/Scratch/HCP/$subject
MNI_ROIS_DIR=/home/zcqsfey/Scratch/MNI/ROIs

# ------------------------ Job Run ---------------------------

cd $TMPDIR

mkdir -p ROIs dwi{1..5}/rois_and_tracks

stamp=$(date +"%m%d%H%M")


# Create mrtrix usable xform
mrconvert $SUBJECT_DIR/MNI/standard2acpc_dc.nii.gz tmp-[].nii
mv tmp-0.nii x.nii
mrcalc x.nii -neg tmp-0.nii -force
warpconvert tmp-[].nii displacement2deformation mni2acpc_mrtrix.nii.gz
rm x.nii tmp-?.nii
# Transform MNI ROIs to subject space in TMPDIR
for roi in $MNI_ROIS_DIR/*; do
  mrtransform -warp mni2acpc_mrtrix.nii.gz $roi -template $SUBJECT_DIR/t1.nii.gz - | \
    mrthreshold - -abs 0.01 ROIs/$(basename $roi)
done

## Tractography PARTY

for d in dwi{1,2}; do
for h in l r; do
tckgen $SUBJECT_DIR/$d/modelling/csd_wm.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_PLIC.mif \
      -include ROIs/"$h"_CST_PONS.mif -include ROIs/"$h"_CP.mif \
      -exclude ROIs/ML.mif -exclude ROIs/CPs.mif -exclude ROIs/SFOF.mif \
      -exclude ROIs/CST_post.mif -exclude ROIs/midline.mif \
      $d/rois_and_tracks/"$h"cst"$stamp".tck -quiet &

tckgen $SUBJECT_DIR/$d/modelling/csd_wm.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_LGN.mif -seed_unidirectional \
      -include ROIs/"$h"_SS.mif -exclude ROIs/OR_exclude.mif \
      $d/rois_and_tracks/"$h"or"$stamp".tck -quiet  &

tckgen $SUBJECT_DIR/$d/modelling/csd_wm.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_AF_cor.mif -include ROIs/"$h"_AF_ax.mif \
      -exclude ROIs/midline.mif -exclude ROIs/"$h"_CP.mif -exclude ROIs/SFOF.mif \
      -exclude ROIs/CR.mif -exclude ROIs/"$h"_SS.mif -exclude ROIs/EC.mif \
      $d/rois_and_tracks/"$h"af"$stamp".tck -quiet &
done
done

for d in dwi{3..5}; do
# MSMT 1
for h in l r; do
tckgen $SUBJECT_DIR/$d/modelling/msmt/csd_wm1.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_PLIC.mif \
      -include ROIs/"$h"_CST_PONS.mif -include ROIs/"$h"_CP.mif \
      -exclude ROIs/ML.mif -exclude ROIs/CPs.mif -exclude ROIs/SFOF.mif \
      -exclude ROIs/CST_post.mif -exclude ROIs/midline.mif \
      $d/rois_and_tracks/"$h"cst"$stamp"_msmt1.tck -quiet &

tckgen $SUBJECT_DIR/$d/modelling/msmt/csd_wm1.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_LGN.mif -seed_unidirectional \
      -include ROIs/"$h"_SS.mif -exclude ROIs/OR_exclude.mif \
      $d/rois_and_tracks/"$h"or"$stamp"_msmt1.tck -quiet  &

tckgen $SUBJECT_DIR/$d/modelling/msmt/csd_wm1.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_AF_cor.mif -include ROIs/"$h"_AF_ax.mif \
      -exclude ROIs/midline.mif -exclude ROIs/"$h"_CP.mif -exclude ROIs/SFOF.mif \
      -exclude ROIs/CR.mif -exclude ROIs/"$h"_SS.mif -exclude ROIs/EC.mif \
      $d/rois_and_tracks/"$h"af"$stamp"_msmt1.tck -quiet  &
done

# MSMT 2
for h in l r; do
tckgen $SUBJECT_DIR/$d/modelling/msmt/csd_wm2.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_PLIC.mif \
      -include ROIs/"$h"_CST_PONS.mif -include ROIs/"$h"_CP.mif \
      -exclude ROIs/ML.mif -exclude ROIs/CPs.mif -exclude ROIs/SFOF.mif \
      -exclude ROIs/CST_post.mif -exclude ROIs/midline.mif \
      $d/rois_and_tracks/"$h"cst"$stamp"_msmt2.tck -quiet &

tckgen $SUBJECT_DIR/$d/modelling/msmt/csd_wm2.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_LGN.mif -seed_unidirectional \
      -include ROIs/"$h"_SS.mif -exclude ROIs/OR_exclude.mif \
      $d/rois_and_tracks/"$h"or"$stamp"_msmt2.tck -quiet  &

tckgen $SUBJECT_DIR/$d/modelling/msmt/csd_wm2.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_AF_cor.mif -include ROIs/"$h"_AF_ax.mif \
      -exclude ROIs/midline.mif -exclude ROIs/"$h"_CP.mif -exclude ROIs/SFOF.mif \
      -exclude ROIs/CR.mif -exclude ROIs/"$h"_SS.mif -exclude ROIs/EC.mif \
      $d/rois_and_tracks/"$h"af"$stamp"_msmt2.tck -quiet  &
done

# SSST
for h in l r; do
tckgen $SUBJECT_DIR/$d/modelling/ssst/csd_wm.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_PLIC.mif \
      -include ROIs/"$h"_CST_PONS.mif -include ROIs/"$h"_CP.mif \
      -exclude ROIs/ML.mif -exclude ROIs/CPs.mif -exclude ROIs/SFOF.mif \
      -exclude ROIs/CST_post.mif -exclude ROIs/midline.mif \
      $d/rois_and_tracks/"$h"cst"$stamp"_ssst.tck -quiet &

tckgen $SUBJECT_DIR/$d/modelling/ssst/csd_wm.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_LGN.mif -seed_unidirectional \
      -include ROIs/"$h"_SS.mif -exclude ROIs/OR_exclude.mif \
      $d/rois_and_tracks/"$h"or"$stamp"_ssst.tck -quiet  &

tckgen $SUBJECT_DIR/$d/modelling/ssst/csd_wm.mif -mask $SUBJECT_DIR/mask.nii.gz \
      -seed_image ROIs/"$h"_AF_cor.mif -include ROIs/"$h"_AF_ax.mif \
      -exclude ROIs/midline.mif -exclude ROIs/"$h"_CP.mif -exclude ROIs/SFOF.mif \
      -exclude ROIs/CR.mif -exclude ROIs/"$h"_SS.mif -exclude ROIs/EC.mif \
      $d/rois_and_tracks/"$h"af"$stamp"_ssst.tck -quiet  &

done
done

## Tractmaps

wait # Finish all background tckgen processes first
for f in dwi?/rois_and_tracks/*.tck; do
  tckmap -template $SUBJECT_DIR/dwi1/dwi.mif $f ${f/.tck/.nii.gz}
done
