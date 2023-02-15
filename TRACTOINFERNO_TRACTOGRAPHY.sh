#!/bin/bash -l

# Wallclock time (format hours:minutes:seconds).
#$ -l h_rt=0:15:0

# Cores
#$ -pe smp 36

# RAM (must be an integer followed by M, G, or T)
#$ -l mem=16G

# TMPDIR space (default is 10 GB - remove if cluster is diskless)
# avg 5MB per .tck x 3 tracts x 2 hems x 144 sbj = ca. 4.32GB total, 30MB per subject. 6 tractograms per subject
#$ -l tmpfs=1G

# Set up the job array task number range
#$ -t 1-144

# Set the name of the job.
#$ -N tractoinferno-tractography

# Set the working directory to somewhere in scratch space.
#$ -wd /home/zcqsfey/Scratch/output

#Local2Scratch

# ------------------------ Job Setup -------------------------

# Parse parameter file to get variables.
number=$SGE_TASK_ID
paramfile=${1:-/home/zcqsfey/Scratch/input/tractoinferno.txt}
subject=sub-""$(sed -n ${number}p $paramfile )""

SUBJECT_DIR=/home/zcqsfey/Scratch/tractoinferno/$subject
MNI_DIR=/home/zcqsfey/Scratch/MNI
MNI_ROIS_DIR=/home/zcqsfey/Scratch/MNI/ROIs

# ------------------------ Job Run ---------------------------

cd $TMPDIR

module load ants

mkdir -p reg rois tractography-ich

stamp=$(date +"%m%d%H%M")


## ANTs registration

# Nonlinear registration
$ANTSPATH/antsRegistrationSyNQuick.sh \
	-d 3 \
	-n 4 \
	-m $MNI_DIR/MNI152.nii.gz \
	-f $SUBJECT_DIR/anat/${subject}__T1w.nii.gz \
	-o reg/${subject}_transform-MNI152-to-anat_ants_

# Transform ROIs
for roi in $MNI_ROIS_DIR/*.nii.gz; do
	${ANTSPATH}antsApplyTransforms \
	-d 3 \
	-i $roi \
	-r $SUBJECT_DIR/anat/${subject}__T1w.nii.gz \
	-t reg/${subject}_transform-MNI152-to-anat_ants_0GenericAffine.mat \
	-t reg/${subject}_transform-MNI152-to-anat_ants_1Warp.nii.gz \
	-n NearestNeighbor \
	-o rois/${roi##*/}
done

## Tractography PARTY

for h in l r; do

tckgen $SUBJECT_DIR/csd/csd_wm1.mif -mask $SUBJECT_DIR/mask/mask.nii.gz \
      -seed_image rois/"$h"_PLIC.nii.gz \
      -include rois/"$h"_CST_PONS.nii.gz -include rois/"$h"_CP.nii.gz \
      -exclude rois/ML.nii.gz -exclude rois/CPs.nii.gz -exclude rois/SFOF.nii.gz \
      -exclude rois/CST_post.nii.gz -exclude rois/midline.nii.gz \
      tractography-ich/"$h"cst"$stamp".tck -quiet &

tckgen $SUBJECT_DIR/csd/csd_wm1.mif -mask $SUBJECT_DIR/mask/mask.nii.gz \
      -seed_image rois/"$h"_LGN.nii.gz -seed_unidirectional \
      -include rois/"$h"_SS.nii.gz -exclude rois/OR_exclude.nii.gz \
      tractography-ich/"$h"or"$stamp".tck -quiet  &


tckgen $SUBJECT_DIR/csd/csd_wm1.mif -mask $SUBJECT_DIR/mask/mask.nii.gz \
      -seed_image rois/"$h"_AF_cor.nii.gz -include rois/"$h"_AF_ax.nii.gz \
      -exclude rois/midline.nii.gz -exclude rois/"$h"_CP.nii.gz -exclude rois/SFOF.nii.gz \
      -exclude rois/CR.nii.gz -exclude rois/"$h"_SS.nii.gz -exclude rois/EC.nii.gz \
      tractography-ich/"$h"af"$stamp".tck -quiet &
done


## Tractmaps

wait # Finish all background tckgen processes first
for f in tractography-ich/*.tck; do
  tckmap -template $SUBJECT_DIR/mask/mask.nii.gz $f ${f/.tck/.nii.gz}
done
