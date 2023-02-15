# shell processing script for running tractfinder and TractSeg on preprocessed Tractoinferno data (v1.1.1).
#
# Fiona Young
# University College London
# December 2022
#
# This script is not intended to function out-of-the-box
# on a new system and may need adapting
#
# This script utilises MRtrix (v3.0.3) and FSL (6.0.5:9e026117) software packages:
# https://mrtrix.readthedocs.io/en/3.0.3/
# https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/
#
# See this article for information on the data dataset: https://www.nature.com/articles/s41597-022-01833-1
# Data downloaded from OpenNeuro: doi:10.18112/openneuro.ds003900.v1.1.1 (https://openneuro.org/datasets/ds003900/versions/1.1.1)

IP() {
  mrcalc $1 $2 -mult -|mrmath - sum -axis 3 $3 -force
}

if [ -z ${DATAROOT+x} ];
  then echo "DATAROOT is unset"; exit
fi

while read -r s; do

dir=$DATAROOT/images_and_data/tractoinferno/sub-${s}
T1=$dir/anat/sub-${s}__T1w.nii.gz
dwi=$dir/dwi/sub-${s}__dwi

mkdir $dir/{reg,csd,forage,tractseg}

## --- Convert to .mif and downsample

mrconvert $dwi.nii.gz -fslgrad $dwi.bvec $dwi.bval \
  -strides -1,2,3,4 - | \
  mrgrid - regrid $dwi.mif.gz -voxel 2.3

## --- Create brain mask

mask=$dir/mask/sub-${s}__mask_brain.nii.gz
#dwi2mask $dwi.mif.gz $mask
mrcalc $dir/dti/sub-${s}__md.nii.gz 0 -gt - |\
 mrgrid - regrid $mask -template $dwi.mif.gz -interp nearest -datatype uint8 -force
ln -sf $mask $dir/mask/mask.nii.gz

## --- Register MNI and T1

flirt -in $DATAROOT/images_and_data/MNI152.nii.gz -ref $T1 -dof 9 \
  -omat $dir/reg/sub-${s}_transform-MNI152-to-anat_flirt.txt
transformconvert $dir/reg/sub-${s}_transform-MNI152-to-anat_flirt.txt \
  $DATAROOT/images_and_data/MNI152.nii.gz $T1 flirt_import \
  $dir/reg/sub-${s}_transform-MNI152-to-anat_mrtrix.txt -force

## --- CSD

dwi2response dhollander $dwi.mif.gz -mask $mask \
  $dir/csd/response_sfwm.txt $dir/csd/response_gm.txt $dir/csd/response_csf.txt -quiet -force
sed '3d' $dir/csd/response_sfwm.txt > $dir/csd/response_sfwm_ssst.txt

# MSMT WM + GM
dwi2fod msmt_csd -mask $mask $dwi.mif.gz \
  $dir/csd/response_sfwm.txt $dir/csd/sub-${s}_model-msmtcsd_comp-wm_desc-1.mif \
  $dir/csd/response_gm.txt $dir/csd/sub-${s}_model-msmtcsd_comp-gm.mif -force
ln -sf $dir/csd/sub-${s}_model-msmtcsd_comp-wm_desc-1.mif $dir/csd/csd_wm1.mif

# SSST CSD
dwi2fod csd -mask $mask $dwi.mif.gz \
  $dir/csd/response_sfwm_ssst.txt $dir/csd/sub-${s}_model-csd_fod.mif -force
ln -sf $dir/csd/sub-${s}_model-csd_fod.mif $dir/csd/csd_wm.mif

# Extract peaks for TractSeg
sh2peaks $dir/csd/csd_wm.mif $dir/csd/tmp.nii.gz -force
fslreorient2std $dir/csd/tmp.nii.gz $dir/csd/sub-${s}_peaks.nii.gz && \
  rm $dir/csd/tmp.nii.gz
ln -sf $dir/csd/sub-${s}_peaks.nii.gz $dir/csd/peaks.nii.gz

## --- Forage

for t in cst or af; do
  for h in l r; do

  # Transform atlas
  mrtransform -linear $dir/reg/sub-${s}_transform-MNI152-to-anat_mrtrix.txt \
    $DATAROOT/images_and_data/mni${t}/${h}*_filtered_tod_atlas.mif \
    $dir/forage/sub-${s}_atlas-${h}${t}tod.mif.gz \
    -template $dwi.mif.gz -reorient_fod yes
  ln -sf $dir/forage/sub-${s}_atlas-${h}${t}tod.mif.gz $dir/forage/${h}_${t}_tod_atlas.mif.gz

  # Compute inner product with FOD
  IP $dir/csd/csd_wm1.mif $dir/forage/${h}_${t}_tod_atlas.mif.gz $dir/forage/sub-${s}__${h}${t}.nii.gz
  ln -sf $dir/forage/sub-${s}__${h}${t}.nii.gz $dir/forage/${h}_${t}_tractmap.nii.gz

  done
done

# TractSeg
tractseg -i $dir/csd/peaks.nii.gz -o $dir/tractseg \
  --get_probabilities --nr_cpus 4 --tract_definition TractQuerier+

tractseg -i $dir/csd/peaks.nii.gz -o $dir/tractseg \
  --get_probabilities --nr_cpus 4 --tract_definition xtract

done < $1
