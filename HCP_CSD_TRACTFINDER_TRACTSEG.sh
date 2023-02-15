#!/bin/bash -l

# Wallclock time (format hours:minutes:seconds).
#$ -l h_rt=1:30:0

# Cores
#$ -pe smp 36

# RAM (must be an integer followed by M, G, or T)
#$ -l mem=32G

# TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=15G

# Set up the job array task number range
#$ -t 1-18

# Set the name of the job.
#$ -N hcp-forage-tractseg

# Where to send output log files
#$ -wd /home/zcqsfey/Scratch/output

#Local2Scratch

# ------------------------ Job Setup -------------------------

# Parse parameter file to get variables. Default list of 49
# subject numbers is used unless alternative paramfile passed
number=$SGE_TASK_ID
paramfile=${1:-/home/zcqsfey/Scratch/input/49.txt}
subject=""$(sed -n ${number}p $paramfile )""

# Set input files
mnidir=/home/zcqsfey/Scratch/mni
subjectdir=/home/zcqsfey/Scratch/HCP/$subject

# MNI t1 reference
mni=$mnidir/mni152_t1.nii.gz

# Atlases
cstatlas=$mnidir/cst_atlas.mif
oratlas=$mnidir/or_atlas.mif
afatlas=$mnidir/af_atlas.mif

# Subject specific files
t1=$subjectdir/t1.nii.gz
mask=$subjectdir/mask.nii.gz

# Define inner product function
IP() {
  mrcalc $1 $2 -mult -|mrmath - sum -axis 3 $3
}

# ------------------------ Job Run ---------------------------

cd $TMPDIR

# MNI registration
flirt -in $mni -ref $t1 -dof 9 -omat t_mni_2_t1_fsl.txt -out mnireg.nii.gz
transformconvert t_mni_2_t1_fsl.txt $mni $t1 flirt_import t_mni_2_t1_mr.txt

# NB: -reoriend_fod removed because myriad has MRtrix v.3.0_RC3
mrtransform -linear t_mni_2_t1_mr.txt -template $subjectdir/dwi1/dwi.mif $cstatlas cst_atlas.mif.gz
mrtransform -linear t_mni_2_t1_mr.txt -template $subjectdir/dwi1/dwi.mif $oratlas or_atlas.mif.gz
mrtransform -linear t_mni_2_t1_mr.txt -template $subjectdir/dwi1/dwi.mif $afatlas af_atlas.mif.gz

#----------------------#
# Multi shell datasets #
#----------------------#

for dir in dwi1 dwi2; do

  mkdir $TMPDIR/$dir
  cd $TMPDIR/$dir

  mkdir modelling

  # Do MSMT CSD
  dwi2response dhollander -mask $mask $subjectdir/$dir/dwi.mif modelling/response_sfwm.txt modelling/response_gm.txt modelling/response_csf.txt
  dwi2fod msmt_csd -mask $mask $subjectdir/$dir/dwi.mif modelling/response_sfwm.txt modelling/csd_wm.mif modelling/response_gm.txt modelling/csd_gm.mif modelling/response_csf.txt modelling/csd_csf.mif

  # IP with (transformed) atlas
  IP $TMPDIR/cst_atlas.mif.gz  modelling/csd_wm.mif tractmap_cst.mif.gz
  IP $TMPDIR/or_atlas.mif.gz 	modelling/csd_wm.mif tractmap_or.mif.gz
  IP $TMPDIR/af_atlas.mif.gz 	modelling/csd_wm.mif tractmap_af.mif.gz

  ## TractSeg

  # Convert SH FODs to peaks
  sh2peaks modelling/csd_wm.mif modelling/peaks.nii.gz
  # Run for DKFZ (TractQuerier+) tract definitions
  TractSeg -i modelling/peaks.nii.gz -o tractseg_dkfz --tract_definition TractQuerier+ --get_probabilities --nr_cpus 20
  # Run for XTract tract definitions
  TractSeg -i modelling/peaks.nii.gz -o tractseg_xtract --tract_definition xtract --get_probabilities --nr_cpus 20

done

#-----------------------#
# Single shell datasets #
#-----------------------#

for dir in dwi{3..5}; do

  mkdir $TMPDIR/$dir
  cd $TMPDIR/$dir


  mkdir -p modelling/{msmt,ssst}

  cd $TMPDIR/$dir/modelling/msmt

  ## Multi shell, multi (=2) tissue

  # Response function
  dwi2response dhollander -mask $mask $subjectdir/$dir/dwi.mif response_sfwm.txt response_gm.txt response_csf.txt
  # WM and GM
  dwi2fod msmt_csd -mask $mask $subjectdir/$dir/dwi.mif response_sfwm.txt csd_wm1.mif response_gm.txt csd_gm.mif
  # WM and CSF
  dwi2fod msmt_csd -mask $mask $subjectdir/$dir/dwi.mif response_sfwm.txt csd_wm2.mif response_csf.txt csd_csf.mif

  ## Single Tissue
  cd $TMPDIR/$dir/modelling/ssst

  dwi2response tournier -mask $mask $subjectdir/$dir/dwi.mif response_sfwm.txt
  dwi2fod csd -mask $mask $subjectdir/$dir/dwi.mif response_sfwm.txt csd_wm.mif

  cd $TMPDIR/$dir

  ## IP with (transformed) atlas

  IP $TMPDIR/cst_atlas.mif.gz modelling/msmt/csd_wm1.mif tractmap_cst_msmt1.mif.gz
  IP $TMPDIR/cst_atlas.mif.gz modelling/msmt/csd_wm2.mif tractmap_cst_msmt2.mif.gz
  IP $TMPDIR/cst_atlas.mif.gz modelling/ssst/csd_wm.mif tractmap_cst_ssst.mif.gz

  IP $TMPDIR/or_atlas.mif.gz modelling/msmt/csd_wm1.mif tractmap_or_msmt1.mif.gz
  IP $TMPDIR/or_atlas.mif.gz modelling/msmt/csd_wm2.mif tractmap_or_msmt2.mif.gz
  IP $TMPDIR/or_atlas.mif.gz modelling/ssst/csd_wm.mif tractmap_or_ssst.mif.gz

  IP $TMPDIR/af_atlas.mif.gz  modelling/msmt/csd_wm1.mif tractmap_af_msmt1.mif.gz
  IP $TMPDIR/af_atlas.mif.gz  modelling/msmt/csd_wm2.mif tractmap_af_msmt2.mif.gz
  IP $TMPDIR/af_atlas.mif.gz  modelling/ssst/csd_wm.mif tractmap_af_ssst.mif.gz

  ## TractSeg

  # Convert SH FODs to peaks
  sh2peaks modelling/msmt/csd_wm1.mif modelling/msmt/peaks1.nii.gz
  sh2peaks modelling/msmt/csd_wm2.mif modelling/msmt/peaks2.nii.gz
  sh2peaks modelling/ssst/csd_wm.mif modelling/ssst/peaks.nii.gz
  # Run for DKFZ (TractQuerier+) tract definitions
  TractSeg -i modelling/msmt/peaks1.nii.gz -o tractseg_dkfz_msmt1 --tract_definition TractQuerier+ --get_probabilities --nr_cpus 20
  TractSeg -i modelling/msmt/peaks2.nii.gz -o tractseg_dkfz_msmt2 --tract_definition TractQuerier+ --get_probabilities --nr_cpus 20
  TractSeg -i modelling/ssst/peaks.nii.gz -o tractseg_dkfz_ssst --tract_definition TractQuerier+ --get_probabilities --nr_cpus 20
  # Run for XTract tract definitions
  TractSeg -i modelling/msmt/peaks1.nii.gz -o tractseg_xtract_msmt1 --tract_definition xtract --get_probabilities --nr_cpus 20
  TractSeg -i modelling/msmt/peaks2.nii.gz -o tractseg_xtract_msmt2 --tract_definition xtract --get_probabilities --nr_cpus 20
  TractSeg -i modelling/ssst/peaks.nii.gz -o tractseg_xtract_ssst --tract_definition xtract --get_probabilities --nr_cpus 20

done
