## NOTE: this is an old script for creating a bilateral TOD atlas.
# Newer analyses use left/right split atlases

stamp=$(date +"%m%d%H%M")
# Set matching pattern of tck files to manipulate
# NOTE: this is kinda icky? find better way of using wildcards
pre1=$1
pre2=$2
tract=$3

root=$DATAROOT/images_and_data

for s in 32 {35..40} {42..50}
do

	dir=$root/$s
	cd $dir

	# Combine left and right versions of both with and without ML waypoint
	tckedit rois_and_tracks/*"$tract""$pre1".tck rois_and_tracks/*"$tract""$pre2".tck tmp.tck &&\
	#tckedit rois_and_tracks/*"$tract""$pre1".tck tmp.tck &&\

	# Filter spurious streamlines with a 2mm track density map (min 3 streamlines per (2mm)^3)
	tckmap tmp.tck -vox 2 -| mrcalc - 3 -lt -| tckedit tmp.tck -exclude - rois_and_tracks/b"$tract""$stamp".tck &&\
	rm tmp.tck

	# Transform the bilateral tractogram into mni space
	tcktransform rois_and_tracks/b"$tract""$stamp".tck $dir/wi_transf.mif rois_and_tracks/b"$tract""$stamp"_mni.tck

	# Map to TOD
	tckmap -tod 8 -template $root/mni/mni152v2009a/mni_icbm152_t1_tal_nlin_sym_09a.nii rois_and_tracks/b"$tract""$stamp"_mni.tck maps/b"$tract""$stamp"_mni_tod.mif

	# Normalise and remove NaNs
	mrconvert -coord 3 0 $dir/maps/b"$tract""$stamp"_mni_tod.mif - |\
	mrcalc $dir/maps/b"$tract""$stamp"_mni_tod.mif - 4 pi -mult -sqrt -mult -div - |\
	mrcalc - -isnan 0 - -if $dir/maps/b"$tract""$stamp"_mni_tod_norm.mif

done
cd $root

# Create average TOD map
mrmath $root/*/maps/b"$tract""$stamp"_mni_tod_norm.mif mean $root/mni"$tract"/b"$tract""$stamp"_todmean.mif
