camera=$1

cd cmake-build-debug

echo "Correcting images..."
if [ -d "../raw_images/$camera" ]
then
	imageFiles=$(ls "../raw_images/$camera" -1 | grep ".ppm")
	for image in $imageFiles;
	do
		img_name=$(basename "$image" ".ppm")
		./sensor_correction "../raw_images/$camera/$image"
		mv "../corrected_xyz.ppm" "../${img_name}_xyz.ppm"		
		mv "../corrected_srgb.ppm" "../${img_name}_srgb.ppm"
		mv "../corrected_rec2020.ppm" "../${img_name}_rec2020.ppm"
	done;
else
	echo "No images to correct for $camera!"
fi
