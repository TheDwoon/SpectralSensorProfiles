#!/bin/bash

spectral_images=("resources/woman_face.spimg")

if [ $# -lt 1 ]
then
	cameras=./cameras/*.dat
else
	cameras=$1
fi

IFS=$'\n'
for camFile in $cameras;
do
	# generate some paths and names
	camera=$(basename "$camFile" ".dat")
	camFitDir="./cameras/$camera.fit"

	# move into build folder and execute solver for current camera
	cd cmake-build-debug


	# ******************************************************************
	# Solving CFA
	# ******************************************************************
	echo "Running solver for $camera"
	cat "../$camFile" | ./sensor_spectra


	# ******************************************************************
	# Rendering spectral images for camera
	# ******************************************************************
	echo "Rendering spectral images"
	for spectral_image in ${spectral_images[@]};
	do
		./sensor_images "../$spectral_image"
		img_name=$(basename "$spectral_image" ".spimg")
		
		mv "../spectral_camrgb.ppm" "../${img_name}_camrgb.ppm"
		mv "../spectral_xyz.ppm" "../${img_name}_xyz.ppm"
		mv "../spectral_srgb.ppm" "../${img_name}_srgb.ppm"
		mv "../spectral_rec2020.ppm" "../${img_name}_rec2020.ppm"
	done;


	# ******************************************************************
	# Color correction if images are present
	# ******************************************************************
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

	# ******************************************************************
	# Plotting things. Move results into designated folder
	# ******************************************************************

	cd ..

	# Generate plots and delete/create folder to contain fit information
	[ -d "$camFitDir" ] && rm -rf "$camFitDir"
	mkdir "$camFitDir"

	gnuplot cfa_plots.gp

	# move generated ppm images (just dump them all)
	mv *.ppm "$camFitDir"

	# move generated plots
	mv "cfa_fitted.png" "$camFitDir"
	mv "cfa_cost.png" "$camFitDir"
	mv "cfa_iterations.png" "$camFitDir"
	mv "cfa_errors.png" "$camFitDir"
	mv "cfa_sampling.png" "$camFitDir"
	mv "cfa_inverse.png" "$camFitDir"

	# move data files
	mv "cfa_fitting.dat" "$camFitDir"
	mv "cfa_sampling.dat" "$camFitDir"
	mv "cfa_spectra.dat" "$camFitDir"
	mv "cfa_vector.dat" "$camFitDir"
	mv "inverse.dat" "$camFitDir"
done;
