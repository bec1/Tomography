# Tomography

Put slices of atom clouds together for a volumetric display.

## Files:

* tomography.m

usage:

1. Use listfits to creat a list of your fits images:
	sideimages = listfits('side_images_directory');
	topimages = listfits('top_images_directory');

2. Apply tomography(sideimages,topimages) to get an isosurface plot and stuff

* Old files from may: use for smoothing and stuff


## Parameters to play with:

side_a = 170; % width of square crop
side_crop = [66,304,side_a,side_a];

top_w = 60;  % width of crop (radial)
top_l = 400; % length of crop (axial)
top_r = 0;  % degrees of rotation
top_crop = [225,62,top_w,top_l];

surface_level = 0.1; 