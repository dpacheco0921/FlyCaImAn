# CaImProPi-MATLAB
Calcium Imaging Processing Pipeline of 3D time series

# Organization of data

# Pipeline
To process either 3DxT or 2DxT datasets:
- 1) Convert tiffs to mat files (generate image data variable 'Data', and image metadata variable 'iDat').
- 2) Pull extra information from metadata files ('.bin', '_vDat.mat', or '.mat' (LEDcontroler)).
- 3) Do motion correction (using [NoRMCorre](https://github.com/flatironinstitute/NoRMCorre)).
- 4) Do spatial and/or temporal resampling
    - this includes re-slicing for volumetric datasets and aligment relative to stimuli delivery.
    - generates main metadata variable used for ROI segmentation (wDat).
- 5) select brain pixels (generate a binary mask).
- 6) i) format stacks for ROI segmentation.
- 6) ii) stitch (along z axis) and format stacks for ROI segmentation.
- 7) ROI segmentation (in progress)
- 8) Detect stimulus-modulated ROIs (in progress)

To register image segments to local whole brain and to in vivo atlas
- (in progress)

# Dependencies

For processing of imaging data this pipeline requires the following packages:
- [CaImAn](https://github.com/flatironinstitute/CaImAn-MATLAB), see link for dependencies.
- [NoRMCorre](https://github.com/flatironinstitute/NoRMCorre), see link for dependencies.
- [CMTK_matlab_wrapper](https://github.com/dpacheco0921/CMTK_matlab_wrapper), this requires the Computational Morphometry Toolkit [CMTK](https://www.nitrc.org/projects/cmtk)

For processing of behavior videos with fictrac this pipeline requires the following packages:
- [Fictrac](http://rjdmoore.net/fictrac/), see link for installation & dependencies (for windows see [fic-trac-win](https://github.com/murthylab/fic-trac-win)).
- generate calibration-transform.dat, and update that file in /toolbox/fictrac_offline/fictrac_settings

For saving figures install:
- [export_fig](https://github.com/altmany/export_fig)

For stitching Z-stacks:
- [fiji](https://imagej.net/Fiji/Downloads), make sure you have the [stitching plugin](https://imagej.net/Image_Stitching)

# Usage

To use this package
1) Copy and edit fiji_fullpathedit.m (see that file for details), and save as fiji_fullpath.m.

For examples see all the demos: FlyCaImAn_*_demo.m

# Acknowledgements

Special thanks to:
- [Eftychios Pnevmatikakis](https://github.com/epnev) and [Andrea A. Giovannuci](https://github.com/agiovann) for help with [CaImAn-MATLAB](https://github.com/flatironinstitute/CaImAn-MATLAB) and [NoRMCorre](https://github.com/flatironinstitute/NoRMCorre) toolboxes
- [Gregory Jefferis](https://github.com/jefferis) and Torsten Rohlfing for help with [CMTK toolbox](https://www.nitrc.org/projects/cmtk)

# Citation

If you use this code please cite the following corresponding paper:
[Diego Pacheco, Stephan Thiberge, Eftychios Pnevmatikakis, Mala Murthy (2019). Auditory Activity is Diverse and Widespread Throughout the Central Brain of Drosophila](https://doi.org/10.1101/709519)
