# CaImProPi-MATLAB
Calcium Imaging Processing Pipeline of 3D time series

# Organization of data

# Pipeline
To process either 3DxT or 2DxT datasets:
- 1) Convert tiffs to mat files (generate image data variable 'Data', and image metadata variable 'iDat')
- 2) Pull extra information from metadata files ('.bin', '_vDat.mat', or '.mat' (LEDcontroler))
- 3) Do motion correction
- 4) Optional (do spatial and/or temporal resampling (this includes re-slicing for volumetric datasets and aligment relative to stimuli delivery))
- 5) Optional (stitch imaged volumes along the z-axes)
- 6) ROI segmentation
- 7) Detect stimulus-modulated ROIs

To register image segments to local whole brain and to in vivo atlas
- **

# Dependencies

This pipeline requires the following packages:
- [CaImAn](https://github.com/flatironinstitute/CaImAn-MATLAB), see link for dependencies.
- [NoRMCorre](https://github.com/flatironinstitute/NoRMCorre), see link for dependencies.
- [CMTK_matlab_wrapper](https://github.com/dpacheco0921/CMTK_matlab_wrapper), this requires the Computational Morphometry Toolkit [CMTK](https://www.nitrc.org/projects/cmtk)

# Acknowledgements

Special thanks to:
- [Eftychios Pnevmatikakis](https://github.com/epnev) and [Andrea A. Giovannuci](https://github.com/agiovann) for help with CaImAn and NoRMCorre toolboxes
- [Gregory Jefferis](https://github.com/jefferis) and Torsten Rohlfing for help with [CMTK toolbox](https://www.nitrc.org/projects/cmtk)

# Citation

If you use this code please cite the following corresponding paper:
[Diego Pacheco, Stephan Thiberge, Eftychios Pnevmatikakis, Mala Murthy (2019). Auditory Activity is Diverse and Widespread Throughout the Central Brain of Drosophila](https://doi.org/10.1101/709519)