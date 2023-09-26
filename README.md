# Spectral and continuum stacking codes 

## Overview
The spectral and continuum stacking codes are two standalone Python scripts for performing spectral stacking (in 3D cubes) and continuum stacking (images) of ALMA datasets, respectively. 
The codes are presented and described in Guerrero et al. submitted. 

Each code can be used separately within CASA (Common Astronomy Software Applications).

Author: Andrea Guerrero

## Installation 
### Before running the codes
The stacking codes need be executed within CASA. To install CASA, please read the [CASA documentation](https://casa.nrao.edu/casa_obtaining.shtml). In particular, the codes were tested in CASA version 5.6.

Before using the codes, the python modules _astropy_ and _lineid_plot_ must be installed within CASA. These modules can be installed following the [astropy documentation](https://docs.astropy.org/en/stable/install.html\#installing-astropy-into-casa). First, open CASA in your console and install _pip_ using:

> CASA <1> : from setuptools.command import easy_install

> CASA <2> : easy_install.main(['--user', 'pip'])

Quit CASA and re-open it. Now to install _astropy_ and _lineid_plot_ type,

> CASA <1>: import subprocess, sys

> CASA <2>: subprocess.check_call([sys.executable, '-m', 'pip','install','--user', 'astropy', 'lineid_plot'])

Quit CASA and re-open it. Now it should be possible to import either module using, 

> CASA <1> : import astropy

Both codes are now ready to be run within CASA.

### Examples and Tests
Before running the codes with your own data, we encourage you to run the provided examples to make sure the codes are working properly. 
First, download the folder(s) desired (continuum and/or spectral). Inside each folder, you will find the stacking code and a folder named: _stack_example_.
The _stack_example_ folder contains a fits file and a mock catalog. For the spectral code, the fits file is a 3D cube and the catalog is a csv file with three columns: ra (deg), dec (deg) and redshift. 
For the continuum code, the fits file is a continuum image and the catalog is a csv file with two columns: ra(deg) and dec (deg). Each code has its inputs set to be run with the examples.

To run, open your console, go to the directory in which the stacking code is located and run CASA. Then use execfile to run the code desired in the CASA console. For instance:
> CASA <1> : execfile('continuum_stacking.py')

The stacking results will be saved in the _stack_example_ folder. The results include intermediate files and plots. Details on each code can be found further below. 

## Usage 

<details>
  <summary> Details for spectral stacking code </summary>
   
## Spectral stacking code 
The spectral stacking code reads one or more 3D spectral cubes (with axes ra, dec, frequency) and a csv (comma separated value) catalog of galaxy positions (ra, dec, redshift). It extracts the observed spectra of the positions in the catalog from the cubes when available, combines the rest-frame spectra for all sources, and returns the rest-frame stacked spectrum (1D array; mean and median) and the rest-frame frequency (1D array). 

### Inputs 
At the end of each code you will find the following inputs: 

| Name       | Type  | Default       | Description|
|------------|-------|---------------|------------| 
|dirname     |string |               | Name of output folder in which stacking results will be saved. If this folder is not in the same folder as spectral_stacking.py, the full path must be provided. |
|catalog_name|string | Empty string  | Name of the catalog file. The catalog file must be a csv (comma separated values) file with no header and with three required columns: RA, DEC and redshift. RA and DEC coordinates must be J2000 coordinates in degrees. An optional fourth column can be used as a weight. If the file is not in the same folder as spectral_stacking.py, the full path must be provided. In the example, the catalog_name is given by the "incat" variable.|
|stampsize   | int  | 3''           | Size in arcseconds [''] over which each object spectrum is extracted.|
| cubenames  | list  |  Empty list   | List of input 3D cubes that will be used in the stack. Each requires to be a fits file. If the fits file is not in the same folder as spectral_stacking.py, the full path must be provided.  |
|weightstack | bool  |   False       | If True, stacking is done using weighting, i.e. a galaxy can contribute in a different percentage to the combined spectra. If set to True, the fourth column of catalog_name is used as the weight for stacking. If False, the fourth column of catalog_name is ignored and stacking is done without weighting (each catalog object contribute the same amount to the final stack).  | 
| overwrite  | bool  |  False        | On executing the code, the extracted spectra of all target in a given input cube are saved to a file. If overwrite is False, the code will only extract and save these files for new cubes added to cubenames (in case of running the code again with more cubes). If overwrite is True, the code will re-extract and save spectra for all input cubes in cubenames. Note that this is meant to save time in case the code was previously run, the catalog remains the same, but only a few new input cubes are now added to the stack. If the catalog has been updated, then this parameter should be set to False.           | 
|  verbose   | int   |    1          | This value defines the amount of information printed to the message window  when running the code. Verbose values can be 0, 1, 2 or 3, where 0 provides minimal information and 3 provides the most detailed messages. | 


### Outputs 

#### Main results 
The main results of the stack are text files containing the stacked spectrum (mean and/or median), rest frequency, and additional information such as the number of objects which entered the stack, and the standard deviation of the stacked spectrum at each frequency. These files are: 

| File type  | File name                |  Description                           |
|------------|--------------------------|-------------------------------------------------------|
| Text       | stacked-meanspec.txt     | Text file of array with mean stacked spectrum.   |
| Text       | stacked-medianspec.txt   | Text file of array with median stacked spectrum. |
| Text       | stacked-frequency.txt    | Text file of array with rest-frame frequency of each value in the stacked spectra. |
| Text       | stacked-nobs.txt         | Text file of array with the number of objects stacked in each frequency channel. |
| Text       | stacked-stdspec.txt      | Text file of array with the standard deviation of the stacked spectrum in each frequency channel. |

#### Intermediate files
Several intermediate files are generated when the code is executed. These files are generated per cube (or cubename) used in the stack. Some correspond to coordinate files with ra (rad), dec (rad), redshift, weight, ra (deg) and dec (deg) for sources. These files are:

| File type  | File name               |  Description           |
|------------|-------------------------|------------------------|
| Text       | dirname+cubename+extractspec.txt  | Text file with the extracted observed-frame spectra of each source in this input cubename. |
| Text       | dirname+cubename+freqdef.txt      | Text file with the observed frequency of each spectrum in the above file (i.e. for each source extracted in cubename. |
| Text       | dirname+cubename+galaxyinfo.txt   | Coordinate file for all sources found in cubename.    |
| Text       | dirname+cubename+notfound.txt     | Coordinate file for all sources not found in cubename.|
| Text       | dirname+cubename+speccube.txt     | Text file with the rest-frame spectra of each source in this input cubename.|
| Text       | dirname+cubename+spec_corr.txt    | Text file with the rest-frame spectra of each source in this input cubename, but corrected for continuum contribution.|
| Text       | dirname+cubename+restfreqcube.txt | Text file with the rest-frame frequency of each spectrum in the above two files (i.e. rest-frame spectra of each source in this input cubename.  |


#### Figures
The stacking code provides standard figures to visualize the stacking results. If you do not wish to generate a specific figure, simply comment the line in which the figure is generated (which is right below the name of the figure). Some more detailed figures are currently commented; these can be uncommented if you are interested in the figure. 

All figures are generated in the function plot_spectra, which is executed at the end of the stacking code. The stacking and the plotting can thus be done separately by commenting one or the other. 
 
At the beginning of the plot_spectra function the user can decide a unit to show the results in mJy, Jy or μJy (unit = 'mJy', or 'Jy', 'muJy'). This can also be changed in a specific figure (in the unit argument). The table below describes the figures avaliable to plot with the name suggested by the code, however, this name can be changed as pleased (but be aware of the matplotlib bug that won't create a figure of file name is too long). Also, all figures are .png files for size purposes, though other formats like pdf can also be used.  


| Figure type                       |      File name                                               |  Description                               |
|-----------------------------------|--------------------------------------------------------------|--------------------------------------------|
| Full stack                        | full_stacked_spectra_mean.png                                | Two panel figure with stacked spectra (mean). Top panel shows the flux per rest-frame frequency (GHz) in blue, solid red line shows zero level flux and red dashed line shows the mean standard deviation of the spectrum in red (value shown on legend). Bottom panel shows the number of objects contributing to stack per frequency. If a line is within the frequency range, a dashed black line is position if line frequency with the name of the line at the top. |
|                                   | full_stacked_spectra_median.png                              | Same as full_stacked_spectra_mean.png but for median.              |
| For every galaxy in stack         | dirname+cubename+id+number+observed_spec.png                 | Observed sprectrum per frequency for galaxy that went into the stack. The values for id number is the row number from coordinate file (starting at zero). Includes same information as figures above, but it also provides cubename, RA (deg), Dec (deg) and redshift information.  |
|                                   | dirname+cubename+id+number+restframe_spec.png                | Same as dirname+cubename+id+number+observed_spec.png, but for rest-frame spectrum.        |
|                                   | dirname+cubename+id+number+rest-frame-continuum-corrected.png| Same as dirname+cubename+id+number+restframe_spec.png, but with rest-frame corrected by continuum contribution.  |
| Specific range in stacked spectra | stack_mean_+step+GHz+number.png                              | Figures of stacked spectra (mean) every specified GHz range. For example, figures every 50 GHz. User can also define first and last frequency to include in plot. |
|                                   | stack_median_+step+GHz+number.png                            | Same as stack_mean_+step+GHz+number.png but for median.  |
| Every common emission line        | line_stack_mean+linename.png                                 | Stacked spectra (mean) but centered on the most common emission lines (defined in lines_label and lines_freq, more lines can be added/discarded from these variables).   |
|                                   | line_stack_median+linename.png                               | Same as line_stack_mean+linename.png, but for median. |
| Gaussian fit on line              | stacked_mean_fit_line_+linename.png                          | Similar for line figures, but a Gaussian fit is performed using _curve_fit_ from scipy. The name of the line must be provided as input in the function speclinesnr. The legend of the figure also provides information of peak flux, total flux, FWHM, rms and SNR. |
|                                   | stacked_median_fit_line_+linename.png                        | Same as stacked_mean_fit_line_+linename.png, but for median.     |

</details>


<details>
  <summary> Details for continuum stacking code </summary>
   
## Continuum stacking code 
The continuum stacking code reads 2D images of continuum emission (e.g. moment 0 images) and a catalog of galaxy positions (ra, dec), extracts a region of the image for the galaxies in the catalog, combines all images extracted and return a stacked image (2d image) for mean and median. It also returns intermediate results for all images, all detected images and all undetected images.
Note that the stamp extracted from the image depends on the CDELT of the image (degrees per pixel information from the header), therefore, beware that although this code has the option of using images with different CDELTS values (by resizing the maps), this code was mainly tested using images from the same observation (i.e same CDELT).
 
### Inputs 
Towards the end of the code you will find the following inputs: 
 
| Name       | Type  | Default       | Description|
|------------|-------|---------------|------------| 
|dirname     |string |               | Name of output folder in which stacking results will be saved. If this folder is not in the same folder as continuum_stacking.py, full path must be provided. |
|catalog_name|string | Empty string  | Name of catalog file. Catalog file must be csv file with no header and two required columns: RA and DEC. RA and DEC coordinates must be J2000 coordinates in degrees. An optional third column can be used as a weight. If the file is not in the same folder as continuum_stacking.py, full path must be provided. In the example, its value is given by "incat" variable.|
|stampsize   | list  | 10''           | Size in arcseconds [''] of region extracted from each image to be stacked. |
| images  | list  |  Empty list   | List of images that will be used in the stack. They must be fits files. If the file is not in the same folder as continuum_stacking.py, full path must be provided.  |
|weightstack | bool  |   False       | If True, stacking is done using weighting, i.e. a galaxy can contribute in a different percentage to the combined emission. If this is the case, the fourth column of catalog is used as weight for stacking. If false, fourth column is ignored and stacking is done without weighting (all galaxies contribute the same amount).  | 
|threshold   | int   | 5             | The stacking code returns results for all sources but also it divides the stacking for all detected and all undetected. To define what is considered a detection, the standard deviation (σ) of each galaxy image is measured and if the center of the image has a total flux higher than threshold*σ, the image (or that galaxy) is considered detected. For instance, if threshold = 5, all images with center flux > 5σ will go into a stack of detected galaxies. |
| overwrite  | bool  |  False        | When this code is executed, files of the continuum emission from galaxies in each image are saved. If False, the code will only save these files for new files added to images (in case of running the code again with more images). If True, the code will save these files for all files in images list. Note that this is meant to save time in case the extracted expectra for a image was already done, but it is not helpful if you want to change the catalog between runs (in this case, set to false).          | 
|  verbose   | int   |    1          | Its value defines the amount of prints when running the code. Verbose values can be 0, 1, 2 or 3, where 0 provides minimal amount of prints and 3 provides more detailed messages. | 

### Outputs 

#### Main results 
The main results of the stack are the stacked images (fits files). These images can be analized in CASA or other softwares to measure the flux in case of detection. The headers of the files are come from one of the images that went into the stack, therefore, the information from the headers are not intended for use. Also, some text files are provided as coordinate files, that include ra (rad), dec (rad), weight, ra (deg) and dec (deg) for sources.

The following files are generated everytime the code is executed:

| File type  | File name             |  Description                                              |
|------------|-----------------------|-----------------------------------------------------------|
| FITS       | mean_stack.fits       | Stacked image in mean                                     |
| FITS       | median_stack.fits     | Stacked image in median                                   |
| FITS       | all_sources_cube.fits | Cube with all images that went into the stack             |
| Text       | all_sources.txt       | Coordinate file for all sources that went into the stack  |

The following files are generated if detections are found:

| File type  | File name             |  Description                                       |
|------------|-----------------------|----------------------------------------------------|
| FITS       | sources_cube_det.fits |  Cube with all images that have a detection        |
| Text       | sources_det.txt       |  Coordinate file for sources with a detection      |
| FITS       | mean_stack_det.fits   |  Stacked image of detections in mean               |
| FITS       | median_stack_det.fits |  Stacked image of detections in median             |
| FITS       | std_stack_det.fits    |  Standard deviation map for images with detections |

The following files are generated if non-detections are found:

| File type  | File name                |  Description                                          |
|------------|--------------------------|-------------------------------------------------------|
| FITS       | sources_cube_nondet.fits | Cube with all images that don't have a detection      |
| Text       | sources_non_det.txt      | Coordinate file for sources without a detection       |
| FITS       | mean_stack_nondet.fits   |  Stacked image of non-detections in mean              |
| FITS       | median_stack_nondet.fits |  Stacked image of non-detections in median            |
| FITS       | std_stack_nondet.fits    |  Standard deviation map for images without detections |

The following files are generated if edge detections (within the outer 1/3rd of the image) are found:

| File type  | File name                  |  Description                                             |
|------------|----------------------------|----------------------------------------------------------|
| FITS       | sources_cube_edge_det.fits | Cube with all images that have a detection on the edge   |
| Text       | sources_edge_det.txt       | Coordinate file for sources with a detection on the edge |

#### Intermediate files
Some files are generated as intermediate files when the code is executed. These files are generated per image (or imagename) used in the stack: 

| File type  | File name                                     |  Description                                       |
|------------|-----------------------------------------------|----------------------------------------------------|
| FITS       | dirname+imagename+extracted_sources_cube.fits |  Cube with all images found in imagename           |
| Text       | dirname+imagename+galaxyinfo.txt              | Coordinate file for sources found in imagename     |
| Text       | dirname+imagename+notfound.txt                | Coordinate file for sources not found in imagename |


#### Figure files
The stacking code provides standard figures to visualize the stacking results. If you don't want to generate a specific figure, simply comment the line in which the figure is generated (which is right below the name of the figure). All figures are included in the function plots_cont, which is executed after the stacking code finished, so both the stacking and the plotting can be done at separate times (commenting one or the other). 

At the beginning of the function the user can decide a unit to show the results in mJy, Jy or μJy (unit = 'mJy', or 'Jy', 'muJy'). This can also be changed in a specific figure (in the unit argument). The table below describes the figures avaliable to plot with the name suggested by the code, however, this name can be changed as pleased (but be aware of the matplotlib bug that won't create a figure of file name is too long). Also, all figures are .png files for size purposes, but another format like pdf can also be used.  If you want to change the number of decimals shown for the variables, this must be changed within the function that creates the figure, in the lines using plt.title.
 

| File name                                 |  Description                                        |
|-------------------------------------------|-----------------------------------------------------|
| mean_stack.png - median_stack.png         | Figures of mean/median stack images. It provides a colorbar for flux, method of stacking (mean or median), number of sources that went into the stack (N), mean standard deviation σ of every image that went into the stack, and the standard deviation σ of the stack map. Black contours are shown for -3σ,-2σ,2σma,3σ,4σ,5σ levels of flux, the 1/3rd center is denoted by black square, the center of the image is denoted by black cross and the beam size is shown in the left bottom corner in black.|
| mean_stack_fit.png - median_stack_fit.png |  The same map (mean/median) as above but with the fitted function _fitcomponents_ from the iatool from CASA within python, using one_panel_fit function. If the fit is successful, the title of the image includes the total flux and SNR of the detection. If not, the title will include 'Fit unsucessfull'. Also, if the upperlimit argument of the function is set to True, the title will say 'Upper limit'.|
| mean_stack_smooth.png - median_stack_smooth.png               | Similar to me(di)an_stack.png, but three panel figures (mean/median) with different smoothing levels (0, 3 and 6 pixels) applied to map.|    
| mean_stackdet.png - median_stackdet.png                       | Same as me(di)an_stack.png but for stack of only detections.            |  
| mean_stackdet_smooth.png - median_stackdet_smooth.png         | Same as me(di)an_stack_smooth.png but for stack of only detections.     |  
| mean_stackdet_fit.png - median_stackdet_fit.png               | Same as me(di)an_stack_fit.png but for stack of only detections.        |  
| mean_stacknondet.png - median_stacknondet.png                 | Same as me(di)an_stack.png but for stack of only non-detections.        |    
| mean_stacknondet_smooth.png - median_stacknondet_smooth.png   | Same as me(di)an_stack_smooth.png but for stack of only non-detections. |    
| mean_stacknondet_fit.png - median_stacknondet_fit.png         | Same as me(di)an_stack_fit.png but for stack of only non-detections.    |    

</details>


