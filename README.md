# Spectral and continuum stacking codes 

## Overview
The spectral and continuum stacking codes are two standalone Python scripts for performing spectral stacking (in 3d cubes) and continuum stacking (images), respectively. 
These codes were presented and described in Guerrero et al. submitted, in which stacking was applied in ALMA cubes and images.

Each code can be used separately and are meant to be run within CASA (Common Astronomy Software Applications).

Authors: Andrea Guerrero and Neil Nagar.

## Preparation 
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

Then both codes are ready to be run within CASA.

### Examples
Before running the codes with your own data, we encorage you to run the examples provided to make sure the codes are working properly. 
First, download the folder(s) desired (continuum and/or spectral). Inside the folder, you will find the stacking code and a folder named: _stack_example_.
Inside _stack_example_ there's a fits file and a catalog. For the spectral code, the fits file is a 3d cube and the catalog is a csv with three columns: ra (deg), dec (deg) and redshift. 
For the continuum code, the fits file is a continuum image and the catalog is a csv file with two columns: ra(deg) and dec (deg). Each code has the inputs ready to be run with the examples.

To run, open your console, go to the directory in which the stacking code is located and run CASA. Then use execfile to run the code desired on CASA console. For instance:
> CASA <1> : execfile('continuum_stacking.py')

The stacking results will be saved on the _stack_example_. The results include intermediate files and plots. Details can be found further below. 

## Spectral stacking code 
** Documentation under construction ** 

## Continuum stacking code 
** Documentation under construction ** 





