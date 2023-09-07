# Documentation under construction

## Instructions

Before using the codes, the python modules _astropy_ and _lineid_plot_ must be installed within CASA. These modules can be installed following the [astropy documentation](https://docs.astropy.org/en/stable/install.html\#installing-astropy-into-casa). First, open CASA in your console and install _pip_ using:

CASA <1> : from setuptools.command import easy_install
CASA <2> : easy_install.main(['--user', 'pip'])

Quit CASA and re-open it. Now to install _astropy_ and _lineid_plot_ type,

CASA <1>: import subprocess, sys
CASA <2>: subprocess.check_call([sys.executable, '-m', 'pip','install','--user', 'astropy', 'lineid_plot'])

Quit CASA and re-open it. Now it should be possible to import either module using, 

CASA <1> : import astropy

Then, go to the directory in which the stacking code is located and use execfile to run on CASA console:
CASA <1> : execfile('continuum_stacking.py')
