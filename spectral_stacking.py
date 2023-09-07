# This Python file uses the following encoding: utf-8

import os
import numpy as np
import matplotlib.pyplot as plt
import lineid_plot
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.table import Table
from astropy.modeling import models,fitting
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
import csv
import taskinit
import warnings

class CoordList(list):
    """
        Extended list to contain list of coordinates.

    """
    def __init__(self, imagenames=[], coord_type='physical', unit='rad'):
        """
        Requires an image list in case of pixel coordinates to work properly.
        """

        super(CoordList, self).__init__()

        if isinstance(imagenames, str):
            imagenames = [imagenames]

        self.coords = []
        self.imagenames = imagenames
        self.coord_type = coord_type
        self.unit = unit

    def __getitem__(self, i):
        return self.coords[i]

    def __setitem__(self, i, x):
        self.coords[i] = x

    def append(self, x):
        self.coords.append(x)

    def __len__(self):
        return len(self.coords)

    def __iter__(self):
        for x in self.coords:
            yield x

    def __getslice__(self, i, j):
        new_a = CoordList(self.imagenames, self.coord_type, self.unit)
        new_a.coords = self.coords.__getslice__(i, j)
        return new_a

    def __repr__(self):
        ret = []
        for x in self.coords:
            ret.append('(' + x.__str__() + ')')
        return '\n'.join(ret)

    def __str__(self):
        ret = []
        for x in self.coords:
            ret.append('(' + x.__str__() + ')')
        return '\n'.join(ret)


class Coord:
    """
        Describes a stacking position.
        Class used internally to represent coordinates. May describe a
        physical coordinate or a pixel coordinate.
    """


    def __init__(self, x, y, redshift, weight=1., ra=0., dec=0.,image=0):
        """
            Create a coordinate. A pixel coordinate should always specify
            to which image it belongs. Physical coordinates should be in
            J2000 radians.
        """
        self.x = x
        self.y = y
        self.redshift= redshift
        self.weight = weight
        self.ra = ra
        self.dec = dec
        self.image = image

    def __str__(self):
        return '{0},{1},{2},{3},{4},{5},{6}'.format(self.x, self.y ,self.redshift, self.weight, self.ra, self.dec,self.image)


def readCoords(coordfile,weightstack, unit='deg'):
    """
        Reads a coordinate file from disk and produces a list of coordinates in radians.

        coordfile:
            Path to coordinate file. A file in csv format. x and y should
            be in J2000 . A weight may be added in third column
            to weight positions for stacking. If no weight are wanted
            put no third column in coordinate file.
        unit:
            Unit of input coordinates. Allows two values, 'deg' and 'rad'.
    """
    coordreader = csv.reader(open(dirname+'/'+coordfile, 'rb'), delimiter=',')
    coords = CoordList()
    for row in coordreader:
        x = float(row[0])
        y = float(row[1])
        ra = x    #maintain a copy of ra and dec in deg
        dec= y
        redshift= float(row[2])
        if unit == 'deg':  #to radians
            x = np.pi/180.*x
            y = np.pi/180.*y

        if weightstack == 'True':
            weight = float(row[3])
        else:
            weight = 1.

        if x > 2*np.pi:
            x -= 2*np.pi
        if y > np.pi:
            y -= 2*np.pi

        coords.append(Coord( x, y, redshift, weight, ra, dec))

    return coords


def _getPixelCoords1Im(coords, imagename):
    """
    Called by _getPixelCoords.
    Makes a pixel coord list of the coordinates that are inside an image.
    """
    ia = taskinit.iatool()
    ia.open(imagename)
    cs = ia.coordsys()
    Nx = ia.shape()[0]
    Ny = ia.shape()[1]
    ia.done()
    x0 = cs.referencevalue()['numeric'][0]
    y0 = cs.referencevalue()['numeric'][1]
    x_pix_ref = cs.referencepixel()['numeric'][0]
    y_pix_ref = cs.referencepixel()['numeric'][1]
    x_pix_inc = cs.increment()['numeric'][0]       #in rad/pixel
    y_pix_inc = cs.increment()['numeric'][1]

    cenimagera = (Nx/2. - x_pix_ref)* x_pix_inc + x0  #radians
    cenimagedec= (Ny/2. - y_pix_ref)* y_pix_inc + y0
    sizeimagex = (Nx*x_pix_inc)   #radians
    sizeimagey = (Ny*y_pix_inc)  


    pixcoords = []
    for coord in coords:
        if (abs(coord.x - cenimagera) < 1.5*abs(sizeimagex)) and (abs(coord.y - cenimagedec) < 1.5*abs(sizeimagey)): 
            dx = (coord.x - x0)*np.cos(coord.y)
            dy = np.arcsin(np.sin(coord.y)/np.cos(dx)) - y0
            x = dx/x_pix_inc+x_pix_ref
            y = dy/y_pix_inc+y_pix_ref
            if (x > 0 and x < Nx-1) and (y > 0 and y < Ny-1):
                c = Coord(x, y, coord.redshift, coord.weight, coord.ra, coord.dec, coord.image)
                try:
                    c.index = coord.index
                except AttributeError:
                    pass
            
                pixcoords.append(c)
            
    return pixcoords


def getPixelCoords(coords, imagenames):
    """
        Creates pixel coordinate list from a physical coordinate
        list and a list of images.
    """

    pixcoords = CoordList(imagenames, 'pixel', unit='pix')
    for (i, imagename) in enumerate(pixcoords.imagenames):
        for coord in _getPixelCoords1Im(coords, imagename):
            coord.image = i
            pixcoords.append(coord)
    return pixcoords


def getCubes(coords, overwrite):
    """
    Gets list of cubes names to stack.
    If overwrite is True, cubenames will be all cubes, if not,
    it will be only the cubes missing from the output directory

    """
    cubenames = []
    if overwrite == False:
        for cubename in coords.imagenames:
            outspecfile = cubename.replace('.fits', '-extractspec.txt')  
            if not os.path.exists(dirname + '/' + outspecfile):
                cubenames.append(cubename)

    if overwrite == True:
        cubenames = coords.imagenames

    return cubenames


def stack(catalog_name='', stampsize = 3, cubenames= [], weightstack = False, overwrite= False, verbose=1):
    """
        Performs spectral stacking in the image domain.

        Inputs:
            catalog_name  -- Name of csv catalog which has the coordinates of the targets (ra[deg], dec[deg] and redshift columns required).
            stampsize     -- Size of target region in arcseconds that will be extracted for stacking.
            cubenames     -- Name of cubes to extract spectra from.
            weightstack   -- If True, takes fourth column of csv catalog and uses it as weight column for stacking.
            overwrite     -- When this code is executed, files of the continuum emission from each input image are saved. If overwrite = False,
                             the code will only save these files for new images added. 
                             If overwrite = True, the code will save these files for all input images. 
            verbose       -- Verbose level from 0 to 3, where 0 provides minimal amount of prints and 3 provides more detailed messages.
        Outputs:
            Stacked spectrum  (and writes out other text and data files)
    """

    coords = readCoords(catalog_name,weightstack)
    
    if verbose >= 1:
        print('Total number of objects in catalog: %s'%len(coords))
        print('Cubes to be used for spectral stacking: %s'%cubenames)
        print('Stamp size [arcsec]: %s'%stampsize)
    #transforms coordinates from radians to pixels
    if coords.coord_type == 'physical':
        coords = getPixelCoords(coords, cubenames)
        
    imagenames = getCubes(coords, overwrite)

    if verbose >= 1:
        print('\n Starting spectral stacking extraction')

    index=0
    for imagename in imagenames:
        startline=100000
        endline=0
        count=0
        for coord in coords:
            if coord.image==index:
                if startline==100000:
                    startline=count
                endline=count
            count=count+1
                
        imcoords=coords[startline:endline+1]

        ia = taskinit.iatool()
        ia.open(str(imagename))
        cs = ia.coordsys()
        x_pix_inc = cs.increment()['numeric'][0]   #CDELT1   #in rad/pixel
        y_pix_inc = cs.increment()['numeric'][1]   #CDELT
        ia.done()

        new_stampsize = int((stampsize * np.pi/(180.*3600.))/y_pix_inc)
        _load_stackspec(imagename,new_stampsize,imcoords,verbose,weightstack)
        index=index+1
        
    if verbose >= 1:
        print('\n Spectral extraction finished')
        print('\n Starting spectral stacking')

    stacked_spec  =  _stack_stackspec(coords,catalog_name)    

    if verbose >= 1:
        print('\n Spectral stacking finished')

    return coords,stacked_spec


def _load_stackspec(imagename,new_stampsize,coords,verbose,weightstack):                      
    """
    Extracts spectral and frequency information from a subcube. Each subcube is the size of the stampsize, 
    with the center coordinate from the coordinate file.  

    3 files are generated: -spectra file that ends with -extractspec.txt
                           -frequency per channel that ends with -freqdef.txt
                           -file that contains (coord x, coord y , coord redshift, weight, coord ra, coord dec) for positions that will be use in stack
    """
    
    if verbose >= 2:
        print('\n Reading image: %s'%imagename)
    if verbose >= 3:
        print('The coordinates found in this cube are (ra(rad),dec(rad),z,weight,ra(dec),dec(dec)): %s'%coords)

    imagesize = []
    ia = taskinit.iatool()
    ia.open(str(imagename))
    imagesize.append((ia.shape()[0], ia.shape()[1]))
    ia.done()

    hdu = fits.open(imagename)
    extractcubedata = hdu[0].data                                             #shape is (polarisations, number of channels, size in y, size in x)    
    extractcubeheader = hdu[0].header

    if len(extractcubedata.shape) == 3:
        dataspec = np.zeros((len(coords), extractcubedata.shape[0]))*np.nan

    if len(extractcubedata.shape) == 4:
        dataspec = np.zeros((len(coords), extractcubedata.shape[1]))*np.nan

    if verbose >= 3:
        print('Reading header of cube')
   
    f0 = float(extractcubeheader['CRVAL3'])                                # reference freq in Hz
    df = float(extractcubeheader['CDELT3'])                                # channel width in Hz                                                                                      
    i0 = extractcubeheader['CRPIX3'] - 1                                   # reference pixel (note that first channel is '1')

    if len(extractcubedata.shape) == 3:
        freqspec = ((np.arange(extractcubedata.shape[0])-(i0))*df + f0)/1.e9

    if len(extractcubedata.shape) == 4:
        freqspec = ((np.arange(extractcubedata.shape[1])-(i0))*df + f0)/1.e9

    extractcubedata=np.squeeze(extractcubedata)                            #shape:(number of channels, size in y, size in x) 

    galaxyinfo=[]   #every position which has a spectrum 
    notfound=[]
    count=0
    
    listblcx=[]
    listblcy=[]
    listtrcx=[]
    listtrcy=[]
    for (i,coord) in enumerate(coords): 
        blcx = int(coord.x - new_stampsize/2 +0.5)
        blcy = int(coord.y - new_stampsize/2 +0.5)
        trcx = blcx + int(new_stampsize)
        trcy = blcy + int(new_stampsize)
        
        listblcx.append(blcx)  
        listblcy.append(blcy) 
        listtrcx.append(trcx) 
        listtrcy.append(trcy) 
    
    beamtable = Table.read(hdu, hdu=1)
    
    cube_beammajor = np.array(beamtable['BMAJ'] /3600.)     #degree
    cube_beamminor = np.array(beamtable['BMIN'] /3600.)     #degree
    cube_cdelt = extractcubeheader['CDELT2']                #degree

    beam_area = ((np.pi/(4.*np.log(2)))*cube_beammajor*cube_beamminor)/ (cube_cdelt**2) #pix

    if verbose >= 2:
        print('Extracting spectral stamp per position')

    for (i,coord) in enumerate(coords): 
        if (listblcx[i]>0 and listblcx[i]<imagesize[0][0]-1) and (listtrcx[i]>0 and listtrcx[i]<imagesize[0][0]-1) and (listblcy[i]>0 and listblcy[i]<imagesize[0][1]-1) and (listtrcy[i]>0 and listtrcy[i]<imagesize[0][1]-1):
            valspec  = extractcubedata[:,listblcy[i]:listtrcy[i],listblcx[i]:listtrcx[i]]

            valspec=valspec.sum(axis=1)  #sum in RA
            valspec=valspec.sum(axis=1)  #sum in DEC   # now a single spectrum
            
            valspec_c = np.divide(valspec, beam_area)
            
            if df < 0:
                valspec_c=valspec_c[::-1]
                       
            if weightstack == True:
               valspec_w = valspec_c*coord.weight
               dataspec[count,:]=valspec_w.tolist() 
            else:
               dataspec[count,:]=valspec_c.tolist()

            galaxyinfo.append([coord.x,coord.y,coord.redshift,coord.weight,coord.ra,coord.dec])        
            count=count+1
            
        else:   
            notfound.append([coord.x,coord.y,coord.redshift,coord.weight,coord.ra,coord.dec])
    
    if verbose >= 3:
        print('Number of regions not found in the image: %s'%len(notfound))

    galaxyinfo = np.array(galaxyinfo)
    dataspec=dataspec[0:count,:]

    ind = []
    countnonnan=0    
    for i in range(len(dataspec)):
        arraynan=np.argwhere(np.isnan(dataspec[i,:]))
        if len(arraynan) != len(dataspec[i,:]):
            countnonnan=countnonnan+1
        else:
            ind.append(i)

    if verbose >= 3:
        print('Deleting empty stamps (nan flux)')

    dataspec   = np.delete(dataspec  , np.array(ind), axis=0)
    galaxyinfo = np.delete(galaxyinfo, np.array(ind), axis=0)

    cube_name = imagename.replace("/","-")

    outspecfile=cube_name.replace('.fits','-extractspec.txt')
    outfreqdeffile=cube_name.replace('.fits','-freqdef.txt')
    outgalfile=cube_name.replace('.fits','-galaxyinfo.txt')
    outnotfound=cube_name.replace('.fits','-notfound.txt')

    nanvalues = np.isnan(dataspec)
    dataspec[nanvalues] =  -99.0 
    np.savetxt(dirname+'/'+str(outspecfile),dataspec)          #saves file with spectra                   #2D array, shape is (number of pos, number of channels)
    
    nanvalues = np.where(dataspec == -99.0)
    dataspec[nanvalues] = np.nan

    if df < 0:
        freqspec=freqspec[::-1]

    filed=np.savetxt(dirname+'/'+str(outfreqdeffile),freqspec)  #saves file with frecuency per channel    #1D array, len is number of pos
    np.savetxt(dirname+'/'+str(outgalfile),galaxyinfo)
    np.savetxt(dirname+'/'+str(outnotfound),notfound)
    
    if verbose >= 2:
        print('Number of objects with non nan flux in image: %s'%countnonnan)
    
    if verbose >= 2:
        print('Saving files of spectra, observed frequency and galaxy information')

    return
    

def _stack_stackspec(coords, catalog_name):   
    """
    Takes observed spectra, converts it into rest-frame spectra and inserts it into a large array to be stack.
    """
    restfreqmin=[]                            #list of minimum frequency per cube
    restfreqmax=[]                            #list of maximum frequency per cube 
    df=[]                                     #list of frequency width per channel per cube
    totnumpos=0
    maxnumchan=0
    
    if verbose >= 2:
        print('\n Reading each cube to save frequency range')
    
    for imagename in coords.imagenames:
        cube_name = imagename.replace("/","-")
        
        outfreqdeffile=cube_name.replace('.fits','-freqdef.txt')
        outgalfile=cube_name.replace('.fits','-galaxyinfo.txt')
                
        freqfile= np.loadtxt(dirname+'/'+outfreqdeffile)                     #reads file of frequency per channel (of entire cube)
        galfile=np.loadtxt(dirname+'/'+outgalfile)                           #reads file of redshift of galaxies used in stacking
    
        if len(galfile.shape) == 1:
            galfile = galfile.reshape(1, len(galfile))
        
        if galfile.shape[1] == 0:
            continue
    
        z_vector = galfile[:,2]

        trestfreqmin = int(np.min(freqfile)*(np.min(z_vector)+1))            
        trestfreqmax = int((np.max(freqfile)*(np.max(z_vector)+1)))+1 
        
        if verbose >= 3:
            print('Creating rest-frame frequency array')
            print('Frequency range for cube %s: %.1f - %.1f GHz'%(imagename,trestfreqmin,trestfreqmax))
            
        tdf= freqfile[2]-freqfile[1]                                   
        
        restfreqmin.append(trestfreqmin)
        restfreqmax.append(trestfreqmax)
        df.append(abs(tdf))

        totnumpos = totnumpos + galfile.shape[0]                                
        if freqfile.shape[0] > maxnumchan:                                      
            maxnumchan=freqfile.shape[0]
        
    restdf=np.min(df)/2.                           
    restdf= 1./(int(1./restdf))      
    
    restfreqspec=np.arange(np.min(restfreqmin),np.max(restfreqmax)+restdf,restdf)   
    if verbose >= 3:
        print('Final frequency array is between %.1f and %.1f GHz'%(restfreqspec[0],restfreqspec[-1]))

    tottalspectra=np.zeros((5,(len(restfreqspec))))
    tottalspectra[0,:]=np.copy(restfreqspec)

    imagecounter=0    
    cubenamedat=[]
    cubenamefreq=[]
    indiceslist=[]
    
    if verbose >= 2:
        print('Converting observed spectra to rest-frame spectra')

    for imagename in coords.imagenames:
        cube_name = imagename.replace("/","-")
        outspecfile   =cube_name.replace('.fits','-extractspec.txt')
        outfreqdeffile=cube_name.replace('.fits','-freqdef.txt')
        outgalfile    =cube_name.replace('.fits','-galaxyinfo.txt')
        
        obsspec= np.loadtxt(dirname+'/'+outspecfile)                   #array of observed spectra in cube;  shape is (number of coords used, number of channels)
        obsfreq= np.loadtxt(dirname+'/'+outfreqdeffile)                #array of observed frequency in the cube; shape is (number of channels)
        galfile=np.loadtxt(dirname+'/'+outgalfile)                     #array with galaxy information. Redshifts of used coords are in column 2

        nanvalues = np.where(obsspec == -99.0)
        obsspec[nanvalues] = np.nan
        
        if len(galfile.shape) == 1:
            galfile = galfile.reshape(1, len(galfile))
            obsspec = obsspec.reshape(1, len(obsspec))
            
        if galfile.shape[1] == 0:
            continue
        
        z_vector = galfile[:,2]

        restfreqspecimage=np.arange(restfreqmin[imagecounter],restfreqmax[imagecounter]+restdf,restdf)  
        sizerestfreq=len(restfreqspecimage)

        imagecounter = imagecounter+1
        numpos = obsspec.shape[0]                                
        
        restspec = np.zeros((numpos,len(restfreqspecimage)))    

        for poscounter in range(numpos):
            freq_with_z = obsfreq*(1+z_vector[poscounter])        
            restfreqspecpos = np.arange(int(np.min(freq_with_z)),np.max(freq_with_z)+restdf,restdf)         
            restspec[poscounter,:]=griddata(freq_with_z,(obsspec[poscounter,:])/(1+z_vector[poscounter]),restfreqspecimage,method='linear')                                                             
            tposspec=griddata(freq_with_z,obsspec[poscounter,:],restfreqspecpos,method='linear') 
 
        #####################################################################
        restspec_corr = np.copy(restspec)

        for i in range(restspec_corr.shape[0]):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sigclip = SigmaClip(sigma=3, cenfunc = np.nanmedian, stdfunc=np.nanstd)
                spec_sigclip = sigclip(restspec_corr[i,:], axis=0, copy=True)
                spec_sigclip=spec_sigclip.filled(np.nan)

            arraynan=np.argwhere(np.isnan(spec_sigclip))
            if len(arraynan) != len(spec_sigclip):
                medianspec = np.nanmedian(spec_sigclip)
                restspec_corr[i,:] = restspec_corr[i,:] - medianspec
        #####################################################################
        cubecounter='cubefreq' + str(imagecounter)
        exec(cubecounter  + " = restfreqspecimage ")
        cubenamefreq.append(cubecounter)
        
        cubecont='cubedat' + str(imagecounter) 
        exec(cubecont  + " = restspec_corr ")
        cubenamedat.append(cubecont)
        cube_name = imagename.replace("/","-")
        
        outspeccube=cube_name.replace('.fits','-speccube.txt')
        outresfreqcube= cube_name.replace('.fits','-restfreqcube.txt')
        outspec_corr = cube_name.replace('.fits','-spec_corr.txt')
        
        nanvalues = np.isnan(restspec)
        restspec[nanvalues] = -99.0
        
        nanvalues = np.isnan(restspec_corr)
        restspec_corr[nanvalues] = -99.0
        
        np.savetxt(dirname+'/'+outspeccube,restspec)
        np.savetxt(dirname+'/'+outresfreqcube,restfreqspecimage)
        np.savetxt(dirname+'/'+outspec_corr,restspec_corr)
        
        nanvalues = np.where(restspec == -99.0)
        restspec[nanvalues] = np.nan
        
        nanvalues = np.where(restspec_corr == -99.0)
        restspec_corr[nanvalues] = np.nan
        
    if verbose >= 2:
        print('Performing spectral stack')
       
    for j in range(len(restfreqspec)):
        indices=[]
        tempfreq=restfreqspec[j]
        datapoints=[]
        for i in range(len(cubenamedat)):
            tfile1=cubenamedat[i]
            cubedat0=eval(tfile1)
            tfile2=cubenamefreq[i]
            cubefreq0=eval(tfile2)
            indicesarray=np.where(np.isclose(cubefreq0,tempfreq,atol=restdf/(2+2*tempfreq),rtol=restdf/(2+2*tempfreq)))[0]
            if len(indicesarray) > 0:   
                index= indicesarray[0]
                datapoints.append(cubedat0[:,index])
                indices.append(index)
           
        datapoints=[x for y in datapoints for x in y]
        datapoints=np.array(datapoints)
        
        countnonan=np.sum(np.isnan(datapoints))
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tottalspectra[1,j]=np.nanmedian(datapoints)
            tottalspectra[2,j]=np.nanmean(datapoints)
            tottalspectra[3,j]=np.nanstd(datapoints)
            tottalspectra[4,j]=np.count_nonzero(~np.isnan(datapoints))

    if verbose >= 1:
        print('Final number of objects stacked: %i'%len(tottalspectra[4,:]))

    nanvalues1 = np.isnan(tottalspectra[1,:])
    tottalspectra[1,nanvalues1] = -99.0
    
    nanvalues2 = np.isnan(tottalspectra[2,:])
    tottalspectra[2,nanvalues2]= -99.0

    nanvalues3 = np.isnan(tottalspectra[3,:])
    tottalspectra[3,nanvalues2]= -99.0

    tottalspectra[4,:] = np.nan_to_num(tottalspectra[4,:])

    if verbose >= 2:
        print('Saving stack results in text files')
    
    np.savetxt(dirname+'/stacked-frequency.txt',tottalspectra[0,:])
    np.savetxt(dirname+'/stacked-medianspec.txt',tottalspectra[1,:])
    np.savetxt(dirname+'/stacked-meanspec.txt',tottalspectra[2,:])
    np.savetxt(dirname+'/stacked-stdspec.txt',tottalspectra[3,:])
    np.savetxt(dirname+'/stacked-nobs.txt',tottalspectra[4,:])


    

    return tottalspectra
    
def plots_style():
    ### Defines the plots style, fonts sizes, tick sizes ####
    fsize = 18
    lsize = 18
    tsize =  18
    tlsize = 18
    msize = 8
    tdir = 'in'
    major = 6.0
    minor = 4.0
    lwidth = 0.9
    lhandle = 2.0
    plt.rcParams['font.size'] = fsize
    plt.rcParams['lines.markersize'] = msize
    plt.rcParams['legend.fontsize'] = tsize
    plt.rcParams['xtick.direction'] = tdir
    plt.rcParams['ytick.direction'] = tdir
    plt.rcParams['xtick.major.size'] = major
    plt.rcParams['xtick.minor.size'] = minor
    plt.rcParams['ytick.major.size'] = 5.0
    plt.rcParams['ytick.minor.size'] = 3.0
    plt.rcParams['axes.linewidth'] = lwidth
    plt.rcParams['axes.labelsize'] = lsize
    plt.rcParams['legend.handlelength'] = lhandle
    plt.rcParams['font.family'] = "serif"
    plt.rcParams['xtick.labelsize'] = tlsize
    plt.rcParams['ytick.labelsize'] = tlsize
    return

def define_unit(unit):
    if unit == 'mJy':
        c = 1.e3
        unit_text='m'
    if unit == 'Jy':
        c = 1.
        unit_text=' '
    if unit== 'muJy':
        c = 1.e6
        unit_text = '$\mu$'
    return c, unit_text 

def everyplot(spec,std,n,restfreq,title,figname, unit='mJy'):
    global lines_freq
    global lines_label
    
    plots_style()
    c, unit_text = define_unit(unit)
    
    s = 2.
    gauss_width = 10.  
    w = gauss_width
    t = (((w-1.)/2.)-0.5)/s
    
    spec = gaussian_filter1d(spec, sigma = s, truncate = t)
    
    fig=plt.figure(figsize=(15,10))
    ax1 = fig.add_axes([0.08, 0.3, 0.88, 0.6],xticklabels=[])
    ax2 = fig.add_axes([0.08, 0.1, 0.88, 0.2],sharex=ax1)
    minspec=np.nanmin(np.vstack([d for d in spec]))
    maxspec=np.nanmax(np.vstack([d for d in spec]))
    
    ax1.plot(restfreq,c*spec)
    ax1.plot([restfreq[0],restfreq[-1]],[0,0],'r-')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ax1.fill_between(restfreq,c*spec, where=spec>0,alpha=0.6, interpolate= True, label=None)
        ax1.fill_between(restfreq,c*spec, where=spec<= 0,alpha=0.6, interpolate= True, label=None)
    
    spec_no_nan_ind= np.argwhere(~np.isnan(spec))
    ax1.hlines(np.nanmedian(spec),restfreq[spec_no_nan_ind[0]], restfreq[spec_no_nan_ind[-1]] )
    
    ax1.set_xlim(restfreq[spec_no_nan_ind[0]], restfreq[spec_no_nan_ind[-1]])
    ax1.set_ylim(1.1*c*minspec,c*maxspec*2.5 )
    ax1.set_ylabel('Flux [%sJy/beam]'%unit_text)
    substd = np.nanstd(spec)
    ax1.plot([restfreq[0],restfreq[-1]], [c*substd,c*substd],'r:', lw = 5,  label= r'Mean $\sigma$= %.3f %sJy/beam'%(c*substd,unit_text))
    ax1.legend()

    ax2.plot(restfreq,n)
    ax2.set_xlabel('Rest-frame frequency [GHz]')
    ax2.set_ylabel(u'N')
    ax2.set_ylim(0,np.nanmax(n)*1.2)
    
    plt.locator_params(axis='y', nbins=5)
    lines_freq =[115271.202,230538.0,345795.99,461040.768,576267.931,691473.076,806651.801,921799.704,1036912.385,1151985.443,492160.651,809341.97,1900536.9,52033.143,987926.759,1097364.790,1162911.60200,1207638.73,1228788.7190,1410618.069,88631.60,265886.43,354505.47,443116.14,531716.35,708877.00,797433.26,885970.69,110201.35,220398.86,330587.97,440765.17,550926.29,661067.28,771184.13,881272.81,89188.53,178375.07,267557.63,356734.24,445908.91,535061.64,624208.46,713341.37,2060068.86,3393006.24,5785879.59,48990.96,97980.95,146969.03,195954.21,244935.56,293912.09,342882.86,391846.89,440803.24,489750.93,538689.00]        
    lines_freq= list(np.array(lines_freq)/1000.)
    
    lines_label=['CO1-0','CO2-1','CO3-2','CO4-3','CO5-4','CO6-5','CO7-6','CO8-7','CO9-8','CO10-9','CI 3P1-3P0','CI 3P2-3P1','CII','H2O 2(1,1)-2(0,2)','H2O 2(0,1)-1(1,1)','H2O 3(1,2)-3(0,3)','H2O 3(2,1)-3(1,2)','H20 4(2,2)-4(1,3)','H2O 2(2,0)-2(1,1)','H20 5(2,3)-5(1,4)','HCN 1-0','HCN 3-2','HCN 4-3','HCN 5-4','HCN 6-5','HCN 8-7','HCN 9-8','HCN 10-9','13CO 1-0','13CO 2-1','13CO3-2','13CO 4-3','13CO 5-4','13CO 6-5','13CO 7-6','13CO 8-7','HCO+ 1-0','HCO+ 2-1','HCO+ 3-2','HCO+ 4-3','HCO+ 5-4','HCO+ 6-5','HCO+ 7-6','HCO+ 8-7','OI','OIII 3P1-3P0','OIII 3P2-3P1','CS 1-0','CS 2-1','CS 3-2','CS 4-3','CS 5-4','CS 6-5','CS 7-6','CS 8-7','CS 9-8','CS 10-9','CS 11-10']
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lineid_plot.plot_line_ids(restfreq,c*spec,lines_freq,lines_label,arrow_tip=c*maxspec*1.0, ax=ax1, label1_size = np.array([20]*len(lines_label)))

    plt.savefig(str(figname), dpi = 150,bbox_inches='tight')
    plt.close(fig)
    
    return 

def gaussian(x,a,mu,sigma):
    g=a*np.exp(((-(x-mu)**2.)/(2.*(sigma**2))))
    return g

def line_flux(peak, width, peak_err, width_err):
    flux = peak * width * np.sqrt(2*np.pi)
    flux_err = flux*np.sqrt((peak_err/peak)**2 + (width_err/width)**2)
    return flux, flux_err

def speclinesnr(line,stacked_frequency,stacked_spec,stacked_spec_median,stacked_spec_std,stacked_galaxies, figname, title, makeplot, unit='mJy'):
    plots_style()
    c, unit_text = define_unit(unit)
 
    i = np.where(np.array(lines_label) == str(line))[0][0]
    temp_freq = lines_freq[i]
    temp_label= lines_label[i]

    dv = (2000./3e5) * temp_freq
    startfreq= temp_freq - (dv)
    endfreq= temp_freq + dv
    chans=np.where((stacked_frequency>startfreq)&(stacked_frequency<=endfreq))[0]
    if len(chans) < 5:
        if verbose >=2:
            print('Emission line %s not found in frequency range or not enough channels'%line)
    else:
        subspec=stacked_spec[chans]
        subfreq=stacked_frequency[chans]
        substd=stacked_spec_std[chans]
        subnobs=stacked_galaxies[chans]

        subspec_median = stacked_spec_median[chans]
        #subfreq_median=stacked_frequency[chans]
        #substd_median=stacked_spec_std[chans]
        #subnobs_median=stacked_galaxies[chans]

        ##################
        subspec_wo_line=np.copy(subspec)
        ######################
        #parameter to use a smoothing by
        s = 2.
        gauss_width = 10.   #channels
        w = gauss_width
        t = (((w-1.)/2.)-0.5)/s

        subspec = gaussian_filter1d(subspec, sigma = s, truncate = t)
        subspec_median =   gaussian_filter1d(subspec_median, sigma = s, truncate = t)

        arraynan=np.argwhere(np.isnan(subspec))
        ind = np.where(np.isnan(subspec))[0]

        minspec=np.nanmin(np.vstack([d for d in subspec]))
        maxspec=np.nanmax(np.vstack([d for d in subspec]))

        #gaussian fit to line
        ind = np.where(~np.isnan(subspec))[0]
        subspec_red =  subspec[ind]
        subfreq_red =  subfreq[ind]
        
        subspec_wo_line=np.copy(subspec_red)
        subfreqkms = (subfreq_red - temp_freq)*3e5 / temp_freq

        g1 = models.Gaussian1D(amplitude=np.nanmax(c*subspec_red), mean=0, stddev=80)
        fit_g = fitting.LevMarLSQFitter()
        g = fit_g(g1, subfreqkms, c*subspec_red, maxiter = 1000)
        x_g = np.linspace(np.nanmin(subfreqkms), np.nanmax(subfreqkms), 10000)

        peak, mean, width = g.parameters
        fit_errs = np.sqrt(np.diag(fit_g.fit_info['param_cov']))
        peak_err, mean_err, width_err = fit_errs
        total_flux, total_flux_err = line_flux(peak, width, peak_err, width_err)

        popt,pcov = curve_fit(gaussian, subfreq_red ,c*subspec_red, p0=[np.nanmax(c*subspec_red), temp_freq, 0.2], maxfev=360000)
        chans_line=np.where((subfreq_red<=(popt[1] + (5*popt[2])) ) & (subfreq_red>(popt[1] - (5*popt[2]))  ))[0]
        subspec_wo_line[chans_line] = np.nan
        peak_line = popt[0]
        p_errors = np.sqrt(np.diag(pcov))
        substd_woline= np.nanstd(subspec_wo_line)
        snr_line = peak /  (c*substd_woline)

        ################################
        if makeplot == 1:
            plt.clf()
            fig1=plt.figure(figsize=(15,8))
            ax1 = fig1.add_axes([0.08, 0.2, 0.88, 0.8],xticklabels=[])
            ax2 = fig1.add_axes([0.08, 0.0, 0.88, 0.2],sharex=ax1)

            ax1.plot(subfreq,c*subspec, color='k', label=None)

            ax1.plot([startfreq,endfreq],[0,0],'r-', label=None)
            ax1.plot(subfreq_red, gaussian(subfreq_red, *popt), 'r--',linewidth =3, label= 'Peak flux = %.1f $\pm$ %.1f %sJybeam$^{-1}$ \nTotal flux = %.2f $\pm$ %.2f Jy km s$^{-1}$ \nFWHM = %i $\pm$ %i km s$^{-1}$ '%(peak,peak_err,unit_text,1e-3*total_flux,1e-3*total_flux_err,width*2.2,width_err*2.2))

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ax1.fill_between(subfreq,c*subspec, where=subspec>0,alpha=0.6, interpolate= True, label=None)
                ax1.fill_between(subfreq,c*subspec, where=subspec<= 0,alpha=0.6, interpolate= True, label=None)

            ax1.plot([startfreq, endfreq], [c*substd_woline,c*substd_woline],'r:', lw = 5, label= 'rms= %.2f %sJy \nSNR = %i'%(c*substd_woline, unit_text, snr_line))
            ax1.plot([startfreq, endfreq], [c*substd_woline,c*substd_woline],'r:', lw = 5 )
            ax1.plot(subfreq,c*subspec_median,'k-', linewidth =1 , label=None)

            ax1.set_xlim(startfreq,endfreq)
            ax1.set_ylim(1.1*c*minspec,c*maxspec*1.5)
            ax1.set_ylabel('Flux [%sJy beam$^{-1}$]'%unit_text)
            ax1.tick_params(axis='x',  which='both',  bottom='on',  top='on', labelbottom='off' )

            plt.locator_params(axis='y', nbins=3)

            ax2.plot(subfreq,subnobs, label=None)
            ax2.set_ylim(0,np.nanmax(subnobs)*1.2)
            ax2.set_xlabel('Rest-frame frequency [GHz]')
            ax2.set_ylabel(u'N')
            ax1.legend(handlelength=2)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                lineid_plot.plot_line_ids(subfreq,c*subspec,lines_freq,lines_label,arrow_tip=c*maxspec*1.1,ax=ax1, label=None, label1_size = np.array([20]*len(lines_label)))

            plt.savefig(figname+str(temp_label)+'.png', dpi=150,bbox_inches='tight')
            plt.close(fig1)
    return  


def plot_range_freq(spec,stacked_frequency,stacked_spec_std,stacked_galaxies,figname, title,startfreq, endfreq,freqstep=50):
    cont=1
    while startfreq < endfreq:
        chans=np.where((stacked_frequency<startfreq+freqstep)&(stacked_frequency>startfreq))[0]
        subspec=spec[chans]
        subfreq=stacked_frequency[chans]
        substd=stacked_spec_std[chans]
        subnobs=stacked_galaxies[chans]

        arraynan=np.argwhere(np.isnan(subspec))

        if len(arraynan) < 0.95*len(subspec):    #if array is full of nan, dont plot
            plot=everyplot(subspec,substd,subnobs,subfreq,title,figname+str(cont)+'.png')

        startfreq=startfreq+freqstep
        cont=cont+1
    return

def plot_every_common_line_fit(spec,stacked_frequency,stacked_spec_std,stacked_galaxies,figname, title, velocity_width= 2000., unit = 'mJy'):
    c, unit_text = define_unit(unit)
    
    for i in range(len(lines_freq)):
        temp_freq = lines_freq[i]
        temp_label= lines_label[i]

        dv = (velocity_width / 3e5) * temp_freq
        startfreq= temp_freq - ( dv)
        endfreq= temp_freq + dv

        chans=np.where((stacked_frequency<=endfreq)&(stacked_frequency>startfreq))[0]
        subspec=spec[chans]
        subfreq=stacked_frequency[chans]
        substd=stacked_spec_std[chans]
        subnobs=stacked_galaxies[chans]

        s = 2.
        gauss_width = 10.
        t = (((gauss_width-1.)/2.)-0.5)/s

        subspec = gaussian_filter1d(subspec, sigma = s, truncate = t)
        arraynan=np.argwhere(np.isnan(subspec))

        if len(arraynan) < 0.8*len(subspec):
            minspec=np.nanmin(np.vstack([d for d in subspec]))
            maxspec=np.nanmax(np.vstack([d for d in subspec]))
            fig1=plt.figure(figsize=(20,10))
            ax1 = fig1.add_axes([0.08, 0.3, 0.88, 0.6],xticklabels=[])
            ax2 = fig1.add_axes([0.08, 0.1, 0.88, 0.2],sharex=ax1)

            ax1.plot(subfreq,c*subspec)
            ax1.plot([startfreq,endfreq],[0,0],'r-')

            ind = np.where(~np.isnan(subspec))[0]
            subspec_red =  subspec[ind]
            subfreq_red =  subfreq[ind]
            popt,pcov = curve_fit(gaussian, subfreq_red,c*subspec_red, p0=[np.nanmax(c*subspec_red), temp_freq, np.nanstd(c*subspec_red)],maxfev=360000)
            freq_gaussian = np.linspace(np.nanmin(subfreq_red), np.nanmax(subfreq_red), 500)
            ax1.plot(subfreq_red, gaussian(subfreq_red, *popt), 'r--')

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ax1.fill_between(subfreq,c*subspec, where=subspec>0,alpha=0.6, interpolate= True, label=None)
                ax1.fill_between(subfreq,c*subspec, where=subspec<= 0,alpha=0.6, interpolate= True, label=None)

            ax1.set_xlim(startfreq,endfreq)
            ax1.set_ylim(1.1*c*minspec,c*maxspec*1.8)
            ax1.set_ylabel('Flux [%sJy/beam]'%unit_text)
            ax2.plot(subfreq,subnobs)
            ax2.set_ylim(0,np.nanmax(subnobs)*1.2)
            ax2.set_xlabel('Rest-frame frequency [GHz]')
            ax2.set_ylabel(u'N of objects')
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                lineid_plot.plot_line_ids(subfreq,c*subspec,lines_freq,lines_label,arrow_tip=c*maxspec*1.0,ax=ax1)
           
            plt.savefig(figname+str(temp_label)+'.png', dpi=150)
            plt.close(fig1)
    return

def plot_every_galaxy_observed(lines_freq, lines_label, spec_type = 'observed',figname = '.png', unit='mJy'):
    plots_style()

    c, unit_text = define_unit(unit)

    for imagename in incubes:
        cube_name = imagename.replace("/","-")
        outgalfile=cube_name.replace('.fits','-galaxyinfo.txt')
        galfile=np.loadtxt(dirname+'/'+outgalfile)                     #array with galaxy information. Redshifts of used coords are in column 2

        if len(galfile.shape) == 1:
            galfile = galfile.reshape(1, len(galfile))

        if galfile.shape[1] == 0:
            continue

        z_vector= galfile[:,2]                           #redshift column
        ra   = galfile[:,4]        #ra in deg
        dec = galfile[:,5]        #dec in deg

        if spec_type == 'observed':
            outspecfile=cube_name.replace('.fits','-extractspec.txt')
            outfreqdeffile=cube_name.replace('.fits','-freqdef.txt')
            spec= np.loadtxt(dirname+'/'+outspecfile)                   #array of observed spectra in cube;  shape is (number of coords used, number of channels)
            freq= np.loadtxt(dirname+'/'+outfreqdeffile)                #array of observed frequency in the cube; shape is (number of channels)
            nanvalues = np.where(spec == -99.0)
            spec[nanvalues] = np.nan
            fig_title= 'Observed spectrum of galaxy'
            fig_y = 'Observed '

        elif spec_type == 'rest-frame':
            outspeccube=cube_name.replace('.fits','-speccube.txt')
            outresfreqcube= cube_name.replace('.fits','-restfreqcube.txt')
            spec= np.loadtxt(dirname+'/'+outspeccube)                   #array of rest-frame regridded spectra in cube;  shape is (number of coords used, number of channels)
            freq= np.loadtxt(dirname+'/'+outresfreqcube)                #array of rest-frame frequency in the cube; shape is (number of channels)
            nanvalues = np.where(spec == -99.0)
            spec[nanvalues] = np.nan
            fig_title= 'Rest-frame spectrum of galaxy'
            fig_y = 'Rest-frame '

        elif spec_type == 'rest-frame-continuum-corrected':
            outspec_corr = cube_name.replace('.fits','-spec_corr.txt')
            outresfreqcube= cube_name.replace('.fits','-restfreqcube.txt')
            spec= np.loadtxt(dirname+'/'+outspec_corr)
            freq= np.loadtxt(dirname+'/'+outresfreqcube)
            nanvalues = np.where(spec == -99.0)
            spec[nanvalues] = np.nan
            fig_title= 'Rest-frame (corrected) spectrum of galaxy'
            fig_y = 'Rest-frame '

        if len(spec.shape) == 1:
            spec = spec.reshape(1, len(spec))

        minspec=c *np.nanmin(np.vstack([d for d in spec]))
        maxspec=c *np.nanmax(np.vstack([d for d in spec]))
        for i in range(len(galfile)):
            if spec_type == 'observed':
                lines_freq_final =list(np.array(lines_freq)/(z_vector[i]+1))

            elif spec_type == 'rest-frame' or spec_type == 'rest-frame-continuum-corrected':
                lines_freq_final = list(np.copy(np.array(lines_freq)))

            spec[i,:] = c * spec[i,:]
            arraynan=np.argwhere(np.isnan(spec[i,:]))
            if len(arraynan) != len(spec[i,:]):
                plt.clf()
                fig=plt.figure()
                ax=fig.add_subplot(111)
                ax.plot(freq,spec[i,:],label='Cube: %s \nRA (deg) = %.6f , DEC (deg) = %.6f , z = %0.3f' %(imagename,ra[i],dec[i],z_vector[i]))
                ax.plot([freq[0],freq[-1]],[0,0],'r-')

                x_axis = np.linspace(freq[0], freq[-1]+(freq[2]-freq[1]), len(freq)*1000)
                y_axis = np.interp(x_axis, freq, spec[i,:])

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    ax.fill_between(x_axis,y_axis , where=y_axis>0,alpha=0.6)
                    ax.fill_between(x_axis,y_axis , where=y_axis<= 0,alpha=0.6)

                spec_no_nan_ind= np.argwhere(~np.isnan(spec[i,:]))
                ax.set_ylabel('Flux [%sJy/beam]'%unit_text)
                ax.set_xlabel(fig_y+'Frequency [GHz]')
                ax.set_title(str(fig_title)+' %i ' %i)

                ax.set_ylim(1.1*minspec,maxspec*2)
                plt.gcf().subplots_adjust(top=0.85,bottom=0.25)
                plt.legend(bbox_to_anchor=(0,-0.57, 1, 0.4 ), loc="upper left", fancybox=True,mode='expand',borderaxespad=0, prop={'size': 12})
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    lineid_plot.plot_line_ids(freq,spec[i,:],lines_freq_final,lines_label,arrow_tip=maxspec*1.5,ax=ax)
                
                fig.savefig(dirname+'/' +str(cube_name)+'id'+str(i)+str(figname))
                plt.close(fig)
    return


def plot_spectra():
    global lines_freq
    global lines_label

    stacked_frequency    = np.loadtxt(dirname+'/stacked-frequency.txt')
    stacked_spec_median  = np.loadtxt(dirname+'/stacked-medianspec.txt')
    stacked_spec_mean    = np.loadtxt(dirname+'/stacked-meanspec.txt')
    stacked_spec_std     = np.loadtxt(dirname+'/stacked-stdspec.txt')
    stacked_galaxies     = np.loadtxt(dirname+'/stacked-nobs.txt')
    
    stacked_spec_median[stacked_spec_median==-99.0] = np.nan      
    stacked_spec_mean[stacked_spec_mean==-99.0] = np.nan
    stacked_spec_std[stacked_spec_std==-99] = np.nan
    stacked_galaxies[stacked_galaxies==0] = np.nan

    lines_freq =[115271.202,230538.0,345795.99,461040.768,576267.931,691473.076,806651.801,921799.704,1036912.385,1151985.443,492160.651,809341.97,1900536.9,52033.143,987926.759,1097364.790,1162911.60200,1207638.73,1228788.7190,1410618.069,88631.60,265886.43,354505.47,443116.14,531716.35,708877.00,797433.26,885970.69,110201.35,220398.86,330587.97,440765.17,550926.29,661067.28,771184.13,881272.81,89188.53,178375.07,267557.63,356734.24,445908.91,535061.64,624208.46,713341.37,2060068.86,3393006.24,5785879.59,48990.96,97980.95,146969.03,195954.21,244935.56,293912.09,342882.86,391846.89,440803.24,489750.93,538689.00]        
    lines_freq= list(np.array(lines_freq)/1000.)
    lines_label=['CO1-0','CO2-1','CO3-2','CO4-3','CO5-4','CO6-5','CO7-6','CO8-7','CO9-8','CO10-9','CI 3P1-3P0','CI 3P2-3P1','CII','H2O 2(1,1)-2(0,2)','H2O 2(0,1)-1(1,1)','H2O 3(1,2)-3(0,3)','H2O 3(2,1)-3(1,2)','H20 4(2,2)-4(1,3)','H2O 2(2,0)-2(1,1)','H20 5(2,3)-5(1,4)','HCN 1-0','HCN 3-2','HCN 4-3','HCN 5-4','HCN 6-5','HCN 8-7','HCN 9-8','HCN 10-9','13CO 1-0','13CO 2-1','13CO 3-2','13CO 4-3','13CO 5-4','13CO 6-5','13CO 7-6','13CO 8-7','HCO+ 1-0','HCO+ 2-1','HCO+ 3-2','HCO+ 4-3','HCO+ 5-4','HCO+ 6-5','HCO+ 7-6','HCO+ 8-7','OI','OIII 3P1-3P0','OIII 3P2-3P1','CS 1-0','CS 2-1','CS 3-2','CS 4-3','CS 5-4','CS 6-5','CS 7-6','CS 8-7','CS 9-8','CS 10-9','CS 11-10']
    
    plots_style()
    unit = 'mJy'

    if verbose >= 1:
        print('\n Plotting spectral stacking results ')

    ###################################################
    #Full stacked spectra (mean and median)
    if verbose >= 2:
        print('Plot: full stacked spectra')

    figname = dirname+'/full_stacked_spectra_mean.png'
    title='\nStacked spectra (mean)' 
    everyplot(stacked_spec_mean,stacked_spec_std,stacked_galaxies,stacked_frequency,title,figname,unit=unit)

    # figname = dirname+'/full_stacked_spectra_median.png'
    # title='Stacked spectra (median)'
    # everyplot(stacked_spec_median,stacked_spec_std,stacked_galaxies,stacked_frequency,title,figname,unit=unit)
    ####################################################
    #Individual spectra of stacked galaxies (observed, restframe or corrected restframe)
    if verbose >= 2:
        print('Plot: individual spectra of stacked galaxies')

    figname = 'observed_spec.png'
    plot_every_galaxy_observed(lines_freq, lines_label,spec_type = 'observed', figname = figname,unit=unit)
    # figname = 'restframe_spec.png'
    # plot_every_galaxy_observed(lines_freq, lines_label,spec_type = 'rest-frame', figname = figname,unit=unit)
    # figname = 'rest-frame-continuum-corrected.png'
    # plot_every_galaxy_observed(lines_freq, lines_label,spec_type = 'rest-frame-continuum-corrected', figname = figname,unit=unit)
    ##################################################### 
    #Stacked spectrum every x GHz (mean and median)
    freq1 = np.nanmin(stacked_frequency)
    freq2 = np.nanmax(stacked_frequency)
    #freq1, freq2= 327, 347
    step = 50

    if verbose >= 2:
        print('Plot: stacked spectrum every %s GHz'%str(step))

    figname=dirname+'/stack_mean_'+str(step)+'GHz'
    title=' \nStacked spectra (mean)'
    plot_range_freq(stacked_spec_mean,stacked_frequency,stacked_spec_std,stacked_galaxies,figname, title,startfreq = freq1, endfreq = freq2,freqstep=step)
 
    # figname=dirname+'/stack_median_'+str(step)+'GHz'
    # title=' \nStacked spectra (median)'
    # plot_range_freq(stacked_spec_median,stacked_frequency,stacked_spec_std,stacked_galaxies,figname, title,startfreq = freq1, endfreq = freq2,freqstep=step)
    ####################################################
    #Stacked spectrum centered in common emission lines (mean and median)
    if verbose >= 2:
        print('Plot: common emission lines')

    figname=dirname+'/line_stack_mean'
    title=' \nStacked spectra (mean)'
    plot_every_common_line_fit(stacked_spec_mean,stacked_frequency,stacked_spec_std,stacked_galaxies,figname, title, velocity_width= 2000.,unit=unit)
    
    # figname=dirname+'/line_stack_median'
    # title=' \nStacked spectra (median)'
    # plot_every_common_line_fit(stacked_spec_median,stacked_frequency,stacked_spec_std,stacked_galaxies,figname, title, velocity_width= 2000.,unit=unit)
    ####################################
    #Emission line plot with Gaussian fit 
    if verbose >= 2:
        print('Plot: emission line fit')

    figname=dirname+'/stacked_mean_fit_line_'
    title=' \nStacked spectra (mean)'
    speclinesnr('CO3-2',stacked_frequency,stacked_spec_mean,stacked_spec_median,stacked_spec_std,stacked_galaxies, figname, title, makeplot = 1,unit=unit)
    #speclinesnr('CO4-3',stacked_frequency,stacked_spec_mean,stacked_spec_median,stacked_spec_std,stacked_galaxies, figname, title, makeplot = 1,unit=unit)

    # figname=dirname+'/stacked_median_fit_line_'
    # title=' \nStacked spectra (median)'
    # speclinesnr('CO3-2',stacked_frequency,stacked_spec_median,stacked_spec_mean,stacked_spec_std,stacked_galaxies, figname, title, makeplot = 1,unit=unit)
    # speclinesnr('CO4-3',stacked_frequency,stacked_spec_median,stacked_spec_mean,stacked_spec_std,stacked_galaxies, figname, title, makeplot = 1,unit=unit)
    return 
    
plt.ioff()
#################################################################

dirname='Subsamples/detections-spec-all/new'
incubes = [
'Abell370/spec/cube1_60km_dirty_nat.fits','Abell370/spec/cube2_60km_dirty_nat.fits',
'Abell2744/spec/cube1_60km_dirty_nat.fits','Abell2744/spec/cube2_60km_dirty_nat.fits',
'MACS1931.8-2635/spec/cube1_60km_dirty_nat.fits','MACS1931.8-2635/spec/cube2_60km_dirty_nat.fits']
incat='catalog.csv'
stampsize = 2
verbose = 1

if os.path.isdir(dirname) == False:
    print('\n Output folder not found. Enter correct path.')
else:
    if os.path.isfile(dirname+'/'+incat) == False:
        print('\n Catalogue not found. Enter correct path.')
    else:
        cubes=[]
        for file in incubes:
            cubes.append(os.path.isfile(file))
            if os.path.isfile(file) == False:
                print('\n File {} not found. Enter correct path.'.format(file))

        if all(cubes):
            print('Working in folder : {}'.format(dirname))
            coords,stacked_spec = stack(catalog_name = incat, stampsize = stampsize , cubenames = incubes, overwrite = True, weightstack = False, verbose = verbose)
            plot_spectra()
            print('Done')
