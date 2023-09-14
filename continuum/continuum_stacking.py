# This Python file uses the following encoding: utf-8
import os
import csv
import taskinit
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Ellipse
import warnings
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

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
        return

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


    def __init__(self, x, y, weight=1., ra=0., dec=0.,image=0):
        """
            Create a coordinate. A pixel coordinate should always specify
            to which image it belongs. Physical coordinates should be in
            J2000 radians.
        """
        self.x = x
        self.y = y
        self.weight = weight
        self.ra = ra
        self.dec = dec
        self.image = image

    def __str__(self):
        return '{0}, {1}, {2},{3},{4},{5}'.format(self.x, self.y , self.weight, self.ra, self.dec,self.image)


def readCoords(coordfile,weightstack,unit='deg'):
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
    coordreader = csv.reader(open(coordfile, 'rb'), delimiter=',')
    coords = CoordList()
    for row in coordreader:
        x = float(row[0])
        y = float(row[1])
        ra = x    
        dec= y

        if unit == 'deg':  #to radians
            x = np.pi/180.*x
            y = np.pi/180.*y

        if weightstack == 'True':
            weight = float(row[2])
        else:
            weight = 1.

        if x > 2*np.pi:
            x -= 2*np.pi
        if y > np.pi:
            y -= 2*np.pi

        coords.append(Coord( x, y, weight, ra, dec))

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
    sizeimagex = (Nx*x_pix_inc)   # radians
    sizeimagey = (Ny*y_pix_inc)

    pixcoords = []
    for coord in coords:   
        if (abs(coord.x - cenimagera) < 1.2*abs(sizeimagex)) and (abs(coord.y - cenimagedec) < 1.2*abs(sizeimagey)):
            dx = (coord.x - x0)*np.cos(coord.y)
            dy = np.arcsin(np.sin(coord.y)/np.cos(dx)) - y0
            x = dx/x_pix_inc+x_pix_ref
            y = dy/y_pix_inc+y_pix_ref
            if (x > 0 and x < Nx-1) and (y > 0 and y < Ny-1):
                c = Coord(x, y, coord.weight, coord.ra, coord.dec, coord.image)
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


def getImages(coords, overwrite):
    """
    Gets list of images names to stack.
    If overwrite is True, imagenames will be all images, if not,
    it will be only the images missing from the output directory

    """
    imagenames = []
    if overwrite == False:
        for imagename in coords.imagenames:
            outspecfile = imagename.replace('.fits', '-extractspec.txt')
            if not os.path.exists(dirname + '/' + outspecfile):
                imagenames.append(imagename)

    if overwrite == True:
        imagenames = coords.imagenames

    return imagenames


def stackcont(catalog_name='',stampsize=10,images=[],weightstack=False,threshold=5,overwrite=False,verbose=1):
    """
        Performs continuum stacking in the image domain.
        Inputs:
            catalog_name  -- Name of the csv catalog which has the coordinates of the targets (ra[deg] and dec[deg] columns required).
            stampsize     -- Size of target region in arcseconds that will be extracted for stacking.
            images        -- Name of images to extract continuum from.
            weightstack   -- If True, takes third column of csv catalog and uses it as weight column for stacking. 
            threshold     -- Threshold of signal-to-noise ratio to consider stack image as detected or undetected. 
            overwrite     -- When this code is executed, files of the continuum emission from each input image are saved. If overwrite = False,
                             the code will only save these files for new images added. 
                             If overwrite = True, the code will save these files for all input images. 
            verbose       -- Verbose level from 0 to 3, where 0 provides minimal amount of prints and 3 provides more detailed messages.
        Outputs:
            Stacked continuum  (and writes out other text and data files)
    """

    coords = readCoords(catalog_name,weightstack)

    if verbose >= 1:
        print('Total number of objects in catalog: %s'%len(coords))
        print('Images to be used for continuum stacking: %s'%images)
        print('Stamp size [arcsec]: %s'%stampsize)
    
    if coords.coord_type == 'physical':
        coords = getPixelCoords(coords, images)

    imagenames = getImages(coords, overwrite)

    if verbose >= 1:
        print('\n Starting continuum stacking extraction')

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

        new_stampsize = int((stampsize * np.pi/(180.*3600.))/y_pix_inc) #arcsec to rad to pix
        _load_stackcont(imagename,new_stampsize,imcoords,verbose,weightstack)
        index=index+1

    if verbose >= 1:
        print('\n Continuum extraction finished')
        print('\n Starting continuum stacking')

    stacked_cont= _stack_stackcont(coords,catalog_name,threshold, verbose)
    
    if verbose >= 1:
        print('\n Stacking continuum completed')

    return coords,stacked_cont


def _load_stackcont(imagename,new_stampsize,coords,verbose,weightstack): 
    """
    Generates subimages that will be stacked for each image. Each subimage is the size of the stampsize, 
    with the center coordinate from the coordinate file. This step discards subimages with more than 25% 
    of pixels with nan values. 
    Then, each subimage is saved in a fits file with name imagename-extracted_sources_cube.fits.
    Two text files are saved:
    1. imagename-galaxyinfo.txt, with ra(rad), dec(rad), weight, ra(deg) and dec(deg) for sources found in the image
    2. imagename-notfound.txt, that saves the same information but sources from the catalog not found in the image
    For this, a source is considered found, if the coordinate is inside the spatial range of the image and if that range
    includes more than 75% pixels non nan.  
    """
    if verbose >= 2:
        print('Reading image: %s'%imagename)
        
    if verbose >= 3:
        print('The coordinates found in this image are: %s'%coords)

    hdu = fits.open(imagename)[0]
    mom0 = hdu.data                          #shape is (number of polarizations, number of channels, size in y, size in x)
    mom0=np.squeeze(mom0)                    #shape is (size in y, size in x)

    datacont = np.zeros((len(coords),new_stampsize, new_stampsize))
    imagesize = [(mom0.shape[1], mom0.shape[0])]

    if verbose >= 2:
        print('Extracting continuum stamp per position')

    galaxyinfo=[]
    notfound=[]
    count = 0
    for (i,coord) in enumerate(coords):
        if (coord.x>0 and coord.x<imagesize[0][0]-1) and (coord.y>0 and coord.y<imagesize[0][1]-1):
            blcx = int(coord.x - new_stampsize/2. + 0.5)
            blcy = int(coord.y - new_stampsize/2. + 0.5)
            trcx = blcx + int(new_stampsize)
            trcy = blcy + int(new_stampsize)
            if (blcx>0 and blcx<imagesize[0][0]-1) and (trcx>0 and trcx<imagesize[0][0]-1) and (blcy>0 and blcy<imagesize[0][1]-1) and (trcy>0 and trcy<imagesize[0][1]-1):
                if weightstack == True:
                      datacont[count,:,:] = mom0[blcy:trcy,blcx:trcx]*coord.weight
                else:
                    datacont[count,:,:] = mom0[blcy:trcy,blcx:trcx]
                count = count + 1
                galaxyinfo.append([coord.x,coord.y,coord.weight,coord.ra,coord.dec])
            
            else:
                data_image = np.pad(mom0, (new_stampsize, new_stampsize), 'constant', constant_values=0)
                
                new_x, new_y = coord.x+new_stampsize, coord.y+new_stampsize
                blcx = int(new_x - new_stampsize/2. + 0.5)
                blcy = int(new_y - new_stampsize/2. + 0.5)
                trcx = blcx + int(new_stampsize)
                trcy = blcy + int(new_stampsize)

                if weightstack == True:
                    datacont[count,:,:] = data_image[blcy:trcy,blcx:trcx]*coord.weight
                else:
                    datacont[count,:,:] = data_image[blcy:trcy,blcx:trcx]
                count = count + 1
                galaxyinfo.append([coord.x,coord.y,coord.weight,coord.ra,coord.dec])
        else:
            notfound.append([coord.x,coord.y,coord.weight,coord.ra,coord.dec])

    if verbose >= 3:
        print('Number of regions not found in the image: %s'%len(notfound))

    datacont=datacont[0:count,:,:]

    nan=[]
    for i in range(datacont.shape[0]):
        image_galaxy = datacont[i]
        total_pixels = len(image_galaxy[0])* len(image_galaxy[1])
        nan_pixels = np.count_nonzero(np.isnan(image_galaxy))
        ratio = float(nan_pixels)/float(total_pixels)
        if ratio  > 0.25:
            nan.append(i)
    
    if verbose >= 3:
        print('Number of regions discarded (>25% of nan values): %s'%len(nan))

    datacont =  np.delete(datacont, np.array(nan), axis=0)
    galaxyinfo = np.delete(galaxyinfo, np.array(nan), axis =  0)

    countnonnan=0
    for i in range(datacont.shape[0]):
        if not np.isnan(np.nansum(datacont[i,:,:])):
            countnonnan=countnonnan+1

    cube_name = imagename.replace("/","-")

    outcontname=cube_name.replace('.fits','-extracted_sources_cube.fits')
    outgalfile=cube_name.replace('.fits','-galaxyinfo.txt')
    outnotfound=cube_name.replace('.fits','-notfound.txt')

    np.savetxt(dirname+'/'+outgalfile,galaxyinfo)
    np.savetxt(dirname+'/'+outnotfound,notfound)

    with fits.open(imagename) as hdus:
        headermain = hdus[0].header

    hduexp = fits.PrimaryHDU(datacont, header = headermain)
    hdulist = fits.HDUList([hduexp])
    hdulist.writeto(dirname+'/'+outcontname,overwrite = True)

    if verbose >= 2:
        print('Number of objects with non nan flux in image: %s'%countnonnan)

    if verbose >= 2:
        print('Saving files of continuum emission, galaxy information and galaxies not found')
        
    return

def _stack_stackcont(coords, catalog_name,threshold, verbose):        
    """
    Reads extracted continuum saved in files and performs the stack.
    """

    imagenames_with_sources = []
    for imagename in coords.imagenames:
        cube_name = imagename.replace("/","-")
        outcontname=cube_name.replace('.fits','-extracted_sources_cube.fits')
        outgalfile=cube_name.replace('.fits','-galaxyinfo.txt')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            galaxyinfo = np.loadtxt(dirname+'/'+outgalfile)

        if len(galaxyinfo) > 0:
            imagenames_with_sources.append(imagename)

    if verbose >= 2:
        print('\n Reading each image to save image size')

    larger_cube_size = 0
    larger_cube = 0
    for imagename in imagenames_with_sources:
        cube_name = imagename.replace("/","-")
        outcontname=cube_name.replace('.fits','-extracted_sources_cube.fits')
        hdu = fits.open(dirname+'/'+outcontname)[0]
        sourcescube = hdu.data                                  

        if larger_cube_size < sourcescube.shape[1]:
            larger_cube_size = sourcescube.shape[1]
            larger_cube = sourcescube

    all_sources_cube = np.empty(shape=(0,larger_cube_size,larger_cube_size))
    all_sources_galaxyinfo = np.empty(shape=(0, 5))
    pos_tot = 0
    pos_files = 0
    for imagename in imagenames_with_sources:
        cube_name = imagename.replace("/","-")
        outcontname=cube_name.replace('.fits','-extracted_sources_cube.fits')
        outgalfile =cube_name.replace('.fits','-galaxyinfo.txt')
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            galaxyinfo = np.loadtxt(dirname+'/'+outgalfile)

        if galaxyinfo.shape == 1:
            galaxyinfo = galaxyinfo.reshape(1, len(galaxyinfo))

        hdu = fits.open(dirname+'/'+outcontname)[0]
        sourcescube = hdu.data
        npos= sourcescube.shape[0]
        larger_cube_points = np.meshgrid(larger_cube.shape[2],larger_cube.shape[1])

        regridded_sources_cube=np.zeros((npos, larger_cube_size, larger_cube_size))

        if sourcescube.shape[1] < larger_cube_size :
            for pos in range(sourcescube.shape[0]):                                     
                points = np.meshgrid(sourcescube[2],sourcescube[1])
                regridded_sources_cube[pos,:,:]=griddata(points,sourcescube[pos,:,:], larger_cube_points,method='linear' )
        else:
            regridded_sources_cube=np.copy(sourcescube)

        all_sources_cube= np.append(all_sources_cube, regridded_sources_cube, axis=0)
        all_sources_galaxyinfo = np.vstack([all_sources_galaxyinfo, galaxyinfo])

        if verbose >= 3:
            print('Adding %i sources from image %s to the stack'%(sourcescube.shape[0],imagename))

        pos_tot = pos_tot + npos
        pos_files = pos_files + len(galaxyinfo)

    if verbose >= 1:
        print('\n The final cube for stacking has %i regions'%all_sources_cube.shape[0])

    if verbose >= 2:
        print('Performing continuum stack')

    mean_stack = np.nanmean(all_sources_cube, axis=0)
    median_stack = np.nanmedian(all_sources_cube, axis=0)

    if verbose >= 2:
        print('\n Saving continuum stacking results')

    with fits.open(imagename) as hdus:
        headermain = hdus[0].header

    hdu_all_sources_cube = fits.PrimaryHDU(all_sources_cube,  header = headermain )
    hdulist_all_sources_cube = fits.HDUList([hdu_all_sources_cube])
    hdulist_all_sources_cube.writeto(dirname+'/all_sources_cube.fits',overwrite = True)

    hdu_mean_stack = fits.PrimaryHDU(mean_stack, header = headermain )
    hdulist_mean_stack = fits.HDUList([hdu_mean_stack])
    hdulist_mean_stack.writeto(dirname+'/mean_stack.fits',overwrite = True)

    hdu_median_stack = fits.PrimaryHDU(median_stack, header = headermain )
    hdulist_median_stack = fits.HDUList([hdu_median_stack])
    hdulist_median_stack.writeto(dirname+'/median_stack.fits',overwrite = True)

    np.savetxt(dirname + '/all_sources.txt', all_sources_galaxyinfo)
    
    sources_cube_edge= np.copy(all_sources_cube)
    sources_file_edge = np.copy(all_sources_galaxyinfo)

    sources_cube_edge[:,int(0.3*larger_cube_size):int(0.7*larger_cube_size),int(0.3*larger_cube_size):int(0.7*larger_cube_size)]  = np.nan
    nondetindex=[]
    detindex=[]
    edgedetindex= []

    for i in range(all_sources_cube.shape[0]):
        sources_cube_center = all_sources_cube[i,int(0.35*larger_cube_size):int(0.65*larger_cube_size),int(0.35*larger_cube_size):int(0.65*larger_cube_size)]

        centerpeak = np.nanmax(sources_cube_center)
        centerflux = np.nansum(sources_cube_center)     #totalflux
        centerstd  = np.nanstd(sources_cube_center)

        tsources_cube_edge= sources_cube_edge[i,:]
        edgeflux  = np.nansum(tsources_cube_edge)
        edgestd   = np.nanstd(tsources_cube_edge)
        edgepeak = np.nanmax(tsources_cube_edge)

        if ((centerpeak > threshold*edgestd)):
            detindex.append(i)
        elif (edgepeak > threshold*edgestd):
            edgedetindex.append(i)
        else:
            nondetindex.append(i)

    sources_cube_edge_det = np.delete(all_sources_cube, np.array(detindex[:]+nondetindex[:]), axis=0)
    sources_cube_det    = np.delete(all_sources_cube, np.array(nondetindex[:]+edgedetindex[:]), axis=0)
    sources_cube_nondet = np.delete(all_sources_cube, np.array(detindex[:]+edgedetindex[:]), axis=0)

    sources_file_edge_det = np.delete(all_sources_galaxyinfo, np.array(detindex[:]+nondetindex[:]), axis=0)
    sources_file_det    = np.delete(all_sources_galaxyinfo, np.array(nondetindex[:]+edgedetindex[:]), axis=0)
    sources_file_nondet = np.delete(all_sources_galaxyinfo, np.array(detindex[:]+edgedetindex[:]), axis=0)

    step= 0

    for i in range(sources_cube_nondet.shape[0]):
        sources_cube_nondet[i] = np.rot90(sources_cube_nondet[i], k= step)
        step = step + 1
        if step == 4:
            step = 0

    with fits.open(imagename) as hdus:
        headermain = hdus[0].header

    if sources_cube_det.shape[0] > 0 :
        hdu_sources_cube_det = fits.PrimaryHDU(sources_cube_det,  header = headermain )
        hdulist_sources_cube_det = fits.HDUList([hdu_sources_cube_det])
        hdulist_sources_cube_det.writeto(dirname+'/'+'sources_cube_det.fits',overwrite = True)

        mean_stack_det = np.nanmean(sources_cube_det, axis=0)
        median_stack_det = np.nanmedian(sources_cube_det.astype(np.float64), axis=0)
        std_stack_det = np.nanstd(sources_cube_det, axis=0)

        np.savetxt(dirname + '/' + 'sources_det.txt', sources_file_det)

        ############Continuum stacking of images with detections - Mean
        hdu_mean_stackdet = fits.PrimaryHDU(mean_stack_det, header = headermain)
        hdulist_mean_stackdet = fits.HDUList([hdu_mean_stackdet])
        hdulist_mean_stackdet.writeto(dirname+'/'+'mean_stack_det.fits',overwrite = True)

        ############Continuum stacking of images with detections - median
        hdu_median_stackdet = fits.PrimaryHDU(median_stack_det, header = headermain)
        hdulist_median_stackdet = fits.HDUList([hdu_median_stackdet])
        hdulist_median_stackdet.writeto(dirname+'/'+'median_stack_det.fits',overwrite = True)

        ############Continuum stacking of images with detections - STD
        hdu_std_stackdet = fits.PrimaryHDU(std_stack_det, header = headermain)
        hdulist_std_stackdet = fits.HDUList([hdu_std_stackdet])
        hdulist_std_stackdet.writeto(dirname+'/'+'std_stack_det.fits',overwrite = True)

    #If non-detections are found, they are stored in a cube and then stacked.
    if sources_cube_nondet.shape[0] > 0 :
        hdu_sources_cube_nondet = fits.PrimaryHDU(sources_cube_nondet,  header = headermain )
        hdulist_sources_cube_nondet = fits.HDUList([hdu_sources_cube_nondet])
        hdulist_sources_cube_nondet.writeto(dirname+'/'+'sources_cube_nondet.fits',overwrite = True)

        np.savetxt(dirname + '/' + 'sources_non_det.txt', sources_file_nondet)

        mean_stack_nondet= np.nanmean(sources_cube_nondet, axis=0)
        median_stack_nondet= np.nanmedian(sources_cube_nondet, axis=0)
        std_stack_nondet= np.nanstd(sources_cube_nondet, axis=0)

        ############Continuum stacking of images without detections - Mean
        hdu_mean_stacknondet= fits.PrimaryHDU(mean_stack_nondet, header = headermain)
        hdulist_mean_stacknondet = fits.HDUList([hdu_mean_stacknondet])
        hdulist_mean_stacknondet.writeto(dirname+'/'+'mean_stack_nondet.fits',overwrite = True)

        ############Continuum stacking of images without detections - median
        hdu_median_stacknondet = fits.PrimaryHDU(median_stack_nondet, header = headermain)
        hdulist_median_stacknondet = fits.HDUList([hdu_median_stacknondet])
        hdulist_median_stacknondet.writeto(dirname+'/'+'median_stack_nondet.fits',overwrite = True)

        ############Continuum stacking of images without detections - STD
        hdu_std_stacknondet = fits.PrimaryHDU(std_stack_nondet, header = headermain)
        hdulist_std_stacknondet = fits.HDUList([hdu_std_stacknondet])
        hdulist_std_stacknondet.writeto(dirname+'/'+'std_stack_nondet.fits',overwrite = True)

    if sources_cube_edge_det.shape[0] > 0 :
        hdu_sources_cube_edge_det = fits.PrimaryHDU(sources_cube_edge_det,  header = headermain )
        hdulist_sources_cube_edge_det = fits.HDUList([hdu_sources_cube_edge_det])
        hdulist_sources_cube_edge_det.writeto(dirname+'/'+'sources_cube_edge_det.fits',overwrite = True)

        np.savetxt(dirname + '/' + 'sources_edge_det.txt', sources_file_edge_det)

    if verbose >= 1:
        print('Number of individual detections : %i'%len(detindex))
        print('Number of non detection : %i'%len(nondetindex))
        print('Number of individual offset detections : %i'%len(edgedetindex))
    return

def three_panel_fig(stack_image_file, stack_cube_file, figname , unit='mJy'):
    stack_cube      = fits.open(dirname+'/'+stack_cube_file)
    stack_cube_data = stack_cube[0].data

    stack_image = fits.open(dirname+'/'+stack_image_file)
    stack_image_data   = stack_image[0].data
    stack_image_header = stack_image[0].header

    image_beammajor = stack_image_header['BMAJ']                  #degree/pix
    image_beamminor = stack_image_header['BMIN']                  #degree/pix
    image_beampa    = stack_image_header['BPA']
    image_cdelt     = stack_image_header['CDELT2']                #degree/pix

    if unit == 'mJy':
        c = 1.e3
        unit_text='m'
    if unit == 'Jy':
        c = 1.
        unit_text=' '
    if unit== 'muJy':
        c = 1.e6
        unit_text = '$\mu$'

    cube_width = stack_cube_data.shape[1]
    cube_edge= np.copy(stack_cube_data)
    cube_edge[:,int(0.3*cube_width):int(0.7*cube_width),int(0.3*cube_width):int(0.7*cube_width)]  = np.nan

    if stack_cube_data.shape[0] == 1:
        cube_std  = np.nanmean(np.nanstd(cube_edge))
    else:
        cube_std  = np.nanmean(np.nanstd(cube_edge, axis = (1,2)))

    cube_std = c*cube_std
    stack_image_data = c*stack_image_data

    im_width = len(stack_image_data[0])
    center_im = stack_image_data[int(0.3*im_width):int(0.7*im_width),int(0.3*im_width):int(0.7*im_width)]

    edge_im= np.copy(stack_image_data)
    edge_im[int(0.3*im_width):int(0.7*im_width),int(0.3*im_width):int(0.7*im_width)]  = np.nan

    fig, axes = plt.subplots(1, 3, figsize=(15, 6), dpi=80)
    fig.subplots_adjust(left=0.1, bottom=0.04, right=0.93, top=0.98, wspace=0.39)
    s = 0
    i = 0
    for ax in axes.flat:
        t = 3*s
        image_data_f = gaussian_filter(stack_image_data, sigma = s, truncate = t)
        center_im_filt = image_data_f[int(0.3*im_width):int(0.7*im_width),int(0.3*im_width):int(0.7*im_width)]
        edge_im_filt = np.copy(image_data_f)
        edge_im_filt[int(0.3*im_width):int(0.7*im_width),int(0.3*im_width):int(0.7*im_width)] = np.nan
        arcsecdata = (image_data_f.shape[1] * image_cdelt ) *3600.
        im = ax.imshow(image_data_f,origin='lower',extent=[-arcsecdata/2., arcsecdata/2., -arcsecdata/2., arcsecdata/2. ])
        cb=fig.colorbar(im, ax = ax, shrink = 0.45)
        cb.set_label('Flux [%sJy/beam]' %(unit_text))
        ax.set_title('Smoothing by %i pixels' '\n' '$\sigma_{stack}$ = %.3f %sJy ' %(s, np.nanstd(edge_im_filt),unit_text), fontsize = 12)
        ellipse = Ellipse(xy=(-arcsecdata/2.*0.8,-arcsecdata/2.*0.8), width=image_beamminor*3600., height=image_beammajor*3600., edgecolor='w', fc='None', lw=1,angle = image_beampa)
        ax.add_patch(ellipse)
        s = s + 3
        i = i+1
    #######################################################################
    axes[0].set_ylabel(r'$\Delta$Dec [arcsec]')
    axes[0].set_xlabel(r'$\Delta$RA [arcsec]')
    axes[0].set_title('No smoothing' '\n' '$\sigma_{stack}$ = %.3f %sJy' %(np.nanstd(edge_im), unit_text), fontsize = 12)
    axes[1].set_xlabel(r'$\Delta$RA [arcsec]')
    axes[2].set_xlabel(r'$\Delta$RA [arcsec]')

    figtitle = fig.suptitle(r'N = %i, $\bar{\sigma}_{all}$ = %.1f %sJy'  %(stack_cube_data.shape[0],cube_std,unit_text), fontsize=20, y=1.05)
    figtitle.set_position([.5,0.95])
    plt.savefig(dirname+'/'+figname)
    #plt.close()
    return

def stack_fig(stack_image_file, stack_cube_file, figname , unit='mJy', title = ''):
    lsize = 24
    tlsize = 20
    plt.rcParams['xtick.labelsize'] = tlsize
    plt.rcParams['ytick.labelsize'] = tlsize
    plt.rcParams['axes.labelsize'] = lsize

    stack_cube      = fits.open(dirname+'/'+stack_cube_file)
    stack_cube_data = stack_cube[0].data

    stack_image = fits.open(dirname+'/'+stack_image_file)
    stack_image_data   = stack_image[0].data
    stack_image_header = stack_image[0].header

    image_beammajor = stack_image_header['BMAJ']                  #degree/pix
    image_beamminor = stack_image_header['BMIN']                  #degree/pix
    image_beampa    = stack_image_header['BPA']
    image_cdelt     = stack_image_header['CDELT2']                #degree/pix

    if unit == 'mJy':
        c = 1.e3
        unit_text='m'
    if unit == 'Jy':
        c = 1.
        unit_text=' '
    if unit== 'muJy':
        c = 1.e6
        unit_text = '$\mu$'

    stack_image_data = c*stack_image_data
    cube_width = stack_cube_data.shape[1]
    cube_edge= np.copy(stack_cube_data)
    cube_edge[:,int(0.3*cube_width):int(0.7*cube_width),int(0.3*cube_width):int(0.7*cube_width)]  = np.nan

    if stack_cube_data.shape[0] == 1:
        cube_std  = np.nanmean(np.nanstd(cube_edge))
    else:
        cube_std  = np.nanmean(np.nanstd(cube_edge, axis = (1,2)))

    
    cube_std = c*cube_std

    im_width = len(stack_image_data[0])
    center_im = stack_image_data[int(0.3*im_width):int(0.7*im_width),int(0.3*im_width):int(0.7*im_width)]
    
    edge_im= np.copy(stack_image_data)
    edge_im[int(0.3*im_width):int(0.7*im_width),int(0.3*im_width):int(0.7*im_width)]  = np.nan

    fig = plt.figure(figsize = (8,6))
    arcsecdata = (stack_image_data.shape[1] * image_cdelt ) *3600.
    plt.imshow(stack_image_data,origin='lower', extent=[-arcsecdata/2. , arcsecdata/2., -arcsecdata/2., arcsecdata/2. ])
    ellipse = Ellipse(xy=(-arcsecdata/2.*0.8,-arcsecdata/2.*0.8), width=image_beamminor*3600., height=image_beammajor*3600., edgecolor='k', fc='None', lw=4,angle = image_beampa)
    plt.plot([0],[0], 'k', marker = 'x' , scalex= False, scaley=False, markersize = 14, markeredgewidth = 5)
    ax = plt.gca()
    
    ax.add_patch(ellipse)
    cb=plt.colorbar()
    cb.set_label('Flux [%sJy/beam]' %(unit_text), fontsize = lsize)
    plt.xlabel(r'$\Delta$RA [arcsec]')
    plt.ylabel(r'$\Delta$Dec [arcsec]')
    sigma = np.nanstd(edge_im)
    plt.plot([-2.5,-2.5,2.5,2.5,-2.5],[-2.5,2.5,2.5,-2.5,-2.5],'-k')
    figtitle = plt.title('%s - N = %i' '\n' r' $\bar{\sigma}_{all}$ = %.1f %sJy, $\sigma_{stack}$ = %.1f %sJy'  %(title, stack_cube_data.shape[0],cube_std,unit_text,np.nanstd(edge_im),unit_text), fontsize=32)
    figtitle.set_position([.5,1.03])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.contour(stack_image_data, levels=[-3*sigma,-2*sigma,2*sigma,3*sigma,4*sigma,5*sigma], colors='black',origin='lower',extent=[-arcsecdata/2. , arcsecdata/2., -arcsecdata/2., arcsecdata/2. ])
        plt.savefig(dirname+'/'+figname,bbox_inches='tight')

    #plt.close()

def one_panel_fit(stack_image_file, stack_cube_file, figname ,unit='mJy', upperlimit=False):
    lsize = 22
    tlsize = 20
    plt.rcParams['xtick.labelsize'] = tlsize
    plt.rcParams['ytick.labelsize'] = tlsize
    plt.rcParams['axes.labelsize'] = lsize

    stack_cube      = fits.open(dirname+'/'+stack_cube_file)
    stack_cube_data = stack_cube[0].data

    stack_image = fits.open(dirname+'/'+stack_image_file)
    stack_image_data   = stack_image[0].data
    stack_image_header = stack_image[0].header

    image_beammajor = stack_image_header['BMAJ']                  #degree/pix
    image_beamminor = stack_image_header['BMIN']                  #degree/pix
    image_beampa    = stack_image_header['BPA']                 #degree
    image_cdelt     = stack_image_header['CDELT2']                #degree/pix

    if unit == 'mJy':
        c = 1.e3
        unit_text='m'
    if unit == 'Jy':
        c = 1.
        unit_text=' '
    if unit== 'muJy':
        c = 1.e6
        unit_text = '$\mu$'

    stack_image_data = c*stack_image_data
    cube_width = stack_cube_data.shape[1]
    cube_edge= np.copy(stack_cube_data)
    cube_edge[:,int(0.3*cube_width):int(0.7*cube_width),int(0.3*cube_width):int(0.7*cube_width)]  = np.nan

    if stack_cube_data.shape[0] == 1:
        cube_std  = np.nanmean(np.nanstd(cube_edge))
    else:
        cube_std  = np.nanmean(np.nanstd(cube_edge, axis = (1,2)))

    cube_std = c*cube_std

    im_width = len(stack_image_data[0])
    center_im = stack_image_data[int(0.3*im_width):int(0.7*im_width),int(0.3*im_width):int(0.7*im_width)]

    edge_im= np.copy(stack_image_data)
    edge_im[int(0.3*im_width):int(0.7*im_width),int(0.3*im_width):int(0.7*im_width)]  = np.nan

    fig = plt.figure(figsize = (8,6))
    arcsecdata = (stack_image_data.shape[1] * image_cdelt ) *3600.
    plt.imshow(stack_image_data,vmin=np.nanmin(stack_image_data), vmax=np.nanmax(stack_image_data),origin='lower', extent=[-arcsecdata/2. , arcsecdata/2., -arcsecdata/2., arcsecdata/2. ])
    ellipse = Ellipse(xy=(-arcsecdata/2.*0.8,-arcsecdata/2.*0.8), width=image_beamminor*3600., height=image_beammajor*3600., edgecolor='k', fc='None', lw=4,angle = image_beampa)
    plt.plot([0],[0], 'k', marker = 'x' , scalex= False, scaley=False, markersize = 14, markeredgewidth = 5)
    ax = plt.gca()
    ax.add_patch(ellipse)
    cb=plt.colorbar()
    cb.set_label('Flux [%sJy/beam]' %(unit_text), fontsize = lsize)
    plt.xlabel(r'$\Delta$RA [arcsec]')
    plt.ylabel(r'$\Delta$Dec [arcsec]')
    #######################################################################
    myia = iatool()
    myia.open(dirname+'/'+stack_image_file)

    blc1,trc1, blc2, trc2 = int(im_width/3.), int(im_width/3.), int(im_width/3.)*2, int(im_width/3.)*2
    box = str(blc1)+","+str(trc1)+","+str(blc2)+","+str(trc2)
    region = ""
    res = myia.fitcomponents(box=box, region =region)
    clist = res['results']
    total_flux = 0
    try:
       total_flux = c*clist['component0']['flux']['value'][0]
       total_flux_error = c*clist['component0']['flux']["error"][0]
       total_flux_error_mjy = 1e3*clist['component0']['flux']["error"][0]
       peak = c*clist['component0']["peak"]['value']
       peakerror = c*clist['component0']["peak"]['error']
       snr = total_flux / np.nanstd(edge_im)
       fit= 1
    except:
       fit = 0
       print('fit unsucessfull')

    if fit == 1:
        if upperlimit == True:
            figtitle = plt.title(r'N = %i, $\bar{\sigma}_{all}$ = %.1f %sJy, $\sigma_{stack}$ = %.1f %sJy' '\n' ' Upperlimit'  %(stack_cube_data.shape[0],cube_std,unit_text,np.nanstd(edge_im),unit_text), fontsize=24, y=1.05)
        if upperlimit == False:
            figtitle = plt.title( r'N = %i, $\bar{\sigma}_{all}$ = %.1f %sJy, $\sigma_{stack}$ = %.1f %sJy' '\n' 'Flux = $%.1f \pm %.1f$ %sJy, SNR = %.1f '  %(stack_cube_data.shape[0],cube_std,unit_text,np.nanstd(edge_im),unit_text, total_flux,total_flux_error, unit_text, snr), fontsize=24, y=1.05)

    if fit == 0:
        if upperlimit == True:
            figtitle = plt.title(r'N = %i, $\bar{\sigma}_{all}$ = %.1f %sJy, $\sigma_{stack}$ = %.1f %sJy' '\n' ' Upperlimit'  %(stack_cube_data.shape[0],cube_std,unit_text,np.nanstd(edge_im),unit_text), fontsize=24, y=1.05)
        if upperlimit == False:
            figtitle = plt.title( r'N = %i, $\bar{\sigma}_{all}$ = %.1f %sJy, $\sigma_{stack}$ = %.1f %sJy' '\n' 'Fit unsucessfull '  %(stack_cube_data.shape[0],cube_std,unit_text,np.nanstd(edge_im),unit_text), fontsize=24, y=1.05)

    sigma = np.nanstd(edge_im)
    plt.plot([-2.5,-2.5,2.5,2.5,-2.5],[-2.5,2.5,2.5,-2.5,-2.5],'-k')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.contour(stack_image_data, levels=[-3*sigma,-2*sigma,2*sigma,3*sigma,4*sigma,5*sigma], colors='black',origin='lower',extent=[-arcsecdata/2. , arcsecdata/2., -arcsecdata/2., arcsecdata/2. ])

    figtitle.set_position([.5,1.01])
    plt.savefig(dirname+'/'+figname,bbox_inches='tight')
    return

def plots_cont(verbose):
    """
    Plots of continuum stacking results 
    """
    if verbose >= 1:
        print('\n Plotting continuum stacking results ')

    #valid units: mJy, Jy, muJy
    unit = 'mJy'
    cube_file = 'all_sources_cube.fits'
    
    if verbose >= 2:
        print('Plot: all sources - mean')

    image_file = 'mean_stack.fits'
    figname = 'mean_stack.png'
    stack_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit, title = 'All - Mean')
    #figname = 'mean_stack_fit.png'
    #one_panel_fit(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit)
    figname = 'mean_stack_smooth.png'
    three_panel_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit)
    if verbose >= 2:
        print('Plot: all sources - median')

    cube_file = 'all_sources_cube.fits'
    image_file = 'median_stack.fits'
    figname = 'median_stack.png'
    stack_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit, title = 'All - Median')
    figname = 'median_stack_smooth.png'
    three_panel_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit)

    if os.path.isfile(dirname+'/'+'sources_cube_det.fits'):
        cube_file = 'sources_cube_det.fits'
        if verbose >= 2:
            print('Plot: only detections - mean')

        image_file = 'mean_stack_det.fits'
        figname = 'mean_stackdet.png'
        stack_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit, title = 'Detections - Mean')
        figname = 'mean_stackdet_smooth.png'
        three_panel_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit)
        figname = 'mean_stackdet_fit.png'
        one_panel_fit(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit,upperlimit=False)

        if verbose >= 2:
            print('Plot: only detections - median')

        image_file = 'median_stack_det.fits'
        figname = 'median_stackdet.png'
        stack_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit, title = 'Detections - Median')
        figname = 'median_stackdet_smooth.png'
        three_panel_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit)

    if os.path.isfile(dirname+'/'+'sources_cube_nondet.fits'):
        cube_file = 'sources_cube_nondet.fits'
        if verbose >= 2:
            print('Plot: only non-detections - mean')

        cube_file = 'sources_cube_nondet.fits'
        image_file = 'mean_stacknondet.fits'
        figname = 'mean_stacknondet.png'
        stack_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit, title = 'Non detections - Mean')
        figname = 'mean_stacknondet_fit.png'
        one_panel_fit(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit,upperlimit=False)
        #one_panel_fit(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit,upperlimit=True)
        figname = 'mean_stacknondet_smooth.png'
        three_panel_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit)

        if verbose >= 2:
            print('Plot: only non-detections - median')

        image_file = 'median_stack_nondet.fits'
        figname = 'median_stack_nondet.png'
        stack_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit, title = 'Non detections - Median')
        figname = 'median_stack_nondet_fit.png'
        one_panel_fit(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit)
        figname = 'median_stack_nondet_smoothing.png'
        three_panel_fig(stack_image_file = image_file, stack_cube_file = cube_file, figname = figname , unit=unit)
    return

#disables interactive matplotlib to avoid opening multiple figure windows
plt.ioff()
##
################################################################# 
dirname='stack_example'
images=['stack_example/cont_example.fits']
incat='stack_example/lines.csv'
stampsize = 15
verbose = 1

if os.path.isdir(dirname) == False:
    print('\n Output folder not found. Enter correct path.')
else:
    if os.path.isfile(incat) == False:
        print('\n Catalog not found. Enter correct path.')
    else:
        im=[]
        for file in images:
            im.append(os.path.isfile(file))
            if os.path.isfile(file) == False:
                print('\n File %s not found. Enter correct path'%file)

        if all(im):
            print('Working in folder : %s'%dirname)
            coords,stacked_cont = stackcont(catalog_name = incat,stampsize = stampsize ,images = images,weightstack = False,threshold = 5,overwrite = True,verbose = verbose)
            plots_cont(verbose)
            print('Done')