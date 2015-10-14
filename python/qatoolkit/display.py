#!/usr/bin/env python

import os,sys
import pyfits
import fitsio
import time
from despyastro import wcsutil
import despydb.desdbi
import numpy
# From here
import ds9DES as ds9

class DisplayDES:

    """
    A Class to display DECam detections using ds9' xpa
    Felipe Menanteau, NCSA (felipe@illinois.edu)
    """
    def __init__(self,**keys):
        
        self.ImageName  = keys.get('ImageName')
        self.CatName    = keys.get('CatName') 
        self.scalemode  = keys.get('scalemode','zscale')
        self.verb       = keys.get('verb',True) 

        self.SExColor   = keys.get('SExColor','blue')
        self.StarsColor = keys.get('StarsColor','yellow')


    def showSCI(self,fr=1):

        """ Display the Science Image """
        
        print "# Will try to display %s" % self.ImageName
        
        # Get the sci extension
        self.sci_hdu,self.sci_hdr = get_hdu_number(self.ImageName,extname='IMAGE')
        
        # Check env variables and set to proper values
        ds9.check_env(verb=self.verb);

        # Check if DS9 is already running -- if not start
        status = ds9.check_ds9(verb=self.verb);
        
        # Will display images in new frame if existing session
        if status == 'old':
            os.system("xpaset -p ds9 frame new")

        # Science on frame 1
        sci = "%s[%d]" % (self.ImageName,self.sci_hdu)
        print "# Displaying %s" % sci
        ds9.put(sci,fr=fr, zoom=True,
                scale     = 'linear', 
                scalemode = self.scalemode,
                advance   = None,
                limits    = None,
                cmap      = 'Gray') 


    def showWGT(self,fr=2):

        """ Display the Weight Image """
        
        print "# Will try to display %s" % self.ImageName

        # Get the wgt extension
        self.wgt_hdu,self.wgt_hdr = get_hdu_number(self.ImageName,extname='WEIGHT')

        # Check env variables and set to proper values
        ds9.check_env(verb=self.verb);

        # Check if DS9 is already running -- if not start
        status = ds9.check_ds9(verb=self.verb);
        
        # Will display images in new frame if existing session
        if status == 'old':
            os.system("xpaset -p ds9 frame new")

        # Weight on frame 2
        sci = "%s[%d]" % (self.ImageName,self.wgt_hdu)
        print "# Displaying %s" % sci
        ds9.put(sci,fr=fr, zoom=True,
                scale     = 'linear', 
                scalemode = self.scalemode,
                advance   = None,
                limits    = None,
                cmap      = 'Gray') 



    def whichCat(self):

        """ Figure out of if this is a flat SExtractor or a concatenated Scamp catalog"""

        print "# Reading: %s" % self.CatName
        hdulist = pyfits.open(self.CatName)

        Nhdu = len(hdulist)
        if Nhdu > 4: # Large Scamp catalog
            CCDNUM = self.sci_hdr['CCDNUM']
            print "# Will search for catalog for CCDNUM:%s" % CCDNUM
            hdu = (CCDNUM-1)*3 + 2
            cols  = hdulist[hdu].columns
            self.CatType = "multi"
        else: # Flat SEx catalog
            hdu = Nhdu-1
            self.CatType = "flat"

        self.hduCat     = hdu
        self.hdulistCat = hdulist 

    def readCat(self):

        # Figure out which type of Catalog, returns self.CatType
        self.whichCat()

        """ Read a SExtractor Single Epoch/Coadd or Scamp Catalog using pyfits"""
        
        tbdata = self.hdulistCat[self.hduCat].data
        cols   = self.hdulistCat[self.hduCat].columns
        
        self.NUMBER      = tbdata['NUMBER']
        self.X_IMAGE     = tbdata['X_IMAGE']
        self.Y_IMAGE     = tbdata['Y_IMAGE']
        self.FLAGS       = tbdata['FLAGS']
        
        if self.CatType == 'flat':
            # Store the Relevant positional information existan in Flat SEx cats
            self.KRON_RADIUS = tbdata['KRON_RADIUS']
            self.A_IMAGE     = tbdata['A_IMAGE']*self.KRON_RADIUS
            self.B_IMAGE     = tbdata['B_IMAGE']*self.KRON_RADIUS
            self.THETA_IMAGE = tbdata['THETA_IMAGE']
            try:
                self.IMAFLAGS_ISO = tbdata['IMAFLAGS_ISO']
            except:
                print "# WARNING: Could not lead IMAFLAGS_ISO"


    def displayCatalog(self,imaflags=False):

        """ Make the call to display the relevant catalog"""

        if self.CatType == 'flat':
            if imaflags:
                self.displaySExDetectionsFlags()
            else:
                self.displaySExDetections()
        elif self.CatType == 'multi':
            self.displaySCampDetections()
        else:
            sys.exit("ERROR: Catalog Type not defined")
            
    def displaySExDetections(self):

        """ Display detection for SExtractor flat catalog"""
        shapes = zip(self.X_IMAGE,self.Y_IMAGE,self.A_IMAGE,self.B_IMAGE,self.THETA_IMAGE)
        ds9.ellipses(shapes,color=self.SExColor,text='',options='',units='deg',out=None)      

    def displaySExDetectionsFlags(self):

        """ Display detection for SExtractor flat catalog"""
        idx = numpy.where(self.IMAFLAGS_ISO > 0)
        shapes = zip(self.X_IMAGE[idx],self.Y_IMAGE[idx],self.A_IMAGE[idx],self.B_IMAGE[idx],self.THETA_IMAGE[idx])

        ds9.ellipses(shapes,color='blue',text='',options='',units='deg',out=None)      

        idx  = numpy.where(self.FLAGS < 4)
        shapes = zip(self.X_IMAGE[idx],self.Y_IMAGE[idx],self.A_IMAGE[idx],self.B_IMAGE[idx],self.THETA_IMAGE[idx])
        ds9.ellipses(shapes,color=self.SExColor,text='',options='',units='deg',out=None)      


    def displaySCampDetections(self):

        """ Display detection for SCamp concatenate """
        ds9.circles(self.X_IMAGE,self.Y_IMAGE,radius=20,color=self.SExColor,system='physical')
        
    def getCCDlimits(self):

        """
        Get the corners of a DECam CCD using Erin Sheldon wcsutil
        """

        nx = self.sci_hdr['NAXIS1']
        ny = self.sci_hdr['NAXIS2']

        # Create the object call
        wcs = wcsutil.WCS(self.sci_hdr)

        ra1,dec1 = wcs.image2sky(1, 1)
        ra2,dec2 = wcs.image2sky(nx,1)
        ra3,dec3 = wcs.image2sky(nx,ny)
        ra4,dec4 = wcs.image2sky(1, ny)

        self.ra_min  = min(ra1,ra2,ra3,ra4)
        self.ra_max  = max(ra1,ra2,ra3,ra4)
        self.dec_min = min(dec1,dec2,dec3,dec4)
        self.dec_max = max(dec1,dec2,dec3,dec4)

        self.wcs = wcs

        
    def getCatNOMAD(self):

        """
        Get the Star Catalogs from NOMAD using desar
        """

        # Setup desar queries here for later
        section = "db-desoper"
        try:
            desdmfile = os.environ["des_services"]
        except KeyError:
            desdmfile = None
        dbh = despydb.desdbi.DesDbi(desdmfile,section)
        cur = dbh.cursor()

        query = """
        select ra, dec, B,V from NOMAD
        where ra  > %s and
              ra  < %s and
              dec > %s and
              dec < %s and
              ( B < 19 or V < 19)""" % (self.ra_min,self.ra_max,self.dec_min,self.dec_max)

        print "# Will execute the SQL query:\n********\n** %s\n********" % query
        cur.execute(query)
        (x,y,m1,m2) = zip(*cur.fetchall())
        self.ra_NOMAD  = numpy.array(x)
        self.dec_NOMAD = numpy.array(y)
        self.B_NOMAD   = numpy.array(m1) 
        self.V_NOMAD   = numpy.array(m2) 

    def displayNOMAD(self):

        self.getCatNOMAD()
        # Get from (ra,dec) --> (x,y)
        x,y = self.wcs.sky2image(self.ra_NOMAD,self.dec_NOMAD)
        ds9.crosses(x,y,color=self.StarsColor,system='physical')
        


###########################
# General useful functions
###########################
def get_hdu_numbers(fitsFile):
    """
    Simple function to figure the HDU extensions for DESDM fits files,
    and better than relying on the .fz vs .fits ending of the files
    """
    sci_ext = None
    msk_ext = None
    wgt_ext = None
    # Loop trough each HDU on the fits file
    FITS = pyfits.open(fitsFile)
    for i in range(len(FITS)):
        h = FITS[i].header       # Get the header

        if ('DES_EXT' in h.keys()) :
            extname = h['DES_EXT'].strip()
            if (extname == 'IMAGE') :
                sci_ext = i
            elif (extname == 'MASK') :
                msk_ext = i
            elif (extname == 'WEIGHT'):
                wgt_ext = i

    if (sci_ext is None or msk_ext is None or wgt_ext is None):
        sys.exit("Cannot find IMAGE, MASK, WEIGHT extensions via DES_EXT keyword in header")

    FITS.close()
    return sci_ext,msk_ext,wgt_ext

def get_sci_hdu(fitsFile):

    """
    Simple function to figure the HDU extensions for DESDM fits files,
    and extract the science header.
    Better than relying on the .fz vs .fits ending of the files
    """ 
     
    sci_hdu = None
    # Loop trough each HDU on the fits file
    FITS = pyfits.open(fitsFile)
    for i in range(len(FITS)):
        h = FITS[i].header       # Get the header

        if ('DES_EXT' in h.keys()) :
            extname = h['DES_EXT'].strip()
            if (extname == 'IMAGE') :
                sci_hdu = i
                sci_hdr = h.copy()

    if (sci_hdu is None):
        sci_hdu = 0
        sci_hdr = FITS[0].header.copy()
        print "# WARNING: Cannot find IMAGE extensions via DES_EXT keyword in header, defaulting to 0"

    FITS.close()
    return sci_hdu,sci_hdr


def get_hdu_number(fitsFile,extname='IMAGE'):

    """
    Simple function to figure the HDU extensions for DESDM fits files,
    and extract the science header.
    Better than relying on the .fz vs .fits ending of the files
    """ 
     
    hdu_num = None
    # Loop trough each HDU on the fits file
    FITS = pyfits.open(fitsFile)
    for i in range(len(FITS)):
        h = FITS[i].header       # Get the header

        if ('DES_EXT' in h.keys()) :
            des_ext = h['DES_EXT'].strip()
            if (des_ext == extname) :
                hdu_num = i
                hdu_hdr = h.copy()

    if (hdu_num is None):
        hdu_num = 0
        hdu_hdr = FITS[0].header.copy()
        print "# WARNING: Cannot find extensions via DES_EXT keyword in header, defaulting to 0"

    FITS.close()
    return hdu_num,hdu_hdr



# The main fuction to fill up all of the options
def cmdline():

    import argparse
    parser = argparse.ArgumentParser(description="Display SExtractor detections of DECam images")

    # The positional arguments
    parser.add_argument("ImageName", help="Fits file Image [fits.fz or .fits]")
    parser.add_argument("CatName",   help="SExtractor Catalog [fits format]")
    
    # The optional arguments for display
    parser.add_argument("--scalemode", default='zscale',
                        help="scalemode to use [default=zscale]")

    parser.add_argument("--SExColor", default='yellow',
                        help="Color for SExtractor detections [default=blue]")

    parser.add_argument("--StarsColor", default='yellow',
                        help="Color for Star's catalog [default=yellow]")
    
    args = parser.parse_args()

    print "# Will run:"
    print "# %s " % parser.prog
    for key in vars(args):
        print "# \t--%-10s\t%s" % (key,vars(args)[key])

    return args


# The main fuction to fill up all of the options
def cmdline_coadd():

    import argparse
    parser = argparse.ArgumentParser(description="Display SExtractor detections of DECam images")

    # The positional arguments
    parser.add_argument("ImageName", help="Fits file Image [fits.fz or .fits]")
    parser.add_argument("CatName",   help="SExtractor Catalog [fits format]")
    
    # The optional arguments for display
    parser.add_argument("--scalemode", default='zscale',
                        help="scalemode to use [default=zscale]")

    parser.add_argument("--SExColor", default='yellow',
                        help="Color for SExtractor detections [default=blue]")

    parser.add_argument("--StarsColor", default='red',
                        help="Color for Star's catalog [default=yellow]")

    parser.add_argument("--ShowWeight", action='store_true',default=False,
                        help="Show the weight plane [default=False]")

    parser.add_argument("--imaflags", default=False, action='store_true',
                        help="Color for SExtractor detections [default=blue]")

    args = parser.parse_args()

    print "# Will run:"
    print "# %s " % parser.prog
    for key in vars(args):
        print "# \t--%-10s\t%s" % (key,vars(args)[key])

    return args

