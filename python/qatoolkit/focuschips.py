#!/usr/bin/env python

import os,sys
import fitsio
import numpy
from PIL import Image

def parse_section(section):

    """ Hack to parse sections """
    
    section = section.strip("[")
    section = section.strip("]")

    j1,j2 = map(int,section.split(",")[0].split(":"))
    i1,i2 = map(int,section.split(",")[1].split(":"))

    return i1,i2, j1,j2

def extract_and_png(file,graylevel=0.95,outdir='.' ,TNsize=300):

    """ Extract focus chips and create the PNGs """

    # Current location
    curdir = os.getcwd()
    root = os.path.basename(file).split(".")[0]
    fits = fitsio.FITS(file)

    # Check if outpath exists
    if not os.path.exists(outdir):
        print "# Creating output dir: %s" % outdir
        os.mkdir(outdir)

    # Go to the place we want to work
    os.chdir(outdir)

    # CHANGE: This only work if we use eups, might want more general solution
    if os.uname()[0] == 'Darwin':
        print "# Setting local Dawing Environment"
        stiff_conf = '%s/DESDM-Code/devel/qatoolkit/trunk/etc/default.stiff' % os.environ['HOME']
    else:
        stiff_conf = '%s/etc/default.stiff' % os.environ['QATOOLKIT_DIR']

    # Got trough each HDU
    for k in range(len(fits)):

        # Skip the first HDU in fz files
        if k==0:
            continue

        # Get the header information
        h = fits[k].read_header()

        # Skip CCDNUM below 62 too
        if h['CCDNUM'] <= 62:
            continue

        # Read in the hdu
        newdata = fits[k].read()

        # Full A+B
        newname = '%s_%s.fits' % (root,h['CCDNUM'])
        newfits = fitsio.FITS(newname,'rw')
        newfits.write(newdata, header=h)
        newfits.close()
        print "# Extracted CCDNUM: %s --> %s" % (h['CCDNUM'],newname)
        
        # Amp A
        i1,i2, j1,j2 = parse_section(h['DATASECA'])
        newnameA = '%s_%sA.fits' % (root,h['CCDNUM'])
        newtiffA = '%s_%sA.tif' % (root,h['CCDNUM'])
        newfitsA = fitsio.FITS(newnameA,'rw')
        newfitsA.write(newdata[i1:i2,j1:j2], header=h)
        newfitsA.close()
        iA = j1

        # Amp B
        i1,i2, j1,j2 = parse_section(h['DATASECB'])
        newnameB = '%s_%sB.fits' % (root,h['CCDNUM'])
        newtiffB = '%s_%sB.tif' % (root,h['CCDNUM'])
        newfitsB = fitsio.FITS(newnameB,'rw')
        newfitsB.write(newdata[i1:i2,j1:j2], header=h)
        newfitsB.close()
        iB = j1

        # Figure out the order
        if iA<iB:
            order = 'AB'
        else:
            order = 'BA'

        # Create the PNG for each amplifier 
        cmd1 = 'stiff -c %s %s -OUTFILE_NAME %s -MAX_LEVEL %s > /dev/null 2>&1' % (stiff_conf,newnameA,newtiffA,graylevel)
        cmd2 = 'stiff -c %s %s -OUTFILE_NAME %s -MAX_LEVEL %s > /dev/null 2>&1' % (stiff_conf,newnameB,newtiffB,graylevel)
        os.system(cmd1)
        os.system(cmd2)

        imA = Image.open(newtiffA)
        imB = Image.open(newtiffB)

        # Get the size of the new image
        (nx,ny) = imA.size
        newsize = (nx*2,ny)
        newim = Image.new(mode='L',size=newsize,color=None)

        # Paste them together
        if order == 'AB':
            newim.paste(imA,box=(0,0))
            newim.paste(imB,box=(nx,0))
        elif order == 'BA':
            newim.paste(imB,box=(0,0))
            newim.paste(imA,box=(nx,0))
        else:
            print "Order not defined"

        # Save the png file now
        newpng     = '%s_%s.png' % (root,h['CCDNUM'])
        newim.save(newpng,format='png')

        # Now let's create the thumbnails, using PIL/Image
        newpng_TN  = '%s_%s_TN.png' % (root,h['CCDNUM'])
        make_png_thumbnail(newpng,newpng_TN,TNsize=TNsize)

        # Now we should clean up the A,B files and tif
        files = [newtiffA,newtiffB, newnameA,  newnameB]
        for f in files:
            os.remove(f)

    return


def cmdline_bypath():

    from optparse import OptionParser

    """ Read in the command line options """

    USAGE = "\n"
    USAGE = USAGE + "\t %prog /path/to/file/in/desar/DECam_00XXYYXX.fitz.fz --outdir <some-place>\n"
    USAGE = USAGE + "\t i.e.: \n"
    USAGE = USAGE + "\t %prog /archive_data/Archive/DTS/src/20131014/src/DECam_00244099.fits.fz --outdir ./focus-files \n"

    parser = OptionParser(usage=USAGE)

    parser.add_option("--outdir",
                      dest="outdir", default='.',
                      help="Output Directory to put files in")

    parser.add_option("--grayscale",
                      type='float',dest="grayscale", default=0.95,
                      help="grayscale for png creation [i.e. 0.950]")

    parser.add_option("--TNsize",
                      type='int',dest="TNsize", default=300,
                      help="Size of thumbnail [pixels]")


    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("\n\tERROR:incorrect number of arguments")        
        
    return options,args


def cmdline_SQL():

    from optparse import OptionParser


    """ Read in the command line options """
    USAGE = "\n"
    USAGE = USAGE + "\t %prog DECam_00XXYYXX  --outdir <some-place>\n"
    USAGE = USAGE + "\t i.e.: \n"
    USAGE = USAGE + "\t %prog DECam_00244099.fits  --outdir ./focus-files \n"

    parser = OptionParser(usage=USAGE)

    parser.add_option("--outdir",
                      dest="outdir", default='.',
                      help="Output Directory to put files in")

    parser.add_option("--grayscale",
                      type='float',dest="grayscale", default=0.95,
                      help="grayscale for png creation [i.e. 0.950]")

    parser.add_option("--TNsize",
                      type='int',dest="TNsize", default=300,
                      help="Size of thumbnail [pixels]")


    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("\n\tERROR:incorrect number of arguments")        
        
    return options,args


def make_png_thumbnail(pngName,TN_pngName,TNsize=800):

    """ Make the thumbnail of the PNG for the webpages"""        
    
    print "# Creating TN %s" % TN_pngName
    im = Image.open(pngName)
    (w,h) = im.size # Get width (w) and height (h)
    TNshape = TNsize,int(TNsize*h/w)
    im.thumbnail(TNshape, Image.ANTIALIAS)
    im.save(TN_pngName,"png",options='optimize')
    return
