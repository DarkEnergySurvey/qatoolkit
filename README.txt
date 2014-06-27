To set up for testing:
----------------------
 setup -r ~/DESDM-Code/devel/qatoolkit/trunk
   or
 setup -r ~/DESDM-Code/devel/qatoolkit/tags/XX.YY.ZZ

To install:
-----------
 
 python setup install --home=$HOME/Python
  or
 python setup.py install --prefix=$PRODUCT_DIR --install-lib=$PRODUCT_DIR/python 


Scripts: (inside bin/)
----------------------

 - color_tile : generate RGB color image for a DES Tile

 - projectDECamPNG : Projects a DECam exposure using SWarp and creates
 grayscale PNGs using Python's native PIL/Image module. It reads the
 input files from lists provided command-line. To be used in the
 Refact system.

 - projectDECamPNG_SQL : Projects a DECam exposure using SWarp and
 creates grayscale PNGs using Python's native PIL/Image module. Gets
 the file's location using SQL

 - projectDECamPNG_bypath : Projects a DECam exposure using SWarp and
 creates grayscale PNGs using Python's native PIL/Image module. Gets
 the file's location using the paths and the structure in the old system

 - compare_corners_NewFramework.py : compares the CCD corners values
 stored in the DB (in db-destest for now) against the values computed from
 the WCS information in the image header.

 - make_png_focus_chips_bypath : simple call to the focuschips functions
 inside focuschips.py to create png of the raw focus chips for the
 exposure images. Uses the input filename in arguments

 - make_png_focus_chips_SQL : simple call to the focuschips functions
 inside focuschips.py to create png of the raw focus chips for the
 exposure images. Uses SQL to figure the location of the file in
 desar.

 - display_DECam_MEFs: A very, very simple stand-alone script to
 display DECam CCD mef fits files. It uses a lot of functions from 
 ds9DES.py module, which are available via despyutils.

 - display_DECam_detection_NOMAD: Another simple method/script to
 display DECam detections using ds9' xpa and overlay NOMAD catalog for
 that section in order the check astrometry aligment.
