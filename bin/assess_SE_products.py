#! /usr/bin/env python
"""
Perfoms an assessent of exposures from a first/final cut run.  The quality 
of the exposures is based upon the seeing (FWHM), background, and extinction
due to clouds.


Syntax:
    quick_assess_RUN.py  -r run 
Arguments:
    Requires a run number for which the assessemnt for each exposure will be based.
     
"""

########################
def QAplt_cloud( PlotFile, Data, cat,  band, astrom_good, verbose=False):
    """
    Function to produce a QA plot showing comparison between data and a reference
    catalog.  Plots show the diagnosed amount of atmospheric opacity (extinction) 
    at the time of the observation.

    inputs:
        PlotFile:   filename where the output plot will be written (can include 
                        relative path)
        Data:       Custom dictionary that holds data for the plots.  These include:
                        nmatch:  (integer) number of matches to the reference catalog

                        mag_des: list of rough measured DES magnitude
                        mag_cat: list of associated reference catalog data
                        mag_diff: the difference between the two
                        mag_des_reject: points in mag_des that were outlier rejects
                        mag_cat_reject: analogous points in mag_cat
                        mag_diff_reject: analogous points in mag_diff

                        mindes,maxdes: (float) range of data being plotted
                        med_magdiff:   (float) median magnitide offset/difference 
                                       between objects and refernce catalog

                    For failure cases (number of stars < 100) the following are not needed:
                        mag_des_reject,mag_cat_reject, mag_diff_reject
                        mindes,maxdes,med_magdiff

        cat:        Name of the reference catalog (string for plotting axis)
        band:       Band of observations
        astrom_good: Flag that causes overplotting of text indicating an astrometric failure
        verbose:    Provide verbose output (curently there is none).

    OUTPUT: 
        Nothing is returned to the calling program (a PNG file is written where directed)
    """

    plt.figure()
    plt.subplot(2,1,1)
    if (Data['nmatch']<100):
#       Case no fit
        plt.scatter(Data['mag_des'],Data['mag_cat'],marker='.',color='blue')
        plt.plot([10,18],[10,18],color='red',linewidth=1)
    else:
#       Case normal
        plt.scatter(Data['mag_des'],Data['mag_cat'],marker='.',color='blue')
        if (len(Data['mag_des_reject']) > 0):
            plt.scatter(Data['mag_des_reject'],Data['mag_cat_reject'],marker='.',color='magenta')
        plt.plot([Data['mindes'],Data['maxdes']],[Data['mindes'],Data['maxdes']],color='red',linewidth=1)

    plt.xlabel('DES MAG_AUTO({:s})'.format(band))
    if (band in ["u","g"]):
        plt.ylabel('%s g\'' % cat.upper() )
    elif (band in ["r","VR"]):
        plt.ylabel('%s r\'' % cat.upper() )
    elif (band in ["i","z","Y"]):
        plt.ylabel('%s i\'' % cat.upper() )
    else:
        plt.ylabel('%s [unknown]' % cat.upper() )
    plt.text(10.5,17,'Number of %s matches: %d' % (cat.upper(),Data['nmatch']))
    if (not(astrom_good)):
        plt.text(10.5,16,'Probable failure to find an astrometric solution')
#
#   Second panel
#
    plt.subplot(2,1,2)
    if (Data['nmatch']<100):
#       Case no fit
        plt.plot([10,18],[0,0],color='red',linewidth=3)
    else:
#       Case normal
        plt.scatter(Data['mag_des'],Data['mag_diff'],marker='.',color='blue')
        if (len(Data['mag_des_reject']) > 0):
            plt.scatter(Data['mag_des_reject'],Data['mag_diff_reject'],marker='.',color='magenta')
        plt.plot([Data['mindes'],Data['maxdes']],[Data['med_magdiff'],Data['med_magdiff']],color='red',linewidth=3)

    plt.xlabel('DES MAG_AUTO(%s)' % band)
    if (band in ["u","g"]):
        plt.ylabel('DES - %s g\'' % cat.upper() )
    elif (band in ["r","VR"]):
        plt.ylabel('DES - %s r\'' % cat.upper() )
    elif (band in ["i","z","Y"]):
        plt.ylabel('DES - %s i\'' % cat.upper() )
    else:
        plt.ylabel('DES - %s [unknown]' % cat.upper() )
    plt.savefig(PlotFile)

    return 0


########################
def QAplt_maghist(PlotFileName,Data,band,astrom_good,verbose=False):
    """
    Function to produce a QA plot showing luminosity function of objects in the
    exposure.

    inputs:
        PlotFile:   filename where the output plot will be written (can include 
                        relative path)
        Data:       Custom dictionary that holds data for the plots.
                        nobj:        (integer) total number of objects 
                                         contributing to the histogram
                        bins:        list of magnitude bins
                        mag_hist:    list of number of objects per bin
                        magerr_hist: list of median value of magerr for 
                                        objects in each magnitude bin
                    For failure cases (nobj < 100) no data from lists are plotted

        band:       Band of observations
        astrom_good: Flag that causes overplotting of text indicating an astrometric failure
        verbose:    Provide verbose output (curently there is none).

    OUTPUT: 
        Nothing is returned to the calling program (a PNG file is written where directed)
    """


    plt.figure()
    plt.subplot(2,1,1)
#    plt.scatter(xplt,yplt1,marker='.',color='blue')
    plt.axis([8,26,0.5,10000])
    if (Data['nobj']>10):
        plt.semilogy(Data['bins'],Data['mag_hist'],marker='.',ls='None',color='blue')
    else:
        plt.text(9.,9000.,'Number of objects ({:d}) insufficient for histogram.'.format(Data['nobj']))
    if (not(astrom_good)):
        plt.text(9.,7800,'Probable failure to find an astrometric solution')
    plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
    plt.ylabel('# objects')
#
#   Second panel showing median MAGERR per bin
#    
    plt.subplot(2,1,2)
    if (Data['nobj']>10):
        plt.scatter(Data['bins'],Data['magerr_hist'],marker='.',color='blue')
    plt.axis([8,26,-0.01,0.5])
    plt.xlabel('DES MAG_AUTO({:s})'.format(band))
    plt.ylabel('median(MAGERR_AUTO({:s})'.format(band))
    plt.savefig(PlotFileName)

    return 0

#######################################################################
def MkCatDict(CATcol,cols_retrieve):
    """
    For use with FITSIO.  Make a dictionary that relates the order 
    of the requested columns to the order they are retreived by FITSIO.

    INPUT: CATcol is the list of columns from the FITS table
                  e.g. CATcol=cat_fits[hdu].get_colnames() 
           cols_retrieve is the list of column names requested from the FITS table.
          e.g. CAT=cat_fits[hdu].read(columns=cols_retrieve)

    OUTPUT: catdict

            for row in CAT:
                row[catdict["COLUMN_NAME"]]
            where "COLUMN_NAME" was an element in cols_retrieve and a column in CAT
 
    """
#
#   Figure out which column in CAT corresponds to each requested column
#   Get all columns in table (form dictionary showing their order)
#   Form a reverse dictionary based on columns retrieved
#   Then sort to get the order of the columns retieved in CAT and 
#   form a dictionary/legend that contains their order
#
    full_catdict={}
    for index, item in enumerate(CATcol):
        full_catdict[item]=index
    rev_catdict={}
    catdict_list=[]
    for item in cols_retrieve:
        rev_catdict[full_catdict[item]]=item
        catdict_list.append(full_catdict[item])
    icnt=0
    catdict={}
    for item in sorted(list(set(catdict_list))):
        catdict[rev_catdict[item]]=icnt
        icnt=icnt+1
#
    return catdict
####################################################################
def binsrch_crosscorr_cat(ra,dec,catalog,match_rad):
    """
    Search through a catalog for match to a set of coordinates (only nearest match is returned).

    ra        = RA used for comparison
    dec       = Dec used for comparison
    catalog   = is a list containing records(dictionaries) that have an "ra" and "dec" entries 
                along with associated information.  NOTE: This routine assumes that the catalog
                is RA-ordered so that a binary search can be used to find the small range of 
                entries that need to be checked
    match_rad = minimum distance for an acceptable match (in seconds of arc).
    """

    pi=3.141592654
    halfpi=pi/2.0
    deg2rad=pi/180.0

#
#   Search radius value
#
#    match_rad2=4.*(match_rad/3600.)*(match_rad/3600.)
    match_rad2=4.*(match_rad/3600.)

    ra1_deg=ra
    dec1_deg=dec
    sinl2=numpy.sin(halfpi-(dec1_deg*deg2rad))
    cosl2=numpy.cos(halfpi-(dec1_deg*deg2rad))
    cosdec1=numpy.cos(dec1_deg*deg2rad)
    num_rec=len(catalog)

#
#   Work out maximum RA range that will need to be searched
#
#    print "DRA range calculation: ",(match_rad2/cosdec1),(match_rad2/cosdec1*3600.0)
    if (abs(dec1_deg)>89.0):
        ra_range1=0.0
        ra_range2=360.0
    else:
        ra_range1=ra1_deg-(match_rad2/cosdec1)
        if (ra_range1 > 360.):
            ra_range1=ra_range1-360.0
        if (ra_range1 < 0.0): 
            ra_range1=ra_range1+360.0
        ra_range2=ra1_deg+(match_rad2/cosdec1)
        if (ra_range2 > 360.):
            ra_range2=ra_range2-360.0
        if (ra_range2 < 0.0): 
            ra_range2=ra_range2+360.0
#
#   Create a list of RA ranges that will be searched
#   (two ranges are necessary for initial binary search)
#
    ra_range_list=[]
    if (ra_range1 > ra_range2):
#       print "Need two ranges for search around 0h-24h branch ",ra1_deg,dec1_deg
        ra_range_list.append([0.0,ra_range2])
        ra_range_list.append([ra_range1,360.0])
#       print ra_range_list
    else:
#       print "Single range will suffice for ",ra1_deg,dec1_deg
        ra_range_list.append([ra_range1,ra_range2])
#
#   Use each range to search for a counterpart.
#
    match_record=[]
    rad_dist=[]
    for ra_range in ra_range_list:
#
#   Initiate a binary search to find the appropriate entries in the catalog
#
        klo=0
        khi=num_rec-1
        while (khi-klo > 1):
            knew=(khi+klo)/2
#           print klo,khi,knew
            if (catalog[knew]["ra"] > ra_range[0]):
                khi=knew
            else:
                klo=knew
#
#       Should have initial spot to begin a serious comparison
#       Loop through catalog values until RA exceeds upper limit of range
#           or the end of the catalog has been reached
#
#       print klo,khi,catalog[klo]["ra"]
        ncheck=0
        while ((klo < num_rec-1)and(catalog[klo]["ra"]<ra_range[1])):
            ncheck=ncheck+1
            record=catalog[klo]
            ra2_deg=record["ra"]
            dec2_deg=record["dec"]
            dra=abs(ra2_deg-ra1_deg)/cosdec1
            ddec=abs(dec2_deg-dec1_deg)
            if ((dra < match_rad2)and(ddec < match_rad2)):
#               print ra2_deg,dec2_deg
                ang1=halfpi-(dec2_deg*deg2rad)
                ang2=deg2rad*(ra2_deg-ra1_deg)
                if (ang2 < 0.0):
                    ang2=-ang2
                val=numpy.cos(ang1)*cosl2+numpy.sin(ang1)*sinl2*numpy.cos(ang2)
                if (val > 1.):
                    dist=0.0
                else:
                    val2=numpy.sqrt(1.0-val*val)
                    dist=3600.0*(halfpi-numpy.arctan2(val,val2))/deg2rad
                if (dist <= match_rad):
                    if (len(rad_dist) > 0):
                        if (dist < rad_dist[0]):
                            rad_dist=[]
                            match_record=[]
                            match_record.append(record)
                            rad_dist.append(dist)
                    else:
                        match_record.append(record)
                        rad_dist.append(dist)
            klo=klo+1
        if (ncheck > 1000):
            print "Warning possible brute force use detected (RA,DEC,ncheck,range): ",ra1_deg,dec1_deg,ncheck,ra_range
#       if (ncheck == 0):
#           print "Warning no checks made (RA,DEC,ncheck,range): ",ra1_deg,dec1_deg,ncheck,ra_range

    return match_record
####################################################################

if __name__ == "__main__":

    import argparse
    import os
    import despydb.desdbi
    import stat
    import time
    import math
    import re
    import csv
    import sys
    import datetime
    import numpy
    import fitsio
    import matplotlib 
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    svnid="$Id$"
    svnrev=svnid.split(" ")[2]
    db_table="firstcut_eval"

    parser = argparse.ArgumentParser(description='Assess whether the pipeline products of a FIRSTCUT/FINALCUT processing meet survey quality metrics.')
    parser.add_argument('-i', '--attid',   action='store', type=str, default=None, help='Processing attempt id', required=True)
#    parser.add_argument('-R', '--ReqNum',   action='store', type=str, default=None, help='Processing request number.')
#    parser.add_argument('-U', '--UnitName', action='store', type=str, default=None, help='Unit Name.')
#    parser.add_argument('-A', '--AttNum',   action='store', type=str, default=None, help='Attempt Number.')
    parser.add_argument('-C', '--CamSym',   action='store', type=str, default='D', help='CamSym (default=\"D\").')
    parser.add_argument('-l', '--list',     action='store', type=str, default=None, 
                        help='Optional list of catalogs for the assessment (otherwise query for "red_finalcut" catalogs associated with the Request Unit Attempt)')

    parser.add_argument('-a', '--analyst',  action='store', type=str, default='assess_RUN_v%s'%(svnrev), 
                        help='Provides value for analyst (default: assess_RUN_v%s, \"None\" will use os.getlogin())'%(svnrev))
    parser.add_argument('-c', '--csv',      action='store_true', default=False, help='Flag for optional output of CSV')
    parser.add_argument('-u', '--updateDB', action='store_true', default=False, help='Flag for program to DIRECTLY update DB (%s).' % (db_table) )
    parser.add_argument('--over_table',     action='store', type=str, default=None, help='Override output DB table with specified table')
    parser.add_argument('-D', '--DB_file',  action='store_true', default=False, help='Flag for optional output of DB update file')
    parser.add_argument('-f', '--froot',    action='store',      default=None, help='Root for output file names')
#    parser.add_argument('-M', '--Magout',   action='store_true', default=False,
#                         help='Flag to dump file containing raw data summarizing number of objects per magnitude bin.')
    parser.add_argument('-q', '--qaplot',   action='store_true', default=False, help='Flag to generate QA plots.')
    parser.add_argument('-d', '--delimiter', action='store',  default=',', help='Optional delimiter for CSV (default=,)')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Print extra (debug) messages to stdout')
    parser.add_argument('--format_query',  action='store_true', default=False, help='Print queries with formatting (default=False)')
    parser.add_argument('-s', '--section', action='store', type=str, default=None, help='section of .desservices file with connection info')
    parser.add_argument('-S', '--Schema', action='store', type=str, default=None, help='Schema')

    args = parser.parse_args()
    if (args.verbose):
        print "##### Initial arguments #####"
        print "Args: ",args

#   Some old flags/arguemnts almost ready to remove completely...
##    parser.add_option("-m", "--maglimit", action="store_true", dest="maglimit", default=Fals,
##                         help="Flag to include magnitude limit calculation.")
##############################################################################
#   Arguments are defined... some extra checks before we begin.
#    if ((args.ReqNum is None)or(args.UnitName is None)or(args.AttNum is None)):
#    if (args.PFW_Attempt_ID is None):
#        parser.print_help()
#        exit(1)
    PFW_Attempt_ID=args.attid
#    ReqNum=args.ReqNum
#    UnitName=args.UnitName
#    AttNum=args.AttNum

    print("####################################################")
    print("# Note this execution requests:")
    print("# ")
    print("# Evaluation of PFW_attempt_id={:s}  ".format(args.attid))
    print("# ")
#
#   Obtain Schema (if user specified).
#
    if (args.Schema is None):
        db_Schema=""
    else:
        db_Schema="%s." % (args.Schema)

    if (args.over_table is None):
        db_table='%s%s' % (db_Schema,db_table)
    else:
        if (len(args.over_table.split('.')) > 1):
            db_table=args.over_table
        else:
            db_table='%s%s' % (db_Schema,args.over_table)

    if (args.updateDB):
        print("# DB WRITE evaluation to: {:s}".format(db_table))
    if (args.DB_file):
        print("# File WRITE evaluation to {:s}.dbupdate as though for DB table {:s}".format(args.froot,db_table))
    else:
        print("# WRITE summary evaluation to STDOUT")

#
#   Automatically get analyst name (from os.getlogin()) unless an override value has been specified
#
    if ((args.analyst is None)or(args.analyst == "None")):
        analyst=os.getlogin()
    else:
        analyst=args.analyst
    print("# Using Analyst: {:s}".format(analyst))
    print("# ")
    if (args.qaplot):
        print("# QA Plots requested (with root filename: {:s}).".format(args.froot))
        print("# ")

#
#   Setup for database interactions (through despydb)
#
    dbh = despydb.desdbi.DesDbi(None,args.section,retry=True)
    cur = dbh.cursor()

#################################################
#
#   Define pi, halfpi, deg2rad conversion
#
    pi=3.141592654
    halfpi=pi/2.0
    deg2rad=pi/180.0
#
#   Define pixel size 
#   Define a magnitude bin range to search for completeness/depth
#   Define a dictionary to translate filter/band to an ordinate
#
    pixsize=0.263
    fp_rad=1.2
#   Below (fwhm_DMtoQC_offset was an empirical offset when using FWHM_WORLD)
    fwhm_DMtoQC_offset_world=1.10
#   New veresion (an additive offset needed when comparing FWHM_MEAN (from PSFex) with respect to QC)
    fwhm_DMtoQC_offset_psfex=+0.04
    
    magerr_thresh=0.1
    magnum_thresh=20
    magbin_min=10.
    magbin_max=25.
    magbin_step=0.25
    mbin=numpy.arange(magbin_min,magbin_max,magbin_step)
    band2i={"u":0,"g":1,"r":2,"i":3,"z":4,"Y":5,"VR":6}
#
#   Old ellipticity limits (currently not used)
#
#    ellip_lim=0.13
#    ellip_good=0.07
#
    kolmogorov={"u":1.2,"g":1.103,"r":1.041,"i":1.00,"z":0.965,"Y":0.95,"VR":1.04}
    teff_lim={  "u":0.2,"g":0.2,  "r":0.3,  "i":0.3, "z":0.3,  "Y":0.2,"VR":0.3}
    seeing_lim={}
    seeing_fid={}
#
#   Set seeing cutoff to be 1.6 times Kolmogov except at "Y" which should
#   be forced to match that at g-band
#
    for band in ["u","g","r","i","z","Y","VR"]:
        if (band == "Y"):
            seeing_lim[band]=1.6*kolmogorov["g"]
        else:
            seeing_lim[band]=1.6*kolmogorov[band]
#       Commented version below was needed when using FWHM_WORLD
#        seeing_fid[band]=fwhm_DMtoQC_offset*0.9*kolmogorov[band]
#       Now fiducial value is additive (and applied to the FWHM_MEAN value coming from PSFex)
        seeing_fid[band]=0.9*kolmogorov[band]
#    print kolmogorov
#    print seeing_lim
#    print seeing_fid
#
#   Surface brightness limits from Eric Nielson which were derived "...from a few 
#   exposures from a photometric night in SV with little moon (20121215)"
#
    sbrite_good={"u":0.2,"g":1.05,"r":2.66,"i":7.87,"z":16.51,"Y":14.56,"VR":3.71}
    sbrite_lim={"u":0.8,"g":4.0,"r":9.58,"i":21.9,"z":50.2,"Y":27.6,"VR":13.58}
#
#   These (the above) were originally based on the following estimate by Annis
#   sbrite_good={"u":2.0,"g":1.2,"r":3.8,"i":8.7,"z":20.0,"Y":11.0}
#   roughly equivalent to grizY=22.09,21.21,20.12,18.95,18.00 mag/sq-arcsec
#
#   APASS and NOMAD magnitude limits
#
    glimit=90.0
    rlimit=90.0
    ilimit=90.0
#
    jlimit=16.0
    blimit=18.0
#
#   APASS and NOMAD convergence criteria (stop performing cross_correltions when
#   percentage of last 300 attempts are below the limit
#
    a100_lim=1
    n100_lim=3
#
#   APASS and NOMAD magnitude correction factors to roughly match DES
#   These are given in the sense that they would be added to the APASS/NOMAD catalog value
#   Values *_mag_corr (assume no extinction correction)
#   Values *_kmag_corr (assume a k*airmass correction)
#
    cat_mag_corr={"apass":{"u":3.5,"g":0.205,"r":0.128,"i":0.112,"z":0.0,"Y":0.0,"VR":0.0},
                  "nomad":{"u":3.65,"g":0.341,"r":0.235,"i":1.398,"z":1.201,"Y":1.083,"VR":0.0}}
    cat_kmag_corr={"apass":{"u":0.0,"g":0.000,"r":0.000,"i":0.000,"z":0.0,"Y":0.0,"VR":0.0},
                  "nomad":{"u":0.0,"g":0.111,"r":0.109,"i":1.289,"z":1.139,"Y":1.022,"VR":0.0}}
#
#   NOMAD B-mag correction
#
    tmp_b_corr=[-0.24,-0.24,-0.25,-0.25,-0.32,-0.37,-0.37,-0.35,-0.28,-0.33,-0.33,-0.34,-0.30,-0.28,-0.32,
        -0.32,-0.31,-0.27,-0.25,-0.25,-0.25,-0.22,-0.27,-0.27,-0.28,-0.30,-0.28,-0.23,-0.25,-0.23,
        -0.20,-0.19,-0.20,-0.23,-0.23,-0.21,-0.19,-0.21,-0.28,-0.25,-0.22,-0.21,-0.25,-0.27,-0.22,
        -0.23,-0.18,-0.25,-0.30,-0.26,-0.26,-0.23,-0.29,-0.31,-0.24,-0.20,-0.19,-0.28,-0.32,-0.26,
        -0.21,-0.22,-0.22,-0.27,-0.24,-0.22,-0.19,-0.09,0.15,0.12,0.14,0.19,0.23,0.24,0.17,
        0.12,0.18,0.25,0.27,0.21,0.18,0.23,0.27,0.26,0.19,0.18,0.21,0.22,0.06,0.01,
        0.05,0.10,0.17,0.15,0.05,0.11,0.19,0.20,0.11,0.06,0.11,0.19,0.18,0.10,0.05,
        0.06,0.14,0.17,0.10,0.06,0.09,0.16,0.20,0.12,0.13,0.11,0.15,0.19,0.02,-0.03,
        0.02,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
        0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
        0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
        0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00]
    nomad_bmag_corr=numpy.array(tmp_b_corr,dtype=numpy.float32)
#
#   (IN case for some reason I want to go back to reading from a file)
#
#    bmedian_file='/work/devel/rgruendl/cal_test/APASS/NOMAD.bmedian.txt'
#
#   nomad_bmag_corr=numpy.zeros(180,dtype=numpy.float32)
#   if (os.path.isfile(bmedian_file)):
#       f1=open(bmedian_file,'r')
#       for line in f1:
#           line=line.strip()
#           columns=line.split(' ')
#           if (columns[0] != "#"):
#               spd=90+int(columns[0])
#               nomad_bmag_corr[spd]=float(columns[1])
##              print columns[0],columns[1],spd,nomad_bmag_corr[spd]
#       f1.close()

#
#   A set of aterm(s) that represent an estimate for the nominal zeropoint on a clear, photometric night.
#
    aterm=numpy.zeros(7,dtype=numpy.float32)
    aterm[0]=-25.0
    aterm[1]=-25.428
    aterm[2]=-25.532
    aterm[3]=-25.413
    aterm[4]=-25.086
    aterm[5]=-24.000
    aterm[6]=-25.47
#
#   A set of kterm(s) that represent an estimate for the nominal extinction/airmass correction on a clear, photometric night.
#
    kterm=numpy.zeros(7,dtype=numpy.float32)
    kterm[0]=0.489
    kterm[1]=0.181
    kterm[2]=0.095
    kterm[3]=0.089
    kterm[4]=0.053
    kterm[5]=0.05
    kterm[6]=0.13

    print("####################################################")
    print("# Initialized. ")
##############################################################################
##############################################################################
##############################################################################
#   Initial Query(s) to determine information about the specific processing attempt for the requested item.
#   Starts with obtaining a list of catalogs (either by query or from an input list)
#
    t0=time.time()
    exp_rec={}
    cat_fname=[]
    cat_fpath=[]
    if (args.list is None):
#
#       If no list of catalogs was provided (--list) then go figure out whether they exist based on the PFW_Attempt_ID
#
        query = """
            select c.filename as filename,
                c.filetype as filetype,
                fai.path as path,
                oa.root as root,
                c.ccdnum as ccdnum
            from {schema:s}file_archive_info fai, {schema:s}ops_archive oa, {schema:s}catalog c 
            where c.pfw_attempt_id={aid:} 
                and (c.filetype='cat_finalcut' or c.filetype='cat_firstcut')
                and fai.filename=c.filename 
                and oa.name=fai.archive_name 
            """.format(schema=db_Schema,aid=PFW_Attempt_ID)

        if (args.verbose):
            print "# Executing initial query for catalogs "
            if (args.format_query):
                print("# example query = {:s}".format(query))
            else:
                print("# example query ={:s} ".format(" ".join([d.strip() for d in query.split('\n')])))

        cur.arraysize = 1000 # get 1000 at a time when fetching
        cur.execute(query)
        desc = [d[0].lower() for d in cur.description]
        for row in cur:
            rowd = dict(zip(desc, row))
            cat_fname.append(rowd['filename'])
            cat_fpath.append(os.path.join(rowd['root'],rowd['path'],rowd['filename']))
    else:
#
#       Open list file and proceed.
#
        if (args.verbose):
            print "# Using list of files/catalogs to obtain objects from exposure"
        if (os.path.isfile(args.list)):
            f_cat=open(args.list,'r')
            for line in f_cat:
                line=line.strip()
                columns=line.split(' ')
                if (columns[0] != "#"):
                    tmp_fname=columns[0].split('/')
                    cat_fname.append(tmp_fname[-1])
                    cat_fpath.append(columns[0])
#                   print columns[0]
            f_cat.close()
    if (args.verbose):
        print "# Total number of catalogs/CCDs identified: ",len(cat_fpath)
        print "# "

##############################################################################
#   Now that a starting point has been established...
#   Obtain associated DB information for each catalog/CCD.
# 
    ccd_info={}
    num_ccd_verbose=0
    max_ccd_verbose=1
    if (args.verbose):
        print "# Executing query(s) to obtain image level metadata and QA values "

    for index, tmp_fname in enumerate(cat_fname):
        query = """
            select i.filename as img_fname,
                i.expnum,
                i.band,
                i.ccdnum,
                i.airmass,
                i.exptime,
                i.fwhm as fwhm_old,
                i.elliptic as elli_old,
                i.skybrite,
                i.skysigma,
                i.gaina,i.gainb,
                i.ra_cent,i.dec_cent,
                i.rac1,i.rac2,i.rac3,i.rac4,
                i.decc1,i.decc2,i.decc3,i.decc4
            from {schema:s}image i, {schema:s}catalog c 
            where i.filetype='red_immask' 
                and c.filename='{fname:s}' 
                and c.pfw_attempt_id=i.pfw_attempt_id 
                and c.expnum=i.expnum 
                and c.ccdnum=i.ccdnum 
            """.format(schema=db_Schema,fname=cat_fname[index])

        if ((args.verbose)and(num_ccd_verbose < max_ccd_verbose)):
            if (args.format_query):
                print("# example query = {:s}".format(query))
            else:
                print("# example query ={:s} ".format(" ".join([d.strip() for d in query.split('\n')])))
            num_ccd_verbose=num_ccd_verbose+1

        cur.arraysize = 1000 # get 1000 at a time when fetching
        cur.execute(query)
        desc = [d[0].lower() for d in cur.description]

        for row in cur:
            rowd = dict(zip(desc, row))

            ccdnum=int(rowd['ccdnum'])
            ccd_info[ccdnum]=rowd
#
#           Add entries for filename, filepath, iband 
#           Check entries that are sometimes malformed or null
#
            ccd_info[ccdnum]["cat_fname"]=cat_fname[index]
            ccd_info[ccdnum]["cat_fpath"]=cat_fpath[index]
            if (ccd_info[ccdnum]['band'] is None):
                ccd_info[ccdnum]['band']='Unknown'
                ccd_info[ccdnum]['iband']=-1
            else:
                ccd_info[ccdnum]['iband']=band2i[ccd_info[ccdnum]['band']]

            if (ccd_info[ccdnum]['gaina'] is None):
                ccd_info[ccdnum]['gaina']=-1.
                print "# Warning: Null found for GAINA"

            if (ccd_info[ccdnum]['gainb'] is None):
                ccd_info[ccdnum]['gainb']=-1.
                print "# Warning: Null found for GAINB"

            if (ccd_info[ccdnum]['airmass'] is None):
                ccd_info[ccdnum]['airmass']=1.2
                print "# Warning: Null found for AIRMASS (default to 1.2)"

            if ((ccd_info[ccdnum]['rac1'] is None)or(ccd_info[ccdnum]['rac2'] is None)or
                (ccd_info[ccdnum]['rac3'] is None)or(ccd_info[ccdnum]['rac4'] is None)):
                print "# Warning: Null found among RAC1-4"
                ccd_info[ccdnum]['img_ra']=[]
            else:
                ccd_info[ccdnum]['img_ra']=[ccd_info[ccdnum]['rac1'],ccd_info[ccdnum]['rac2'],ccd_info[ccdnum]['rac3'],ccd_info[ccdnum]['rac4']]
#
            if ((ccd_info[ccdnum]['decc1'] is None)or(ccd_info[ccdnum]['decc2'] is None)or
                (ccd_info[ccdnum]['decc3'] is None)or(ccd_info[ccdnum]['decc4'] is None)):
                print "# Warning: Null found among DECC1-4"
                ccd_info[ccdnum]['img_dec']=[]
            else:
                ccd_info[ccdnum]['img_dec']=[ccd_info[ccdnum]['decc1'],ccd_info[ccdnum]['decc2'],ccd_info[ccdnum]['decc3'],ccd_info[ccdnum]['decc4']]

#            ccd_info[ccdnum]["exptime"]=float(item[coldict["i.exptime"]])
            ccd_info[ccdnum]["fwhm_old"]=pixsize*float(ccd_info[ccdnum]["fwhm_old"])
#            ccd_info[ccdnum]["elli_old"]=float(item[coldict["i.elliptic"]])
#
#           sbrite_good and sbrite_lim are defined in cts/sec
#           Therefore:
#              If exposure time is present then use to normalize brightness to a per second quantity
#              If gains are present then use to express brightness in counts (rather than electrons)
#
            if (ccd_info[ccdnum]['exptime']>0.01):
                efactor=ccd_info[ccdnum]['exptime']
            else:
                efactor=1.0

            gtesta=ccd_info[ccdnum]['gaina']-1.
            gtestb=ccd_info[ccdnum]['gainb']-1.
            if ((abs(gtesta)<0.5)and(abs(gtestb)<0.5)):
#               The case where gains are 1... therefore units are electrons
                gfactor=4.0
                ccd_info[ccdnum]['bunit']='e-'
                if (args.verbose):
                    print("# GAINA/B for ccdnum={:d} are consistent with units of electrons (GAINA:{:.3f},GAINB:{:.3f})".format(
                        ccdnum,ccd_info[ccdnum]['gaina'],ccd_info[ccdnum]['gainb']))
            else:
#               The case where gains are not 1... therefore units are already in counts
                gfactor=1.0
                ccd_info[ccdnum]['bunit']='DN'
                if (args.verbose):
                    print("# GAINA/B for ccdnum={:d} are consistent with units of DN (GAINA:{:.3f},GAINB:{:.3f})".format(
                        ccdnum,ccd_info[ccdnum]['gaina'],ccd_info[ccdnum]['gainb']))
            ccd_info[ccdnum]['skyb']=ccd_info[ccdnum]['skybrite']/efactor/gfactor
            ccd_info[ccdnum]['skys']=ccd_info[ccdnum]['skysigma']/efactor/gfactor

#            if (ccd_info[ccdnum]['scampflg'] is None):
#                ccd_info[ccdnum]['sflag']=2
#            else:
#                ccd_info[ccdnum]['sflag']=ccd_info[ccdnum]['scampflg']

##############################################################################
#   Fill in entries for missing CCDs (if there are any)
#   Make consistency checks across CCDs to check that once fragile bits are consistent:
#       specific checks are: expnum, band, exptime, airmass, bunit 
#   Form 1-d arrays of RA, DEC (corners) for use when cataloging.
#   From 1-d arrays of skyb and skys (for eventual b_eff calculation)
#   Historically 1-d arrays of fwhm and ellipticity [fwhm_old, elli_old]
#
    bunit_chk=[ccd_info[k]['bunit'] for k in ccd_info]
    band_chk=[ccd_info[k]['band'] for k in ccd_info]
    expnum_chk=[ccd_info[k]['expnum'] for k in ccd_info]
    exptime_chk=[ccd_info[k]['exptime'] for k in ccd_info]
    airmass_chk=[ccd_info[k]['airmass'] for k in ccd_info]

    img_ra=numpy.reshape(numpy.array([ccd_info[k]['img_ra'] for k in ccd_info]),-1)
    img_dec=numpy.reshape(numpy.array([ccd_info[k]['img_dec'] for k in ccd_info]),-1)
    skyb=numpy.array([ccd_info[k]['skyb'] for k in ccd_info])
    skys=numpy.array([ccd_info[k]['skys'] for k in ccd_info])
    fwhm_old=numpy.array([ccd_info[k]['fwhm_old'] for k in ccd_info])
    elli_old=numpy.array([ccd_info[k]['elli_old'] for k in ccd_info])

    exp_rec['numccd']=len(expnum_chk)

    print len(exptime_chk)

##############################################################################
#   Perform the actual checks (and set exposure level values when these pass muster.
#   First check exptime.
#
    uniq_exptime_chk=list(set(exptime_chk))
    uniq_band_chk=list(set(band_chk))
    uniq_bunit_chk=list(set(bunit_chk))
    uniq_expnum_chk=list(set(expnum_chk))
    uniq_airmass_chk=list(set(airmass_chk))
    print len(uniq_exptime_chk)
    if (len(uniq_exptime_chk) != 1):
        if (len(uniq_exptime_chk) > 1):
            print "WARNING: Other than one exptime?: ",uniq_exptime_chk
            print "WARNING: Using exptime: ",uniq_exptime_chk[0]
            exp_rec['exptime']=uniq_exptime_chk[0]
        else:
            print "No exptime? Setting to NoneType will try to obtain value from EXPOSURE"
            exp_rec['exptime']=None
    else:
        exp_rec['exptime']=uniq_exptime_chk[0]
#
#   Now check airmass
#
    if (len(uniq_airmass_chk) != 1):
        if (len(uniq_airmass_chk) > 1):
            print "WARNING: Other than one airmass?: ",uniq_airmass_chk
            print "WARNING: Using airmass: ",uniq_airmass_chk[0]
            exp_rec['airmass']=uniq_airmass_chk[0]
        else:
            print "WARNING: No airmass?: "
            print "WARNING: Using defualt setting airmass: 1.2"
            exp_rec['airmass']=1.2
    else:
        exp_rec['airmass']=uniq_airmass_chk[0]
#
#   Now check bunit for consistency
#
    if (len(uniq_bunit_chk) != 1):
        if (len(uniq_bunit_chk) > 1):
            print "WARNING: Other than one bunit?: ",uniq_bunit_chk
            print "WARNING: Using bunit: ",uniq_bunit_chk[0]
            exp_rec['bunit']=uniq_bunit_chk[0]
        else:
            print "WARNING: Assuming bunit = DN"
            exp_rec['bunit']='DN'
    else:
        exp_rec['bunit']=uniq_bunit_chk[0]
#
#   Now check band (also check that it is a sanctioned value)
#
    if (len(uniq_band_chk) != 1):
        if (len(uniq_band_chk)==0):
            print "No band? Really?  Setting to NoneType will try to obtain value from EXPOSURE"
            exp_rec['band']=None
        else:
            print "Abort: Other than one band identified?: ",uniq_band_chk
            exit(1)
    else:
#        if (uniq_band_chk[0] in ['u','g','r','i','z','Y','VR']):
        if (uniq_band_chk[0] in band2i):
            exp_rec['band']=uniq_band_chk[0]
        else:
            print "Abort: Unsupported value for band?: ",uniq_band_chk[0]
            exit(1)
#
#   Now check expnum
#
    if (len(list(set(expnum_chk))) != 1):
        print "Abort: Other than one expnum?: ",uniq_expnum_chk
        exit(1)
    else:
        exp_rec['expnum']=uniq_expnum_chk[0]
##############################################################################
#   Calculated exposure level values (old FWHM, old Ellipticity, Sky Brightness and Sky Sigma)
#
    if (len(fwhm_old) > 0):
        exp_rec["fwhm_img"]=fwhm_old.mean()
    else:
        exp_rec["fwhm_img"]=-1.0
    if (len(elli_old) > 0):
        exp_rec["ellip_avg"]=elli_old.mean()
        exp_rec["ellip_rms"]=elli_old.std()
    else:
        exp_rec["ellip_avg"]=-1.0
        exp_rec["ellip_rms"]=-1.0
    if (len(skyb) > 0):
        exp_rec["skyb_avg"]=skyb.mean()
        exp_rec["skyb_rms"]=skyb.std()
    else:
        exp_rec["skyb_avg"]=-1.0
        exp_rec["skyb_rms"]=-1.0
    if (len(skys) > 0):
        exp_rec["skys_avg"]=skys.mean()
        exp_rec["skys_rms"]=skys.std()
    else:
        exp_rec["skys_avg"]=-1.0
        exp_rec["skys_rms"]=-1.0

###############################################################################
#   Get more exposure level information
#
    query = """
        select e.filename as exp_fname,
            e.expnum,
            e.nite,
            e.date_obs, e.time_obs, e.mjd_obs,
            e.obstype, 
            e.band,
            e.camsym,
            e.object,
            e.telra, e.teldec,
            e.tradeg, e.tdecdeg,
            e.airmass,
            e.exptime,
            e.program,
            e.field
        from {schema:s}exposure e
        where e.expnum={expnum:}
        """.format(schema=db_Schema,expnum=exp_rec['expnum'])

    if (args.verbose):
        print "# Executing query to obtain exposure level metadata"
        if (args.format_query):
            print("# query = {:s}".format(query))
        else:
            print("# query ={:s} ".format(" ".join([d.strip() for d in query.split('\n')])))

    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)
    desc = [d[0].lower() for d in cur.description]

#    print desc
#
#    for key in exp_rec:
#        print key,exp_rec[key]

    num_exp_recs=0
    for row in cur:
        num_exp_recs=num_exp_recs+1
        if (num_exp_recs > 1):
            print "# WARNING: multiple exposure records found for this exposure? (num_exp_recs=",num_exp_recs,")"
        rowd = dict(zip(desc, row))
#       The following seem to be OK to grab wholesale
        for key in ['nite','exp_fname','date_obs','obstype','telra','teldec','tradeg','tdecdeg','object','mjd_obs']:
            exp_rec[key]=rowd[key]
#
#       And then there are some that need more sanity checks.
#
        if (rowd['band'] is None):
            if (exp_rec['band'] is None):
                print("All attempts to acquire BAND have failed.")
                print("Abort!")
                exit(0)
            if (exp_rec["band"] != 'Unknown'):
                print("WARNING: BAND miss-match between exposure-level (Unknown) and image/catalog-level ({:s}) queries.  Using image/cat-result.".format(exp_rec['band']))
        else:
            if (exp_rec['band'] is None):
                if (rowd['band'] in band2i):
                    exp_rec['band']=rowd['band']
                else:
                    print "Abort: Unsupported value for band?: ",rowd['band']
                    exit(1)
            if (rowd['band'] != exp_rec['band']):
                print("WARNING: BAND miss-match between exposure-level ({:s}) and image/catalog-level ({:s}) queries.  Using image/cat-result.".format(rowd['band'],exp_rec['band']))
#
        if (exp_rec['exptime'] is None):
            print("Fixing exptime using value from EXPOSURE")
            exp_rec['exptime']=rowd['exptime']
        else:
            if (rowd['exptime'] != exp_rec["exptime"]):
                print("WARNING: EXPTIME miss-match between exposure-level ({:.1f}) and image/catalog-level ({:.1f}) queries.  Using image/cat-result.".format(rowd['exptime'],exp_rec['exptime']))

        if (rowd['airmass'] is None):
           print("WARNING: AIRMASS miss-match between exposure-level (Unknown) and image/catalog-level ({:.3f}) queries.  Using image/cat-result.".format(exp_rec['airmass']))
        else:
            if (rowd['airmass'] != exp_rec['airmass']):
                print("WARNING: AIRMASS miss-match between exposure-level ({.3f}) and image/catalog-level ({:.3f}) queries.  Using image/cat-result.".format(rowd['airmass'],exp_rec['airmass']))
#
#       Work out whether the exposure is part of the general survey, SN, or other.
#
        if (rowd['program']=="survey"):
            exp_rec["program"]='survey'
            exp_rec["sn_field"]=False
            exp_rec["survey_field"]=True
        elif (rowd['program']=="supernova"):
            exp_rec["program"]='SN'
            exp_rec["sn_field"]=True
            exp_rec["survey_field"]=False
        elif (rowd['program']=="photom-std-field"):
            exp_rec["program"]='phot-std'
            exp_rec["sn_field"]=False
            exp_rec["survey_field"]=False
        else:
            if (exp_rec['obstype'] in ['zero','dark','dome flat','sky flat']):
                exp_rec["program"]='cal'
                exp_rec["sn_field"]=False
                exp_rec["survey_field"]=False
            else:
                exp_rec["program"]='unknown'
                exp_rec["sn_field"]=False
                exp_rec["survey_field"]=False

###############################################################################
#   Query to get SCAMP_QA information
#
    query = """
        select q.astromndets_ref_highsn as astrom_ndets, 
            q.astromchi2_ref_highsn as astrom_chi2, 
            q.astromsigma_ref_highsn_1 as astrom_sig1,
            q.astromsigma_ref_highsn_2 as astrom_sig2,
            q.astromoffset_ref_highsn_1 as astrom_off1, 
            q.astromoffset_ref_highsn_2 as astrom_off2
        from {schema:s}scamp_qa q, {schema:s}miscfile m 
        where m.pfw_attempt_id={aid:s} 
            and m.filetype='xml_scamp' 
            and m.filename=q.filename 
        """.format(schema=db_Schema,aid=PFW_Attempt_ID)

    if (args.verbose):
        print "# Executing query to obtain astrometry QA for this exposure"
        if (args.format_query):
            print("# query = {:s}".format(query))
        else:
            print("# query ={:s} ".format(" ".join([d.strip() for d in query.split('\n')])))

    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)
    desc = [d[0].lower() for d in cur.description]

    num_ast_recs=0
    for row in cur:
        num_ast_recs=num_ast_recs+1
        if (num_ast_recs > 1):
            print("# WARNING: multiple exposure records found for this exposure? (num_ast_recs={:d})".format(num_ast_recs))
        rowd = dict(zip(desc, row))
        for fld in ['astrom_ndets','astrom_chi2','astrom_sig1','astrom_sig2','astrom_off1','astrom_off2']:
            exp_rec[fld]=rowd[fld]

#
#   In case no record was found.
#
    if (num_ast_recs == 0):
        exp_rec["astrom_sig1"]=None
        exp_rec["astrom_sig2"]=None
        exp_rec["astrom_off1"]=None
        exp_rec["astrom_off2"]=None
        exp_rec["astrom_rms2"]=None
        exp_rec["astrom_ndets"]=None
        exp_rec["astrom_chi2"]=None

#
#   Check to see whether this looks like an astrometric failure
#
    astrom_good=True
    for fld in ['astrom_ndets','astrom_chi2','astrom_sig1','astrom_sig2','astrom_off1','astrom_off2']:
        if (exp_rec[fld] is None):
            print("# WARNING: Probable astrometric solution failure: No/null value for {:s} in scamp_qa".format(fld))
            astrom_good=False
            if (fld == 'astrom_ndets'):
                exp_rec[fld]=-1
            else:
                exp_rec[fld]=-1.0

    sigx=exp_rec['astrom_sig1']
    sigy=exp_rec['astrom_sig2']
    exp_rec["astrom_rms2"]=numpy.sqrt((sigx*sigx)+(sigy*sigy))

    if ((exp_rec["astrom_sig1"] < 0.001)or(exp_rec["astrom_sig2"] < 0.001)):
        astrom_good=False
        print("# WARNING: Probable astrometric solution failure: astrom_sig1,2 = {:7.4f},{:7.4f}".format(
            exp_rec["astrom_sig1"],exp_rec["astrom_sig2"]))
    if (exp_rec["astrom_rms2"] > 0.500):
        astrom_good=False
        print("# WARNING: Probable astrometric solution failure: astrom_rms2 = {:.3f}".format(
            exp_rec["astrom_rms2"]))

###############################################################################
#   Query to get PSF_QA information
#
    query = """
        select m.expnum as expnum,
            count(m.ccdnum) as count,
            median(q.fwhm_mean) as med_fwhm
        from {schema:s}psf_qa q, {schema:s}miscfile m 
        where m.pfw_attempt_id={aid:} 
            and m.filetype='xml_psfex' 
            and m.filename=q.filename 
        group by m.expnum 
        """.format(schema=db_Schema,aid=PFW_Attempt_ID)

    if (args.verbose):
        print "# Executing query to obtain PSF QA for this exposure"
        if (args.format_query):
            print("# query = {:s}".format(query))
        else:
            print("# query ={:s} ".format(" ".join([d.strip() for d in query.split('\n')])))

    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)
    desc = [d[0].lower() for d in cur.description]

    psf_good=True
    num_psf_recs=0
    for row in cur:
        num_psf_recs=num_psf_recs+1
        if (num_psf_recs > 1):
            print "# WARNING: multiple exposure records found for PSF QA for this exposure/pfw_attempt_id? (num_psf_recs=",num_psf_recs,")"
        rowd = dict(zip(desc, row))
        exp_rec['num_psfex']=rowd['count']
        exp_rec['psfex_fwhm']=pixsize*rowd['med_fwhm']
    
    if ('psfex_fwhm' not in exp_rec):
        exp_rec['psfex_fwhm']=-1.0

################################################
#   Open generic output file (changed to STDOUT)
#    ftxt=open("%s.txt"% (args.froot), 'w')
    sys.stdout.flush()
    ftxt = sys.stdout 

    t1=time.time()
    print("# TIMING (general queries): {:.2f}".format((t1-t0)))
###############################################################################
#
#    Prepare for optional CSV output if args.csv is TRUE
#
    if (args.csv):
        fcsv = open("%s.csv"% (opts.froot), 'w')
        writer = csv.writer(fcsv, delimiter=opts.delimiter,  quotechar='"', quoting=csv.QUOTE_MINIMAL)
        out_row=[]
        out_row.append("expid")
        out_row.append("dm_done")
        out_row.append("dm_accept")
        out_row.append("comments")
        writer.writerow(out_row)
#
#    Prepare for optional DB update file output if args.DB_file is TRUE
#
    if (args.DB_file):
        fdbout=open("%s.dbupdate"% (args.froot), 'w')

###############################################################################
#
#   Find nominal center and range spanned by focal plane for catalog search/comparison
#
    if ((len(img_ra) < 2)or(len(img_dec) < 2)):
        print "Warning: No CCD corners present in database for these items"
        exp_rec['ra_cen']=exp_rec['tradeg']
        exp_rec['dec_cen']=exp_rec['tdecdeg']
    else:
        if ((len(img_ra) < (4*exp_rec["numccd"]))or(len(img_dec) < (4*exp_rec["numccd"]))):
            print "Warning: Some CCDs may be missing RA/DEC corners)! Working with what was present..."
        ra_min=numpy.min(img_ra)
        ra_max=numpy.max(img_ra)
        dec_min=numpy.min(img_dec)
        dec_max=numpy.max(img_dec)
        if ((ra_max-ra_min)>180.0):
#           print "IT HAPPENS ",ra_min,ra_max
            new_img_ra=[]
            for ra_val in img_ra:
                if (ra_val > 180.0):
                    new_img_ra.append(ra_val-360.0)
                else:
                    new_img_ra.append(ra_val)
            img_ra=numpy.array(new_img_ra)
        ra_cen=img_ra.mean()
        dec_cen=img_dec.mean()
        if (ra_cen < 0.0):
            ra_cen=ra_cen+360.0
        exp_rec['ra_cen']=ra_cen
        exp_rec['dec_cen']=dec_cen

#
#   Denote whether or not this is a special case where (where exposure straddles RA=0h
#
    near_branch=False
    if ((ra_cen < fp_rad)or((ra_cen > 180.0)and(ra_cen > (360.-fp_rad)))):
       near_branch=True

#
#   Report exposure information for current 
#
    print "########################################"
    print "#   date_obs: ",exp_rec['date_obs']
    print "#       Nite: ",exp_rec['nite']
    print "#   Exposure: ",exp_rec['expnum']
    print "#       Band: ",exp_rec['band']
    print("#    Exptime: {:.1f} ".format(exp_rec['exptime']))
    print "#      BUNIT: ",exp_rec['bunit']
    print "# "
    print "#    Obstype: ",exp_rec['obstype']
    print "#     Object: ",exp_rec['object']
    print "#    Program: ",exp_rec['program']
    print "# "
    print("#      Telescope(Ra,Dec): {:9.5f} {:9.5f} ".format(exp_rec['tradeg'],exp_rec['tdecdeg']))
    print "#"
    print("# Image Centroid(Ra,Dec): {:9.5f} {:9.5f} ".format(exp_rec['ra_cen'],exp_rec['dec_cen']))
    print "#"
    print "# "
    if (astrom_good):
        print "# Astrometric solution appears OK "
    else:
        print "# WARNING: Probable astrometric solution failure"
    print "# Astrometry summary:"
    print("#                    high S/N      ")
    print("#        ndets: {:d} ".format(exp_rec['astrom_ndets']))
    print("#         chi2: {:.2f} ".format(exp_rec['astrom_chi2']))
    print("#   astrom_sig: {:7.4f},{:7.4f} ".format(exp_rec['astrom_sig1'],exp_rec['astrom_sig2']))
    print("#   astrom_off: {:7.4f},{:7.4f} ".format(exp_rec['astrom_off1'],exp_rec['astrom_off2']))
    print "#"
    print "# PSF summary:"
    print("# median(FWHM): {:.3f} ".format(exp_rec['psfex_fwhm']))
    print("#      num CCD: {:d} ".format(exp_rec['num_psfex']))
    print "#"
    print "########################################"


###############################################################################
#   Perform search on APASS and/or NOMAD catalog to get entries for this exposure
#
#   Steve Kent's wisdom
#   The comparison catalogs I settled on are as follows:
#    1. For izY, I use NOMAD J magnitudes (from 2MASS) exclusively
#    2. For g and r, I try APASS on 6 CCDs.  If I find enough stars, I use those.
#    3. If not, I use NOMAD as follows:
#       a) For g, I use NOMAD B mags with a declination-dependent correction to put them close to the APASS system.
#       b) For r, I use NOMAD (2*B+J)/3 but do not include the declination-dependent corrections.
#   The offsets from DECam to the various catalogs are determined by medianing the offsets from 16 nights of data.
#
#   I spent a fair amount of time making the matching algorithms robust to various types of errors and failures
#   (e.g., center command failing, junk stars in the catalogs)  but there always seem to be a few edge cases.
#
    print "# Preparing queries for NOMAD and APASS in vicinity of the exposure."
#
#   Setup rules for parsing NOMAD/APASS 
#
    nomad_parse={"ra":"n.ra","dec":"n.dec"}
    if (exp_rec["band"] in ["u","g"]):
        nomad_parse["mag"]="n.b"
        nomad_parse["mlimit"]=blimit
        nomad_parse["simple"]=False
        nomad_parse['keys']=['ra','dec','mag']
    elif (exp_rec["band"] == "r"):
        nomad_parse["mag"]="n.b"
        nomad_parse["mag2"]="n.j"
        nomad_parse["mlimit"]=blimit
        nomad_parse["mlimit2"]=jlimit
        nomad_parse["simple"]=False
        nomad_parse['keys']=['ra','dec','mag','mag2']
    elif (exp_rec["band"] in ["i","z","Y"]):
        nomad_parse["mag"]="n.j"
        nomad_parse["mlimit"]=jlimit
        nomad_parse["simple"]=True
        nomad_parse['keys']=['ra','dec','mag']

    apass_parse={"ra":"a.ra","dec":"a.dec"}
    if (exp_rec["band"] in ["u","g"]):
        apass_parse["mag"]="a.g"
        apass_parse["dmag"]="a.g_err"
        apass_parse["mlimit"]=glimit
        apass_parse['keys']=['ra','dec','mag','dmag']
    elif (exp_rec["band"] in ["r","VR"]):
        apass_parse["mag"]="a.r"
        apass_parse["dmag"]="a.r_err"
        apass_parse["mlimit"]=rlimit
        apass_parse['keys']=['ra','dec','mag','dmag']
    elif (exp_rec["band"] in ["i","z","Y"]):
        apass_parse["mag"]="a.i"
        apass_parse["dmag"]="a.i_err"
        apass_parse["mlimit"]=ilimit
        apass_parse['keys']=['ra','dec','mag','dmag']

#
#   Form query lists (with as to facilitate roughly homogenous catalogs
#
    qlist_nomad=[]
    for key in nomad_parse['keys']:
        qlist_nomad.append("%s as %s" % (nomad_parse[key],key))
    querylist_nomad = ",".join(qlist_nomad)
#    print qlist_nomad

    qlist_apass=[]
    for key in apass_parse['keys']:
        qlist_apass.append("%s as %s" % (apass_parse[key],key))
    querylist_apass = ",".join(qlist_apass)
#    print qlist_apass

#
#   From magnitude constraints for queries
#
    if ('mag2' in nomad_parse):
        nomad_mag_constraint=' and {:s}<{:.1f} and {:s}<{:.1f}'.format(
            nomad_parse['mag'],nomad_parse['mlimit'],nomad_parse['mag2'],nomad_parse['mlimit2'])
    else:
        nomad_mag_constraint=' and {:s}<{:.1f}'.format(nomad_parse['mag'],nomad_parse['mlimit'])

    apass_mag_constraint=' and {:s}<{:.1f}'.format(apass_parse['mag'],apass_parse['mlimit'])
    if (args.verbose):
        print("# Applying NOMAD MAG CONSTRAINT ({:s}) for work with {:s}-band data".format(nomad_mag_constraint,exp_rec["band"]))
        print("# Applying APASS MAG CONSTRAINT ({:s}) for work with {:s}-band data".format(apass_mag_constraint,exp_rec["band"]))
#
#   Form queries for APASS/NOMAD objects 
#
    if (near_branch):
#
#       NOTE need something better if we get closer than 2 degrees from pole
#
        cosdec_cen=numpy.cos(dec_cen*deg2rad)
        ra1=ra_cen+(fp_rad/cosdec_cen)
        if (ra1 > 360.0):
            ra1=ra1-360.0
        ra2=ra_cen-(fp_rad/cosdec_cen)
        if (ra2 < 0.0):
            ra2=ra2+360.0
        dec1=dec_cen-fp_rad
        dec2=dec_cen+fp_rad
        query_nomad = """
            select {qlist:s} 
            from nomad n 
            where n.dec between {d1:.7f} and {d2:.7f} 
                and (n.ra < {r1:.7f} or n.ra > {r2:.7f})
                {mcon:s} 
            order by n.ra 
            """.format(qlist=querylist_nomad,d1=dec1,d2=dec2,r1=ra1,r2=ra2,mcon=nomad_mag_constraint)

        query_apass = """
            select {qlist:s} 
            from apass_dr7 a 
            where a.dec between {d1:.7f} and {d2:.7f} 
                and (a.ra < {r1:.7f} or a.ra > {r2:.7f})
                {mcon:s} 
            order by a.ra 
            """.format(qlist=querylist_apass,d1=dec1,d2=dec2,r1=ra1,r2=ra2,mcon=apass_mag_constraint)
    else:
#
#       NOTE need something better if we get closer than 2 degrees from pole
#
        cosdec_cen=numpy.cos(dec_cen*deg2rad)
        ra1=ra_cen-(fp_rad/cosdec_cen)
        ra2=ra_cen+(fp_rad/cosdec_cen)
        dec1=dec_cen-fp_rad
        dec2=dec_cen+fp_rad
        query_nomad = """
            select {qlist:s} 
            from nomad n 
            where n.dec between {d1:.7f} and {d2:.7f} 
                and n.ra between {r1:.7f} and {r2:.7f} 
                {mcon:s}
            order by n.ra
            """.format(qlist=querylist_nomad,d1=dec1,d2=dec2,r1=ra1,r2=ra2,mcon=nomad_mag_constraint)
        query_apass = """
            select {qlist:s} 
            from apass_dr7 a 
            where a.dec between {d1:.7f} and {d2:.7f}
                and a.ra between {r1:.7f} and {r2:.7f} 
                {mcon:s} 
            order by a.ra
             """.format(qlist=querylist_apass,d1=dec1,d2=dec2,r1=ra1,r2=ra2,mcon=apass_mag_constraint)
    if args.verbose:
        print "# Executing queries for NOMAD and APASS in vicinity of the exposure."
        print("# {:s}".format(query_nomad))
        print("# {:s}".format(query_apass))
    cur.arraysize = 1000 # get 1000 at a time when fetching

    t2=time.time()
    print("# TIMING (pre-NOMAD/APASS): {:.2f}".format((t2-t1)))
#
#   Acquire APASS data (currently always)
#
    cur.execute(query_apass)
    desc = [d[0].lower() for d in cur.description]

    apass_cat=[]
    for row in cur:
        rowd = dict(zip(desc, row))
        if (rowd['mag'] > apass_parse['mlimit']):
            rowd['mag']=99.0
            rowd['dmag']=99.0
        if (rowd['mag']<98.0):
            apass_cat.append(rowd)

    print "# "
    print "# Number of APASS entries found (after magnitude cuts) is: ",len(apass_cat)
    t2a=time.time()
#
#   Acquire NOMAD data (currently always)
#
    cur.execute(query_nomad)
    desc = [d[0].lower() for d in cur.description]

    nomad_cat=[]
    for row in cur:
        rowd = dict(zip(desc, row))

        if (nomad_parse["simple"]):
            if (rowd['mag'] > nomad_parse['mlimit']):
                rowd['mag']=99.0
        else:
#
#           When B magnitudes are involved then corrections (excissions are needed)
#           When red magnitudes are involved then nasty extra extrapolations also occur
#
            if (exp_rec["band"] in ["u","g"]):
                bmag=rowd['mag']+nomad_bmag_corr[int(rowd['dec']+90.0)]
                if (bmag > nomad_parse['mlimit']):
                    rowd['mag']=99.0
                else:
                    if ((bmag > 18.37)and(bmag < 18.39)):
                        rowd['mag']=99.0
                    else:
                        rowd['mag']=bmag
            elif (exp_rec["band"] in ["r","VR"]):
                bmag=rowd['mag']+nomad_bmag_corr[int(rowd['dec']+90.0)]
                jmag=rowd['mag2']
                if ((jmag > nomad_parse['mlimit2'])or(bmag> nomad_parse['mlimit'])):
                    rowd['mag']=99.0
                else:
                    if ((bmag > 18.37)and(bmag < 18.39)):
                        rowd['mag']=99.0
                    else:
                        rowd['mag']=0.333333333*((2.0*bmag)+jmag)
        if (rowd['mag']<98.0):
            nomad_cat.append(rowd)
    print "# Number of NOMAD entries found (after magnitude cuts) is: ",len(nomad_cat)
#
    t3=time.time()
    print("# TIMING (NOMAD query): {:.2f}".format((t2a-t2)))
    print("# TIMING (APASS query): {:.2f}".format((t3-t2a)))

###############################################################################
#   Time to read in the finalcut_cat (was red_cat) products.
#   Setup to read catalogs
#
    mtime=2.5*numpy.log10(exp_rec["exptime"])
    aval=aterm[band2i[exp_rec["band"]]]
    if (exp_rec["bunit"]=="e-"):
        bval=-2.5*numpy.log10(4.0)
    else:   
        bval=-2.5*numpy.log10(1.0)
#    print "RAG: bval=",bval
#
    mcorr=mtime-aval-25.0-bval

    nstar_found=0

###############################################################################
#   Columns to be retrieved from FITS tables.
    cols_retrieve=["ALPHAWIN_J2000","DELTAWIN_J2000","MAG_AUTO","MAGERR_AUTO","SPREAD_MODEL","FWHM_WORLD",
                   "A_IMAGE","B_IMAGE","FLUX_RADIUS","KRON_RADIUS","CLASS_STAR","FLAGS"]
    cols_datacheck={"ALPHAWIN_J2000":{"min":-0.1,"max":360.1},
                    "DELTAWIN_J2000":{"min":-90.0,"max":90.0},
                    "FWHM_WORLD":{"min":(0.1/3600.),"max":(10.0/3600.)},
                    "A_IMAGE":{"min":(0.1/pixsize),"max":(10.0/pixsize)},
                    "B_IMAGE":{"min":(0.1/pixsize),"max":(10.0/pixsize)},
                    "FLUX_RADIUS":{"min":(0.1/pixsize),"max":(5.0/pixsize)},
                    "KRON_RADIUS":{"min":(0.1/pixsize),"max":(5.0/pixsize)}
                    }
###############################################################################
#   Current columns avaliable in finalcut_cat objects.
#
#   ['NUMBER', 'FLAGS', 'DURATION_ANALYSIS', 'X_IMAGE', 'Y_IMAGE', 'XMIN_IMAGE', 'XMAX_IMAGE', 'YMIN_IMAGE', 
#    'YMAX_IMAGE', 'X2_IMAGE', 'Y2_IMAGE', 'XY_IMAGE', 'X2_WORLD', 'Y2_WORLD', 'XY_WORLD', 'ERRX2_IMAGE', 
#    'ERRY2_IMAGE', 'ERRXY_IMAGE', 'ERRX2_WORLD', 'ERRY2_WORLD', 'ERRXY_WORLD', 'X2WIN_IMAGE', 'Y2WIN_IMAGE', 
#    'XYWIN_IMAGE', 'ERRX2WIN_IMAGE', 'ERRY2WIN_IMAGE', 'ERRXYWIN_IMAGE', 'X2WIN_WORLD', 'Y2WIN_WORLD', 
#    'XYWIN_WORLD', 'ERRX2WIN_WORLD', 'ERRY2WIN_WORLD', 'ERRXYWIN_WORLD', 'ELLIPTICITY', 'FWHM_WORLD', 
#    'XWIN_IMAGE', 'YWIN_IMAGE', 'ERRAWIN_IMAGE', 'ERRBWIN_IMAGE', 'ERRTHETAWIN_IMAGE', 'ALPHA_J2000', 
#    'DELTA_J2000', 'ALPHAWIN_J2000', 'DELTAWIN_J2000', 'ERRAWIN_WORLD', 'ERRBWIN_WORLD', 'ERRTHETAWIN_J2000', 
#    'XPEAK_IMAGE', 'YPEAK_IMAGE', 'ALPHAPEAK_J2000', 'DELTAPEAK_J2000', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 
#    'THETA_J2000', 'ELONGATION', 'FLUX_RADIUS', 'MAG_ISO', 'MAGERR_ISO', 'ISOAREA_IMAGE', 'ISOAREAF_IMAGE', 
#    'ISOAREA_WORLD', 'ISOAREAF_WORLD', 'FLUX_APER', 'FLUXERR_APER', 'MAG_APER', 'MAGERR_APER', 'FLUX_AUTO', 
#    'FLUXERR_AUTO', 'MAG_AUTO', 'MAGERR_AUTO', 'KRON_RADIUS', 'BACKGROUND', 'THRESHOLD', 'MU_THRESHOLD', 
#    'FLUX_MAX', 'MU_MAX', 'XPSF_IMAGE', 'YPSF_IMAGE', 'ERRAPSF_IMAGE', 'ERRBPSF_IMAGE', 'ERRTHETAPSF_IMAGE', 
#    'ALPHAPSF_J2000', 'DELTAPSF_J2000', 'ERRAPSF_WORLD', 'ERRBPSF_WORLD', 'ERRTHETAPSF_J2000', 'FLUX_PSF', 
#    'FLUXERR_PSF', 'MAG_PSF', 'MAGERR_PSF', 'NITER_PSF', 'CHI2_PSF', 'SPREAD_MODEL', 'SPREADERR_MODEL', 
#    'CLASS_STAR', 'FWHMPSF_IMAGE', 'FWHMPSF_WORLD', 'MAG_POINTSOURCE']
#
#   Read each catalog and build list of dictionaries containing contents.
#   Provisionally check data while reading against limits set in cols_datacheck.
#

    if (args.verbose):
        print("# Obtaining catalog information from cat_finalcut catalog files")

    exp_cat=[]
    for iccd in range(1,63):
        if (iccd in ccd_info):
            if (os.path.isfile(ccd_info[iccd]['cat_fpath'])):
                cat_fits = fitsio.FITS(ccd_info[iccd]['cat_fpath'],'r')
                cat_hdu=2
#               OPTIONAL line to insert and dump out list of parameters/columns in catalog.
#               CATcol=cat_fits[cat_hdu].get_colnames()
#               print CATcol
                CAT=cat_fits[cat_hdu].read(columns=cols_retrieve)
                CATcol=cat_fits[cat_hdu].get_colnames()
                cdict=MkCatDict(CATcol,cols_retrieve)
#
                for row in CAT:
                    if (row[cdict["FLAGS"]] == 0):
                        check_entry=True
                        tmp_dict={}
                        for colnam in cols_retrieve:
                            tmp_dict[colnam]=row[cdict[colnam]]
                            if (colnam in cols_datacheck):
                                if (tmp_dict[colnam]<cols_datacheck[colnam]["min"]): check_entry=False
                                if (tmp_dict[colnam]>cols_datacheck[colnam]["max"]): check_entry=False
                        if (check_entry):
                            exp_cat.append(tmp_dict)
                cat_fits.close()
#   Re-order sorted on magnitude (allows for speed up when comparing with APASS/NOMAD below)
    new_expcat=sorted(exp_cat,key=lambda k: k['MAG_AUTO'])

    if (args.verbose):
        print("# Total number of cat_finalcut objects: {:d}".format(len(new_expcat)))
###################################################################
#
#   Form some simple sets for analysis
#

    mag=numpy.array([new_expcat[k]['MAG_AUTO']+mcorr for k, object in enumerate(new_expcat)])
    magerr=numpy.array([new_expcat[k]['MAGERR_AUTO']    for k, object in enumerate(new_expcat)])

    rstruct=[{'a_image':obj['A_IMAGE'],'b_image':obj['B_IMAGE'],'fwld':obj['FWHM_WORLD'],
              'frad':obj['FLUX_RADIUS'],'krad':obj['KRON_RADIUS']} 
              for obj in new_expcat if ((abs(obj['SPREAD_MODEL'])<0.0015)and(obj['MAGERR_AUTO']<0.1))]
    rrad1=numpy.array([2.*obj['a_image'] for obj in rstruct])
    rrad2=numpy.array([obj['a_image']+obj['b_image'] for obj in rstruct])
    fwld=numpy.array([obj['fwld'] for obj in rstruct])
    frad=numpy.array([obj['frad'] for obj in rstruct])
    krad=numpy.array([obj['krad'] for obj in rstruct])

###################################################################
#
#   Perform matching between cat_finalcut and APASS/NOMAD
#
    if (args.verbose):
        print("# Entering matching phase between cat_finalcut with APASS/NOMAD")

    cat_match={"apass":{"cat":apass_cat,"continue_ccor":True,"ccorr_mlimit":23.,"rate100":[],"rate100_lim":1,"mag100":[],"ncnt":0,"nfnd":0,"nfnd_sum":0},
               "nomad":{"cat":nomad_cat,"continue_ccor":True,"ccorr_mlimit":23.,"rate100":[],"rate100_lim":3,"mag100":[],"ncnt":0,"nfnd":0,"nfnd_sum":0}}

    nstar_found=0
    match_starset=[]
    for object in new_expcat:
#       Work on objects that are constrained to be stars (by SPREAD_MODEL).
        smval=abs(object["SPREAD_MODEL"])
        if ((smval < 0.0015)and(object["MAGERR_AUTO"]<0.1)):
            nstar_found=nstar_found+1
#
#           Check to see whether continued comparisons to catalogs APASS/NOMAD are warranted.
#
            continue_ccorr=False
            for cat in cat_match:
                if (cat_match[cat]["continue_ccor"]): 
                    continue_ccorr=True
            if (continue_ccorr):
                tmp_starset={}
                tmp_starset["ra"]=object["ALPHAWIN_J2000"]
                tmp_starset["dec"]=object["DELTAWIN_J2000"]
                tmp_starset["mag"]=object["MAG_AUTO"]+mcorr
                tmp_starset["magerr"]=object["MAGERR_AUTO"]
#
#               Cross correlation vs. NOMAD/APASS
#
                add_tmp_starset_record=False
                for cat in cat_match:
                    if (cat_match[cat]["continue_ccor"]):
                        match_rec=binsrch_crosscorr_cat(tmp_starset["ra"],tmp_starset["dec"],cat_match[cat]["cat"],1.0)
                        cat_match[cat]["ncnt"]=cat_match[cat]["ncnt"]+1
                        if (len(match_rec)>0):
                            cat_match[cat]["nfnd"]=cat_match[cat]["nfnd"]+1
                            tmp_starset[cat]=match_rec[0]
                            tmp_starset[cat]["magcorr"]=tmp_starset[cat]["mag"]+cat_mag_corr[cat][exp_rec["band"]]
                            tmp_starset[cat]["magdiff"]=tmp_starset["mag"]-tmp_starset[cat]["mag"]-cat_mag_corr[cat][exp_rec["band"]]
                            add_tmp_starset_record=True
#
#                       Look for the point of diminishing returns (very few new matches found)
#
                        if (cat_match[cat]["ncnt"]%100 == 0):
                            cat_match[cat]["rate100"].append(cat_match[cat]["nfnd"])  
                            cat_match[cat]["nfnd_sum"]=cat_match[cat]["nfnd_sum"]+cat_match[cat]["nfnd"]
                            cat_match[cat]["nfnd"]=0
                            if ((len(cat_match[cat]["rate100"])>10)and(cat_match[cat]["nfnd_sum"]>500)):
                                if ((cat_match[cat]["rate100"][-3]<cat_match[cat]["rate100_lim"])and
                                    (cat_match[cat]["rate100"][-2]<cat_match[cat]["rate100_lim"])and
                                    (cat_match[cat]["rate100"][-1]<cat_match[cat]["rate100_lim"])):
                                    cat_match[cat]["continue_corr"]=False
#               If add_tmp_starset_record is True then add the record to the matched star set.
                if (add_tmp_starset_record):
                    match_starset.append(tmp_starset)

    t4=time.time()
    exp_rec["numobj"]=len(mag)
    print "# "
    print("# Number of objects from exposure catalogs: {:d}".format(exp_rec["numobj"]))
    print("# Number of stellar objects from exposure catalogs: {:d}".format(nstar_found))

    sat_check={}
    num_matches={}
    for cat in ['apass','nomad']:
        sat_check[cat]=[]
        num_matches[cat]=0
    for record in match_starset:
        for cat in ['apass','nomad']:
            if (cat in record):
                sat_check[cat].append((record["mag"],record[cat]["magcorr"],record[cat]["magdiff"]))
                num_matches[cat]=num_matches[cat]+1
    for cat in ['apass','nomad']:
        print("# Number of {:s} matches: {:d}".format(cat,num_matches[cat]))
    print("# TIMING (APASS/NOMAD crosscorr): {:.2f}".format((t4-t3)))

###################################################################
#
#   Now that matching is complete...
#   Check for significant #'s of saturated stars.
#
    sat_limit=5.0
    sat_sig_thresh=2.0
    print("###################################################################")
    print("# Evaluating saturation limit")
    for cat in ['apass','nomad']:
        print("# ")
        adjust_satlimit=False
        sat_check[cat+'_sorted']=sorted(sat_check[cat],key=lambda x: x[1])
        if (len(sat_check[cat+'_sorted'])>0):
            magdes,magcorr,magdiff=zip(*sat_check[cat+'_sorted'])
        else:
            magdes=[]
            magcorr=[]
            magdiff=[]
        a_magcorr=numpy.array(magcorr)
        a_magdiff=numpy.array(magdiff)
        if (len(magcorr)>100):
            sample_end=len(magcorr)/5
            if (sample_end > 100):
                sample_end=100
            b_med_magdiff=numpy.median(a_magdiff[:sample_end])
            f_med_magdiff=numpy.median(a_magdiff[-sample_end:])
            b_std_magdiff=numpy.std(a_magdiff[:sample_end])
            f_std_magdiff=numpy.std(a_magdiff[-sample_end:])

            f_diff_cut=f_med_magdiff+(sat_sig_thresh*f_std_magdiff)
            print("# For {:s} median mag_diff for FAINT  stars {:7.3f} +/- {:7.3f}".format(cat,f_med_magdiff,f_std_magdiff))
            print("# For {:s} median mag_diff for BRIGHT stars {:7.3f} +/- {:7.3f}".format(cat,b_med_magdiff,b_std_magdiff))
            if (b_med_magdiff>f_diff_cut):
                print("# Magdiff for faint stars in {:s} more than {:.1f}-sigma offset from bright stars.".format(cat,sat_sig_thresh))
                adjust_satlimit=True
            nstar=0
            dnstar=10
            while((nstar<(len(magcorr)-dnstar))and(adjust_satlimit)):
                num_sat=len(numpy.where(a_magdiff[nstar:(nstar+dnstar)] > f_diff_cut)[0])
                if (num_sat >= dnstar/2):
                    b_avg=numpy.average(a_magcorr[nstar:(nstar+dnstar)])
                    if (b_avg > sat_limit):
                        sat_limit=b_avg
                else:
                    adjust_satlimit=False
                nstar=nstar+dnstar 
            print("# After examining {:s}, saturation limit now: {:.3f} ".format(cat,sat_limit))
#
#   Find out how many stars are left.
#
    print("# ")
    num_matches={}
    for cat in ['apass','nomad']:
        if (len(sat_check[cat])>0):
            magdes,magcorr,magdiff=zip(*sat_check[cat])
        else:
            magdes=[]
            magcorr=[]
            magdiff=[]
        a_magcorr=numpy.array(magcorr)
        num_matches[cat]=len(numpy.where(a_magcorr > sat_limit)[0])
    for cat in ['apass','nomad']:
        print("# Number of {:s} matches remaining: {:d}".format(cat,num_matches[cat]))

########################################################################################
#   Proceed to assessing the atmospheric opacity.
#
    cat_use={"apass":['g','r'],
             "nomad":['u','g','r','i','z','Y','VR']}
#
    transparency_calc=False
#   if (exp_rec["band"] in ["g","r","i"]):
    for cat in ['apass','nomad']:
        if (num_matches[cat]<100):
#           For a limited number of stars do not make comparison... output will be sentinel value.
            exp_rec[cat+"_magdiff"]=99.0
            exp_rec[cat+"_kmagdiff"]=99.0
            exp_rec[cat+"_num"]=num_matches[cat]
        else:
#
#           Now the workhorse case (for APASS)
#           Workhorse case... find offset
#
            if (len(sat_check[cat])>0):
                mag_des,mag_cat,mag_diff=zip(*sat_check[cat])
            else:
                mag_des=[]
                mag_cat=[]
                mag_diff=[]
            amag_cat=numpy.array(mag_cat)
            amag_diff=numpy.array(mag_diff)
#           numpy.where(numpy.logical_and((numpy.logical_and((r2_acolor>=cntcut),(i2_acolor>=cntcut))),(numpy.logical_and((z2_acolor>=cntcut),(g1_acolor>75.)))
            cut_val=numpy.where(amag_cat > sat_limit)
            med_magdiff=numpy.median(amag_diff[cut_val])
            std_magdiff=numpy.std(amag_diff[cut_val])
            print("# ")
            print("# (PRELIM {:s}) Median magnitude offset: {:8.3f} {:8.3f} ".format(cat,med_magdiff,std_magdiff))
#
#           Perform a rejection iteration.
#
            mag_fit=[]
            mag_reject=[]
            for record in sat_check[cat]:
                junkval=record[0]-record[1]
                reject_val=False
                if (abs((junkval-med_magdiff)/std_magdiff)>3.0): reject_val=True
                if (record[1] < sat_limit): reject_val=True
                if (reject_val):
                    mag_reject.append(record)
                else:
                    mag_fit.append(record)
            if (len(mag_fit)>0):
                mag_des,mag_cat,mag_diff=zip(*mag_fit)
            else:
                mag_des=[]
                mag_cat=[]
                mag_diff=[]
            if (len(mag_reject)>0):
                mag_des_reject,mag_cat_reject,mag_diff_reject=zip(*mag_reject)
            else:
                mag_des_reject=[]
                mag_cat_reject=[]
                mag_diff_reject=[]
#
#           Go on and calculate magnitude(s) of opacity
#
            amag_diff=numpy.array(mag_diff)
            mindes=min(mag_des)
            maxdes=max(mag_des)
            med_magdiff=numpy.median(amag_diff)
            std_magdiff=numpy.std(amag_diff)
            exp_rec[cat+"_magdiff"]=med_magdiff
            exp_rec[cat+"_kmagdiff"]=med_magdiff-cat_mag_corr[cat][exp_rec['band']]+cat_kmag_corr[cat][exp_rec['band']]-(kterm[band2i[exp_rec['band']]]*exp_rec["airmass"])
            exp_rec[cat+"_num"]=num_matches[cat]
#
#           LOGIC HERE DECIDES RULES FOR WHICH MEASURE FROM WHICH CATALOG TO USE
#              -- order they go through is important (APASS, then NOMAD cleans up)
#              -- once a choice is made it simple reports the value (with NOUSE)
#
            if(not(transparency_calc)):
                if (exp_rec["band"] in cat_use[cat]):
                    transparency_calc=True
                    exp_rec["magdiff"]=med_magdiff
                    exp_rec["cloud_cat"]=cat.upper()
                    print("# (USING {:s}) Median magnitude offset: {:8.3f} {:8.3f} ".format(cat.upper(),med_magdiff,std_magdiff))
                else:
                    print("# (NOUSE {:s}) Median magnitude offset: {:8.3f} {:8.3f} ".format(cat.upper(),med_magdiff,std_magdiff))
            else:
                print("# (NOUSE {:s}) Median magnitude offset: {:8.3f} {:8.3f} ".format(cat.upper(),med_magdiff,std_magdiff))
#
#       Finished analysis
#
        if (args.qaplot):
#
#           If --qaplot option is active then generate QA plots for atmospheric exinctiont.
#
            PlotFileName='{:s}_{:s}.png'.format(args.froot,cat)
            PlotData={}
            PlotData['nmatch']=num_matches[cat]
            if (num_matches[cat]<100):
#               If QA plots have been requested then output an empty QA plot that gives 
#               information about the failure.
                if (len(sat_check[cat])>0):
                    PlotData['mag_des'],PlotData['mag_cat'],PlotData['mag_diff']=zip(*sat_check[cat])
                else:
                    PlotData['mag_des']=[]
                    PlotData['mag_cat']=[]
                    PlotData['mag_diff']=[]
            else:
                PlotData['mag_des']=mag_des
                PlotData['mag_cat']=mag_cat
                PlotData['mag_diff']=mag_diff
                PlotData['mag_des_reject']=mag_des_reject
                PlotData['mag_cat_reject']=mag_cat_reject
                PlotData['mag_diff_reject']=mag_diff_reject
                PlotData['mindes']=mindes
                PlotData['maxdes']=maxdes
                PlotData['med_magdiff']=med_magdiff

            QAplt_cloud(PlotFileName,PlotData,cat,exp_rec['band'],astrom_good,args.verbose)
#
#   Now the case where all attempts to make a transparency calcultion have failed... give a sentinel value
#
    if (not(transparency_calc)):
        exp_rec["magdiff"]=99.
        exp_rec["cloud_cat"]="Failed"
    t5=time.time()
    print("# TIMING (Transparency calculation): {:.2f}".format((t5-t4)))
    print(" ")
##############################################################################
#   Load all the various FWHM measuments into the dictionary
#
    if (nstar_found > 2):                   
        exp_rec["nstar_found"]=nstar_found
        exp_rec["fwhm_rrad1"]=pixsize*rrad1.mean()
        exp_rec["fwhm_rrad2"]=pixsize*rrad2.mean()
        exp_rec["fwhm_world"]=3600.*fwld.mean()
        exp_rec["fwhm_krad"]=2.0*pixsize*krad.mean()
        exp_rec["fwhm_frad"]=2.0*pixsize*frad.mean()
    else:
        exp_rec["nstar_found"]=nstar_found
        exp_rec["fwhm_rrad1"]=-1.
        exp_rec["fwhm_rrad2"]=-1.
        exp_rec["fwhm_world"]=-1.
        exp_rec["fwhm_krad"]=-1.
        exp_rec["fwhm_frad"]=-1.
    t5a=time.time()

#########################################################################################
#   In the following section calculate some QA on the magnitude depth (needs work to make it 
#   robust (or remotely useful).
#
    magerr_hist=numpy.zeros((len(mbin)),dtype=numpy.float32)
    magnum_hist=numpy.zeros((len(mbin)),dtype=numpy.float32)

    if (len(mag)<10):
#
#       If insufficient data is present to even begin to try then simply report failure,
#       put out a super-unintersting plot, and move along.  
#
        print("#")
        print("# Insufficient data to estimate magnitude depth")
        mag_thresh=99.9
        if (args.qaplot):   
#           If QA plots were requested then write a dummy plot that indicates failure
            PlotFileName='{:s}_magplot.png'.format(args.froot)
            PlotData={}
            PlotData['nobj']=len(mag)
            QAplt_maghist(PlotFileName,PlotData,exp_rec['band'],astrom_good,args.verbose)

    else:
#
#       Sort the resuling mag and magerr based on mag... (sort indexes then apply)
#
        indmag=numpy.argsort(mag)
        emag_sort=numpy.take(mag,indmag)
        emagerr_sort=numpy.take(magerr,indmag)
#
#       Form histogram of objects
#
        icnt=0
        m0=mbin[0]-magbin_step
        for ibin, m1 in enumerate(mbin):
            i0=icnt
            while((emag_sort[icnt]<mbin[ibin])and(icnt<(len(emag_sort)-1))):
                icnt=icnt+1

            ncnt=icnt-1-i0
            if (ncnt > 1):
                x=numpy.median(emagerr_sort[i0:icnt-1])
            elif (ncnt == 1):
                x=emagerr_sort[i0]
            else:
                x=-1.
        
            magerr_hist[ibin]=x
            magnum_hist[ibin]=ncnt
#
#       Calculate (simple interpolation) the magnitude where the median magerr became greater than the threshhold
#
        m0=mbin[0]-magbin_step
        ibin_thresh=-1
        ifirst_bin=-1
        ilast_bin=-1
        ilast_bin2=-1
#       Find first and last good data point (bin had more than magnum_thresh)
        for ibin, m1 in enumerate(mbin):
            if (magnum_hist[ibin] > magnum_thresh):
                if (ifirst_bin < 0):
                    ifirst_bin=ibin
                ilast_bin2=ilast_bin
                ilast_bin=ibin
                if ((magerr_hist[ibin] > magerr_thresh)and(ibin_thresh < 0)):
                    ibin_thresh=ibin
#       Check that the magerr threshold has been exceeded
        if (ibin_thresh > 0):
            mag0=mbin[ibin_thresh]-(1.5*magbin_step)
            mag1=mbin[ibin_thresh]-(0.5*magbin_step)
            me0=magerr_hist[ibin_thresh-1]
            me1=magerr_hist[ibin_thresh]
            mag_thresh=mag0+((magerr_thresh-me0)*(mag1-mag0)/(me1-me0))
            mag_method="Interpolated"
#       If not exceeded but there are at least two data points then extrapolate using last two points
        elif ((ilast_bin > 0)and(ilast_bin2 > 0)):
#           If slope is not increasing then there is going to be a problem
            if (magerr_hist[ilast_bin] > magerr_hist[ilast_bin2]):
                mag0=mbin[ilast_bin2]-(0.5*magbin_step)
                mag1=mbin[ilast_bin]-(0.5*magbin_step)
                me0=magerr_hist[ilast_bin2]
                me1=magerr_hist[ilast_bin]
                mag_thresh=mag0+((magerr_thresh-me0)*(mag1-mag0)/(me1-me0))
                mag_method="Extrapolated"
#           Try taking one step back in the histogram and recheck
            elif ((magnum_hist[ilast_bin2-1] > magnum_thresh)and(magerr_hist[ilast_bin2] > magerr_hist[ilast_bin2-1])):
                mag0=mbin[ilast_bin2-1]-(0.5*magbin_step)
                mag1=mbin[ilast_bin2]-(0.5*magbin_step)
                me0=magerr_hist[ilast_bin2-1]
                me1=magerr_hist[ilast_bin2]
                mag_thresh=mag0+((magerr_thresh-me0)*(mag1-mag0)/(me1-me0))
                mag_method="Extrapolated_v2"
#           Failure is indeed an option (so return lowest bin)
            else:
                mag_thresh=mbin[0]
                mag_method="Extrapolate Failed"
#           print("{:d} {:d} {:8.3f} {:8.3f} {:8.3f} {:8.3f} ".format(ilast_bin,ilast_bin2,m0,m1,me0,me1))
#       Else return lowest/first magnitude bin
        else:
            mag_thresh=mbin[0]
            mag_method="Failed"
#
#       Optional QA file giving histograrm information.
#
#       if (args.Magout):
#           fmag = open("%s.%s.mag"% (args.froot,exp_rec["expnum"]), 'w')
#           fmag.write("####################################\n")
#           fmag.write("# {:12d} mag_thresh={:8.3f} \n".format(exp_rec["expnum"],mag_thresh))
#           fmag.write("# Type={:s} \n".format(mag_method))
#           for ibin, m1 in enumerate(mbin):
#               fmag.write("  {:8.3f} {:8.0f} {:8.3f}  \n".format(mbin[ibin],magnum_hist[ibin],magerr_hist[ibin]))
#           fmag.close()
    
        if (args.qaplot):   
            PlotFileName='{:s}_magplot.png'.format(args.froot)
            PlotData={}
            PlotData['nobj']=len(mag)
            PlotData['bins']=mbin-(0.5*magbin_step)
            PlotData['mag_hist']=magnum_hist
            PlotData['magerr_hist']=magerr_hist
            QAplt_maghist(PlotFileName,PlotData,exp_rec['band'],astrom_good,args.verbose)

    exp_rec["mag_thresh"]=mag_thresh
    t6=time.time()
    print("# TIMING (Object number counts): {:.2f}".format((t6-t5a)))

###############################################################################
#   Now the calculations for the Teff (and of course the individual components)
#
#   Calculate F_eff
#   Note code is now updated to use the psfex_fwhm (with fwhm_world used as a fallback)
#
    use_fwhm=-1.0
#   Uncomment the following line if you want to force runs to use FWHM_WORLD (i.e. for tests)
#    exp_rec['psfex_fwhm']=-1.0
    if (exp_rec['psfex_fwhm']>0.0):
#       RAG: Addative, empirical, offset now applied here prior to calculating T_eff
        use_fwhm=exp_rec['psfex_fwhm']+fwhm_DMtoQC_offset_psfex
    else:
#       RAG: Multiplicative, empirical, offset (for FWHM_WORLD) now applied here prior to calculating T_eff
        if (exp_rec['fwhm_world']>0.0):
            print("# WARNING:  No FWHM from PSFex... attempting to use FWHM_WORLD")
            use_fwhm=exp_rec['fwhm_world']/fwhm_DMtoQC_offset_world
#
#   OK so I lied above... calculate F_eff (NOW!)
#
    if (use_fwhm > 0.0):
        exp_rec["teff_f"]=(seeing_fid[exp_rec["band"]]*seeing_fid[exp_rec["band"]]/(use_fwhm*use_fwhm))
    else:
        print("# WARNING:  No FWHM measure available. F_EFF set to -1.0")
        exp_rec["teff_f"]=-1.0
#
#   Calculate B_eff
#
    if (exp_rec["skyb_avg"]>0.0):
        exp_rec["teff_b"]=sbrite_good[exp_rec["band"]]/exp_rec["skyb_avg"]
    else:
        print("# WARNING:  No SKY BRIGHTNESS measure available. B_EFF set to -1.0")
        exp_rec["teff_b"]=-1.0
#
#   Calculate C_eff
#
    if ((exp_rec["magdiff"]>-95.)and(exp_rec["magdiff"]<95.0)):
        if (exp_rec["magdiff"]<0.2):
            exp_rec["teff_c"]=1.0
        else:
            exp_rec["teff_c"]=math.pow(10.0,(-2.0*(exp_rec["magdiff"]-0.2)/2.5))
    else:
        print("# WARNING:  No CLOUD measure available. C_EFF set to -1.0")
        exp_rec["teff_c"]=-1.0
#
#   Calculate T_eff
#
    value_teff=1.0
    if (exp_rec["teff_f"]>=0):
        value_teff=value_teff*exp_rec["teff_f"]
    if (exp_rec["teff_b"]>=0):
        value_teff=value_teff*exp_rec["teff_b"]
    if (exp_rec["teff_c"]>=0):
        value_teff=value_teff*exp_rec["teff_c"]
    if ((exp_rec["teff_f"]<0)or(exp_rec["teff_b"]<0)):
        exp_rec["teff"]=-1.
    else:
        exp_rec["teff"]=value_teff
    t7=time.time()

#############################################################################
#   Summary of all timing information
#   Summary of all FWHM measurements

    print("# ")
    print("# ")
    print("#         TIMING Summary  ")
    print("#------------------------------")
    print("#  general query(s): {:7.2f} ".format((t1-t0)))
    print("#             APASS: {:7.2f} ".format((t3-t2a)))
    print("#             NOMAD: {:7.2f} ".format((t2a-t2)))
    print("# Cross-correlation: {:7.2f} ".format((t4-t3)))
    print("#    Cat comparison: {:7.2f} ".format((t5-t4)))
    print("#      Mag Analysis: {:7.2f} ".format((t6-t5)))
    print("# ")
    print("#               all: {:.2f} ".format((t7-t0)))
    print("# ")

    print("# ")
    print("# ")
    print("#              FWHM Summary   ")
    print("#------------------------------------")
    print("#          NSTAR_FOUND = {:d} ".format(exp_rec['nstar_found']))
    print("#        FWHM_IMG(OLD) = {:7.3f} ".format(exp_rec['fwhm_img']))
    print("#           FWHM_WORLD = {:7.3f} ".format(exp_rec['fwhm_world']))
    print("#          2 * A_IMAGE = {:7.3f} ".format(exp_rec['fwhm_rrad1']))
    print("#    A_IMAGE + B_IMAGE = {:7.3f} ".format(exp_rec['fwhm_rrad2']))
    print("#          Kron Radius = {:7.3f} ".format(exp_rec['fwhm_krad']))
    print("#               F(RAD) = {:7.3f} ".format(exp_rec['fwhm_frad']))
    print("# ")
    print("#     PSFex(FWHM_MEAN) = {:7.3f} ".format(exp_rec['psfex_fwhm']))
    print("# ")

###############################################################################
#   Everything is now ready to make an assessment and to output it to the 
#   proper locations.
#
    new_decide="none"
    if ((exp_rec["teff"]<0.)or(exp_rec["teff_c"]<0.0)):
        new_decide="unkn"
    elif (exp_rec["teff"]>teff_lim[exp_rec["band"]]):
        new_decide="good"
        if (exp_rec["fwhm_world"]>seeing_lim[exp_rec["band"]]):
            new_decide="badF"
    else:
        new_decide="badT"
#
    dm_process="True"
    if (new_decide == "good"):
        dm_accept="True"
    elif (new_decide == "unkn"):
        dm_accept="Unknown"
    else:
        dm_accept="False"
#
#       Write CSV if command line option present 
#
    if (args.csv):
        out_row=[]
        out_row.append(exp_rec["expnum"])
        out_row.append(dm_process)
        out_row.append(dm_accept)
        out_row.append("no comment")
        writer.writerow(out_row)  
#
#   Write DB_file if command line option present 
#
#    if((args.updateDB)or(args.DB_file)):
    prog_name='Unknown'
    if (exp_rec["sn_field"]):
        prog_name='supernova'
    if (exp_rec["survey_field"]):
        prog_name='survey'
#
#       First the columns that should always be present
#
    db_cname="INSERT INTO {dbtab:s}(EXPOSURENAME,EXPNUM,PFW_ATTEMPT_ID,CAMSYM,PROCESSED,ACCEPTED,PROGRAM,ANALYST,LASTCHANGED_TIME,T_EFF,F_EFF,B_EFF,C_EFF,FWHM_ASEC,ELLIPTICITY,SKYBRIGHTNESS".format(dbtab=db_table)
    db_value=""" VALUES ('{expname:s}', {expnum:d}, {pfw_att_id:s}, '{camsym:s}', '{flag_proc:s}', '{flag_accept:s}', '{prog:s}', '{analyst:s}', sysdate, {teff:.3f}, {feff:.3f}, {beff:.3f}, {ceff:.3f}, {fwhm:.3f}, {ellip:.3f}, {skyb:.2f}""".format(
        expname=exp_rec["exp_fname"],
        expnum=exp_rec["expnum"], 
        pfw_att_id=args.attid, 
        camsym=args.CamSym, 
        flag_proc=dm_process,
        flag_accept=dm_accept,
        prog=prog_name,
        analyst=analyst,
        teff=exp_rec["teff"],
        feff=exp_rec["teff_f"],
        beff=exp_rec["teff_b"],
        ceff=exp_rec["teff_c"],
        fwhm=exp_rec["psfex_fwhm"],
        ellip=exp_rec["ellip_avg"],
        skyb=exp_rec["skyb_avg"])

#
#       Now the columns that may or may not be present
#
    if ((exp_rec["apass_magdiff"]>-95.)and(exp_rec["apass_magdiff"]<95.)):
        db_cname=db_cname+','+'CLOUD_APASS'
        db_value=db_value+','+"%.3f" % exp_rec["apass_magdiff"]

    if ((exp_rec["nomad_magdiff"]>-95.)and(exp_rec["nomad_magdiff"]<95.)):
        db_cname=db_cname+','+'CLOUD_NOMAD'
        db_value=db_value+','+"%.3f" % exp_rec["nomad_magdiff"]

    db_cname=db_cname+','+'CLOUD_CATALOG'
    db_value=db_value+','+"'%s'" % exp_rec["cloud_cat"]

    if ((exp_rec["apass_kmagdiff"]>-95.)and(exp_rec["apass_kmagdiff"]<95.)):
        db_cname=db_cname+','+'KLOUD_APASS'
        db_value=db_value+','+"%.3f" % exp_rec["apass_kmagdiff"]

    if ((exp_rec["nomad_kmagdiff"]>-95.)and(exp_rec["nomad_kmagdiff"]<95.)):
        db_cname=db_cname+','+'KLOUD_NOMAD'
        db_value=db_value+','+"%.3f" % exp_rec["nomad_kmagdiff"]

    if (exp_rec["apass_num"]>-1):
        db_cname=db_cname+','+'N_APASS'
        db_value=db_value+','+"%d" % exp_rec["apass_num"]

    if (exp_rec["nomad_num"]>-1):
        db_cname=db_cname+','+'N_NOMAD'
        db_value=db_value+','+"%d" % exp_rec["nomad_num"]

#       Finish off the strings that contain columns and value and combine to form INSERT command
#       Then write and/or perform the insert
#
    db_cname=db_cname+")"
    db_value=db_value+")"
    insert_command=db_cname+db_value

    if (args.verbose):
        print(" {:s}".format(insert_command))
    if(args.updateDB):
        print("# Executing insert command")
        cur.execute(insert_command)
    if (args.DB_file):
        fdbout.write("{:s};\n".format(insert_command))

#
#   Write the summary of the assessment to STDOU
#
    ftxt.write("#                                                                                                         Astrometry \n")     
    ftxt.write("# Exposure                                     FWHM                             #        AIR    EXP  SCMP  HighSN     Depth       APASS           NOMAD      FWHM            \n")
    ftxt.write("#  Num   STATE    t_eff F_eff  B_eff   C_eff   PSFex  Ellip    Sky_B   N_src   CCD BAND  MASS   TIME Flag sig1  sig2     Assess     dmag    #       dmag    #     Wld    Object   \n")
    ftxt.write("#                                               [\"]           [DN/s]                             [s]      [\"]    [\"]     [mag]    [mag]            [mag]         \n")
#
    ftxt.write(" {expnum:9d} {decide:4s}  {teff:6.2f} {feff:6.2f} {beff:6.2f} {ceff:6.2f}   {fwhm:5.2f} {ellip:6.3f} {skyb:8.2f} {numobj:7d}   {numccd:3d} {band:1s} {iband:1d} {airmass:6.3f} {exptime:6.1f} {dummy:2d} {asig1:6.3f} {asig2:6.3f}  {mdiff:8.3f} {amdiff:8.3f} {anum:6d} {nmdiff:8.3f} {nnum:6d} {fwhm_wld:6.3f} {object:s} \n".format(
        expnum=exp_rec["expnum"],
        decide=new_decide,
        teff=exp_rec["teff"],
        feff=exp_rec["teff_f"],
        beff=exp_rec["teff_b"],
        ceff=exp_rec["teff_c"],
        fwhm=exp_rec["psfex_fwhm"],
        ellip=exp_rec["ellip_avg"],
        skyb=exp_rec["skyb_avg"],
        numobj=exp_rec["numobj"],
        numccd=exp_rec["numccd"],
        band=exp_rec["band"],
        iband=band2i[exp_rec["band"]],
        airmass=exp_rec["airmass"],
        exptime=exp_rec["exptime"],
        dummy=99,
        asig1=exp_rec["astrom_sig1"],
        asig2=exp_rec["astrom_sig2"],
        mdiff=exp_rec["magdiff"],
        amdiff=exp_rec["apass_magdiff"],
        anum=exp_rec["apass_num"],
        nmdiff=exp_rec["nomad_magdiff"],
        nnum=exp_rec["nomad_num"],
        fwhm_wld=exp_rec["fwhm_world"],
        object=exp_rec["object"]))

    if(args.updateDB):
        dbh.commit()
        print "DB update complete and committed"
    if (args.csv):
        fcsv.close()
    if (args.DB_file):
        fdbout.close()

