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
                is RA-ordered so that a binary search can be used to find the small range of entries
        that need to be checked
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
    parser.add_argument('-R', '--ReqNum',   action='store', type=str, default=None, help='Processing request number.')
    parser.add_argument('-U', '--UnitName', action='store', type=str, default=None, help='Unit Name.')
    parser.add_argument('-A', '--AttNum',   action='store', type=str, default=None, help='Attempt Number.')
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
    if ((args.ReqNum is None)or(args.UnitName is None)or(args.AttNum is None)):
        parser.print_help()
        exit(1)
    ReqNum=args.ReqNum
    UnitName=args.UnitName
    AttNum=args.AttNum

    print("####################################################")
    print("# Note this execution requests:")
    print("# ")
    print("# Evaluation of REQNUM={:s}  ".format(ReqNum))
    print("#             UNITNAME={:s}  ".format(UnitName))
    print("#               ATTNUM={:s}  ".format(AttNum))
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

    try:
        desdmfile = os.environ["des_services"]
    except KeyError:
        desdmfile = None
    dbh = despydb.desdbi.DesDbi(desdmfile,args.section)
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
    fwhm_DMtoQC_offset=1.10
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
        seeing_fid[band]=fwhm_DMtoQC_offset*0.9*kolmogorov[band]
#    print kolmogorov
#    print seeing_lim
#    print seeing_fid
#
#   Surface brightness limits from Eric Nielson which were derived "...from a few 
#   exposures from a photometric night in SV with little moon (20121215)"
#
    sbrite_good={"u":2.0,"g":1.05,"r":2.66,"i":7.87,"z":16.51,"Y":14.56,"VR":3.71}
    sbrite_lim={"u":8.0,"g":4.0,"r":9.58,"i":21.9,"z":50.2,"Y":27.6,"VR":13.58}
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
    apass_mag_corr={"u":0.0,"g":0.205,"r":0.128,"i":0.112,"z":0.0,"Y":0.0,"VR":0.0}
    apass_kmag_corr={"u":0.0,"g":0.000,"r":0.000,"i":0.000,"z":0.0,"Y":0.0,"VR":0.0}
    nomad_mag_corr={"u":0.0,"g":0.341,"r":0.235,"i":1.398,"z":1.201,"Y":1.083,"VR":0.0}
    nomad_kmag_corr={"u":0.0,"g":0.111,"r":0.109,"i":1.289,"z":1.139,"Y":1.022,"VR":0.0}
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
#       If no list of catalogs was provided (--list) then go figure out whether they exist.
#
        queryitems = ["wgb.filename","wgb.filetype","wgb.exec_name","task.status", "fai.path", "oa.root","c.ccdnum"]
        coldict={}
        for index, item in enumerate(queryitems):
            coldict[item]=index
        querylist = ",".join(queryitems)
        query = """select %s from %swgb wgb, %stask task, %sfile_archive_info fai, %sops_archive oa, %scatalog c where wgb.exec_task_id=task.id and wgb.reqnum=%s and wgb.unitname='%s' and wgb.attnum=%s and wgb.filetype='cat_finalcut' and fai.filename=wgb.filename and c.filename=wgb.filename and oa.name=fai.archive_name order by task.start_time """ % (querylist,db_Schema,db_Schema,db_Schema,db_Schema,db_Schema,ReqNum,UnitName,AttNum)

        print "# Executing initial query for catalogs "
        if args.verbose:
            print("# {:s}".format(query))
        cur.arraysize = 1000 # get 1000 at a time when fetching
        cur.execute(query)
#        desc = [d[0].lower() for d in cur.description]
#        print desc
 
        for item in cur:
            if (item[coldict["task.status"]]==0):
                cat_fname.append(item[coldict["wgb.filename"]])
                cat_fpath.append(os.path.join(item[coldict["oa.root"]],item[coldict["fai.path"]],item[coldict["wgb.filename"]]))
    else:
#
#       Open list file and proceed.
#
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
    print "# Total number of catalogs/CCDs identified: ",len(cat_fpath)
##############################################################################
#   Now that a starting point has been established...
#   Obtain associated DB information for each catalog/CCD.
# 
    queryitems = ["i.filename","i.expnum","i.band","i.ccdnum","i.airmass","i.exptime","i.fwhm","i.elliptic","i.skybrite","i.skysigma","i.scampflg","i.gaina","i.gainb","i.ra_cent","i.dec_cent","i.rac1","i.rac2","i.rac3","i.rac4","i.decc1","i.decc2","i.decc3","i.decc4"]
    coldict={}
    for index, item in enumerate(queryitems):
        coldict[item]=index
    querylist = ",".join(queryitems)

    ccd_info={}
    num_ccd_verbose=0
    max_ccd_verbose=1
    print "# Executing query(s) to obtain image level metadata and QA values "
    for index,tmp_fname in enumerate(cat_fname):
        query = """select %s from %simage i, %swdf wdf where i.filename=wdf.parent_name and wdf.child_name='%s' """ % (querylist,db_Schema,db_Schema,cat_fname[index])
        if ((args.verbose)and(num_ccd_verbose < max_ccd_verbose)):
            print("# {:s}".format(query))
            num_ccd_verbose=num_ccd_verbose+1
        cur.arraysize = 1000 # get 1000 at a time when fetching
        cur.execute(query)

        for item in cur:
            tmp_dict={}
            tmp_dict["cat_fname"]=cat_fname[index]
            tmp_dict["cat_fpath"]=cat_fpath[index]
            tmp_dict["img_fname"]=item[coldict["i.filename"]]
            tmp_dict["expnum"]=item[coldict["i.expnum"]]
            if (item[coldict["i.band"]] is None):
                tmp_dict["band"]='Unknown'
                tmp_dict["iband"]=-1
            else:
                tmp_dict["band"]=item[coldict["i.band"]]
                tmp_dict["iband"]=band2i[item[coldict["i.band"]]]
            if (item[coldict["i.gaina"]] is None):
                tmp_dict["gaina"]=-1.
                print "# Warning: Null found for GAINA"
            else:
                tmp_dict["gaina"]=float(item[coldict["i.gaina"]])
            if (item[coldict["i.gainb"]] is None):
                tmp_dict["gainb"]=-1.
                print "# Warning: Null found for GAINB"
            else:
                tmp_dict["gainb"]=float(item[coldict["i.gainb"]])
            tmp_dict["ccdnum"]=int(item[coldict["i.ccdnum"]])
            tmp_dict["airmass"]=float(item[coldict["i.airmass"]])
            tmp_dict["exptime"]=float(item[coldict["i.exptime"]])
            tmp_dict["fwhm_old"]=pixsize*float(item[coldict["i.fwhm"]])
            tmp_dict["elli_old"]=float(item[coldict["i.elliptic"]])
#
            if ((item[coldict["i.rac1"]] is None)or(item[coldict["i.rac2"]] is None)or(item[coldict["i.rac3"]] is None)or(item[coldict["i.rac4"]] is None)):
                print "# Warning: Null found among RAC1-4"
                tmp_dict["img_ra"]=[]
            else:
                tmp_dict["img_ra"]=[float(item[coldict["i.rac1"]]),float(item[coldict["i.rac2"]]),float(item[coldict["i.rac3"]]),float(item[coldict["i.rac4"]])]
#
            if ((item[coldict["i.decc1"]] is None)or(item[coldict["i.decc2"]] is None)or(item[coldict["i.decc3"]] is None)or(item[coldict["i.decc4"]] is None)):
                print "# Warning: Null found among DECC1-4"
                tmp_dict["img_dec"]=[]
            else:
                tmp_dict["img_dec"]=[float(item[coldict["i.decc1"]]),float(item[coldict["i.decc2"]]),float(item[coldict["i.decc3"]]),float(item[coldict["i.decc4"]])]
#
#           sbrite_good and sbrite_lim are defined in cts/sec
#           Therefore:
#              If exposure time is present then use to normalize brightness to a per second quantity
#              If gains are present then use to express brightness in counts (rather than electrons)
#
            if (tmp_dict["exptime"]>0.01):
                efactor=tmp_dict["exptime"]
            else:
                efactor=1.0
#
            gtesta=tmp_dict["gaina"]-1.
            gtestb=tmp_dict["gainb"]-1.
            if ((abs(gtesta)<0.25)and(abs(gtestb)<0.25)):
#               The case where gains are 1... therefore units are electrons
                if (args.verbose):
                    print("# GAINA/B are consistent with units of electrons (GAINA:{:.3f},GAINB:{:.3f})".format(tmp_dict["gaina"],tmp_dict["gainb"]))
                gfactor=4.0
                tmp_dict["bunit"]='e-'
            else:
#               The case where gains are not 1... therefore units are already in counts
                if (args.verbose):
                    print("# GAINA/B are consistent with units of DN (GAINA:{:.3f},GAINB:{:.3f})".format(tmp_dict["gaina"],tmp_dict["gainb"]))
                gfactor=1.0
                tmp_dict["bunit"]='DN'
            tmp_dict["skyb"]=float(item[coldict["i.skybrite"]])/efactor/gfactor
            tmp_dict["skys"]=float(item[coldict["i.skysigma"]])/efactor/gfactor
#
#
            if (item[coldict["i.scampflg"]] is None):
                tmp_dict["sflag"]=2
            else:
                tmp_dict["sflag"]=int(item[coldict["i.scampflg"]])
            ccd_info[tmp_dict["ccdnum"]]=tmp_dict


##############################################################################
#   Fill in entries for missing CCDs (if there are any)
#   Also gather information for the first consistency checks and exposure averaged parameters.
#
    fwhm_old=[]
    elli_old=[]
    skyb=[]
    skys=[]
    img_ra=[]
    img_dec=[]
    bunit_chk=[]
    band_chk=[]
    expnum_chk=[]
    exptime_chk=[]
    airmass_chk=[]
    scamp_chk=[]
    numccd=0
    
    for iccd in range(1,63):
        if iccd in ccd_info:
#            print iccd,ccd_info[iccd]
            ccd_info[iccd]['data']=True
            numccd=numccd+1
            bunit_chk.append(ccd_info[iccd]['bunit'])
            band_chk.append(ccd_info[iccd]['band'])
            expnum_chk.append(ccd_info[iccd]['expnum'])
            exptime_chk.append(ccd_info[iccd]['exptime'])
            airmass_chk.append(ccd_info[iccd]['airmass'])
            scamp_chk.append(ccd_info[iccd]['sflag'])
            fwhm_old.append(ccd_info[iccd]['fwhm_old'])
            elli_old.append(ccd_info[iccd]['elli_old'])
            skyb.append(ccd_info[iccd]['skyb'])
            skys.append(ccd_info[iccd]['skys'])
            for ra_val in ccd_info[iccd]['img_ra']:
                img_ra.append(ra_val)
            for dec_val in ccd_info[iccd]['img_dec']:
                img_dec.append(dec_val)
        else:
            tmp_dict={}
            tmp_dict['data']=False
            ccd_info[iccd]=tmp_dict
    exp_rec["numccd"]=numccd
    exp_rec["scamp_sum"]=sum(scamp_chk)
##############################################################################
#   Perform the actual checks (and set exposure level values when these pass muster.
#   First check exptime.
#
    uniq_exptime_chk=list(set(exptime_chk))
    uniq_band_chk=list(set(band_chk))
    uniq_bunit_chk=list(set(bunit_chk))
    uniq_expnum_chk=list(set(expnum_chk))
    uniq_airmass_chk=list(set(airmass_chk))
    if (len(uniq_exptime_chk) != 1):
        if (len(uniq_exptime_chk) > 1):
            print "WARNING: Other than one exptime?: ",uniq_exptime_chk
            print "WARNING: Using exptime: ",uniq_exptime_chk[0]
            exp_rec['exptime']=uniq_exptime_chk[0]
        else:
            print "Aborting: No exptime?: "
            exit(1)
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
            print "Aborting: No airmass?: "
            exit(1)
    exp_rec['airmass']=uniq_airmass_chk[0]
#
#   Now check bunit for consistency
#
    if (len(uniq_bunit_chk) != 1):
        if (len(uniq_bunit_chk) > 1):
            print "WARNING: Other than one bunit?: ",uniq_bunit_chk
            print "WARNING: Using bunit: ",uniq_bunit_chk[0]
            exp_rec['airmass']=uniq_bunit_chk[0]
        else:
            print "WARNING: Assuming bunit = DN"
            exp_rec['bunit']='DN'
    else:
        exp_rec['bunit']=uniq_bunit_chk[0]
#
#   Now check band (also check that it is a sanctioned value)
#
    if (len(uniq_band_chk) != 1):
        print "Abort: Other than one band identified?: ",uniq_band_chk
        exit(1)
    else:
        if (uniq_band_chk[0] in ['u','g','r','i','z','Y','VR']):
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
    afwhm_old=numpy.array(fwhm_old) 
    aelli_old=numpy.array(elli_old) 
    askyb=numpy.array(skyb)
    askys=numpy.array(skys)
    if (len(afwhm_old) > 0):
        exp_rec["fwhm_img"]=afwhm_old.mean()
    else:
        exp_rec["fwhm_img"]=-1.0
    if (len(aelli_old) > 0):
        exp_rec["ellip_avg"]=aelli_old.mean()
        exp_rec["ellip_rms"]=aelli_old.std()
    else:
        exp_rec["ellip_avg"]=-1.0
        exp_rec["ellip_rms"]=-1.0
    if (len(askyb) > 0):
        exp_rec["skyb_avg"]=askyb.mean()
        exp_rec["skyb_rms"]=askyb.std()
    else:
        exp_rec["skyb_avg"]=-1.0
        exp_rec["skyb_rms"]=-1.0
    if (len(askys) > 0):
        exp_rec["skys_avg"]=askys.mean()
        exp_rec["skys_rms"]=askys.std()
    else:
        exp_rec["skys_avg"]=-1.0
        exp_rec["skys_rms"]=-1.0
#
#   Summarize sanity check of CCD level information
#
#    if (args.verbose):
#        for iccd in range(1,63):
#            if ccd_info[iccd]:
#                print iccd,ccd_info[iccd]['data']
#            else:
#                print iccd,' Still no entry (this should not happen!!!!!!!)'
###############################################################################
#   Get more exposure level information
#
    queryitems = ["e.filename", "e.expnum", "e.nite", "e.date_obs", "e.time_obs", "e.mjd_obs", "e.obstype", "e.band", "e.camsym", "e.object", "e.telra", "e.teldec", "e.tradeg", "e.tdecdeg", "e.airmass", "e.exptime", "e.propid", "e.program", "e.field" ]
    querylist = ",".join(queryitems)
    coldict={}
    for index, item in enumerate(queryitems):
         coldict[item]=index

    query = """select %s from %sexposure e where e.expnum=%s """ % ( querylist, db_Schema, exp_rec['expnum'] )

    print "# Executing query to obtain exposure level metadata"
    if args.verbose:
        print("# {:s}".format(query))
    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)

    num_exp_recs=0
    for item in cur:
        num_exp_recs=num_exp_recs+1
        if (num_exp_recs > 1):
            print "# WARNING: multiple exposure records found for this exposure? (num_exp_recs=",num_exp_recs,")"
        exp_rec["exp_fname"]=item[coldict["e.filename"]]
        exp_rec["dateobs"]=item[coldict["e.date_obs"]]
        exp_rec["obstype"]=item[coldict["e.obstype"]]
        exp_rec["telra"]=item[coldict["e.telra"]]
        exp_rec["teldec"]=item[coldict["e.teldec"]]
        exp_rec["ra"]=float(item[coldict["e.tradeg"]])
        exp_rec["dec"]=float(item[coldict["e.tdecdeg"]])
        exp_rec["object"]=item[coldict["e.object"]]
        exp_rec["mjd_obs"]=float(item[coldict["e.mjd_obs"]])
#
#       And then there are some more sanity checks.
#
        if (item[coldict["e.band"]] is None):
            if (exp_rec["band"] != 'Unknown'):
                print "WARNING: BAND miss-match between exposure-level (Unknown) and image/catalog-level (",exp_rec["band"],") queries.  Using image/cat-result."
        else:
            if (item[coldict["e.band"]] != exp_rec["band"]):
                print "WARNING: BAND miss-match between exposure-level (",item[coldict["e.band"]],") and image/catalog-level (",exp_rec["band"],") queries.  Using image/cat-result."
#
        if (float(item[coldict["e.exptime"]]) != exp_rec["exptime"]):
                print "WARNING: EXPTIME miss-match between exposure-level (",item[coldict["e.exptime"]],") and image/catalog-level (",exp_rec["band"],") queries.  Using image/cat-result."

        if (item[coldict["e.airmass"]] is None):
           print "WARNING: BAND miss-match between exposure-level (Unknown) and image/catalog-level (",exp_rec["airmass"],") queries.  Using image/cat-result."
        else:
            if (item[coldict["e.airmass"]] != exp_rec["airmass"]):
                print "WARNING: AIRMASS miss-match between exposure-level (",item[coldict["e.airmass"]],") and image/catalog-level (",exp_rec["airmass"],") queries.  Using image/cat-result."
#
#       Work out whether the exposure is part of the general survey, SN, or other.
#
        if (item[coldict["e.program"]]=="survey"):
            exp_rec["program"]='survey'
            exp_rec["sn_field"]=False
            exp_rec["survey_field"]=True
        elif (item[coldict["e.program"]]=="supernova"):
            exp_rec["program"]='SN'
            exp_rec["sn_field"]=True
            exp_rec["survey_field"]=False
        elif (item[coldict["e.program"]]=="photom-std-field"):
            exp_rec["program"]='phot-std'
            exp_rec["sn_field"]=False
            exp_rec["survey_field"]=False
        else:
            if (exp_rec["obstype"] in ['zero','dark','dome flat','sky flat']):
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
    queryitems = ["q.astromndets_ref_highsn", "q.astromchi2_ref_highsn", "q.astromsigma_ref_highsn_1", "q.astromsigma_ref_highsn_2",
    "q.astromoffset_ref_highsn_1", "q.astromoffset_ref_highsn_2"]
    coldict={}
    for index, item in enumerate(queryitems):
        coldict[item]=index
    querylist = ",".join(queryitems)
    query = """select %s from %sscamp_qa q, %swgb wgb where wgb.reqnum=%s and wgb.unitname='%s' and wgb.attnum=%s and wgb.filetype='xml_scamp' and wgb.filename=q.filename """ % (querylist,db_Schema,db_Schema,ReqNum,UnitName,AttNum)
    print "# Executing query to obtain astrometry QA for this exposure"
    if args.verbose:
        print("# {:s}".format(query))
    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)

    astrom_good=True
    num_ast_recs=0
    for item in cur:
        num_ast_recs=num_ast_recs+1
        if (num_ast_recs > 1):
            print "# WARNING: multiple exposure records found for this exposure? (num_ast_recs=",num_ast_recs,")"
        if (item[coldict["q.astromndets_ref_highsn"]] is None):
            print"# WARNING: Probable astrometric solution failure: No/null value for scamp_qa.astromndets_ref_highsn"
            astrom_good=False
            exp_rec["astrom_ndets"]=-1
        else:
            exp_rec["astrom_ndets"]=int(item[coldict["q.astromndets_ref_highsn"]])
        if (item[coldict["q.astromchi2_ref_highsn"]] is None):
            print"# WARNING: Probable astrometric solution failure: No/null value for scamp_qa.astromchi2_ref_highsn"
            astrom_good=False
            exp_rec["astrom_chi2"]=-1.0
        else:
            exp_rec["astrom_chi2"]=float(item[coldict["q.astromchi2_ref_highsn"]])
        if (item[coldict["q.astromsigma_ref_highsn_1"]] is None):
            print"# WARNING: Probable astrometric solution failure: No/null value for scamp_qa.astromsigma_ref_highsn_1"
            astrom_good=False
            exp_rec["astrom_sig1"]=-1.0
            sigx=-1.0
        else:
            exp_rec["astrom_sig1"]=float(item[coldict["q.astromsigma_ref_highsn_1"]])
            sigx=float(item[coldict["q.astromsigma_ref_highsn_1"]])
        if (item[coldict["q.astromsigma_ref_highsn_2"]] is None):
            print"# WARNING: Probable astrometric solution failure: No/null value for scamp_qa.astromsigma_ref_highsn_2"
            astrom_good=False
            exp_rec["astrom_sig2"]=-1.0
            sigy=-1.0
        else:
            exp_rec["astrom_sig2"]=float(item[coldict["q.astromsigma_ref_highsn_2"]])
            sigy=float(item[coldict["q.astromsigma_ref_highsn_2"]])
        if (item[coldict["q.astromoffset_ref_highsn_1"]] is None):
            print"# WARNING: Probable astrometric solution failure: No/null value for scamp_qa.astromoffset_ref_highsn_1"
            astrom_good=False
            exp_rec["astrom_off1"]=-1.0
        else:
            exp_rec["astrom_off1"]=float(item[coldict["q.astromoffset_ref_highsn_1"]])
        if (item[coldict["q.astromoffset_ref_highsn_2"]] is None):
            print"# WARNING: Probable astrometric solution failure: No/null value for scamp_qa.astromoffset_ref_highsn_2"
            astrom_good=False
            exp_rec["astrom_off2"]=-1.0
        else:
            exp_rec["astrom_off2"]=float(item[coldict["q.astromoffset_ref_highsn_2"]])
        exp_rec["astrom_rms2"]=numpy.sqrt((sigx*sigx)+(sigy*sigy))

    if ((exp_rec["astrom_sig1"] < 0.0001)or(exp_rec["astrom_sig2"] < 0.0001)):
        print("# WARNING: Probable astrometric solution failure: astrom_sig1,2 = {:7.4f},{:7.4f}".format(exp_rec["astrom_sig1"],exp_rec["astrom_sig2"]))
        astrom_good=False
    if (exp_rec["astrom_rms2"] > 0.500):
        print("# WARNING: Probable astrometric solution failure: astrom_rms2 = {:.3f}".format(exp_rec["astrom_rms2"]))
        astrom_good=False

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
    aimg_ra=numpy.array(img_ra) 
    aimg_dec=numpy.array(img_dec)
    if ((len(aimg_ra) < 2)or(len(aimg_dec) < 2)):
        print "Warning: No CCD corners present in database for these items"
        exp_rec["ra_cen"]=exp_rec["ra"]
        exp_rec["dec_cen"]=exp_rec["dec"]
    else:
        if ((len(aimg_ra) < (4*exp_rec["numccd"]))or(len(aimg_dec) < (4*exp_rec["numccd"]))):
            print "Warning: Some CCDs may be missing RA/DEC corners)! Working with what was present..."
        ra_min=numpy.min(aimg_ra)
        ra_max=numpy.max(aimg_ra)
        dec_min=numpy.min(aimg_dec)
        dec_max=numpy.max(aimg_dec)
        if ((ra_max-ra_min)>180.0):
#           print "IT HAPPENS ",ra_min,ra_max
            new_img_ra=[]
            for ra_val in img_ra:
                if (ra_val > 180.0):
                    new_img_ra.append(ra_val-360.0)
                else:
                    new_img_ra.append(ra_val)
            aimg_ra=numpy.array(new_img_ra)
        ra_cen=aimg_ra.mean()
        dec_cen=aimg_dec.mean()
        if (ra_cen < 0.0):
            ra_cen=ra_cen+360.0
        exp_rec["ra_cen"]=ra_cen
        exp_rec["dec_cen"]=dec_cen
    near_branch=False
    if ((ra_cen < fp_rad)or((ra_cen > 180.0)and(ra_cen > (360.-fp_rad)))):
       near_branch=True


#
#   Report exposure information for current 
#
    print "########################################"
    print "#   date_obs: ",exp_rec["dateobs"]
    print "#   Exposure: ",exp_rec["expnum"]
    print "#       Band: ",exp_rec["band"]
    print("#    Exptime: {:.1f} ".format(exp_rec["exptime"]))
    print "#      BUNIT: ",exp_rec["bunit"]
    print "# "
    print "#    Obstype: ",exp_rec["obstype"]
    print "#     Object: ",exp_rec["object"]
    print "#    Program: ",exp_rec["program"]
    print "# "
    print("#      Telescope(Ra,Dec): {:9.5f} {:9.5f} ".format(exp_rec["ra"],exp_rec["dec"]))
    print "#"
    print("# Image Centroid(Ra,Dec): {:9.5f} {:9.5f} ".format(exp_rec["ra_cen"],exp_rec["dec_cen"]))
    print "#"
    print "# "
    if (astrom_good):
        print "# Astrometric solution appears OK "
    else:
        print "# WARNING: Probable astrometric solution failure"
    print "# Astrometry summary:"
    print("#                    high S/N      ")
    print("#        ndets: {:d} ".format(exp_rec["astrom_ndets"]))
    print("#         chi2: {:.2f} ".format(exp_rec["astrom_chi2"]))
    print("#   astrom_sig: {:7.4f},{:7.4f} ".format(exp_rec["astrom_sig1"],exp_rec["astrom_sig2"]))
    print("#   astrom_off: {:7.4f},{:7.4f} ".format(exp_rec["astrom_off1"],exp_rec["astrom_off2"]))
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
    print "# Executing queries for NOMAD and APASS in vicinity of the exposure."
    queryitems_nomad=["n.ra", "n.dec", "n.sra", "n.sde", "n.mura", "n.mudec", "n.b", "n.v", "n.r", "n.j", "n.flags"] 
    queryitems_apass=["a.ra", "a.dec", "a.dra", "a.ddec", "a.b", "a.b_err", "a.v", "a.v_err", "a.g", "a.g_err", "a.r", "a.r_err", "a.i", "a.i_err" ] 
    coldict_nomad={}
    for index, item in enumerate(queryitems_nomad):
        coldict_nomad[item]=index
    querylist_nomad = ",".join(queryitems_nomad)
    coldict_apass={}
    for index, item in enumerate(queryitems_apass):
        coldict_apass[item]=index
    querylist_apass = ",".join(queryitems_apass)
#
#   mag_constraints for NOMAD, APASS, queries
#
    nomad_mag_constraint=''
    apass_mag_constraint=''
    if (exp_rec["band"] == "u"):
        nomad_mag_constraint=' and n.b<%.1f ' % (blimit)
        apass_mag_constraint=' and a.g<%.1f ' % (glimit)
    if (exp_rec["band"] == "g"):
        nomad_mag_constraint=' and n.b<%.1f ' % (blimit)
        apass_mag_constraint=' and a.g<%.1f ' % (glimit)
    elif ((exp_rec["band"] == "r")or(exp_rec["band"] == "VR")):
        nomad_mag_constraint=' and n.b<%.1f and n.j<%.1f ' % (blimit,jlimit)
        apass_mag_constraint=' and a.r<%.1f ' % (rlimit)
    elif (exp_rec["band"] == "i"):
        nomad_mag_constraint=' and n.j<%.1f ' % (jlimit)
        apass_mag_constraint=' and a.i<%.1f ' % (ilimit)
    elif ((exp_rec["band"] == "z")or(exp_rec["band"] == "Y")):
        nomad_mag_constraint=' and n.j<%.1f ' % (jlimit)
        apass_mag_constraint=' and a.i<%.1f ' % (ilimit)
    print "# Applying NOMAD MAG CONSTRAINT ",nomad_mag_constraint," for work with ",exp_rec["band"],"-band data."
    print "# Applying APASS MAG CONSTRAINT ",apass_mag_constraint," for work with ",exp_rec["band"],"-band data."

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
        dec1=dec_cen+fp_rad
        dec2=dec_cen-fp_rad
        query_nomad = """select %s from nomad n where n.dec < %.7f and n.dec > %.7f and ( n.ra < %.7f or n.ra > %.7f ) %s order by n.ra """ % ( querylist_nomad, dec1,dec2,ra1,ra2,nomad_mag_constraint )
        query_apass = """select %s from apass_dr7 a where a.dec < %.7f and a.dec > %.7f and ( a.ra < %.7f or a.ra > %.7f ) %s order by a.ra """ % ( querylist_apass, dec1,dec2,ra1,ra2,apass_mag_constraint)
    else:
#
#       NOTE need something better if we get closer than 2 degrees from pole
#
        cosdec_cen=numpy.cos(dec_cen*deg2rad)
        ra1=ra_cen+(fp_rad/cosdec_cen)
        ra2=ra_cen-(fp_rad/cosdec_cen)
        dec1=dec_cen+fp_rad
        dec2=dec_cen-fp_rad
        query_nomad = """select %s from nomad n where n.dec < %.7f and n.dec > %.7f and n.ra < %.7f and n.ra > %.7f %s order by n.ra """ % ( querylist_nomad, dec1,dec2,ra1,ra2,nomad_mag_constraint )
        query_apass = """select %s from apass_dr7 a where a.dec < %.7f and a.dec > %.7f and a.ra < %.7f and a.ra > %.7f %s order by a.ra """ % ( querylist_apass, dec1,dec2,ra1,ra2,apass_mag_constraint )
    if args.verbose:
        print("# {:s}".format(query_nomad))
        print("# {:s}".format(query_apass))
    cur.arraysize = 1000 # get 1000 at a time when fetching
    t2=time.time()
    print("# TIMING (pre-NOMAD/APASS): {:.2f}".format((t2-t1)))
#
#       If appropriate obtain APASS data  (currently always)
#
#   if (exp_rec["band"] in ["g","r","i"]):

#
#   Setup rules for parsing APASS
#
    apass_parse={"ra":"a.ra","dec":"a.dec"}
    if (exp_rec["band"] in ["u","g"]):
        apass_parse["mag"]="a.g"
        apass_parse["dmag"]="a.g_err"
        apass_parse["mlimit"]=glimit
    elif (exp_rec["band"] in ["r","VR"]):
        apass_parse["mag"]="a.r"
        apass_parse["dmag"]="a.r_err"
        apass_parse["mlimit"]=rlimit
    elif (exp_rec["band"] in ["i","z","Y"]):
        apass_parse["mag"]="a.i"
        apass_parse["dmag"]="a.i_err"
        apass_parse["mlimit"]=ilimit

    cur.execute(query_apass)
    apass_cat=[]
    for item in cur:
        tmp_apassdic={}
        tmp_apassdic["ra"]=float(item[coldict_apass[apass_parse["ra"]]])
        tmp_apassdic["dec"]=float(item[coldict_apass[apass_parse["dec"]]])
        mag=float(item[coldict_apass[apass_parse["mag"]]])
        if (mag > apass_parse["mlimit"]):
            tmp_apassdic["mag"]=99.0
            tmp_apassdic["dmag"]=99.0
        else:
            tmp_apassdic["mag"]=mag
            tmp_apassdic["dmag"]=float(item[coldict_apass[apass_parse["dmag"]]])
        if (tmp_apassdic["mag"]<98.0):
            apass_cat.append(tmp_apassdic)
    print "# "
    print "# Number of APASS entries found (after magnitude cuts) is: ",len(apass_cat)
#
#   If appropriate (currently always) acquire NOMAD data
#
    t2a=time.time()
    cur.execute(query_nomad)
    nomad_cat=[]
    for item in cur:
        tmp_nomaddic={}
        tmp_nomaddic["ra"]=float(item[coldict_nomad["n.ra"]])
        tmp_nomaddic["dec"]=float(item[coldict_nomad["n.dec"]])
        if (exp_rec["band"] in ["i","z","Y"]):
            jmag=float(item[coldict_nomad["n.j"]])
            if (jmag > jlimit):
                tmp_nomaddic["mag"]=99.0
            else:
                tmp_nomaddic["mag"]=jmag
            tmp_nomaddic["mag"]=jmag
        elif (exp_rec["band"] in ["u","g"]):
            bmag=float(item[coldict_nomad["n.b"]])+nomad_bmag_corr[int(tmp_nomaddic["dec"]+90.0)]
            if (bmag > blimit):
                tmp_nomaddic["mag"]=99.0
            else:
                if ((bmag > 18.37)and(bmag < 18.39)):
                    tmp_nomaddic["mag"]=99.0
                else:
                    tmp_nomaddic["mag"]=bmag
            tmp_nomaddic["mag"]=bmag
        elif (exp_rec["band"] in ["r","VR"]):
            jmag=float(item[coldict_nomad["n.j"]])
            bmag=float(item[coldict_nomad["n.b"]])+nomad_bmag_corr[int(tmp_nomaddic["dec"]+90.0)]
            if ((jmag > jlimit)or(bmag>blimit)):
                tmp_nomaddic["mag"]=99.0
            else:
                if ((bmag > 18.37)and(bmag < 18.39)):
                    tmp_nomaddic["mag"]=99.0
                else:
                    tmp_nomaddic["mag"]=0.333333333*((2.0*bmag)+jmag)
            tmp_nomaddic["mag"]=0.333333333*((2.0*bmag)+jmag)
        else:
            tmp_nomaddic["mag"]=99.0
        tmp_nomaddic["b"]=float(item[coldict_nomad["n.b"]])
        tmp_nomaddic["v"]=float(item[coldict_nomad["n.v"]])
        tmp_nomaddic["r"]=float(item[coldict_nomad["n.r"]])
        tmp_nomaddic["j"]=float(item[coldict_nomad["n.j"]])
        tmp_nomaddic["flags"]=item[coldict_nomad["n.flags"]]
        if (tmp_nomaddic["mag"]<98.0):
            nomad_cat.append(tmp_nomaddic)
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
    mcorr=mtime-aval-25.0-bval
    mag=[]
    magerr=[]
    emag=[]
    emagerr=[]
    rrad1=[]
    rrad2=[]
    fwld=[]
    frad=[]
    krad=[]
    magerr_hist=numpy.zeros((len(mbin)),dtype=numpy.float32)
    magnum_hist=numpy.zeros((len(mbin)),dtype=numpy.float32)

    nstar_found=0
    nomad_star_set=[]
    continue_nomad_ccorr=True
    continue_apass_ccorr=True
    mag_ccorr_nomad_limit=23.0
    mag_ccorr_apass_limit=23.0
    arate100=[]
    nrate100=[]
    amag100=[]
    nmag100=[]
    ancnt=0
    anfnd=0
    anfnd_sum=0
    nnfnd=0
    nncnt=0
    nnfnd_sum=0

###############################################################################
#   Columns to be retrieved from FITS tables.
#   cols_retrieve=["ALPHAWIN_J2000","DELTAWIN_J2000","MAG_AUTO","MAGERR_AUTO","SPREAD_MODEL","FWHM_WORLD",
#                  "A_IMAGE","B_IMAGE","THETA_IMAGE","ERRAWIN_IMAGE","ERRBWIN_IMAGE","ERRTHETAWIN_IMAGE",
#                  "FLUX_RADIUS","KRON_RADIUS","CLASS_STAR","FLAGS"]
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
    exp_cat=[]
    for iccd in range(1,63):
        if (ccd_info[iccd]['data']):
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
    match_starset=[]
    for object in new_expcat:
#       Make lists that are actually operated on to obtain histogram
        mag.append(object["MAG_AUTO"]+mcorr)
        magerr.append(object["MAGERR_AUTO"])
#       Now work on objects that are constrained to be stars (by SPREAD_MODEL).
        smval=abs(object["SPREAD_MODEL"])
        if ((smval < 0.0015)and(object["MAGERR_AUTO"]<0.1)):
#           These are for finding sizes.
            nstar_found=nstar_found+1
            aval=object["A_IMAGE"]
            bval=object["B_IMAGE"]
#            if ((aval>20)or(bval>20)or(object["FWHM_WORLD"]>20.)or(object["FLUX_RADIUS"]>20)or(object["KRON_RADIUS"]>20)):
#                print "# AHHHHHH: ",object
            rrad1.append(2.0*aval)
            rrad2.append(aval+bval)
            fwld.append(object["FWHM_WORLD"])
            frad.append(object["FLUX_RADIUS"])
            krad.append(object["KRON_RADIUS"])
#           And this structure is used for comparison to NOMAD/APASS.
            if (continue_nomad_ccorr)or(continue_apass_ccorr):
                tmp_starset={}
                tmp_starset["ra"]=object["ALPHAWIN_J2000"]
                tmp_starset["dec"]=object["DELTAWIN_J2000"]
                tmp_starset["mag"]=object["MAG_AUTO"]+mcorr
                tmp_starset["magerr"]=object["MAGERR_AUTO"]
                add_tmp_starset_record=False
#
#               Cross correlation vs. NOMAD
#
                if (continue_nomad_ccorr):
                    match_rec=binsrch_crosscorr_cat(tmp_starset["ra"],tmp_starset["dec"],nomad_cat,1.0)
                    nncnt=nncnt+1
                    if (len(match_rec)>0):
                        nnfnd=nnfnd+1
                        tmp_starset["nomad"]=match_rec[0]
                        tmp_starset["nomad"]["magcorr"]=tmp_starset["nomad"]["mag"]+nomad_mag_corr[exp_rec["band"]]
                        tmp_starset["nomad"]["magdiff"]=tmp_starset["mag"]-tmp_starset["nomad"]["mag"]-nomad_mag_corr[exp_rec["band"]]
                        add_tmp_starset_record=True

#
#                   Look for the point of diminishing returns (very few new matches found)
#
                    if (nncnt%100 == 0):
                        nrate100.append(nnfnd)  
                        nnfnd_sum=nnfnd_sum+nnfnd
                        nnfnd=0
                        if ((len(nrate100)>10)and(nnfnd_sum>500)):
                            if ((nrate100[-3]<n100_lim)and(nrate100[-2]<n100_lim)and(nrate100[-1]<n100_lim)):
                                continue_nomad_ccorr=False
#
#               Cross correlation vs. APASS
#
#               if (exp_rec["band"] in ["g","r","i"]):
                if (continue_apass_ccorr):
                    match_rec=binsrch_crosscorr_cat(tmp_starset["ra"],tmp_starset["dec"],apass_cat,1.0)
                    ancnt=ancnt+1
                    if (len(match_rec)>0):
                        anfnd=anfnd+1
                        tmp_starset["apass"]=match_rec[0]
                        tmp_starset["apass"]["magcorr"]=tmp_starset["apass"]["mag"]+apass_mag_corr[exp_rec["band"]]
                        tmp_starset["apass"]["magdiff"]=tmp_starset["mag"]-tmp_starset["apass"]["mag"]-apass_mag_corr[exp_rec["band"]]
                        add_tmp_starset_record=True
#
#                   Look for the point of diminishing returns (very few new matches found)
#
                    if (ancnt%100 == 0):
                        arate100.append(anfnd)
                        anfnd_sum=anfnd_sum+anfnd
                        anfnd=0
                        if ((len(arate100)>10)and(anfnd_sum>500)):
                            if ((arate100[-3]<a100_lim)and(arate100[-2]<a100_lim)and(arate100[-1]<a100_lim)):
                                continue_apass_ccorr=False
                # If add_tmp_starset_record is True then add the record to the matched star set.
                if (add_tmp_starset_record):
                    match_starset.append(tmp_starset)

    t4=time.time()
    exp_rec["numobj"]=len(mag)
    print "# "
    print("# Number of objects from exposure catalogs: {:d}".format(exp_rec["numobj"]))

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
#
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

#
#   Proceed to assessing the atmospheric opacity.
#
    transparency_calc=False
#   if (exp_rec["band"] in ["g","r","i"]):
    if (num_matches['apass']<100):
#       For a limited number of stars do not make comparison... output will be sentinel value.
        exp_rec["apass_magdiff"]=99.0
        exp_rec["apass_kmagdiff"]=99.0
        exp_rec["apass_num"]=num_matches['apass']
        if (args.qaplot):
#           If QA plots have been requested then output an empty QA plot that gives 
#           information about the failure.
            plt.figure()
            plt.subplot(2,1,1)
            if (len(sat_check['apass'])>0):
                mag_des,mag_apass,mag_diff=zip(*sat_check['apass'])
            else:
                mag_des=[]
                mag_apass=[]
                mag_diff=[]
            plt.scatter(mag_des,mag_apass,marker='.',color='blue')
            plt.plot([10,18],[10,18],color='red',linewidth=1)
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            if (exp_rec["band"] in ["u","g"]):
                plt.ylabel('APASS g\'')
            elif (exp_rec["band"] in ["r","VR"]):
                plt.ylabel('APASS r\'')
            elif (exp_rec["band"] in ["i","z","Y"]):
                plt.ylabel('APASS i\'')
            else:
                plt.ylabel('APASS [unknown]')
            plt.text(10.5,17,'Number of APASS matches: %d' % (num_matches['apass']))
            if (not(astrom_good)):
                plt.text(10.5,16,'Probable failure to find an astrometric solution')
            plt.subplot(2,1,2)
            plt.plot([10,18],[0,0],color='red',linewidth=3)
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            if (exp_rec["band"] in ["u","g"]):
                plt.ylabel('DES - APASS g\'')
            elif (exp_rec["band"] in ["r","VR"]):
                plt.ylabel('DES - APASS r\'')
            elif (exp_rec["band"] in ["i","z","Y"]):
                plt.ylabel('DES - APASS i\'')
            else:
                plt.ylabel('DES - APASS [unknown]')
            plt.savefig("%s_apass.png" % (args.froot))
    else:
#       Workhorse case... find offset
        if (len(sat_check['apass'])>0):
            mag_des,mag_apass,mag_diff=zip(*sat_check['apass'])
        else:
            mag_des=[]
            mag_apass=[]
            mag_diff=[]
        amag_apass=numpy.array(mag_apass)
        amag_diff=numpy.array(mag_diff)

#       numpy.where(numpy.logical_and((numpy.logical_and((r2_acolor>=cntcut),(i2_acolor>=cntcut))),(numpy.logical_and((z2_acolor>=cntcut),(g1_acolor>75.)))
        cut_val=numpy.where(amag_apass > sat_limit)
        med_magdiff=numpy.median(amag_diff[cut_val])
        std_magdiff=numpy.std(amag_diff[cut_val])
        print("# ")
        print("# (PRELIM APASS) Median magnitude offset: {:8.3f} {:8.3f} ".format(med_magdiff,std_magdiff))
#
#       Perform a rejection iteration.
#
        mag_fit=[]
        mag_reject=[]
        for record in sat_check['apass']:
            junkval=record[0]-record[1]
            reject_val=False
            if (abs((junkval-med_magdiff)/std_magdiff)>3.0): reject_val=True
            if (record[1] < sat_limit): reject_val=True
            if (reject_val):
                mag_reject.append(record)
            else:
                mag_fit.append(record)
        if (len(mag_fit)>0):
            mag_des,mag_apass,mag_diff=zip(*mag_fit)
        else:
            mag_des=[]
            mag_apass=[]
            mag_diff=[]
        if (len(mag_reject)>0):
            mag_des_reject,mag_apass_reject,mag_diff_reject=zip(*mag_reject)
        else:
            mag_des_reject=[]
            mag_apass_reject=[]
            mag_diff_reject=[]
#
#       Go on and calculate magnitude(s) of opacity
#
        amag_diff=numpy.array(mag_diff)
        mindes=min(mag_des)
        maxdes=max(mag_des)
        med_magdiff=numpy.median(amag_diff)
        std_magdiff=numpy.std(amag_diff)
        exp_rec["apass_magdiff"]=med_magdiff
        exp_rec["apass_kmagdiff"]=med_magdiff-apass_mag_corr[exp_rec["band"]]+apass_kmag_corr[exp_rec["band"]]-(kterm[band2i[exp_rec["band"]]]*exp_rec["airmass"])
        exp_rec["apass_num"]=num_matches['apass']

        if (exp_rec["band"] in ["g","r"]):
            transparency_calc=True
            exp_rec["magdiff"]=med_magdiff
            exp_rec["cloud_cat"]="APASS"
            print("# (USING APASS) Median magnitude offset: {:8.3f} {:8.3f} ".format(med_magdiff,std_magdiff))
        else:
            print("# (NOUSE APASS) Median magnitude offset: {:8.3f} {:8.3f} ".format(med_magdiff,std_magdiff))
#
#       Generate APASS QA plot.
#
        if (args.qaplot):
            plt.figure()
            plt.subplot(2,1,1)
            plt.scatter(mag_des,mag_apass,marker='.',color='blue')
            if (len(mag_des_reject) > 0):
                plt.scatter(mag_des_reject,mag_apass_reject,marker='.',color='magenta')
            plt.plot([mindes,maxdes],[mindes,maxdes],color='red',linewidth=1)
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            if (exp_rec["band"] in ["u","g"]):
                plt.ylabel('APASS g\'')
            elif (exp_rec["band"] in ["r","VR"]):
                plt.ylabel('APASS r\'')
            elif (exp_rec["band"] in ["i","z","Y"]):
                plt.ylabel('APASS i\'')
            else:
                plt.ylabel('APASS [unknown]')
            plt.subplot(2,1,2)
            plt.scatter(mag_des,mag_diff,marker='.',color='blue')
            if (len(mag_des_reject) > 0):
                plt.scatter(mag_des_reject,mag_diff_reject,marker='.',color='magenta')
            plt.plot([mindes,maxdes],[med_magdiff,med_magdiff],color='red',linewidth=3)
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            if (exp_rec["band"] in ["u","g"]):
                plt.ylabel('DES - APASS g\'')
            elif (exp_rec["band"] in ["r","VR"]):
                plt.ylabel('DES - APASS r\'')
            elif (exp_rec["band"] in ["i","z","Y"]):
                plt.ylabel('DES - APASS i\'')
            else:
                plt.ylabel('DES - APASS [unknown]')
            plt.savefig("%s_apass.png" % (args.froot))
#           plt.show()
#
#   Repeat for NOMAD
#
    if (num_matches['nomad']<100):
#       For a limited number of stars do not make comparison... output will be sentinel value.
        exp_rec["nomad_magdiff"]=99.0
        exp_rec["nomad_kmagdiff"]=99.0
        exp_rec["nomad_num"]=num_matches['nomad']
        if (args.qaplot):
#           If QA plots have been requested then output an empty QA plot that gives 
#           information about the failure.
            plt.figure()
            plt.subplot(2,1,1)
            if (len(sat_check['nomad'])>0):
                mag_des,mag_nomad,mag_diff=zip(*sat_check['nomad'])
            else:
                mag_des=[]
                mag_nomad=[]
                mag_diff=[]
            plt.scatter(mag_des,mag_nomad,marker='.',color='blue')
            plt.plot([10,18],[10,18],color='red',linewidth=1)
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            if (exp_rec["band"] in ["i","z","Y"]):
                plt.ylabel('NOMAD J')
            elif (exp_rec["band"] in ["u","g"]):
                plt.ylabel('NOMAD B')
            elif (exp_rec["band"] in ["r","VR"]):
                plt.ylabel('NOMAD [(2*B+J)/3]')
            else:
                plt.ylabel('NOMAD [unknown]')
            plt.text(10.5,17,'Number of NOMAD matches: %d' % (num_matches['nomad']))
            if (not(astrom_good)):
                plt.text(10.5,16,'Probable failure to find an astrometric solution')
            plt.subplot(2,1,2)
            plt.plot([10,18],[0,0],color='red',linewidth=3)
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            if (exp_rec["band"] in ["i","z","Y"]):
                plt.ylabel('DES - NOMAD J')
            elif (exp_rec["band"] in ["u","g"]):
                plt.ylabel('DES - NOMAD B')
            elif (exp_rec["band"] in ["r","VR"]):
                plt.ylabel('DES - NOMAD [(2*B+J)/3]')
            else:
                plt.ylabel('DES - NOMAD [unknown]')
            plt.savefig("%s_nomad.png" % (args.froot))
    else:
#       Workhorse case... find offset
        if (len(sat_check['nomad'])>0):
            mag_des,mag_nomad,mag_diff=zip(*sat_check['nomad'])
        else:
            mag_des=[]
            mag_nomad=[]
            mag_diff=[]
        amag_nomad=numpy.array(mag_nomad)
        amag_diff=numpy.array(mag_diff)

#       numpy.where(numpy.logical_and((numpy.logical_and((r2_acolor>=cntcut),(i2_acolor>=cntcut))),(numpy.logical_and((z2_acolor>=cntcut),(g1_acolor>75.)))
        cut_val=numpy.where(amag_nomad > sat_limit)
        med_magdiff=numpy.median(amag_diff[cut_val])
        std_magdiff=numpy.std(amag_diff[cut_val])
        print("# (PRELIM NOMAD) Median magnitude offset: {:8.3f} {:8.3f} ".format(med_magdiff,std_magdiff))
#
#       Perform a rejection iteration.
#
        mag_fit=[]
        mag_reject=[]
        for record in sat_check['nomad']:
            junkval=record[0]-record[1]
            reject_val=False
            if (abs((junkval-med_magdiff)/std_magdiff)>3.0): reject_val=True
            if (record[1] < sat_limit): reject_val=True
            if (reject_val):
                mag_reject.append(record)
            else:
                mag_fit.append(record)
        if (len(mag_fit)>0):
            mag_des,mag_nomad,mag_diff=zip(*mag_fit)
        else:
            mag_des=[]
            mag_nomad=[]
            mag_diff=[]
        if (len(mag_reject)>0):
            mag_des_reject,mag_nomad_reject,mag_diff_reject=zip(*mag_reject)
        else:
            mag_des_reject=[]
            mag_nomad_reject=[]
            mag_diff_reject=[]

#
#       Go on and calculate magnitude(s) of opacity
#
        amag_diff=numpy.array(mag_diff)
        mindes=min(mag_des)
        maxdes=max(mag_des)
        med_magdiff=numpy.median(amag_diff)
        std_magdiff=numpy.std(amag_diff)
        exp_rec["nomad_magdiff"]=med_magdiff
        exp_rec["nomad_kmagdiff"]=med_magdiff-nomad_mag_corr[exp_rec["band"]]+nomad_kmag_corr[exp_rec["band"]]-(kterm[band2i[exp_rec["band"]]]*exp_rec["airmass"])
        exp_rec["nomad_num"]=num_matches['nomad']

#
#       If APASS has not provided a result then use NOMAD result otherwise report result 
#
        if(not(transparency_calc)):
            transparency_calc=True
            exp_rec["magdiff"]=med_magdiff
            exp_rec["cloud_cat"]="NOMAD"
            print("# (USING NOMAD) Median magnitude offset: {:8.3f} {:8.3f} ".format(med_magdiff,std_magdiff))
        else:
            print("# (NOUSE NOMAD) Median magnitude offset: {:8.3f} {:8.3f} ".format(med_magdiff,std_magdiff))
#
#       Generate NOMAD QA plot.
#
        if (args.qaplot):
            plt.figure()
            plt.subplot(2,1,1)
            plt.scatter(mag_des,mag_nomad,marker='.',color='blue')
            if (len(mag_des_reject) > 0):
                plt.scatter(mag_des_reject,mag_nomad_reject,marker='.',color='magenta')
            plt.plot([mindes,maxdes],[mindes,maxdes],color='red',linewidth=1)
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            if (exp_rec["band"] in ["i","z","Y"]):
                plt.ylabel('NOMAD J')
            elif (exp_rec["band"] in ["u","g"]):
                plt.ylabel('NOMAD B')
            elif (exp_rec["band"] in ["r","VR"]):
                plt.ylabel('NOMAD [(2*B+J)/3]')
            else:
                plt.ylabel('NOMAD [unknown]')
            plt.subplot(2,1,2)
            plt.scatter(mag_des,mag_diff,marker='.',color='blue')
            if (len(mag_des_reject) > 0):
                plt.scatter(mag_des_reject,mag_diff_reject,marker='.',color='magenta')
            plt.plot([mindes,maxdes],[med_magdiff,med_magdiff],color='red',linewidth=3)
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            if (exp_rec["band"] in ["i","z","Y"]):
                plt.ylabel('DES - NOMAD J')
            elif (exp_rec["band"] in ["u","g"]):
                plt.ylabel('DES - NOMAD B')
            elif (exp_rec["band"] in ["r","VR"]):
                plt.ylabel('DES - NOMAD [(2*B+J)/3]')
            else:
                plt.ylabel('DES - NOMAD [unknown]')
            plt.savefig("%s_nomad.png" % (args.froot))
#           plt.show()

#
#   Now the case where the transparency calcultion has failed... give an
#

    if (not(transparency_calc)):
        exp_rec["magdiff"]=99.
        exp_rec["cloud_cat"]="Failed"
    t5=time.time()
    print("# TIMING (Transparency calculation): {:.2f}".format((t5-t4)))
    print(" ")
##############################################################################
#   FWHM 
#
    if (nstar_found > 2):                   
        arrad1=numpy.array(rrad1) 
        arrad2=numpy.array(rrad2) 
        afwld=numpy.array(fwld) 
        akrad=numpy.array(krad) 
        afrad=numpy.array(frad) 
        exp_rec["nstar_found"]=nstar_found
        exp_rec["fwhm_rrad1"]=pixsize*arrad1.mean()
        exp_rec["fwhm_rrad2"]=pixsize*arrad2.mean()
        exp_rec["fwhm_world"]=3600.*afwld.mean()
        exp_rec["fwhm_krad"]=2.0*pixsize*akrad.mean()
        exp_rec["fwhm_frad"]=2.0*pixsize*afrad.mean()
    else:
        exp_rec["nstar_found"]=nstar_found
        exp_rec["fwhm_rrad1"]=-1.
        exp_rec["fwhm_rrad2"]=-1.
        exp_rec["fwhm_world"]=-1.
        exp_rec["fwhm_krad"]=-1.
        exp_rec["fwhm_frad"]=-1.
    t5a=time.time()

    if (len(mag)<10):
#
#       If insufficient data is present to even begin to try then simply report faulure 
#       move along.  
#
        mag_thresh=99.9
        if (args.qaplot):   
#           If QA plots were requested then write a dummy plot that indicates failure
            plt.figure()
            plt.subplot(2,1,1)
#           plt.scatter(xplt,yplt1,marker='.',color='blue')
            plt.axis([8,26,0.5,10000])
#            plt.semilogy(xplt,yplt1,marker='.',ls='None',color='blue')
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            plt.ylabel('# objects')
            plt.text(9.,9000.,'Number of objects insufficient for histogram: %d' % (len(mag)))
            if (not(astrom_good)):
                plt.text(9.,7800,'Probable failure to find an astrometric solution')
            plt.subplot(2,1,2)
#            plt.scatter(xplt,yplt2,marker='.',color='blue')
            plt.axis([8,26,-0.01,0.5])
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            plt.ylabel('median(magerr_auto)')
            plt.savefig("%s_magplot.png" % (args.froot))
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
            xplt=[]
            yplt1=[]
            yplt2=[]
            for ibin, m1 in enumerate(mbin):
#               print("# {:8.3f} {:8.0f} {:8.3f} ".format(mbin[ibin],magnum_hist[ibin],magerr_hist[ibin]))
                xplt.append(mbin[ibin]-0.5)
                yplt1.append(magnum_hist[ibin])
                yplt2.append(magerr_hist[ibin])
            plt.figure()
            plt.subplot(2,1,1)
#           plt.scatter(xplt,yplt1,marker='.',color='blue')
            plt.axis([8,26,0.5,10000])
            plt.semilogy(xplt,yplt1,marker='.',ls='None',color='blue')
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            plt.ylabel('# objects')
    
            plt.subplot(2,1,2)
            plt.scatter(xplt,yplt2,marker='.',color='blue')
            plt.axis([8,26,-0.01,0.5])
            plt.xlabel('DES MAG_AUTO(%s)'%exp_rec["band"])
            plt.ylabel('median(magerr_auto)')
            plt.savefig("%s_magplot.png" % (args.froot))
#           plt.show()
    exp_rec["mag_thresh"]=mag_thresh
    t6=time.time()
    print("# TIMING (Object number counts): {:.2f}".format((t6-t5a)))

###############################################################################
#   Now the calculations for the Teff (and of course the individual components
#
#   Calculate F_eff
#
    if (exp_rec["fwhm_world"]>0.0):
        exp_rec["teff_f"]=(seeing_fid[exp_rec["band"]]*seeing_fid[exp_rec["band"]]/(exp_rec["fwhm_world"]*exp_rec["fwhm_world"]))
    else:
        exp_rec["teff_f"]=-1.0
#
#   Calculate B_eff
#
    if (exp_rec["skyb_avg"]>0.0):
        exp_rec["teff_b"]=sbrite_good[exp_rec["band"]]/exp_rec["skyb_avg"]
    else:
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
    print("# TIMING for exposure {:7d}: (all={:6.2f}, general_query={:6.2f}, APASS={:6.2f}, NOMAD={:6.2f}, crosscorr={:6.2f}, cat_compare={:6.2f} ,mag_analyze={:6.2f}) ".format(exp_rec["expnum"],(t7-t0),(t1-t0),(t3-t2a),(t2a-t2),(t4-t3),(t5-t4),(t6-t5)))
    print("# FWHM measures for exposure {:7d}: (NSTAR_FOUND,FWHM_IMG,WORLD,2*A,A+B,Kron,R50) {:6d} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} {:7.3f} ".format(exp_rec["expnum"],exp_rec["nstar_found"],exp_rec["fwhm_img"],exp_rec["fwhm_world"],exp_rec["fwhm_rrad1"],exp_rec["fwhm_rrad2"],exp_rec["fwhm_krad"],exp_rec["fwhm_frad"]))
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
    if((args.updateDB)or(args.DB_file)):
        prog_name='Unknown'
        if (exp_rec["sn_field"]):
            prog_name='supernova'
        if (exp_rec["survey_field"]):
            prog_name='survey'
#
#       First the columns that should always be present
#
#   EXPOSURENAME                        VARCHAR2(100)
#   EXPNUM                    NOT NULL NUMBER(10)
#   REQNUM                    NOT NULL NUMBER(38)
#   UNITNAME                  NOT NULL VARCHAR2(20)
#   ATTNUM                    NOT NULL NUMBER(5)
#   CAMSYM                    NOT NULL VARCHAR2(1)
#   PROCESSED                      VARCHAR2(10)
#   ACCEPTED                       VARCHAR2(10)
#   PROGRAM                        VARCHAR2(10)
#   ANALYST                        VARCHAR2(30)
#   ANALYST_COMMENT                    VARCHAR2(30)
#   T_EFF                          NUMBER(8,3)
#   F_EFF                          NUMBER(8,3)
#   B_EFF                          NUMBER(8,3)
#   C_EFF                          NUMBER(8,3)
#   FWHM_ASEC                      NUMBER(8,3)
#   ELLIPTICITY                        NUMBER(8,3)
#   SKYBRIGHTNESS                      BINARY_FLOAT
#   CLOUD_APASS                        NUMBER(8,3)
#   CLOUD_NOMAD                        NUMBER(8,3)
#   CLOUD_CATALOG                      VARCHAR2(10)
#   KLOUD_APASS                        NUMBER(8,3)
#   KLOUD_NOMAD                        NUMBER(8,3)
#   N_APASS                        NUMBER(8)
#   N_NOMAD                        NUMBER(8)

        db_cname="INSERT INTO %s(EXPOSURENAME,EXPNUM,REQNUM,UNITNAME,ATTNUM,CAMSYM,PROCESSED,ACCEPTED,PROGRAM,ANALYST,LASTCHANGED_TIME,T_EFF,F_EFF,B_EFF,C_EFF,FWHM_ASEC,ELLIPTICITY,SKYBRIGHTNESS" % ( db_table )
        db_value=""" VALUES ('%s', %d, %s, '%s', %s, '%s', '%s', '%s', '%s', '%s', sysdate, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.2f""" % (exp_rec["exp_fname"],exp_rec["expnum"],args.ReqNum, args.UnitName, args.AttNum, args.CamSym, dm_process,dm_accept,prog_name,analyst,exp_rec["teff"],exp_rec["teff_f"],exp_rec["teff_b"],exp_rec["teff_c"],exp_rec["fwhm_world"],exp_rec["ellip_avg"],exp_rec["skyb_avg"])

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

        if(args.updateDB):
            cur.execute(insert_command)
        if (args.DB_file):
            fdbout.write("{:s};\n".format(insert_command))

#
#   Write the summary of the assessment to STDOU
#

    ftxt.write("#                                                                                                         Astrometry \n")     
    ftxt.write("# Exposure                                     FWHM                             #        AIR    EXP  SCMP  HighSN     Depth       APASS           NOMAD      FWHM            \n")
    ftxt.write("#  Num   STATE    t_eff F_eff  B_eff   C_eff   World  Ellip    Sky_B   N_src   CCD BAND  MASS   TIME Flag sig1  sig2     Assess     dmag    #       dmag    #     Img    Object   \n")
    ftxt.write("#                                               [\"]           [DN/s]                             [s]      [\"]    [\"]     [mag]    [mag]            [mag]         \n")
#
    ftxt.write(" {:9d} {:4s}  {:6.2f} {:6.2f} {:6.2f} {:6.2f}   {:5.2f} {:6.3f} {:8.2f} {:7d}   {:3d} {:1s} {:1d} {:6.3f} {:6.1f} {:2d} {:6.3f} {:6.3f}  {:8.3f} {:8.3f} {:6d} {:8.3f} {:6d} {:6.3f} {:s} \n".format(
    exp_rec["expnum"],new_decide,exp_rec["teff"],exp_rec["teff_f"],exp_rec["teff_b"],exp_rec["teff_c"],
    exp_rec["fwhm_world"],exp_rec["ellip_avg"],exp_rec["skyb_avg"],exp_rec["numobj"],
    exp_rec["numccd"],exp_rec["band"],band2i[exp_rec["band"]],exp_rec["airmass"],exp_rec["exptime"],
    exp_rec["scamp_sum"],exp_rec["astrom_sig1"],exp_rec["astrom_sig2"],
    exp_rec["magdiff"],exp_rec["apass_magdiff"],exp_rec["apass_num"],exp_rec["nomad_magdiff"],exp_rec["nomad_num"],exp_rec["fwhm_img"],exp_rec["object"]))

    if(args.updateDB):
        dbh.commit()
        print "DB update complete and committed"
    if (args.csv):
        fcsv.close()
    if (args.DB_file):
        fdbout.close()

