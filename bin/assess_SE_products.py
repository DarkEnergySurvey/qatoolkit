#! /usr/bin/env python3
# $Id$
# $Rev: 44620 $:  # Revision of last commit.
# $LastChangedBy: rgruendl $:  # Author of last commit.

"""
Perfoms an assessent of exposures from a first/final cut run.  The quality 
of the exposures is based upon the seeing (FWHM), background, and extinction
due to clouds.


Syntax:
    quick_assess_RUN.py  -r run 
Arguments:
    Requires a run number for which the assessemnt for each exposure will be based.
     
"""

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
    import pandas as pd
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    import fitsio
#    import qatoolkit.assess_SE_legacy as SEleg
    import qatoolkit.assess_SE_plot   as SEplot
#    import matplotlib 
#    matplotlib.use('Agg')
#    import matplotlib.pyplot as plt
    
    version="0.2.11"
    db_table="firstcut_eval"

    parser = argparse.ArgumentParser(description='Assess whether the pipeline products of a FIRSTCUT/FINALCUT processing meet survey quality metrics.')
    parser.add_argument('-i', '--attid',   action='store', type=str, default=None, help='Processing attempt id', required=True)
#    parser.add_argument('-R', '--ReqNum',   action='store', type=str, default=None, help='Processing request number.')
#    parser.add_argument('-U', '--UnitName', action='store', type=str, default=None, help='Unit Name.')
#    parser.add_argument('-A', '--AttNum',   action='store', type=str, default=None, help='Attempt Number.')
    parser.add_argument('-C', '--CamSym',   action='store', type=str, default='D', help='CamSym (default=\"D\").')
    parser.add_argument('-l', '--list',     action='store', type=str, default=None, 
                        help='Optional list of catalogs for the assessment (otherwise query for "red_finalcut" catalogs associated with the Request Unit Attempt)')

    parser.add_argument('-a', '--analyst',  action='store', type=str, default='assess_SE_products.py', 
                        help='Provides value for analyst (default: assess_SE_products.py, \"None\" will use os.getlogin())')
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
    parser.add_argument('--scisection',    action='store', type=str, default=None, help='section of .desservices file with connection info for DESSCI (default=None)')
    parser.add_argument('-S', '--Schema', action='store', type=str, default=None, help='Schema')

    args = parser.parse_args()
    if (args.verbose):
        print("##### Initial arguments #####")
        print("Args: ",args)

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
    try:
        desdmfile = os.environ["des_services"]
    except KeyError:
        desdmfile = None
    dbh = despydb.desdbi.DesDbi(desdmfile, args.section, retry=True)

    if (not(args.scisection is None)):
        dbhsci = despydb.desdbi.DesDbi(desdmfile,args.scisection,retry=True)

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
    band2i={'u':0,'g':1,'r':2,'i':3,'z':4,'Y':5,'VR':6,'N964':7}
#
#   Old ellipticity limits (currently not used)
#
#    ellip_lim=0.13
#    ellip_good=0.07
#
    kolmogorov={'u':1.2,'g':1.103,'r':1.041,'i':1.00,'z':0.965,'Y':0.95,'VR':1.04,'N964':0.965}
    teff_lim={  'u':0.2,'g':0.2,  'r':0.3,  'i':0.3, 'z':0.3,  'Y':0.2,'VR':0.3,'N964':0.3}
    seeing_lim={}
    seeing_fid={}
#
#   Set seeing cutoff to be 1.6 times Kolmogov except at "Y" which should
#   be forced to match that at g-band
#
    for band in ['u','g','r','i','z','Y','VR','N964']:
        if (band == "Y"):
            seeing_lim[band]=1.6*kolmogorov['g']
        else:
            seeing_lim[band]=1.6*kolmogorov[band]
#       Commented version below was needed when using FWHM_WORLD
#        seeing_fid[band]=fwhm_DMtoQC_offset*0.9*kolmogorov[band]
#       Now fiducial value is additive (and applied to the FWHM_MEAN value coming from PSFex)
        seeing_fid[band]=0.9*kolmogorov[band]
#    print(kolmogorov)
#    print(seeing_lim)
#    print(seeing_fid)
#
#   Surface brightness limits from Eric Nielson which were derived "...from a few 
#   exposures from a photometric night in SV with little moon (20121215)"
#
    sbrite_good={"u":0.2,"g":1.05,"r":2.66,"i":7.87,"z":16.51,"Y":14.56,"VR":3.71,"N964":0.4}
    sbrite_lim={"u":0.8,"g":4.0,"r":9.58,"i":21.9,"z":50.2,"Y":27.6,"VR":13.58,"N964":20.2}
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
    use_PSF_mag=True
    if (use_PSF_mag):
        use_mag_type='PSF'
    else:
        use_mag_type='AUTO'

    cat_mag_corr={'apass':{'u':3.5,'g':0.205,'r':0.128,'i':0.112,'z':0.0,'Y':0.0,'VR':0.0,'N964':3.56},
                  'nomad':{'u':3.65,'g':0.341,'r':0.235,'i':1.398,'z':1.201,'Y':1.083,'VR':0.235,'N964':4.67},
                  'des':{'u':0.0,'g':0.248,'r':0.175,'i':0.078,'z':0.08,'Y':0.06,'VR':0.0,'N964':3.56}}
    cat_kmag_corr={'apass':{'u':0.0,'g':0.000,'r':0.000,'i':0.000,'z':0.0,'Y':0.0,'VR':0.0,'N964':0.0},
                  'nomad':{'u':0.0,'g':0.111,'r':0.109,'i':1.289,'z':1.139,'Y':1.022,'VR':0.0,'N964':0.0},
                  'des':{'u':0.0,'g':0.000,'r':0.000,'i':0.000,'z':0.0,'Y':0.0,'VR':0.0,'N964':0.0}}
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
##              print(columns[0],columns[1],spd,nomad_bmag_corr[spd])
#       f1.close()

#
#   A set of aterm(s) that represent an estimate for the nominal zeropoint on a clear, photometric night.
#
    aterm=numpy.zeros(8,dtype=numpy.float32)
    aterm[band2i['u']]=-25.0
    aterm[band2i['g']]=-25.428
    aterm[band2i['r']]=-25.532
    aterm[band2i['i']]=-25.413
    aterm[band2i['z']]=-25.086
    aterm[band2i['Y']]=-24.000
    aterm[band2i['VR']]=-25.47
    aterm[band2i['N964']]=-25.086
#
#   A set of kterm(s) that represent an estimate for the nominal extinction/airmass correction on a clear, photometric night.
#
    kterm=numpy.zeros(8,dtype=numpy.float32)
    kterm[band2i['u']]=0.489
    kterm[band2i['g']]=0.181
    kterm[band2i['r']]=0.095
    kterm[band2i['i']]=0.089
    kterm[band2i['z']]=0.053
    kterm[band2i['Y']]=0.05
    kterm[band2i['VR']]=0.13
    kterm[band2i['N964']]=0.053

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
#    cat_fname=[]
#    cat_fpath=[]
    cat_dict={}
    if (args.list is None):
#
#       If no list of catalogs was provided (--list) then go figure out whether they exist based on the PFW_Attempt_ID
#
        tq0=time.time()
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
            print("# Executing initial query for catalogs ")
            if (args.format_query):
                print("# query = {:s}".format(query))
            else:
                print("# query ={:s} ".format(" ".join([d.strip() for d in query.split('\n')])))

        cur.arraysize = 1000 # get 1000 at a time when fetching
        cur.execute(query)
        desc = [d[0].lower() for d in cur.description]
        for row in cur:
            rowd = dict(zip(desc, row))
            cat_fname=rowd['filename']
            cat_fpath=os.path.join(rowd['root'],rowd['path'],rowd['filename'])
            cat_dict[cat_fname]={'fpath':cat_fpath}
#            cat_fname.append(rowd['filename'])
#            cat_fpath.append(os.path.join(rowd['root'],rowd['path'],rowd['filename']))
        if (args.verbose):
            print("Timing for query execution and data handling was {:.2f} ".format(time.time()-tq0))
    else:
#
#       Open list file and proceed.
#
        if (args.verbose):
            print("# Using list of files/catalogs to obtain objects from exposure")
        if (os.path.isfile(args.list)):
            f_cat=open(args.list,'r')
            for line in f_cat:
                line=line.strip()
                columns=line.split(' ')
                if (columns[0] != "#"):
                    tmp_fname=columns[0].split('/')
                    cat_fname=tmp_fname[-1]
                    cat_dict[cat_fname]={'fpath':columns[0]}
#                    cat_fname.append(tmp_fname[-1])
#                    cat_fpath.append(columns[0])
#                   print(columns[0])
            f_cat.close()
    if (args.verbose):
        print("# Total number of catalogs/CCDs identified: {:d}".format(len(cat_dict)))
        print("# ")

##############################################################################
#   Now that a starting point has been established...
#   Obtain associated DB information for each catalog/CCD.
# 
    ccd_info={}
#    num_ccd_verbose=0
#    max_ccd_verbose=1
    if (args.verbose):
        print("# Executing query(s) to obtain image level metadata and QA values ")

    CatList=[]
    for cat_fname in cat_dict:
        CatList.append([cat_fname])
#
#   RAG: switched from GTT_FILENAME to GTT_STR for current Oracle implementation.
#        To switch back it is next line... the insert_many below and then the query 
#           needs to constrain on g.filename instead of g.str
#
#    tempTable="GTT_FILENAME"
    tempTable="GTT_STR"
    cur.execute('delete from {:s}'.format(tempTable))
    # load filenames into GTT_FILENAME table
    if (args.verbose):
        print("# Loading {:s} table for secondary queries with entries for {:d} images".format(tempTable,len(CatList)))
#    dbh.insert_many(tempTable,['FILENAME'],CatList)
    dbh.insert_many(tempTable,['STR'],CatList)

#    for index, tmp_fname in enumerate(cat_fname):
    tq0=time.time()
    query = """
        select i.filename as img_fname,
            c.filename as cat_fname,
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
        from {schema:s}image i, {schema:s}catalog c, {ttab:s} g
        where c.filename=g.str
            and c.pfw_attempt_id=i.pfw_attempt_id 
            and i.pfw_attempt_id={aid:}
            and i.filetype='red_immask' 
            and c.expnum=i.expnum 
            and c.ccdnum=i.ccdnum 
        """.format(schema=db_Schema,ttab=tempTable,aid=PFW_Attempt_ID)

    if (args.verbose):
        if (args.format_query):
            print("# query = {:s}".format(query))
        else:
            print("# query ={:s} ".format(" ".join([d.strip() for d in query.split('\n')])))

    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)
    desc = [d[0].lower() for d in cur.description]

    for row in cur:
        rowd = dict(zip(desc, row))

        ccdnum=int(rowd['ccdnum'])
        ccd_info[ccdnum]=rowd
#
#       Add entries for filename, filepath, iband 
#       Check entries that are sometimes malformed or null
#
#        ccd_info[ccdnum]['cat_fname']=rowd[at_fname[index]
        ccd_info[ccdnum]['cat_fpath']=cat_dict[rowd['cat_fname']]['fpath']
        if (ccd_info[ccdnum]['band'] is None):
            ccd_info[ccdnum]['band']='Unknown'
            ccd_info[ccdnum]['iband']=-1
        else:
            ccd_info[ccdnum]['iband']=band2i[ccd_info[ccdnum]['band']]

        if (ccd_info[ccdnum]['gaina'] is None):
            ccd_info[ccdnum]['gaina']=-1.
            print("# Warning: Null found for GAINA")

        if (ccd_info[ccdnum]['gainb'] is None):
            ccd_info[ccdnum]['gainb']=-1.
            print("# Warning: Null found for GAINB")

        if (ccd_info[ccdnum]['airmass'] is None):
            ccd_info[ccdnum]['airmass']=1.2
            print("# Warning: Null found for AIRMASS (default to 1.2)")

        if (ccd_info[ccdnum]['skybrite'] is None):
#
#           Note this can occur if the PCA Sky subrtaction is not used..
#
            ccd_info[ccdnum]['skybrite']=sbrite_good[ccd_info[ccdnum]['band']]
            ccd_info[ccdnum]['skysigma']=1.0
            print("# Warning: Null found for SKYBRITE (default SKYBRITE to {:.2f} and SKYSIGMA to 1.0)".format(sbrite_good[ccd_info[ccdnum]['band']]))

        if (ccd_info[ccdnum]['skysigma'] is None):
#
#           Note this can occur if the PCA Sky subrtaction is not used..
#
            ccd_info[ccdnum]['skysigma']=1.0
            print("# Warning: Null found for SKYSIGMA (but apparently not SKYBRITE?... default SKYSIGMA to 1.0)")

        if ((ccd_info[ccdnum]['rac1'] is None)or(ccd_info[ccdnum]['rac2'] is None)or
            (ccd_info[ccdnum]['rac3'] is None)or(ccd_info[ccdnum]['rac4'] is None)):
            print("# Warning: Null found among RAC1-4")
            ccd_info[ccdnum]['img_ra']=[]
        else:
            ccd_info[ccdnum]['img_ra']=[ccd_info[ccdnum]['rac1'],ccd_info[ccdnum]['rac2'],ccd_info[ccdnum]['rac3'],ccd_info[ccdnum]['rac4']]
#
        if ((ccd_info[ccdnum]['decc1'] is None)or(ccd_info[ccdnum]['decc2'] is None)or
            (ccd_info[ccdnum]['decc3'] is None)or(ccd_info[ccdnum]['decc4'] is None)):
            print("# Warning: Null found among DECC1-4")
            ccd_info[ccdnum]['img_dec']=[]
        else:
            ccd_info[ccdnum]['img_dec']=[ccd_info[ccdnum]['decc1'],ccd_info[ccdnum]['decc2'],ccd_info[ccdnum]['decc3'],ccd_info[ccdnum]['decc4']]
        ccd_info[ccdnum]['fwhm_old']=pixsize*float(ccd_info[ccdnum]['fwhm_old'])
#
#       sbrite_good and sbrite_lim are defined in cts/sec
#       Therefore:
#          If exposure time is present then use to normalize brightness to a per second quantity
#          If gains are present then use to express brightness in counts (rather than electrons)
#
        if (ccd_info[ccdnum]['exptime']>0.01):
            efactor=ccd_info[ccdnum]['exptime']
        else:
            efactor=1.0

        gtesta=ccd_info[ccdnum]['gaina']-1.
        gtestb=ccd_info[ccdnum]['gainb']-1.
        if ((abs(gtesta)<0.5)and(abs(gtestb)<0.5)):
#           The case where gains are 1... therefore units are electrons
            gfactor=4.0
            ccd_info[ccdnum]['bunit']='e-'
            if (args.verbose):
                print("# GAINA/B for ccdnum={:d} are consistent with units of electrons (GAINA:{:.3f},GAINB:{:.3f})".format(
                    ccdnum,ccd_info[ccdnum]['gaina'],ccd_info[ccdnum]['gainb']))
        else:
#           The case where gains are not 1... therefore units are already in counts
            gfactor=1.0
            ccd_info[ccdnum]['bunit']='DN'
            if (args.verbose):
                print("# GAINA/B for ccdnum={:d} are consistent with units of DN (GAINA:{:.3f},GAINB:{:.3f})".format(
                    ccdnum,ccd_info[ccdnum]['gaina'],ccd_info[ccdnum]['gainb']))
        ccd_info[ccdnum]['skyb']=ccd_info[ccdnum]['skybrite']/efactor/gfactor
        ccd_info[ccdnum]['skys']=ccd_info[ccdnum]['skysigma']/efactor/gfactor
    if (args.verbose):
         print("Timing for query execution and data handling was {:.2f} ".format(time.time()-tq0))

##############################################################################
#   Fill in entries for missing CCDs (if there are any...) ACTUALLY THIS HAS BEEN REMOVED, RAG THINKS COMPLETELY
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

#    print(len(exptime_chk))

##############################################################################
#   Perform the actual checks (and set exposure level values when these pass muster.
#   First check exptime.
#
    uniq_exptime_chk=list(set(exptime_chk))
    uniq_band_chk=list(set(band_chk))
    uniq_bunit_chk=list(set(bunit_chk))
    uniq_expnum_chk=list(set(expnum_chk))
    uniq_airmass_chk=list(set(airmass_chk))
#    print(len(uniq_exptime_chk))
    if (len(uniq_exptime_chk) != 1):
        if (len(uniq_exptime_chk) > 1):
            print("WARNING: Other than one exptime?: ",uniq_exptime_chk)
            print("WARNING: Using exptime: {:.2f}".format(uniq_exptime_chk[0]))
            exp_rec['exptime']=uniq_exptime_chk[0]
        else:
            print("No exptime? Setting to NoneType will try to obtain value from EXPOSURE")
            exp_rec['exptime']=None
    else:
        exp_rec['exptime']=uniq_exptime_chk[0]
#
#   Now check airmass
#
    if (len(uniq_airmass_chk) != 1):
        if (len(uniq_airmass_chk) > 1):
            print("WARNING: Other than one airmass?: ",uniq_airmass_chk)
            print("WARNING: Using airmass: {:.3f}".format(uniq_airmass_chk[0]))
            exp_rec['airmass']=uniq_airmass_chk[0]
        else:
            print("WARNING: No airmass?: ")
            print("WARNING: Using defualt setting airmass: 1.2")
            exp_rec['airmass']=1.2
    else:
        exp_rec['airmass']=uniq_airmass_chk[0]
#
#   Now check bunit for consistency
#
    if (len(uniq_bunit_chk) != 1):
        if (len(uniq_bunit_chk) > 1):
            print("WARNING: Other than one bunit?: ",uniq_bunit_chk)
            print("WARNING: Using bunit: {:s} ".format(uniq_bunit_chk[0]))
            exp_rec['bunit']=uniq_bunit_chk[0]
        else:
            print("WARNING: Assuming bunit = DN")
            exp_rec['bunit']='DN'
    else:
        exp_rec['bunit']=uniq_bunit_chk[0]
#
#   Now check band (also check that it is a sanctioned value)
#
    if (len(uniq_band_chk) != 1):
        if (len(uniq_band_chk)==0):
            print("No band? Really?  Setting to NoneType will try to obtain value from EXPOSURE")
            exp_rec['band']=None
        else:
            print("Abort: Other than one band identified?: ",uniq_band_chk)
            exit(1)
    else:
#        if (uniq_band_chk[0] in ['u','g','r','i','z','Y','VR']):
        if (uniq_band_chk[0] in band2i):
            exp_rec['band']=uniq_band_chk[0]
        else:
            print("Abort: Unsupported value for band?: ",uniq_band_chk[0])
            exit(1)
#
#   Now check expnum
#
    if (len(list(set(expnum_chk))) != 1):
        print("Abort: Other than one expnum?: ",uniq_expnum_chk)
        exit(1)
    else:
        exp_rec['expnum']=uniq_expnum_chk[0]
##############################################################################
#   Calculated exposure level values (old FWHM, old Ellipticity, Sky Brightness and Sky Sigma)
#
    if (fwhm_old.size > 0):
        exp_rec['fwhm_img']=fwhm_old.mean()
    else:
        exp_rec['fwhm_img']=-1.0
    if (elli_old.size > 0):
        exp_rec['ellip_avg']=elli_old.mean()
        exp_rec['ellip_rms']=elli_old.std()
    else:
        exp_rec['ellip_avg']=-1.0
        exp_rec['ellip_rms']=-1.0
    if (skyb.size > 0):
        exp_rec['skyb_avg']=skyb.mean()
        exp_rec['skyb_rms']=skyb.std()
    else:
        exp_rec['skyb_avg']=-1.0
        exp_rec['skyb_rms']=-1.0
    if (skys.size > 0):
        exp_rec['skys_avg']=skys.mean()
        exp_rec['skys_rms']=skys.std()
    else:
        exp_rec['skys_avg']=-1.0
        exp_rec['skys_rms']=-1.0

###############################################################################
#   Get more exposure level information
#
    tq0=time.time()
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
        print("# Executing query to obtain exposure level metadata")
        if (args.format_query):
            print("# query = {:s}".format(query))
        else:
            print("# query ={:s} ".format(" ".join([d.strip() for d in query.split('\n')])))

    cur.arraysize = 1000 # get 1000 at a time when fetching
    cur.execute(query)
    desc = [d[0].lower() for d in cur.description]

    num_exp_recs=0
    for row in cur:
        num_exp_recs=num_exp_recs+1
        if (num_exp_recs > 1):
            print("# WARNING: multiple exposure records found for this exposure? (num_exp_recs={:d})".format(num_exp_recs))
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
            if (exp_rec['band'] != 'Unknown'):
                print("WARNING: BAND miss-match between exposure-level (Unknown) and image/catalog-level ({:s}) queries.  Using image/cat-result.".format(exp_rec['band']))
        else:
            if (exp_rec['band'] is None):
                if (rowd['band'] in band2i):
                    exp_rec['band']=rowd['band']
                else:
                    print("Abort: Unsupported value for band?: ",rowd['band'])
                    exit(1)
            if (rowd['band'] != exp_rec['band']):
                print("WARNING: BAND miss-match between exposure-level ({:s}) and image/catalog-level ({:s}) queries.  Using image/cat-result.".format(rowd['band'],exp_rec['band']))
#
        if (exp_rec['exptime'] is None):
            print("Fixing exptime using value from EXPOSURE")
            exp_rec['exptime']=rowd['exptime']
        else:
            if (rowd['exptime'] != exp_rec['exptime']):
                print("WARNING: EXPTIME miss-match between exposure-level ({:.1f}) and image/catalog-level ({:.1f}) queries.  Using image/cat-result.".format(rowd['exptime'],exp_rec['exptime']))

        if (rowd['airmass'] is None):
           print("WARNING: AIRMASS miss-match between exposure-level (Unknown) and image/catalog-level ({:.3f}) queries.  Using image/cat-result.".format(exp_rec['airmass']))
        else:
            if (rowd['airmass'] != exp_rec['airmass']):
                print("WARNING: AIRMASS miss-match between exposure-level ({:.3f}) and image/catalog-level ({:.3f}) queries.  Using image/cat-result.".format(rowd['airmass'],exp_rec['airmass']))
#
#       Work out whether the exposure is part of the general survey, SN, or other.
#
        if (rowd['program']=="survey"):
            exp_rec['program']='survey'
            exp_rec['sn_field']=False
            exp_rec['survey_field']=True
        elif (rowd['program']=="supernova"):
            exp_rec['program']='SN'
            exp_rec['sn_field']=True
            exp_rec['survey_field']=False
        elif (rowd['program']=="photom-std-field"):
            exp_rec['program']='phot-std'
            exp_rec['sn_field']=False
            exp_rec['survey_field']=False
        else:
            if (exp_rec['obstype'] in ['zero','dark','dome flat','sky flat']):
                exp_rec['program']='cal'
                exp_rec['sn_field']=False
                exp_rec['survey_field']=False
            else:
                exp_rec['program']='unknown'
                exp_rec['sn_field']=False
                exp_rec['survey_field']=False
    if (args.verbose):
         print("Timing for query execution and data handling was {:.2f} ".format(time.time()-tq0))

###############################################################################
#   Query to get SCAMP_QA information
#
    tq0=time.time()
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
        print("# Executing query to obtain astrometry QA for this exposure")
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

    if (args.verbose):
         print("Timing for query execution and data handling was {:.2f} ".format(time.time()-tq0))
#
#   In case no record was found.
#
    if (num_ast_recs == 0):
        exp_rec['astrom_sig1']=None
        exp_rec['astrom_sig2']=None
        exp_rec['astrom_off1']=None
        exp_rec['astrom_off2']=None
        exp_rec['astrom_rms2']=None
        exp_rec['astrom_ndets']=None
        exp_rec['astrom_chi2']=None

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
    exp_rec['astrom_rms2']=numpy.sqrt((sigx*sigx)+(sigy*sigy))
#
#   RAG: note have lowered the threshold that looks for very small astrom_sig values to accomodate GAIA
#
    if ((exp_rec['astrom_sig1'] < 0.0001)or(exp_rec['astrom_sig2'] < 0.0001)):
        astrom_good=False
        print("# WARNING: Probable astrometric solution failure: astrom_sig1,2 = {:7.4f},{:7.4f}".format(
            exp_rec['astrom_sig1'],exp_rec['astrom_sig2']))
    if (exp_rec['astrom_rms2'] > 0.500):
        astrom_good=False
        print("# WARNING: Probable astrometric solution failure: astrom_rms2 = {:.3f}".format(
            exp_rec['astrom_rms2']))

###############################################################################
#   Query to get PSF_QA information
#
    tq0=time.time()
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
        print("# Executing query to obtain PSF QA for this exposure")
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
            print("# WARNING: multiple exposure records found for PSF QA for this exposure/pfw_attempt_id? (num_psf_recs={:d}".format(num_psf_recs))
        rowd = dict(zip(desc, row))
        exp_rec['num_psfex']=rowd['count']
        exp_rec['psfex_fwhm']=pixsize*rowd['med_fwhm']
    
    if (args.verbose):
         print("Timing for query execution and data handling was {:.2f} ".format(time.time()-tq0))
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
    if ((img_ra.size < 2)or(img_dec.size < 2)):
        print("Warning: No CCD corners present in database for these items")
        exp_rec['ra_cen']=exp_rec['tradeg']
        exp_rec['dec_cen']=exp_rec['tdecdeg']
    else:
        if ((img_ra.size < (4*exp_rec['numccd']))or(img_dec.size < (4*exp_rec['numccd']))):
            print("Warning: Some CCDs may be missing RA/DEC corners)! Working with what was present...")
        ra_min=numpy.min(img_ra)
        ra_max=numpy.max(img_ra)
        dec_min=numpy.min(img_dec)
        dec_max=numpy.max(img_dec)
        if ((ra_max-ra_min)>180.0):
#           print("IT HAPPENS {:13.7f} {:13.7f}".format(ra_min,ra_max))
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
    print("########################################")
    print("#   date_obs: {:}".format(exp_rec['date_obs']))
    print("#       Nite: {:}".format(exp_rec['nite']))
    print("#   Exposure: {:}".format(exp_rec['expnum']))
    print("#       Band: {:}".format(exp_rec['band']))
    print("#    Exptime: {:.1f} ".format(exp_rec['exptime']))
    print("#      BUNIT: {:}".format(exp_rec['bunit']))
    print("# ")
    print("#    Obstype: {:}".format(exp_rec['obstype']))
    print("#     Object: {:}".format(exp_rec['object']))
    print("#    Program: {:}".format(exp_rec['program']))
    print("# ")
    print("#      Telescope(Ra,Dec): {:9.5f} {:9.5f} ".format(exp_rec['tradeg'],exp_rec['tdecdeg']))
    print("#")
    print("# Image Centroid(Ra,Dec): {:9.5f} {:9.5f} ".format(exp_rec['ra_cen'],exp_rec['dec_cen']))
    print("#")
    print("# ")
    if (astrom_good):
        print("# Astrometric solution appears OK ")
    else:
        print("# WARNING: Probable astrometric solution failure")
    print("# Astrometry summary:")
    print("#                    high S/N      ")
    print("#        ndets: {:d} ".format(exp_rec['astrom_ndets']))
    print("#         chi2: {:.2f} ".format(exp_rec['astrom_chi2']))
#   RAG notes this is now in milli-arcseconds (assumes using GAIA)
    print("#   astrom_sig: {:7.1f},{:7.1f} ".format(exp_rec['astrom_sig1']*1000.0,exp_rec['astrom_sig2']*1000.0))
    print("#   astrom_off: {:7.1f},{:7.1f} ".format(exp_rec['astrom_off1']*1000.0,exp_rec['astrom_off2']*1000.0))
    print("#")
    print("# PSF summary:")
    print("# median(FWHM): {:.3f} ".format(exp_rec['psfex_fwhm']))
    print("#      num CCD: {:d} ".format(exp_rec['num_psfex']))
    print("#")
    print("########################################")


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
    print("# Preparing queries for NOMAD, APASS, DES(Y3A2) in vicinity of the exposure.")
#
#   Setup rules for form queries and parsing results for NOMAD/APASS/DES
#
    qparse={}
    qparse['nomad']={'pname':'NOMAD','tab':'NOMAD','tab_abbrev':'n','db':'oper','ra':'n.ra','dec':'n.dec'}
    if (exp_rec['band'] in ['u','g']):
        qparse['nomad']['mag']='n.b'
        qparse['nomad']['mlimit']=blimit
        qparse['nomad']['simple']=False
        qparse['nomad']['keys']=['ra','dec','mag']
    elif (exp_rec['band'] in ['r','VR']):
        qparse['nomad']['mag']='n.b'
        qparse['nomad']['mag2']='n.j'
        qparse['nomad']['mlimit']=blimit
        qparse['nomad']['mlimit2']=jlimit
        qparse['nomad']['simple']=False
        qparse['nomad']['keys']=['ra','dec','mag','mag2']
    elif (exp_rec['band'] in ['i','z','Y','N964']):
        qparse['nomad']['mag']='n.j'
        qparse['nomad']['mlimit']=jlimit
        qparse['nomad']['simple']=True
        qparse['nomad']['keys']=['ra','dec','mag']
    else:
        print("Band={:s} not technically supported for NOMAD calibration.  Default to use n.j".format(exp_rec['band']))
        qparse['nomad']['mag']='n.j'
        qparse['nomad']['mlimit']=jlimit
        qparse['nomad']['simple']=True
        qparse['nomad']['keys']=['ra','dec','mag']

    if (exp_rec['band'] in ['u','g','r','i','z','Y','N964']):
        qparse['apass']={'pname':'APASS','tab':'APASS_DR7','tab_abbrev':'a','db':'oper','ra':'a.ra','dec':'a.dec'}
        if (exp_rec['band'] in ['u','g']):
            qparse['apass']['mag']='a.g'
            qparse['apass']['dmag']='a.g_err'
            qparse['apass']['mlimit']=glimit
        elif (exp_rec['band'] in ['r','VR']):
            qparse['apass']['mag']='a.r'
            qparse['apass']['dmag']='a.r_err'
            qparse['apass']['mlimit']=rlimit
        elif (exp_rec['band'] in ['i','z','Y','N964']):
            qparse['apass']['mag']='a.i'
            qparse['apass']['dmag']='a.i_err'
            qparse['apass']['mlimit']=ilimit
        qparse['apass']['simple']=True
        qparse['apass']['keys']=['ra','dec','mag','dmag']
    else:
        print("Band={:s} not supported for APASS calibration.  Will skip".format(exp_rec['band']))

    if (exp_rec['band'] in ['g','r','i','z','Y']):
        qparse['des']={'pname':'Y3A2','tab':'Y3A2_COS_SUBSET','tab_abbrev':'cos','db':'oper','ra':'cos.alphawin_j2000','dec':'cos.deltawin_j2000'}
        qparse['des']['mag']='cos.wavg_mag_psf_{:s}'.format(exp_rec['band'])
        qparse['des']['dmag']="cos.wavg_magerr_psf_{:s}".format(exp_rec['band'])
        qparse['des']['spread_model']="cos.wavg_spread_model_{:s}".format(exp_rec['band'])
        qparse['des']['spreaderr_model']="cos.wavg_spreaderr_model_{:s}".format(exp_rec['band'])
        qparse['des']['mlimit']=glimit
        qparse['des']['keys']=['ra','dec','mag','dmag','spread_model','spreaderr_model']
    elif (exp_rec['band'] in ['VR']):
        qparse['des']={'pname':'Y3A2','tab':'Y3A2_COS_SUBSET','tab_abbrev':'cos','db':'oper','ra':'cos.alphawin_j2000','dec':'cos.deltawin_j2000'}
        qparse['des']['mag']='cos.wavg_mag_psf_{:s}'.format('r')
        qparse['des']['dmag']="cos.wavg_magerr_psf_{:s}".format('r')
        qparse['des']['spread_model']="cos.wavg_spread_model_{:s}".format('r')
        qparse['des']['spreaderr_model']="cos.wavg_spreaderr_model_{:s}".format('r')
        qparse['des']['mlimit']=glimit
        qparse['des']['keys']=['ra','dec','mag','dmag','spread_model','spreaderr_model']
    elif (exp_rec['band'] in ['N964']):
        qparse['des']={'pname':'Y3A2','tab':'Y3A2_COS_SUBSET','tab_abbrev':'cos','db':'oper','ra':'cos.alphawin_j2000','dec':'cos.deltawin_j2000'}
        qparse['des']['mag']='cos.wavg_mag_psf_{:s}'.format('z')
        qparse['des']['dmag']="cos.wavg_magerr_psf_{:s}".format('z')
        qparse['des']['spread_model']="cos.wavg_spread_model_{:s}".format('z')
        qparse['des']['spreaderr_model']="cos.wavg_spreaderr_model_{:s}".format('z')
        qparse['des']['mlimit']=glimit
        qparse['des']['keys']=['ra','dec','mag','dmag','spread_model','spreaderr_model']
    else:
        print("Band={:s} not supported for DES calibration.  Will skip".format(exp_rec['band']))
#
#   Form query lists (with as to facilitate roughly homogenous catalogs
#
    for cat in qparse:
        qparse[cat]['qlist']=[]
        for key in qparse[cat]['keys']:
            qparse[cat]['qlist'].append("{:s} as {:s}".format(qparse[cat][key],key))
        qparse[cat]['qselect']=",".join(qparse[cat]['qlist'])
#       print(qparse[cat]['qlist']

#
#   Form magnitude constraints for queries
#
    if ('nomad' in qparse):
        if ('mag2' in qparse['nomad']):
            qparse['nomad']['constraint']=' and {m1:s}<{ml1:.1f} and {m2:s}<{ml2:.1f}'.format(
                m1=qparse['nomad']['mag'],
                ml1=qparse['nomad']['mlimit'],
                m2=qparse['nomad']['mag2'],
                ml2=qparse['nomad']['mlimit2'])
        else:
            qparse['nomad']['constraint']=' and {m1:s}<{ml1:.1f}'.format(
                m1=qparse['nomad']['mag'],
                ml1=qparse['nomad']['mlimit'])
    if ('apass' in qparse):
        qparse['apass']['constraint']=' and {m1:s}<{ml1:.1f}'.format(
            m1=qparse['apass']['mag'],
            ml1=qparse['apass']['mlimit'])
    if ('des' in qparse):
        if (exp_rec['band']=='N964'):
            qparse['des']['constraint']=' and {m1:s}<{ml1:.1f} and cos.nepochs_{bval:s}>0'.format(
                m1=qparse['des']['mag'],
                ml1=qparse['des']['mlimit'],
                bval='z')
        elif (exp_rec['band']=='VR'):
            qparse['des']['constraint']=' and {m1:s}<{ml1:.1f} and cos.nepochs_{bval:s}>0'.format(
                m1=qparse['des']['mag'],
                ml1=qparse['des']['mlimit'],
                bval='r')
        else:
            qparse['des']['constraint']=' and {m1:s}<{ml1:.1f} and cos.nepochs_{bval:s}>0'.format(
                m1=qparse['des']['mag'],
                ml1=qparse['des']['mlimit'],
                bval=exp_rec['band'])

    if (args.verbose):
        for cat in qparse:
            print("# Will apply CONSTRAINT ({const:s}) for work with {bval:s}-band {cname:s} data".format(
                const=qparse[cat]['constraint'],
                bval=exp_rec['band'],
                cname=qparse[cat]['tab']))
#
#   Form queries for APASS/NOMAD/DES objects 
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

        for cat in qparse:
            qparse[cat]['query']="""
                select {qlist:s} from {tname:s} {tabbrev:s}
                where {tabbrev:s}.dec between {d1:.7f} and {d2:.7f} 
                    and ({tabbrev:s}.ra < {r1:.7f} or {tabbrev:s}.ra > {r2:.7f})
                    {mcon:s}""".format(
                qlist=qparse[cat]['qselect'],
                tname=qparse[cat]['tab'],tabbrev=qparse[cat]['tab_abbrev'],
                d1=dec1,d2=dec2,r1=ra1,r2=ra2,
                mcon=qparse[cat]['constraint'])

    else:
#
#       NOTE need something better if we get closer than 2 degrees from pole
#
        cosdec_cen=numpy.cos(dec_cen*deg2rad)
        ra1=ra_cen-(fp_rad/cosdec_cen)
        ra2=ra_cen+(fp_rad/cosdec_cen)
        dec1=dec_cen-fp_rad
        dec2=dec_cen+fp_rad

        for cat in qparse:
            qparse[cat]['query']="""
                select {qlist:s} from {tname:s} {tabbrev:s}
                where {tabbrev:s}.dec between {d1:.7f} and {d2:.7f} 
                    and {tabbrev:s}.ra between {r1:.7f} and {r2:.7f}
                    {mcon:s}""".format(
                qlist=qparse[cat]['qselect'],
                tname=qparse[cat]['tab'],tabbrev=qparse[cat]['tab_abbrev'],
                d1=dec1,d2=dec2,r1=ra1,r2=ra2,
                mcon=qparse[cat]['constraint'])


    t2=time.time()
    print("# TIMING (pre-NOMAD/APASS): {:.2f}".format((t2-t1)))
#
#   Use formed queries to obtain catalog data
#
    cat_match={}
    if args.verbose:
        print("# Executing queries for NOMAD,APASS,DES in vicinity of the exposure.")
    for cat in qparse:
        t30=time.time()
        if (args.verbose):
            print("########################################")
            print("# query = {:s}".format(qparse[cat]['query']))

        if (qparse[cat]['db'] == "oper"):
            cur=dbh.cursor()
        else:
            cur=dbhsci.cursor()

        prefetch=100000
        cur.arraysize=int(prefetch)
        cur.execute(qparse[cat]['query'])
        header=[d[0].lower() for d in cur.description]
        cat_data=pd.DataFrame(cur.fetchall())

        cat_match[cat]={}
        cat_match[cat]['cat']={}
        if (cat_data.empty):
            print("# No values returned from query of {tval:s} ".format(tval=qparse[cat]['tab']))
            for val in header:
                cat_match[cat]['cat'][val]=numpy.array([])
        else:
            cat_data.columns=header
            for val in header:
                cat_match[cat]['cat'][val]=numpy.array(cat_data[val])
        cur.close()
        cat_match[cat]['timing']=time.time()-t30
        print("# TIMING ({tval:s} query): {texec:.2f}".format(tval=qparse[cat]['tab'],texec=cat_match[cat]['timing']))
        print("# Number of {tval:s} entries found is {nval:d} ".format(tval=qparse[cat]['tab'],nval=cat_match[cat]['cat']['ra'].size))
#
#   Specialized updates to catalog data 
#
    if ('nomad' in cat_match):
        if (qparse['nomad']['simple']):
            wsm=numpy.where(cat_match['nomad']['cat']['mag'] > qparse['nomad']['mlimit'])
            cat_match['nomad']['cat']['mag'][wsm]=99.0 
        else:
            if (exp_rec['band'] in ['u','g']):
                bmag=cat_match['nomad']['cat']['mag']+nomad_bmag_corr[(cat_match['nomad']['cat']['dec']+90.0).astype(int)]
                wsm=numpy.where(numpy.logical_or((bmag>qparse['nomad']['mlimit']),numpy.logical_and((bmag>18.37),(bmag<18.39))))
                bmag[wsm]=99.0
                cat_match['nomad']['cat']['mag']=bmag
            elif (exp_rec['band'] in ['r','VR']):
                bmag=cat_match['nomad']['cat']['mag']+nomad_bmag_corr[(cat_match['nomad']['cat']['dec']+90.0).astype(int)]
                jmag=cat_match['nomad']['cat']['mag2']
                rmag=0.333333333*((2.0*bmag)+jmag)
                wsm=numpy.where(numpy.logical_or((bmag>qparse['nomad']['mlimit']),numpy.logical_and((bmag>18.37),(bmag<18.39))))
                rmag[wsm]=99.0
                wsm=numpy.where(jmag>qparse['nomad']['mlimit2'])
                rmag[wsm]=99.0
                cat_match['nomad']['cat']['mag']=rmag
#
#   Values that have magnitude <98 should be kept (others should be removed)...
#
    for cat in cat_match:
        wsm=numpy.where(cat_match[cat]['cat']['mag']<98.)
        for key in cat_match[cat]['cat']:
            cat_match[cat]['cat'][key]=cat_match[cat]['cat'][key][wsm]           

    t3=time.time()
    print("# TIMING (catalog NOMAD/APASS/DES acquisition): {:.2f}".format((t3-t2)))
    print("########################################")
    print("# Starting to obtain exposure catalog data")

###############################################################################
#   Time to read in the finalcut_cat (was red_cat) products.
#   Setup to read catalogs
#
    mtime=2.5*numpy.log10(exp_rec['exptime'])
    aval=aterm[band2i[exp_rec['band']]]
    if (exp_rec['bunit']=="e-"):
        bval=-2.5*numpy.log10(4.0)
    else:   
        bval=-2.5*numpy.log10(1.0)
#    print("RAG: bval=",bval)
#
    mcorr=mtime-aval-25.0-bval

    nstar_found=0

###############################################################################
#   Columns to be retrieved from FITS tables.
    cols_retrieve=['ALPHAWIN_J2000','DELTAWIN_J2000','MAG_AUTO','MAGERR_AUTO','MAG_PSF','MAGERR_PSF',
                    'SPREAD_MODEL','FWHM_WORLD','A_IMAGE','B_IMAGE','FLUX_RADIUS','KRON_RADIUS','CLASS_STAR','FLAGS','IMAFLAGS_ISO']
    cols_datacheck={'ALPHAWIN_J2000':{'min':-0.1,'max':360.1},
                    'DELTAWIN_J2000':{'min':-90.0,'max':90.0},
                    'MAG_PSF':{'min':1.0,'max':50.0},
                    'MAG_AUTO':{'min':1.0,'max':50.0},
                    'FWHM_WORLD':{'min':(0.1/3600.),'max':(20.0/3600.)},
                    'A_IMAGE':{'min':(0.1/pixsize),'max':(10.0/pixsize)},
                    'B_IMAGE':{'min':(0.1/pixsize),'max':(10.0/pixsize)},
                    'FLUX_RADIUS':{'min':(0.1/pixsize),'max':(5.0/pixsize)},
                    'KRON_RADIUS':{'min':(0.1/pixsize),'max':(5.0/pixsize)}
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

    exp_cat={}
    for iccd in range(1,63):
        if (iccd in ccd_info):
            if (os.path.isfile(ccd_info[iccd]['cat_fpath'])):
                CAT = fitsio.read(ccd_info[iccd]['cat_fpath'],ext='LDAC_OBJECTS',columns=cols_retrieve)
                if ('ALPHAWIN_J2000' in exp_cat):
                    for col in cols_retrieve:
                        newcol=numpy.concatenate((exp_cat[col],CAT[col]),axis=0)    
                        exp_cat[col]=newcol
                    tmp_ccdnum=numpy.zeros(CAT['ALPHAWIN_J2000'].size,dtype=numpy.int16)
                    tmp_ccdnum+=iccd
                    newcol=numpy.concatenate((exp_cat['CCDNUM'],tmp_ccdnum),axis=0)
                    exp_cat['CCDNUM']=newcol
                else:
                    for col in cols_retrieve:
                        exp_cat[col]=CAT[col]
                    exp_cat['CCDNUM']=numpy.zeros(CAT['ALPHAWIN_J2000'].size,dtype=numpy.int16)
                    exp_cat['CCDNUM']+=iccd
#
#   Perform checking on the resulting catalogs to remove errant values.
#
    if (args.verbose):
        print("# Found {:d} entries.".format(exp_cat['ALPHAWIN_J2000'].size))
    for col in cols_datacheck:
        if (col in exp_cat):
            wsm=numpy.where(numpy.logical_and((exp_cat[col]>cols_datacheck[col]['min']),(exp_cat[col]<cols_datacheck[col]['max'])))
            for acol in exp_cat:
                exp_cat[acol]=exp_cat[acol][wsm]
            if (args.verbose):
                print("#       {:d} entries remain after range check on {:s} [{:f}:{:f}]".format(
                    exp_cat['ALPHAWIN_J2000'].size,col,cols_datacheck[col]['min'],cols_datacheck[col]['max']))

#   Re-order sorted on magnitude (allows for speed up when comparing with APASS/NOMAD below)
#    new_expcat=sorted(exp_cat,key=lambda k: k['MAG_'+use_mag_type])

    if (args.verbose):
        print("# Total number of cat_finalcut objects: {:d}".format(exp_cat['ALPHAWIN_J2000'].size))
###################################################################
#
#   Form some simple sets for analysis
#

#    mag=numpy.array([new_expcat[k]['MAG_'+use_mag_type]+mcorr for k, object in enumerate(new_expcat)])
#    magerr=numpy.array([new_expcat[k]['MAGERR_'+use_mag_type]    for k, object in enumerate(new_expcat)])
#
#    rstruct=[{'a_image':obj['A_IMAGE'],'b_image':obj['B_IMAGE'],'fwld':obj['FWHM_WORLD'],
#              'frad':obj['FLUX_RADIUS'],'krad':obj['KRON_RADIUS']} 
#              for obj in new_expcat if ((abs(obj['SPREAD_MODEL'])<0.0015)and(obj['MAGERR_'+use_mag_type]<0.1))]
#    rrad1=numpy.array([2.*obj['a_image'] for obj in rstruct])
#    rrad2=numpy.array([obj['a_image']+obj['b_image'] for obj in rstruct])
#    fwld=numpy.array([obj['fwld'] for obj in rstruct])
#    frad=numpy.array([obj['frad'] for obj in rstruct])
#    krad=numpy.array([obj['krad'] for obj in rstruct])

    mag=exp_cat['MAG_'+use_mag_type]+mcorr 
    magerr=exp_cat['MAGERR_'+use_mag_type]
#    starcut=numpy.where(numpy.logical_and((abs(exp_cat['SPREAD_MODEL'])<0.0015),(exp_cat['MAGERR_'+use_mag_type]<0.1)))
#    rrad1=2.0*exp_cat['A_IMAGE'][starcut]
#    rrad2=exp_cat['A_IMAGE'][starcut]+exp_cat['B_IMAGE'][starcut]
#    fwld=exp_cat['FWHM_WORLD'][starcut]
#    frad=exp_cat['FLUX_RADIUS'][starcut]
#    krad=exp_cat['KRON_RADIUS'][starcut]

########################################################################################
#   Proceed to assessing the atmospheric opacity.
#
#    cat_use={"des":['g','r','i','z','Y'],
#             "apass":['g','r'],
#             "nomad":['u','g','r','i','z','Y','VR']}
    transparency_calc=False
    if ('apass' in cat_match):
        cat_match['apass']['banduse']=['g','r','i','N964']
    if ('nomad' in cat_match):
        cat_match['nomad']['banduse']=['u','g','r','i','z','Y','VR','N964']
    if ('des' in cat_match):
        cat_match['des']['banduse']=['g','r','i','z','Y','N964']

#
#   New copy of exp_cat so that can make cuts to remove non-stellar objects
#
    expCat={}
    if (exp_cat['ALPHAWIN_J2000'].size < 1):
        expCat['RA']=numpy.array([])
        expCat['DEC']=numpy.array([])
        expCat['MAG']=numpy.array([])
        expCat['DMAG']=numpy.array([])
        expCat['SPREAD_MODEL']=numpy.array([])
        expCat['IMAFLAGS_ISO']=numpy.array([])
        expCat['CCDNUM']=numpy.array([])
        expCat['A_IMAGE']=numpy.array([])
        expCat['B_IMAGE']=numpy.array([])
        expCat['FWHM_WORLD']=numpy.array([])
        expCat['FLUX_RADIUS']=numpy.array([])
        expCat['KRON_RADIUS']=numpy.array([])
    else:
        expCat['RA']=exp_cat['ALPHAWIN_J2000']
        expCat['DEC']=exp_cat['DELTAWIN_J2000']
        expCat['MAG']=exp_cat['MAG_'+use_mag_type]+mcorr 
        expCat['DMAG']=exp_cat['MAGERR_'+use_mag_type]
        expCat['SPREAD_MODEL']=exp_cat['SPREAD_MODEL']
        expCat['IMAFLAGS_ISO']=exp_cat['IMAFLAGS_ISO']
        expCat['CCDNUM']=exp_cat['CCDNUM']
        expCat['A_IMAGE']=exp_cat['A_IMAGE']
        expCat['B_IMAGE']=exp_cat['B_IMAGE']
        expCat['FWHM_WORLD']=exp_cat['FWHM_WORLD']
        expCat['FLUX_RADIUS']=exp_cat['FLUX_RADIUS']
        expCat['KRON_RADIUS']=exp_cat['KRON_RADIUS']

#
#       Reject objects with flags set (IMAFLAGS_ISO)
#       In particular this takes the place of a previous analysis to search for and remove saturated objects
#
#       define BADPIX_BPM          1  
#       define BADPIX_SATURATE     2  
#       define BADPIX_INTERP       4
#       define BADPIX_BADAMP       8  
#       define BADPIX_CRAY        16
#       define BADPIX_STAR        32  
#       define BADPIX_TRAIL       64  
#       define BADPIX_EDGEBLEED  128  
#       define BADPIX_SSXTALK    256  
#       define BADPIX_EDGE       512  
#       define BADPIX_STREAK    1024  
#       define BADPIX_SUSPECT   2048 
#       define BADPIX_FIXED     4096  
#       define BADPIX_NEAREDGE  8192  
#       define BADPIX_TAPEBUMP 16384  
#
        badpix_val=numpy.uint(1+2+8+16+32+64+128+256+512+1024)
        ibit=numpy.bitwise_and(expCat['IMAFLAGS_ISO'],badpix_val)
        wsm=numpy.where(numpy.logical_not(ibit))
        expCatCut={}
        for key in expCat:
            expCatCut[key]=expCat[key][wsm]
        expCat=expCatCut
        print("# After cutting on IMAFLAGS_ISO {:d} objects remain in the exposure catalog".format(expCat['RA'].size))
#
#       Select for stars (using abs(SPREAD_MODEL) < 0.0015) and also cut to objects with S/N>10
#
        wsm=numpy.where(numpy.logical_and((abs(expCat['SPREAD_MODEL'])<0.0015),(expCat['DMAG']<0.1)))
        expCatCut={}
        for key in expCat:
            expCatCut[key]=expCat[key][wsm]
        expCat=expCatCut
    print("# After cutting on SPREAD_MODEL and MAGERR {:d} objects remain in the exposure catalog".format(expCat['RA'].size))
    exp_rec['numobj']=expCat['RA'].size
    exp_c2=SkyCoord(ra=expCat['RA']*u.degree,dec=expCat['DEC']*u.degree)

#
#   Match each catalog to the exposure catalog to get a median extinction
#

    for cat in ['des','apass','nomad']:
        t400=time.time()
        print("########################################")
        doNewComp=True
        if (cat not in cat_match):
            doNewComp=False
        else:
            if (exp_rec['band'] not in cat_match[cat]['banduse']):
                doNewComp=False
            if (cat_match[cat]['cat']['mag'].size < 1):
                doNewComp=False
        if (not(doNewComp)):
            print("# No entries found for cat={:s}.  Skipping".format(cat))
            mag_exp=numpy.array([])
            mag_cat=numpy.array([])
            mag_diff=numpy.array([])
            mindes=10.
            maxdes=18.

        if (doNewComp):
            print("# Attempting matching to cat={:s}".format(cat))
        
            sky_c1=SkyCoord(ra=cat_match[cat]['cat']['ra']*u.degree,dec=cat_match[cat]['cat']['dec']*u.degree)
            if ((sky_c1.shape[0]<1)or(exp_c2.shape[0]<1)):
                if (sky_c1.shape[0]<1):
                    print("# Warning CAT: {:s} had {:d} entries".format(cat,sky_c1.shape[0]))
                if (exp_c2.shape[0]<1):
                    print("# Warning Exposure catalog had {:d} entries".format(exp_c2.shape[0]))
                mag_exp=numpy.array([])
                mag_cat=numpy.array([])
                mag_diff=numpy.array([])
            else:
                idx2, d2d, d3d = sky_c1.match_to_catalog_sky(exp_c2)
                idx1=numpy.arange(cat_match[cat]['cat']['ra'].size)
                wsm=numpy.where(d2d.arcsecond<1.0)
                if (args.verbose):
                    print("# Matched {:d} objects between exposure catalog and {:s}".format(idx1[wsm].size,cat))
                if (idx1[wsm].size > 0):
                    mag_ccd=expCat['CCDNUM'][idx2[wsm]]
                    mag_exp=expCat['MAG'][idx2[wsm]]
                    mag_cat=cat_match[cat]['cat']['mag'][idx1[wsm]]+cat_mag_corr[cat][exp_rec['band']]
                    mag_diff=mag_exp-mag_cat
                else:
                    mag_exp=numpy.array([])
                    mag_cat=numpy.array([])
                    mag_diff=numpy.array([])
#
#           Setting limits for the plots later
#
            if (mag_exp.size < 1):
                mindes=10.
                maxdes=18.
            else:
                mindes=numpy.amin(mag_exp)
                maxdes=numpy.amax(mag_exp)
#
#               Outlier rejection
#
                sig_reject=4.0
                med_magdiff=numpy.median(mag_diff)
                std_magdiff=numpy.std(mag_diff)
                for iter_reject in [0,1,2]:
                    mag_sig_outlier=numpy.absolute(mag_diff-med_magdiff)/std_magdiff
                    wsm=numpy.where( mag_sig_outlier > sig_reject)
                    xsm=numpy.where(numpy.logical_not(mag_sig_outlier > sig_reject))
                    if (args.verbose):
                        print("# Outlier rejection iteration {:2d} found Median(mag_diff)={:.3f} w/ Stddev(mag_diff)={:.3f} rejecting {:d} of {:d}".format(
                            iter_reject,med_magdiff,std_magdiff,mag_exp[wsm].size,mag_exp.size))
                    med_magdiff=numpy.median(mag_diff[xsm])
                    std_magdiff=numpy.std(mag_diff[xsm])

                wsm=numpy.where( mag_sig_outlier > sig_reject)
                xsm=numpy.where(numpy.logical_not(mag_sig_outlier > sig_reject))
                if (args.verbose):
                    print("# Outlier rejection iteration {:2d} found Median(mag_diff)={:.3f} w/ Stddev(mag_diff)={:.3f} rejecting {:d} of {:d}".format(
                            iter_reject+1,med_magdiff,std_magdiff,mag_exp[wsm].size,mag_exp.size))

                mag_exp_keep=mag_exp[xsm]
                mag_cat_keep=mag_cat[xsm]
                mag_diff_keep=mag_diff[xsm]
            
                mag_exp_reject=mag_exp[wsm] 
                mag_cat_reject=mag_cat[wsm] 
                mag_diff_reject=mag_diff[wsm] 

#                if (args.verbose):
#                    print("#  ")
#                    print("#   CCD     keep     reject ")
#                    for iccd in range(1,63):
#                        print("# CCD{:02d}:  {:7d}  {:7d} ".format(
#                                iccd,mag_ccd[numpy.where(mag_ccd[xsm]==iccd)].size,mag_ccd[numpy.where(mag_ccd[wsm]==iccd)].size))
#                    print("#  ")
#
#               Go on and calculate magnitude(s) of opacity
#
                exp_rec[cat+"_magdiff"]=med_magdiff
                exp_rec[cat+"_kmagdiff"]=med_magdiff-cat_mag_corr[cat][exp_rec['band']]+cat_kmag_corr[cat][exp_rec['band']]-(kterm[band2i[exp_rec['band']]]*exp_rec['airmass'])
                exp_rec[cat+"_num"]=mag_diff.size
#
#               LOGIC HERE DECIDES RULES FOR WHICH MEASURE FROM WHICH CATALOG TO USE
#                   -- order they go through is important (APASS, then NOMAD cleans up)
#                   -- once a choice is made it simple reports the value (with NOUSE)
#
                if(not(transparency_calc)):
                    if (exp_rec['band'] in cat_match[cat]['banduse']):
                        transparency_calc=True
                        exp_rec['magdiff']=med_magdiff
                        exp_rec['cloud_cat']=cat.upper()
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
            PlotData['nmatch']=mag_diff.size
            PlotData['band']=exp_rec['band']
            if (cat in qparse):
                PlotData['cat']=qparse[cat]['pname']
            else:
                PlotData['cat']=cat
            PlotData['magtype']=use_mag_type
            if (mag_diff.size < 100):
#               If QA plots have been requested then output an empty QA plot that gives 
#               information about the failure.
                if (mag_diff.size > 0):
                    PlotData['mag_des']=mag_exp
                    PlotData['mag_cat']=mag_cat
                    PlotData['mag_diff']=mag_diff
                    PlotData['mag_des_reject']=numpy.array([])
                    PlotData['mag_cat_reject']=numpy.array([])
                    PlotData['mag_diff_reject']=numpy.array([])
                else:
                    PlotData['mag_des']=numpy.array([])
                    PlotData['mag_cat']=numpy.array([])
                    PlotData['mag_diff']=numpy.array([])
                    PlotData['mag_des_reject']=numpy.array([])
                    PlotData['mag_cat_reject']=numpy.array([])
                    PlotData['mag_diff_reject']=numpy.array([])
            else:
                PlotData['mag_des']=mag_exp_keep
                PlotData['mag_cat']=mag_cat_keep
                PlotData['mag_diff']=mag_diff_keep
                PlotData['mag_des_reject']=mag_exp_reject
                PlotData['mag_cat_reject']=mag_cat_reject
                PlotData['mag_diff_reject']=mag_diff_reject
                PlotData['mindes']=mindes
                PlotData['maxdes']=maxdes
                PlotData['med_magdiff']=med_magdiff

            SEplot.QAplt_cloud(PlotFileName,PlotData,astrom_good,args.verbose)
#
#       Finished... get the timing information
#
        if (cat+"_magdiff" not in exp_rec):
            exp_rec[cat+"_magdiff"]=99.
            exp_rec[cat+"_kmagdiff"]=99.
            exp_rec[cat+"_num"]=-1
        if (cat in cat_match):
            cat_match[cat]['crosstiming']=time.time()-t400
#
#   If all else fails... meaning no transparency calculation succeeded then there is always  
#   the fallback position.
#
    if (not(transparency_calc)):
        exp_rec['magdiff']=99.
        exp_rec['cloud_cat']="Failed"

#   We now return to the old code...
###################################################################

    t5=time.time()
##############################################################################
#   Load all the various FWHM measuments into the dictionary
#
    if (expCat['A_IMAGE'].size > 2):                   
        exp_rec['nstar_found']=expCat['A_IMAGE'].size
        exp_rec['fwhm_rrad1']=pixsize*2.0*expCat['A_IMAGE'].mean()
        rrad2=expCat['A_IMAGE']+expCat['B_IMAGE']
        exp_rec['fwhm_rrad2']=pixsize*(expCat['A_IMAGE']+expCat['B_IMAGE']).mean()
        exp_rec['fwhm_world']=3600.*expCat['FWHM_WORLD'].mean()
        exp_rec['fwhm_krad']=2.0*pixsize*expCat['KRON_RADIUS'].mean()
        exp_rec['fwhm_frad']=2.0*pixsize*expCat['FLUX_RADIUS'].mean()
    else:
        exp_rec['nstar_found']=expCat['A_IMAGE'].size
        exp_rec['fwhm_rrad1']=-1.
        exp_rec['fwhm_rrad2']=-1.
        exp_rec['fwhm_world']=-1.
        exp_rec['fwhm_krad']=-1.
        exp_rec['fwhm_frad']=-1.
    t5a=time.time()

#########################################################################################
#   In the following section calculate some QA on the magnitude depth (needs work to make it 
#   robust (or remotely useful).
#
    magerr_hist=numpy.zeros((mbin.size),dtype=numpy.float32)
    magnum_hist=numpy.zeros((mbin.size),dtype=numpy.float32)

    if (mag.size<10):
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
            PlotData['nobj']=mag.size
            PlotData['band']=exp_rec['band']
            PlotData['magtype']=use_mag_type
            SEplot.QAplt_maghist(PlotFileName,PlotData,astrom_good,args.verbose)

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
            while((emag_sort[icnt]<mbin[ibin])and(icnt<(emag_sort.size-1))):
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
            PlotData['nobj']=mag.size
            PlotData['bins']=mbin-(0.5*magbin_step)
            PlotData['mag_hist']=magnum_hist
            PlotData['magerr_hist']=magerr_hist
            PlotData['mag_raw']=mag
            PlotData['magerr_raw']=magerr
            PlotData['band']=exp_rec['band']
            PlotData['magtype']=use_mag_type
            SEplot.QAplt_maghist(PlotFileName,PlotData,astrom_good,args.verbose)

    exp_rec['mag_thresh']=mag_thresh
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
        if (use_fwhm < 0.3):
            print("# WARNING:  FWHM is unbelievably small ({:.3f} arcsec or {:.2f} pixels)!  F_EFF set to -1.0 arcsec".format(use_fwhm,use_fwhm/pixsize))
            exp_rec['teff_f']=-1.0
        else:
            exp_rec['teff_f']=(seeing_fid[exp_rec['band']]*seeing_fid[exp_rec['band']]/(use_fwhm*use_fwhm))
    else:
        print("# WARNING:  No FWHM measure available. F_EFF set to -1.0")
        exp_rec['teff_f']=-1.0
#
#   Calculate B_eff
#
    if (exp_rec['skyb_avg']>0.0):
        exp_rec['teff_b']=sbrite_good[exp_rec['band']]/exp_rec['skyb_avg']
    else:
        print("# WARNING:  No SKY BRIGHTNESS measure available. B_EFF set to -1.0")
        exp_rec['teff_b']=-1.0
#
#   Calculate C_eff
#
    if ((exp_rec['magdiff']>-95.)and(exp_rec['magdiff']<95.0)):
        if (exp_rec['magdiff']<0.2):
            exp_rec['teff_c']=1.0
        else:
            exp_rec['teff_c']=math.pow(10.0,(-2.0*(exp_rec['magdiff']-0.2)/2.5))
    else:
        print("# WARNING:  No CLOUD measure available. C_EFF set to -1.0")
        exp_rec['teff_c']=-1.0
#
#   Calculate T_eff
#
    value_teff=1.0
    if (exp_rec['teff_f']>=0):
        value_teff=value_teff*exp_rec['teff_f']
    if (exp_rec['teff_b']>=0):
        value_teff=value_teff*exp_rec['teff_b']
    if (exp_rec['teff_c']>=0):
        value_teff=value_teff*exp_rec['teff_c']
    if ((exp_rec['teff_f']<0)or(exp_rec['teff_b']<0)):
        exp_rec['teff']=-1.
    else:
        exp_rec['teff']=value_teff
    t7=time.time()

#############################################################################
#   Summary of all timing information
#   Summary of all FWHM measurements

    print("# ")
    print("# ")
    print("#         TIMING Summary  ")
    print("#------------------------------")
    print("#  general query(s): {:7.2f} ".format((t1-t0)))
    if ('apass' in cat_match): print("#       query APASS: {:7.2f} ".format(cat_match['apass']['timing']))
    if ('nomad' in cat_match): print("#       query NOMAD: {:7.2f} ".format(cat_match['nomad']['timing']))
    if ('des'   in cat_match): print("#         query DES: {:7.2f} ".format(cat_match['des']['timing']))
    if ('apass' in cat_match): print("#  cross-corr APASS: {:7.2f} ".format(cat_match['apass']['crosstiming']))
    if ('nomad' in cat_match): print("#  cross-corr NOMAD: {:7.2f} ".format(cat_match['nomad']['crosstiming']))
    if ('des'   in cat_match): print("#    cross-corr DES: {:7.2f} ".format(cat_match['des']['crosstiming']))
#    print("#    Cat comparison: {:7.2f} ".format((t5-t4)))
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
    if ((exp_rec['teff']<0.)or(exp_rec['teff_c']<0.0)):
        new_decide="unkn"
    elif (exp_rec['teff']>teff_lim[exp_rec['band']]):
        new_decide="good"
#        if (exp_rec['fwhm_world']>seeing_lim[exp_rec['band']]):
        if (use_fwhm>seeing_lim[exp_rec['band']]):
            print(" FWHM used in assessment = {fval:6.3f}... exceeds seeing limit for band ({bval:s}) of {limval:6.3f} ".format(
                fval=use_fwhm,bval=exp_rec['band'],limval=seeing_lim[exp_rec['band']]))
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
        out_row.append(exp_rec['expnum'])
        out_row.append(dm_process)
        out_row.append(dm_accept)
        out_row.append("no comment")
        writer.writerow(out_row)  
#
#   Write DB_file if command line option present 
#
#    if((args.updateDB)or(args.DB_file)):
    prog_name='Unknown'
    if (exp_rec['sn_field']):
        prog_name='supernova'
    if (exp_rec['survey_field']):
        prog_name='survey'
#
#       First the columns that should always be present
#
    db_cname="INSERT INTO {dbtab:s}(EXPOSURENAME,EXPNUM,PFW_ATTEMPT_ID,CAMSYM,PROCESSED,ACCEPTED,PROGRAM,ANALYST,LASTCHANGED_TIME,T_EFF,F_EFF,B_EFF,C_EFF,FWHM_ASEC,ELLIPTICITY,SKYBRIGHTNESS".format(dbtab=db_table)
    db_value=""" VALUES ('{expname:s}', {expnum:d}, {pfw_att_id:s}, '{camsym:s}', '{flag_proc:s}', '{flag_accept:s}', '{prog:s}', '{analyst:s}', sysdate, {teff:.3f}, {feff:.3f}, {beff:.3f}, {ceff:.3f}, {fwhm:.3f}, {ellip:.3f}, {skyb:.2f}""".format(
        expname=exp_rec['exp_fname'],
        expnum=exp_rec['expnum'], 
        pfw_att_id=args.attid, 
        camsym=args.CamSym, 
        flag_proc=dm_process,
        flag_accept=dm_accept,
        prog=prog_name,
        analyst=analyst,
        teff=exp_rec['teff'],
        feff=exp_rec['teff_f'],
        beff=exp_rec['teff_b'],
        ceff=exp_rec['teff_c'],
        fwhm=exp_rec['psfex_fwhm'],
        ellip=exp_rec['ellip_avg'],
        skyb=exp_rec['skyb_avg'])

#
#       Now the columns that may or may not be present
#
    if ((exp_rec['apass_magdiff']>-95.)and(exp_rec['apass_magdiff']<95.)):
        db_cname=db_cname+','+'CLOUD_APASS'
        db_value=db_value+','+"%.3f" % exp_rec['apass_magdiff']

    if ((exp_rec['nomad_magdiff']>-95.)and(exp_rec['nomad_magdiff']<95.)):
        db_cname=db_cname+','+'CLOUD_NOMAD'
        db_value=db_value+','+"%.3f" % exp_rec['nomad_magdiff']

    if ((exp_rec['des_magdiff']>-95.)and(exp_rec['des_magdiff']<95.)):
        db_cname=db_cname+','+'CLOUD_DES'
        db_value=db_value+','+"%.3f" % exp_rec['des_magdiff']

    db_cname=db_cname+','+'CLOUD_CATALOG'
    db_value=db_value+','+"'%s'" % exp_rec['cloud_cat']

    if ((exp_rec['apass_kmagdiff']>-95.)and(exp_rec['apass_kmagdiff']<95.)):
        db_cname=db_cname+','+'KLOUD_APASS'
        db_value=db_value+','+"%.3f" % exp_rec['apass_kmagdiff']

    if ((exp_rec['nomad_kmagdiff']>-95.)and(exp_rec['nomad_kmagdiff']<95.)):
        db_cname=db_cname+','+'KLOUD_NOMAD'
        db_value=db_value+','+"%.3f" % exp_rec['nomad_kmagdiff']

    if ((exp_rec['des_kmagdiff']>-95.)and(exp_rec['des_kmagdiff']<95.)):
        db_cname=db_cname+','+'KLOUD_DES'
        db_value=db_value+','+"%.3f" % exp_rec['des_kmagdiff']

    if (exp_rec['apass_num']>-1):
        db_cname=db_cname+','+'N_APASS'
        db_value=db_value+','+"%d" % exp_rec['apass_num']

    if (exp_rec['nomad_num']>-1):
        db_cname=db_cname+','+'N_NOMAD'
        db_value=db_value+','+"%d" % exp_rec['nomad_num']

    if (exp_rec['des_num']>-1):
        db_cname=db_cname+','+'N_DES'
        db_value=db_value+','+"%d" % exp_rec['des_num']

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
        cur = dbh.cursor()
        cur.execute(insert_command)
    if (args.DB_file):
        fdbout.write("{:s};\n".format(insert_command))

#
#   Write the summary of the assessment to STDOU
#
    ftxt.write("#                                                                                                         Astrometry \n")     
    ftxt.write("# Exposure                                     FWHM                             #        AIR    EXP  SCMP  HighSN     Depth       APASS           NOMAD            DES       FWHM            \n")
    ftxt.write("#  Num   STATE    t_eff F_eff  B_eff   C_eff   PSFex  Ellip    Sky_B   N_src   CCD BAND  MASS   TIME Flag sig1  sig2     Assess     dmag    #       dmag    #     dmag    #     Wld    Object   \n")
    ftxt.write("#                                               [\"]           [DN/s]                             [s]      [\"]    [\"]     [mag]    [mag]            [mag]            [mag]         \n")
#
    ftxt.write(" {expnum:9d} {decide:4s}  {teff:6.2f} {feff:6.2f} {beff:6.2f} {ceff:6.2f}   {fwhm:5.2f} {ellip:6.3f} {skyb:8.2f} {numobj:7d}   {numccd:3d} {band:1s} {iband:1d} {airmass:6.3f} {exptime:6.1f} {dummy:2d} {asig1:6.3f} {asig2:6.3f}  {mdiff:8.3f} {amdiff:8.3f} {anum:6d} {nmdiff:8.3f} {nnum:6d} {dmdiff:8.3f} {dnum:6d} {fwhm_wld:6.3f} {object:s} \n".format(
        expnum=exp_rec['expnum'],
        decide=new_decide,
        teff=exp_rec['teff'],
        feff=exp_rec['teff_f'],
        beff=exp_rec['teff_b'],
        ceff=exp_rec['teff_c'],
        fwhm=exp_rec['psfex_fwhm'],
        ellip=exp_rec['ellip_avg'],
        skyb=exp_rec['skyb_avg'],
        numobj=exp_rec['numobj'],
        numccd=exp_rec['numccd'],
        band=exp_rec['band'],
        iband=band2i[exp_rec['band']],
        airmass=exp_rec['airmass'],
        exptime=exp_rec['exptime'],
        dummy=99,
        asig1=exp_rec['astrom_sig1'],
        asig2=exp_rec['astrom_sig2'],
        mdiff=exp_rec['magdiff'],
        amdiff=exp_rec['apass_magdiff'],
        anum=exp_rec['apass_num'],
        nmdiff=exp_rec['nomad_magdiff'],
        nnum=exp_rec['nomad_num'],
        dmdiff=exp_rec['des_magdiff'],
        dnum=exp_rec['des_num'],
        fwhm_wld=exp_rec['fwhm_world'],
        object=exp_rec['object']))

    if(args.updateDB):
        dbh.commit()
        print("DB update complete and committed")
    if (args.csv):
        fcsv.close()
    if (args.DB_file):
        fdbout.close()

