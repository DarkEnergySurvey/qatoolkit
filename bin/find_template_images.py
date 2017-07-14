#! /usr/bin/env python
# $Id$
# $Rev$:  # Revision of last commit.
# $LastChangedBy$:  # Author of last commit.

"""
Find images that overlap another image.
"""
from __future__ import print_function

import argparse
import os
import despydb.desdbi
import time
import sys
#import esutil
import numpy as np
#import pandas as pd
#import fgcm_y3a1_tools
#import fitsio
import qatoolkit.template as tm
    
############################################################

if __name__ == "__main__":

    t00=time.time()   
    parser = argparse.ArgumentParser(description='Example code that finds correspondence between COADD OBJECTS and their underlying single-epoch data products (through Mangle)') 

#    parser.add_argument('-i','--image',    action='store', type=str, required=True, help='Tile to work on')
    parser.add_argument('-p','--proctag', action='store', type=str, required=True, help='ProcTag Name')
    parser.add_argument('-r','--radec',   action='store', type=str, required=True, help='Ra,Dec coordinate (comma separated)')
    parser.add_argument('-b','--band',    action='store', type=str, default=None,  help='Band constraint (default no constraint)')
    parser.add_argument('-f','--frac',    action='store', type=str, default="0.75", help='Fractional overlap (default=0.75) will take comma separated list to have frac RA,Dec')
    parser.add_argument('--release',      action='store', type=str, default='Y3A2', help='Prefix specifying a set of release table (Use "None" when working with PROD in DESOPER)')
    parser.add_argument('-s','--section', action='store', type=str, default=None,   help='section of .desservices file with connection info')
    parser.add_argument('-S','--Schema',  action='store', type=str, default=None,   help='DB schema (do not include \'.\').')
    parser.add_argument('--use_blacklist', action='store_true', default=False, help='Flag to constrain results to not appear in the blacklist')
    parser.add_argument('--blacklist',     action='store', type=str, default='Y3A2_BLACKLIST', help='Over-ride name of blacklist table to use (default=Y3A2_BLACKLIST)')
    parser.add_argument('--use_eval',     action='store_true', default=False, help='Flag to constrain results to match DES survey quality cuts')
    parser.add_argument('--evaltable',     action='store', type=str, default='Y3A2_QA_SUMMARY', help='Over-ride name of evaluation table to use (default=Y3A2_QA_SUMMARY)')
    parser.add_argument('--use_zpt',      action='store_true', default=False, help='Flag to constrain results to have a zeropoint available')
    parser.add_argument('--zpttable',     action='store', type=str, default='Y3A2_ZEROPOINT,FGCM,v2.0,16', help='Over-ride default zeropoint criteia (default=Y3A2_ZEROPOINT,FGCM,v2.0,16)')
    parser.add_argument('-T','--Timing',  action='store_true', default=False, help='If set timing information accompanies output')
    parser.add_argument('--debug'       , action='store_true', default=False, help='Debug mode resticts code to work on a handful of objects')
    parser.add_argument('-v','--verbose', action='store', type=int, default=0, help='Verbosity (defualt:0; currently values up to 2)')

    args = parser.parse_args()
    if (args.verbose > 0):
        print("Args: {:s}".format(args))

##########################################################
#   Handle simple args (verbose, Schema, bandlist)
#
    verbose=args.verbose

    if (args.Schema is None):
        dbSchema=""
    else:
        dbSchema="{:s}.".format(args.Schema)

#
#   Special Case where args.release = "None"
#
    if (args.release == "None"):
        releasePrefix=""
    else:
        releasePrefix="{:s}_".format(args.release)

    radec=np.array([float(i) for i in args.radec.split(',')])
    if ("," in args.frac):
        frac=np.array([float(i) for i in args.frac.split(',')])
    else:
        frac=np.array([float(args.frac),float(args.frac)])

#
#   Parse the zpttable into a dictionary
#
    zptDict={}
    zptlist=args.zpttable.split(',')
    zptDict['table']=zptlist[0]
    zptDict['source']=zptlist[1]
    zptDict['version']=zptlist[2]
    zptDict['flag']=zptlist[3]

#
#   Eval Dict
#
    evalDict={}
    evalDict['table']=args.evaltable
    evalDict['teff_cut']={'g':0.2,'r':0.3,'i':0.3,'z':0.3,'Y':0.2}
#   Flat seeing cut
#    evalDict['fwhm_cut']={'g':2.0,'r':2.0,'i':2.0,'z':2.0,'Y':2.0}
#   Kolmogorov cut with FWHM_ASEC in legacy versions of FIRST/FINALCUT_EVAL
#    evalDict['fwhm_cut']={'g':1.7648,'r':1.6656,'i':1.6,'z':1.544,'Y':1.7648,'u':1.7648,'VR':1.7648}
#   Kolmogorov cut with PSF_FWHM used in QA_SUMMARY and later versions of FIRST/FINALCUT_EVAL
    evalDict['fwhm_cut']={'g':1.7248,'r':1.6256,'i':1.56,'z':1.504,'Y':1.7248,'u':1.7248,'VR':1.7648}

#    print(radec)
#    print(frac)

########################################################
#
#   Setup a DB connection
#
    try:
        desdmfile = os.environ["des_services"]
    except KeyError:
        desdmfile = None
    dbh = despydb.desdbi.DesDbi(desdmfile,args.section)

    print("Using tag_constraint of {:s}".format(args.proctag))
#
    imageDict={}
    imageDict=tm.get_image_list(imageDict,radec,frac,args.band,args.proctag,releasePrefix,dbh,dbSchema,Timing=args.Timing,verbose=verbose)

    if (args.use_blacklist):
        imageDict=tm.check_blacklist(imageDict,args.blacklist,releasePrefix,dbh,dbSchema,Timing=args.Timing,verbose=verbose)

    if (args.use_zpt):
        imageDict=tm.check_zeropoint(imageDict,zptDict,releasePrefix,dbh,dbSchema,Timing=args.Timing,verbose=verbose)

    if (args.use_eval):
        imageDict=tm.check_eval(imageDict,evalDict,args.band,releasePrefix,dbh,dbSchema,Timing=args.Timing,verbose=verbose)


    for image in imageDict:
        if (args.use_eval):
            print("/archive_data/desarchive/{path:s}/{fname:s}{compress:s} {t_eff:6.3f} {fwhm:6.3f}".format(
                path=imageDict[image]['path'],
                fname=imageDict[image]['filename'],
                compress=imageDict[image]['compression'],
                t_eff=imageDict[image]['t_eff'],
                fwhm=imageDict[image]['fwhm']
            ))
        else:
            print("/archive_data/desarchive/{:s}/{:s}{:s}".format(imageDict[image]['path'],imageDict[image]['filename'],imageDict[image]['compression']))

    t11=time.time()
    print("Total execution time was {:.2f} seconds".format(t11-t00))

    exit(0)
