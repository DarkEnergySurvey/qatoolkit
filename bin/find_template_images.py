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

    for image in imageDict:
        print("/archive_data/desarchive/{:s}/{:s}{:s}".format(imageDict[image]['path'],imageDict[image]['filename'],imageDict[image]['compression']))

    t11=time.time()
    print("Total execution time was {:.2f} seconds".format(t11-t00))

    exit(0)
