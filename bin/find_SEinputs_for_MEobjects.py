#! /usr/bin/env python3
# $Id$
# $Rev$:  # Revision of last commit.
# $LastChangedBy$:  # Author of last commit.

"""
Link COADD objects to their associated single-epoch images.
"""
from __future__ import print_function

import argparse
import os
import despydb.desdbi
import time
import sys
#import esutil
import numpy as np
import pandas as pd
#import fgcm_y3a1_tools
#import fitsio
import qatoolkit.mangle_systematics as ms
    
############################################################

if __name__ == "__main__":

    t00=time.time()   
    parser = argparse.ArgumentParser(description='Example code that finds correspondence between COADD OBJECTS and their underlying single-epoch data products (through Mangle)') 

    parser.add_argument('-t','--tile',    action='store', type=str, required=True, help='Tile to work on')
    parser.add_argument('-p','--proctag', action='store', type=str, required=True, help='ProcTag Name')
    parser.add_argument('-r','--release', action='store', type=str, default='None', help='Prefix of release that identify tables being used (default="None" which uses an empty string)')
    parser.add_argument('--dbTable',      action='store', type=str, default='coadd_object', help='DB table with objects')
    parser.add_argument('-s','--section', action='store', type=str, default=None,   help='section of .desservices file with connection info')
    parser.add_argument('-S','--Schema',  action='store', type=str, default=None,   help='DB schema (do not include \'.\').')

    parser.add_argument('-T','--Timing',  action='store_true', default=False, help='If set timing information accompanies output')
    parser.add_argument('--debug'       , action='store_true', default=False, help='Debug mode resticts code to work on a handful of objects')
    parser.add_argument('-v','--verbose', action='store', type=int, default=0, help='Verbosity (defualt:0; currently values up to 2)')

    args = parser.parse_args()
    if (args.verbose > 0):
        print("Args: {:}".format(args))

#    pi=3.141592654
#    halfpi=pi/2.0
#    deg2rad=pi/180.0
#    depth=[1.0,2.0,3.0,4.0,5.0,6.0,7.0,10.0,15.,20.,40.,80.,160.]

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

    BandList=['g','r','i','z','Y']
#   Setup a standard proxy index for each band (in molyArray)
    BDict={}
    for i in range(len(BandList)):
        BDict[BandList[i]]=i

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
    attemptID=ms.get_tile_attempt(args.tile,args.proctag,dbh,dbSchema,releasePrefix,Timing=args.Timing,verbose=verbose)
    print(attemptID)
    if (attemptID is None):
        print("Failed to identify an attempt for TILE={tile:s} for PROCTAG={ptag:s}.".format(tile=args.tile,ptag=args.proctag))
        print("Aborting")
        exit(1)

    molyArray,IdArray = ms.get_moly_array(attemptID,BDict,args.dbTable,dbh,dbSchema,Timing=args.Timing,verbose=verbose)
#    for i, Id in enumerate(IdArray):
#        print(" {:d} {:d} {:d} {:d} {:d} {:d} {:d}".format(i,Id,molyArray[0,i],molyArray[1,i],molyArray[2,i],molyArray[3,i],molyArray[4,i]))

    molygonDict,ccdgonDict = ms.get_CCDGON_Dict(molyArray,BDict,dbh,dbSchema,releasePrefix,Timing=args.Timing,verbose=verbose)

#    print(ccdgonDict)

    collateDict={
        'id':'i8',
        'band':'i2',
        'expnum':'i4',
        'ccdnum':'i2',
        'inverse_variance_weight':'f8'}

    SystematicDict = ms.collate_object_systematics(collateDict,IdArray,molyArray,BDict,molygonDict,ccdgonDict,Timing=args.Timing,verbose=verbose)



    t11=time.time()
    print("Total execution time was {:.2f} seconds".format(t11-t00))

    exit(0)
