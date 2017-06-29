#! /usr/bin/env python
# $Id$
# $Rev$:  # Revision of last commit.
# $LastChangedBy$:  # Author of last commit.

"""
Routines for identifying images (that overlap another) and obtain QA information.
"""
from __future__ import print_function

import despydb.desdbi
import time
import numpy as np
#import pandas as pd
    
######################################################################################
def get_image_list(imageDict,RaDec,frac,band,ProcTag,releasePrefix,dbh,dbSchema,Timing=False,verbose=0):

    """ Query code to obtain list of images that overlap another

        Inputs:
            imageDict  Dictionary that holds results (passed in so that results from multiple queries can be used.
            RaDec:     Right Ascension, Declination for search
            ProcTag:   Proctag name containing set to be worked on
            releasePrefix: Prefix string (including _'s) to identify a specific set of tables
            dbh:       Database connection to be used
            dbSchema:  Schema over which queries will occur.
            verbose:   Integer setting level of verbosity when running.

        Returns:
            ImageDict: Resulting Image dictionary
    """

    t0=time.time()

    pixsize=0.263/3600.0
    imsize=np.array([4096,2048])
    fsize=frac*pixsize*imsize*np.array([(1.0/np.cos(RaDec[1]*np.pi/180.0)),1.0])
    radec_min=RaDec-fsize
    radec_max=RaDec+fsize
    print(radec_min[0],radec_max[0])
    print(radec_min[1],radec_max[1])

    if ((radec_min[0]<0.0)or(radec_max[0]>360.0)):
        crossRA0=True
        if (radec_min[0]<0.0):
            radec_min[0]=radec_min[0]+360.0
        if (radec_max[0]>360.0):
            radec_max[0]=radec_max[0]-360.0
    else:
        crossRA0=False


    if (band is None):
        bandConstraint=""
    else:
        bandConstraint="and band='{:s}'".format(band)

    if (crossRA0):
        query="""SELECT fai.path,i.filename,fai.compression
            FROM {schema:s}{rpref:s}image i,{schema:s}{rpref:s}proctag t, {schema:s}{rpref:s}file_archive_info fai
            WHERE t.tag='{ptag:s}'
                and t.pfw_attempt_id=i.pfw_attempt_id
                and i.filetype='red_immask'
                and ((i.ra_cent between {ra1:.7f} and 360.0)or(i.ra_cent between 0.0 and {ra2:.7f}))
                and i.dec_cent between {dec1:.7f} and {dec2:.7f} 
                {Bconst:s}
                and i.filename=fai.filename
            """.format(
            schema=dbSchema,ptag=ProcTag,rpref=releasePrefix,Bconst=bandConstraint,
            ra1=radec_min[0],ra2=radec_max[0],dec1=radec_min[1],dec2=radec_max[1])

    else:
        query="""SELECT fai.path,i.filename,fai.compression 
            FROM {schema:s}{rpref:s}image i,{schema:s}{rpref:s}proctag t, {schema:s}{rpref:s}file_archive_info fai 
            WHERE t.tag='{ptag:s}'
                and t.pfw_attempt_id=i.pfw_attempt_id
                and i.filetype='red_immask'
                and i.ra_cent between {ra1:.7f} and {ra2:.7f} 
                and i.dec_cent between {dec1:.7f} and {dec2:.7f} 
                {Bconst:s}
                and i.filename=fai.filename 
            """.format(
            schema=dbSchema,ptag=ProcTag,rpref=releasePrefix,Bconst=bandConstraint,
            ra1=radec_min[0],ra2=radec_max[0],dec1=radec_min[1],dec2=radec_max[1])

    if (verbose > 0):
        if (verbose == 1):
            QueryLines=query.split('\n')
            QueryOneLine='sql = '
            for line in QueryLines:
                QueryOneLine=QueryOneLine+" "+line.strip()
            print("{:s}".format(QueryOneLine))
        if (verbose > 1):
            print("{:s}".format(query))
#
#   Establish a DB connection
#
    curDB = dbh.cursor()
    curDB.execute(query)
    desc = [d[0].lower() for d in curDB.description]

    for row in curDB:
        rowd = dict(zip(desc, row))
        fname=rowd['filename']
        if (fname in imageDict):
            print("Found more than one record for image={:s}. Skipping".format(fname))
        else:
            imageDict[fname]=rowd

    if (verbose>0):
        print(" Query found {:d} images.".format(len(imageDict)))
    if (Timing):
        t1=time.time()
        print(" Query execution time: {:.2f}".format(t1-t0))

    curDB.close()

    return imageDict

