#! /usr/bin/env python3
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
        print(" Initial search query found {:d} images.".format(len(imageDict)))
    if (Timing):
        t1=time.time()
        print(" Initial search query execution time: {:.2f}".format(t1-t0))

    curDB.close()

    return imageDict


######################################################################################
def check_excludelist(imageDict,excludelistTable,releasePrefix,dbh,dbSchema,Timing=False,verbose=0):

    """ Query code to check whether images have been included in the exclude_list.

        Inputs:
            imageDict  Dictionary that holds previous results.
            excludelistTable: Name of excludelist table
            releasePrefix: Prefix string (including _'s) to identify a specific set of tables
            dbh:       Database connection to be used
            dbSchema:  Schema over which queries will occur.
            verbose:   Integer setting level of verbosity when running.

        Returns:
            ImageDict: Image dictionary (with excludelist images removed)
    """

    t0=time.time()

    ImgList=[]
    for img in imageDict:
        ImgList.append([img])

    # Make sure the GTT_STR table is empty
    tempTable="GTT_STR"
#    tempTable="gruendl.my_tmp_filename"

    curDB = dbh.cursor()
    curDB.execute('delete from {:s}'.format(tempTable))
    # load filenames into GTT_STR table
    if (verbose > 0):
        print("# Loading {:s} table for secondary queries with entries for {:d} images".format(tempTable,len(ImgList)))
    dbh.insert_many(tempTable,['STR'],ImgList)

#
#   The main query
#
    query="""SELECT g.str as filename
            FROM {schema:s}{rpref:s}image i, {ttab:s} g
            WHERE g.str=i.filename
                and not exists (select 1 from {schema:s}{elist:s} b where b.expnum=i.expnum and b.ccdnum=i.ccdnum)
            """.format(
            schema=dbSchema, rpref=releasePrefix, elist=excludelistTable, ttab=tempTable)

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
    curDB.execute(query)
    desc = [d[0].lower() for d in curDB.description]

#
#   gather the results
#
    keepFileList=[]
    for row in curDB:
        rowd = dict(zip(desc, row))
        keepFileList.append(rowd['filename'])

#
#   Output fname that will be removed.
#
    if (verbose > 2):
        print("#------------------------------------------")  
        if (len(imageDict)>len(keepFileList)):
            for fname in imageDict:
                if (fname not in keepFileList):
                    print("#  Exclude list detected.  Will remove image: {:s}".format(fname))
        else:
            print("#  No images present in the Exclude list")  


#
#   Construct new image dictionary with only entries that made it through the check.
#
    newImageDict={}
    for fname in keepFileList:
        if (fname in imageDict):
            newImageDict[fname]=imageDict[fname]

    if (verbose>0):
        print("# Exclude ist query results removed {:d} images.  {:d} images remain in dict".format((len(imageDict)-len(newImageDict)),len(newImageDict)))
    if (Timing):
        t1=time.time()
        print("# Exclude list query execution time: {:.2f}".format(t1-t0))
    print("#------------------------------------------")  

    curDB.close()

    return newImageDict


######################################################################################
def check_zeropoint(imageDict,zptDict,releasePrefix,dbh,dbSchema,Timing=False,verbose=0):

    """ Query code to check whether images have been included in the excludelist.

        Inputs:
            imageDict  Dictionary that holds previous results.
            zpttable:  Zeropoint [table,source,version,flag]
            dbh:       Database connection to be used
            dbSchema:  Schema over which queries will occur.
            verbose:   Integer setting level of verbosity when running.

        Returns:
            ImageDict: Image dictionary (with excluded images removed)
    """

    t0=time.time()

    ImgList=[]
    for img in imageDict:
        ImgList.append([img])

    # Make sure the GTT_STR table is empty
    tempTable="GTT_STR"
#    tempTable="gruendl.my_tmp_filename"

    curDB = dbh.cursor()
    curDB.execute('delete from {:s}'.format(tempTable))
    # load filenames into GTT_STR table
    if (verbose > 0):
        print("# Loading {:s} table for secondary queries with entries for {:d} images".format(tempTable,len(ImgList)))
    dbh.insert_many(tempTable,['STR'],ImgList)

#
#   The main query
#
    query="""SELECT g.str as filename, z.mag_zero, z.sigma_mag_zero
            FROM {ttab:s} g, {schema:s}{zptab:s} z 
            WHERE g.str=z.imagename
                and z.source='{zsrc:s}' 
                and z.version='{zver:s}' 
                and z.flag<{zflag:s} 
            """.format(
            schema=dbSchema, zptab=zptDict['table'], zsrc=zptDict['source'], zver=zptDict['version'], zflag=zptDict['flag'], ttab=tempTable)

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
    curDB.execute(query)
    desc = [d[0].lower() for d in curDB.description]

#
#   gather the results
#
    keepFileList=[]
    for row in curDB:
        rowd = dict(zip(desc, row))
        fname=rowd['filename']
        if (fname in imageDict):
            imageDict[fname]['mag_zero']=rowd['mag_zero']
            imageDict[fname]['sigma_mag_zero']=rowd['sigma_mag_zero']
        keepFileList.append(fname)
#
#   Output fname that will be removed.
#
    if (verbose > 2):
        print("#------------------------------------------")  
        if (len(imageDict)>len(keepFileList)):
            for fname in imageDict:
                if (fname not in keepFileList):
                     print("#  Zeropoint not found.  Will remove image: {:s}".format(fname))
        else:
            print("#  All images have a zeropoint")  

#
#   Construct new image dictionary with only entries that made it through the check.
#
    newImageDict={}
    for fname in keepFileList:
        if (fname in imageDict):
            newImageDict[fname]=imageDict[fname]

    if (verbose>0):
        print("# Zeropoint query results removed {:d} images.  {:d} images remain in dict".format((len(imageDict)-len(newImageDict)),len(newImageDict)))
    if (Timing):
        t1=time.time()
        print("# Zeropoint query execution time: {:.2f}".format(t1-t0))
    print("#------------------------------------------")  

    curDB.close()

    return newImageDict


######################################################################################
def check_eval(imageDict,evalDict,band,releasePrefix,dbh,dbSchema,Timing=False,verbose=0):

    """ Query code to check whether images have been excluded.

        Inputs:
            imageDict: Dictionary that holds previous results.
            evalDict:  Dictionary configuration for query [table]
            dbh:       Database connection to be used
            dbSchema:  Schema over which queries will occur.
            verbose:   Integer setting level of verbosity when running.

        Returns:
            ImageDict: Image dictionary (with exclude images removed)
    """

    t0=time.time()

    ImgList=[]
    for img in imageDict:
        ImgList.append([img])

    # Make sure the GTT_STR table is empty
    tempTable="GTT_STR"
#    tempTable="gruendl.my_tmp_filename"

    curDB = dbh.cursor()
    curDB.execute('delete from {:s}'.format(tempTable))
    # load filenames into GTT_STR table
    if (verbose > 0):
        print("# Loading {:s} table for secondary queries with entries for {:d} images".format(tempTable,len(ImgList)))
    dbh.insert_many(tempTable,['STR'],ImgList)

#
#   The main query
#   RAG: Note that work would need to be done here to facilitate using tables other than QA_SUMMARY
#
    query="""SELECT g.str as filename, 
                q.psf_fwhm as fwhm, 
                CASE WHEN q.t_eff < 0.0 THEN q.f_eff*q.b_eff ELSE q.t_eff END as t_eff
            FROM {ttab:s} g, {schema:s}{rpref:s}image i, {schema:s}{etab:s} q 
            WHERE g.str=i.filename
                and i.pfw_attempt_id=q.pfw_attempt_id
            """.format(
            schema=dbSchema, etab=evalDict['table'], ttab=tempTable, rpref=releasePrefix)

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
    curDB.execute(query)
    desc = [d[0].lower() for d in curDB.description]

#
#   gather the results
#
    keepFileList=[]
    for row in curDB:
        rowd = dict(zip(desc, row))
        fname=rowd['filename']
        if (fname in imageDict):
            imageDict[fname]['fwhm']=rowd['fwhm']
            imageDict[fname]['t_eff']=rowd['t_eff']
        keepFileList.append(fname)
#
#   Output fname that will be removed.
#
    if (verbose > 2):
        print("#------------------------------------------")  
        if (len(imageDict)>len(keepFileList)):
            for fname in imageDict:
                if (fname not in keepFileList):
                     print("#  Evaluation (QA_SUMMARY) not found.  Will remove image: {:s}".format(fname))
        else:
            print("#  All images have an evaluation")  

#
#   Construct new image dictionary with only entries that made it through the check.
#
    size_imageDict=len(imageDict)
    newImageDict={}
    for fname in keepFileList:
        if (fname in imageDict):
            newImageDict[fname]=imageDict[fname]

#
#   Now check remaining entries to see whether they meet quality criteria.
#

    imageDict=newImageDict
    keepFileList=[]
    for fname in imageDict:
        keepFile=True
        if (imageDict[fname]['t_eff']<evalDict['teff_cut'][band]):
            keepFile=False
            print("#  T_EFF= {:6.3f}.  Will remove image: {:s}".format(imageDict[fname]['t_eff'],fname))
        if (imageDict[fname]['fwhm']>evalDict['fwhm_cut'][band]):
            keepFile=False
            print("#   FWHM= {:6.3f}.  Will remove image: {:s}".format(imageDict[fname]['fwhm'],fname))
        if (keepFile):
            keepFileList.append(fname)
    
    newImageDict={}
    for fname in keepFileList:
        if (fname in imageDict):
            newImageDict[fname]=imageDict[fname]

    if (verbose>0):
        print("# QA query results removed {:d} images.  {:d} images remain in dict".format((len(imageDict)-len(newImageDict)),len(newImageDict)))
    if (Timing):
        t1=time.time()
        print("# QA query execution time: {:.2f}".format(t1-t0))
    print("#------------------------------------------")  

    curDB.close()

    return newImageDict


######################################################################################







