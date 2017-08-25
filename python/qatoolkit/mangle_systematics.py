#! /usr/bin/env python
# $Id: mangle_systematics.py 44620 2016-11-17 17:22:35Z rgruendl $
# $Rev: 44620 $:  # Revision of last commit.
# $LastChangedBy$:  # Author of last commit.

"""
Utilities to tie COADD object IDs to their associated single-epoch images and 
systematics measurements.
"""
from __future__ import print_function

import despydb.desdbi
import time
import numpy as np
import pandas as pd
    
######################################################################################
def get_tile_attempt(TileName,ProcTag,dbh,dbSchema,releasePrefix,Timing=False,verbose=0):
    """ Query code to obtain COADD tile PFW_ATTEMPT_ID after constraining
        that results are part of a specific PROCTAG.

        Inputs:
            TileName:  Tilename to be search for
            ProcTag:   Proctag name containing set to be worked on
            dbh:       Database connection to be used
            dbSchema:  Schema over which queries will occur.
            releasePrefix: Prefix string (including _'s) to identify a specific set of tables
            verbose:   Integer setting level of verbosity when running.

        Returns:
            AttemptID: Resulting AttemptID
    """

    t0=time.time()
    query="""SELECT
            distinct t.pfw_attempt_id as pfw_attempt_id
        FROM {schema:s}{rpref:s}proctag t, {schema:s}{rpref:s}catalog c
        WHERE t.tag='{ptag:s}' 
            and t.pfw_attempt_id=c.pfw_attempt_id
            and c.filetype='coadd_cat'
            and c.tilename='{tname:s}'
        """.format(
            schema=dbSchema,ptag=ProcTag,tname=TileName,rpref=releasePrefix)

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

    attval=None
    for row in curDB:
        rowd = dict(zip(desc, row))
        if (attval is None):
            attval=rowd['pfw_attempt_id']
        else:
            print("Found more than one attempt for tile={:s} attval={:ld} vs {:ld} ".format(TileName,attval,rowd['pfw_attempt_id']))
            attval=rowd['pfw_attempt_id']

    t1=time.time()
    if (Timing):
        print(" Query to find attempt execution time: {:.2f}".format(t1-t0))
    curDB.close()

    return attval


########################################################
def get_moly_array(attemptId,BDict,tableName,dbh,dbSchema,Timing=False,debug=False,verbose=0,returnRadec=False):
    """
    Query to obtain a set of COADD_OBJECT_IDs and their associated MOLYGON_NUMBERs

    attemptID:  Is the attempt ID for the run that produced the COADD_OBJECTs
    BDict:      a dictionary that contains bands that will be queried (keys) and 
                then the values for their associated  column in molyArray
                    Format is {band1:col1,band2:col2}
    tableName:  DB table from which to draw the coadd objects.
                    At a minimum it must have COADD_OBJECT_ID, ALPHAWIN_J2000, DELTAWIN_J2000, and MOLY_NUMBER_{band}
    dbh:        DB connection
    dbSchema:   DB schema to use
    Timing:     Flag to provide timing information
    debug:      Flag that will return only the first 50 entries for debugging purposes.
    verbose:    integer that controls the level of verbosity

    Output:
        molyArray:  a numpy array (Band, Object)
    """

    prefetch=100000

    debug_constraint=''
    if (debug):
        debug_constraint='and rownum<51'

    molyBand=""
    for band in BDict:
        molyBand=molyBand+",o.moly_number_"+band

    t0=time.time()
    query="""SELECT
        o.coadd_object_id,
        o.ALPHAWIN_J2000 as ra,
        o.DELTAWIN_J2000 as dec
        {mband:s}
        FROM {schema:s}{DBtab:s} o
        WHERE o.pfw_attempt_id={AttID:d} {dconst:s}
        """.format(mband=molyBand,schema=dbSchema,DBtab=tableName,AttID=attemptId,dconst=debug_constraint)

#        FROM {DBtab:s} o, {schema}y3a1_catalog c
#            and c.filetype='coadd_det_cat'
#            and c.filename=o.filename {dconst:s}

#   Show the query as it executes (if requested).
    if (verbose > 0):
        if (verbose == 1):
            QueryLines=query.split('\n')
            QueryOneLine='sql = '
            for line in QueryLines:
                QueryOneLine=QueryOneLine+" "+line.strip()
            print("{:s}".format(QueryOneLine))
        if (verbose > 1):
            print("{:s}".format(query))

#   Establish a DB connection... execute query... assemble the results
    curDB = dbh.cursor()
    curDB.arraysize=int(prefetch)
    curDB.execute(query)
    curDB.arraysize=int(prefetch)
    header=[columns[0] for columns in curDB.description]
    data=pd.DataFrame(curDB.fetchall())
    data.columns=header

    if (Timing):
        t1=time.time()
        print(" Query for objects found {:d} objects. Execution time: {:.2f}".format(len(data),t1-t0))
    curDB.close()

#   now collate into fgcm_y3a1_tools format
    ids=np.array(data['COADD_OBJECT_ID'],dtype=np.int64)
    ras=np.array(data['RA'],dtype=np.float64)
    decs=np.array(data['DEC'],dtype=np.float64)
#    print(max(data['COADD_OBJECT_ID']))
#    exit(0)
#   how many objects are there?
    u,uind=np.unique(ids,return_index=True)
   
    molyArray=np.zeros((len(BDict),u.size),dtype=np.int64)
    for band in BDict:
        molyArray[BDict[band]]=np.nan_to_num(data['MOLY_NUMBER_'+band.upper()])

    if (Timing):
        t2=time.time()
        print(" Formed final Moly_number arrays with {:d} objects. Execution time: {:.2f}".format(molyArray.shape[1],t2-t1))

    if returnRadec:
        return molyArray,ids,ras,decs
    else:
        return molyArray,ids


########################################################
def get_CCDGON_Dict(molyArray,BDict,dbh,dbSchema,releasePrefix,Timing=False,verbose=0):
    """
    Query to take a list of coadd object IDs and obtain dicts that hold:
        1) the CCDGONs associated with a MOLYGON,
        2) the metadata for each CCDGON

    molyArray:  a numpy array (Band, Object)
    BDict:      a dictionary that gives column in molyArray correspondence to a band
                Format is {band1:col1,band2:col2}
    dbh:        DB connection
    dbSchema:   DB schema to use
    releasePrefix: Prefix string (including _'s) to identify a specific set of tables
    Timing:     Flag to provide timing information
    verbose:    integer that controls the level of verbosity

    Output:
        ccdgonDict:     Dict with metadata from each CCDGON (key is CCDGON_NUMBER)
        molygonDict:    Dict containing lists of CCDGON_NUMBER that correspond to each MOLYGON_NUMBER (key)
   
    An ideal upgrade would have the metadata (second query) be configurable based on a Dict.
    """

    t0=time.time()
#   If you do not have access to a global temp table you would need to create a suitbale table
#   in your environment (e.g. "CREATE TABLE GRUENDL.MY_TMP_ID (ID NUMBER(25));"
#   And then replace the variable below with its name (not you also need to look below
#   and then uncomment a line to make commits into this table (which is not the case for 
#   GTT-type tables)
#   tempTable="gruendl.my_tmp_id"
    tempTable="GTT_ID"
#   Setup cursor for DB retrieval
    curDB=dbh.cursor()

#   Create empty output Dictionaries
    molygonDict={}
    ccdgonDict={}

#   Loop over molyArray (one column at a time)

    for band in BDict:
        IdList=[]
        bandind=BDict[band]
        # get the unique set of molygon_number's 
        w=np.unique(molyArray[bandind],return_index=False)
        for i in w:
            # Check to make sure a the MOLYNUMBER is sane
            # This is necessary because null result when an object has no data at a specific band.
            #    these in turn become NaN (I think this is because Pandas (are bad pandas))
            #    that is in turn caught and they are turned into zero's
            if (i>0):
                IdList.append([i])

        # Make sure the GTT_ID table is empty
        curDB.execute('delete from {:s}'.format(tempTable))
        # load molygon ids into gtt_id table
        print("# Loading {:s} table for secondary queries with entries for {:d} attempt IDs".format(tempTable,len(IdList)))
        dbh.insert_many(tempTable,['ID'],IdList)
#       This commit would need to be added is using a temporary table in your personal space
#        dbh.commit()

#
#       The first work horse query.  Obtain lists of CCDGONs and their association to MOLYGONs.
#
        query="""SELECT
            mc.moly_number,
            mc.ccdgon_number
        FROM {schema:s}{rpref:s}molygon_ccdgon mc, {tT:s} g
        WHERE g.id=mc.moly_number
        """.format(schema=dbSchema,rpref=releasePrefix,tT=tempTable)

#       Show the query as it executes (if requested).
        if (verbose > 0):
            if (verbose == 1):
                QueryLines=query.split('\n')
                QueryOneLine='sql = '
                for line in QueryLines:
                    QueryOneLine=QueryOneLine+" "+line.strip()
                print("{:s}".format(QueryOneLine))
            if (verbose > 1):
                print("{:s}".format(query))

        curDB.execute(query)
        desc = [d[0].lower() for d in curDB.description]

        ccdgonList=[]
        for row in curDB:
            rowd = dict(zip(desc, row))
            mnum=rowd['moly_number']
            if (mnum not in molygonDict):
                molygonDict[mnum]=[]
            molygonDict[mnum].append(rowd['ccdgon_number'])
            ccdgonList.append(rowd['ccdgon_number'])

#
#       Now form a unique set of CCDGONs, ingest IDs into temp table and then proceed with second query
#
        uccdgonList=list(set(ccdgonList))
        print(len(ccdgonList),len(uccdgonList))

        IdList=[]
        for i in uccdgonList:
            IdList.append([i])

        # Make sure the GTT_ID table is empty
        curDB.execute('delete from {:s}'.format(tempTable))
        # load molygon ids into gtt_id table
        print("# Loading {:s} table for secondary queries with entries for {:d} attempt IDs".format(tempTable,len(IdList)))
        dbh.insert_many(tempTable,['ID'],IdList)
#       This commit would need to be added is using a temporary table in your personal space
#        dbh.commit()

#
#       The second work horse query.  Obtain metadata information for each CCDGON
#
        query="""SELECT
            c.ccdgon_number,
            c.red_image_filename,
            i.expnum as expnum,
            i.ccdnum as ccdnum,
            c.ccd_amp as amp,
            c.inverse_variance_weight
        FROM {schema:s}{rpref:s}image i, {schema:s}{rpref:s}ccdgon c, {tT:s} g
        WHERE g.id=c.ccdgon_number
            and c.red_image_filename=i.filename
        """.format(schema=dbSchema,tT=tempTable,rpref=releasePrefix)

#       Show the query as it executes (if requested).
        if (verbose > 0):
            if (verbose == 1):
                QueryLines=query.split('\n')
                QueryOneLine='sql = '
                for line in QueryLines:
                    QueryOneLine=QueryOneLine+" "+line.strip()
                print("{:s}".format(QueryOneLine))
            if (verbose > 1):
                print("{:s}".format(query))

        curDB.execute(query)
        desc = [d[0].lower() for d in curDB.description]

#       Append query results into ccdgonDict
        for row in curDB:
            rowd = dict(zip(desc, row))
            cnum=rowd['ccdgon_number']
            if (cnum not in ccdgonDict):
                ccdgonDict[cnum]=rowd
            else:
                print("CCDGON {:d} already existed in ccdgonDict?".format(cnum))

#   Report total query execution statistics
    if (Timing):
        t1=time.time()
        print(" Queries to obtain ccdgon association and ccdgon metadata execution time: {:.2f}".format(t1-t0))

    curDB.execute('delete from {:s}'.format(tempTable))
    curDB.close()

    return molygonDict,ccdgonDict


############################################################
def collate_object_systematics(collateDict,IdArray,molyArray,BDict,molygonDict,ccdgonDict,Timing=False,verbose=0):
    """
    Take a set of Object IDs, use their associated molygon_number, to collated CCD-base metadata 
    to form numpy arrays that contain requested systematics.

    collateDict:    Dict of metadata objects and their associated data type
                    format: {systematic:dtype} currently must contain systematics for 'id' and 'band'
    IdArray:        Numpy array (of object IDs... corresponding to rows in molyArray   
    molyArray:      A numpy array (Band, Object)
    BDict:          A dictionary that gives column in molyArray correspondence to a band
                    format: {band1:col1,band2:col2}
    molygonDict:    Dict containing lists of CCDGON_NUMBER that correspond to each MOLYGON_NUMBER (key)
    ccdgonDict:     Dict with metadata from each CCDGON (key is CCDGON_NUMBER)
  
    Output:
        SystematicDict: A dictionary of numpy arrays that correspond to the keys in collateDict
    """

    t0=time.time()

    # we don't know how many objects/observations there are yet, but we know the maximum... 
    # it is number of IDs * number of CCDGONs
    nccdgon=len(ccdgonDict)

    # create a dictionary that holds numpy Arrays of systematics
    SystematicDict={}
    for key in collateDict:
       SystematicDict[key]=np.zeros(IdArray.size * nccdgon, dtype=collateDict[key])
    
    ctr = 0
    # loop over objects
    for i in xrange(IdArray.size):
        for band in BDict:
            iband=BDict[band]
            if (molyArray[iband,i]>0):
                for key in collateDict:
                    if (key not in ['id','band']):
                        key_array=np.array([ccdgonDict[ccdgon][key] for ccdgon in molygonDict[molyArray[iband,i]] ],dtype=collateDict[key])
                        SystematicDict[key][ctr:ctr+key_array.size]=key_array
                key_array=np.array([ccdgon for ccdgon in molygonDict[molyArray[iband,i]]])
                SystematicDict['id'][ctr:ctr+key_array.size][:]=IdArray[i]
                SystematicDict['band'][ctr:ctr+key_array.size][:]=iband
                ctr+=key_array.size

#   Trim off the excess array
    for key in SystematicDict:
        SystematicDict[key]=SystematicDict[key][0:ctr]

    if (Timing):
        t1=time.time()
        print(" Collate objects execution time: {:.2f}".format(t1-t0))

    return SystematicDict

############################################################
def find_object_images(ObjList,ImgDict,Timing=False,verbose=0):
    """ Associate Objects to a set of images (that presumably contain them)

        Current method uses a simple check that the coordinates of an object
        fall in the RA and DEC ranges of each image.

        Inputs:
            ObjList: List of Dictionaries containing: ID, ra, and dec for each object.
            ImgDict: Dict of Dict of Lists of Dicts
                        ImgDict[Band][CROSSRA0]=list of dicts of individual image metadat
            Timing:  flag to output timing information
            verbose: Integer setting level of verbosity of output when running.

        Returns:
            AttemptID: Resulting AttemptID
    """

#
#   Needed when writing objects (because each Object also has keys for 'ID','ra','dec')
#
    BandList=['u','g','r','i','z','Y']

    t0=time.time()
    NewObjList=[]
#
#   Work one object at a time
#
    for Obj in ObjList:
        NewObj=Obj
        ora=Obj['ra']
        odec=Obj['dec']
        for Band in ImgDict:
#
#           Currently uses simplest determination whether an object is in an image
#           
#           Could be replaced by function to check whether an object is interior to a polygon (from image corners)
#           Could also be replaced by a function using WCS info to detemine whether RA,Dec falls in the pixel range of an image.
#           
            N_ObjImg=[Img for Img in ImgDict[Band]['N'] if ((ora>Img['racmin'])and(ora<Img['racmax'])and(odec>Img['deccmin'])and(odec<Img['deccmax']))]
            Y_ObjImg=[Img for Img in ImgDict[Band]['Y'] if (((ora>Img['racmin'])or(ora<Img['racmax']))and(odec>Img['deccmin'])and(odec<Img['deccmax']))]
            NewObj[Band]=[]
            for Img in N_ObjImg:
                NewObj[Band].append({'filename':Img['filename'],'expnum':Img['expnum'],'ccdnum':Img['ccdnum']})
            for Img in Y_ObjImg:
                NewObj[Band].append({'filename':Img['filename'],'expnum':Img['expnum'],'ccdnum':Img['ccdnum']})
        NewObjList.append(NewObj)

    if (verbose > 2):
        for Obj in NewObjList:
            print("##### ID={:d} #####".format(Obj['id']))
            for Band in Obj:
                if (Band in BandList):
                    print("##### band={:s} #####".format(Band))
                    for Img in Obj[Band]:
                        print(Img)
                        print("  {:8d} {:2d} {:s} ".format(Img['expnum'],Img['ccdnum'],Img['filename']))
    if (Timing):
        t1=time.time()
        print(" Association execution time: {:.2f}".format(t1-t0))

    return NewObjList
