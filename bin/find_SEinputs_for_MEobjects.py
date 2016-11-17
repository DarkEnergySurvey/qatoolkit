#! /usr/bin/env python
# $Id$
# $Rev: 44569                            $:  # Revision of last commit.
# $LastChangedBy: rgruendl               $:  # Author of last commit.

"""
Link COADD objects to their associated single-epoch images.
"""
    
######################################################################################
def get_tile_attempt(TileName,ProcTag,dbh,dbSchema,Timing=False,verbose=0):
    """ Query code to obtain COADD tile PFW_ATTEMPT_ID after constraining
        that results are part of a specific PROCTAG.

        Inputs:
            TileName:  Tilename to be search for
            dbh:       Database connection to be used
            dbSchema:  Schema over which queries will occur.
            verbose:   Integer setting level of verbosity when running.

        Returns:
            AttemptID: Resulting AttemptID
    """

    t0=time.time()
    query="""SELECT
            t.pfw_attempt_id as pfw_attempt_id
        FROM {schema:s}proctag t, {schema:s}pfw_attempt_val av
        WHERE t.tag='{ptag:s}' 
            and t.pfw_attempt_id=av.pfw_attempt_id
            and av.key='tilename' 
            and av.val='{tname:s}'
        """.format(
            schema=dbSchema,ptag=ProcTag,tname=TileName)

    if (verbose > 0):
        if (verbose == 1):
            QueryLines=query.split('\n')
            QueryOneLine='sql = '
            for line in QueryLines:
                QueryOneLine=QueryOneLine+" "+line.strip()
            print QueryOneLine
        if (verbose > 1):
            print query
#
#   Establish a DB connection
#
    curDB = dbh.cursor()
#    curDB.arraysize = 1000 # get 1000 at a time when fetching
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
def get_ObjList(AttemptID,TableName,dbh,dbSchema,Timing=False,debug=False,verbose=0):
    """ Query to get a list of objects for a specific COADD tile
    
        TableName is provided so the code could work with both COADD_OBJECT
        and possibly a release table.
    """

    debug_constraint=''
    if (debug):
        debug_constraint='and rownum<51'

    t0=time.time()
    query="""SELECT
        o.id as id,
        o.ALPHAWIN_J2000 as ra,
        o.DELTAWIN_J2000 as dec
        FROM {DBtab:s} o, {schema}catalog c
        WHERE c.pfw_attempt_id={AttID:d}
            and c.filetype='coadd_det_cat'
            and c.filename=o.filename {dconst:s}
        """.format(schema=dbSchema,DBtab=TableName,AttID=AttemptID,dconst=debug_constraint)

#   Show the query as it executes (if requested).
    if (verbose > 0):
        if (verbose == 1):
            QueryLines=query.split('\n')
            QueryOneLine='sql = '
            for line in QueryLines:
                QueryOneLine=QueryOneLine+" "+line.strip()
            print QueryOneLine
        if (verbose > 1):
            print query

#   Establish a DB connection... execute query... assemble the results
    curDB = dbh.cursor()
    curDB.execute(query)
    desc = [d[0].lower() for d in curDB.description]

#    ObjDict={}
    ObjList=[]
    for row in curDB:
        rowd = dict(zip(desc, row))
#        ObjDict[rowd['id']]={'ra':rowd['ra'],'dec':rowd['dec']}
        ObjList.append(rowd)

    if (Timing):
        t1=time.time()
        print(" Query for objects found {:d} objects. Execution time: {:.2f}".format(len(ObjList),t1-t0))
    curDB.close()

    return ObjList


########################################################
def get_ImgDict(AttemptID,dbh,dbSchema,Timing=False,verbose=0):
    """ Query to get a dict of dicts of lists of dicts (yikes)
        that carry metadata for each images that went into a 
        COADD tile's production.

        Output dictionary, ImgDict, has the format:
            ImgDict[band][crossra0] = list of dicts where each dict carries metadata for an image.
    
        Separate dicts for the cases crossra0='Y/N' to facilitate logical search of results. 
    """

    t0=time.time()
#
#   The following query is simpler and would work if IMAGE table was 
#   properly filled by COADD pipeline.
#
#    query="""SELECT
#        i.filename as filename, 
#        i.expnum as expnum,
#        i.ccdnum as ccdnum,
#        i.band as band,
#        i.rac1,i.rac2,i.rac3,i.rac4,
#        i.decc1,i.decc2,i.decc3,i.decc4,
#        i.racmin,i.racmax,i.deccmin,i.deccmax,i.crossra0
#        FROM {schema:s}image i
#        WHERE i.pfw_attempt_id={AttID:d}
#            and i.filetype='coadd_nwgint'
#        """.format(schema=dbSchema,AttID=AttemptID)

#
#   This version of the query is necessary because the Declination ranges and corners
#   is IMAGE have been found to have a problem for the filetype='coadd_nwgint' (i.e. the
#   DB values populated by COADD pipeline have an issue... this circumvents that problem
#   by going back to the filetype='red_immask' images (i.e. the single-epoch pipeline 
#   products).
#
    query="""SELECT
        j.filename as filename, 
        j.expnum as expnum,
        j.ccdnum as ccdnum,
        j.band as band,
        j.rac1,j.rac2,j.rac3,j.rac4,
        j.decc1,j.decc2,j.decc3,j.decc4,
        j.racmin,j.racmax,j.deccmin,j.deccmax,j.crossra0
        FROM {schema:s}image i, {schema:s}image j, proctag t
        WHERE i.pfw_attempt_id={AttID:d}
            and i.filetype='coadd_nwgint'
            and i.expnum=j.expnum
            and i.ccdnum=j.ccdnum
            and j.filetype='red_immask'
            and j.pfw_attempt_id=t.pfw_attempt_id
            and t.tag='Y3A1_FINALCUT'
        """.format(schema=dbSchema,AttID=AttemptID)

#   Show the query as it executes (if requested).
    if (verbose > 0):
        if (verbose == 1):
            QueryLines=query.split('\n')
            QueryOneLine='sql = '
            for line in QueryLines:
                QueryOneLine=QueryOneLine+" "+line.strip()
            print QueryOneLine
        if (verbose > 1):
            print query

#   Establish a DB connection... execute query... assemble the results
    curDB = dbh.cursor()
    curDB.execute(query)
    desc = [d[0].lower() for d in curDB.description]

    ImgDict={}
    for row in curDB:
        rowd = dict(zip(desc, row))
        FileName=rowd['filename']
        BandVal=rowd['band']
        if (BandVal not in ImgDict):
            ImgDict[BandVal]={}
            ImgDict[BandVal]['Y']=[]
            ImgDict[BandVal]['N']=[]
        ImgDict[BandVal][rowd['crossra0']].append(rowd)

    if (Timing):
        t1=time.time()
        print(" Query for image metadata execution time: {:.2f}".format(t1-t0))
    curDB.close()

    return ImgDict


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
                        print Img
                        print("  {:8d} {:2d} {:s} ".format(Img['expnum'],Img['ccdnum'],Img['filename']))
    if (Timing):
        t1=time.time()
        print(" Association execution time: {:.2f}".format(t1-t0))

    return NewObjList


############################################################

if __name__ == "__main__":

    import argparse
    import os
    import despydb.desdbi
#    import stat
    import time
    import csv
    import sys
#    import math
#    import numpy 
#    import scipy
    import datetime
   
    parser = argparse.ArgumentParser(description='Code to find single-epoch images that contributed to a set of COADD objects') 

    parser.add_argument('-t','--tile',    action='store', type=str, required=True, help='Tile to work on')
    parser.add_argument('-p','--proctag', action='store', type=str, required=True, help='ProcTag Name')
    parser.add_argument('--dbTable',      action='store', type=str, default='coadd_object', help='DB table with objects')
    parser.add_argument('-s','--section', action='store', type=str, default=None,   help='section of .desservices file with connection info')
    parser.add_argument('-S','--Schema',  action='store', type=str, default=None,   help='DB schema (do not include \'.\').')

    parser.add_argument('-T','--Timing',  action='store_true', default=False, help='If set timing information accompanies output')
    parser.add_argument('--debug'       , action='store_true', default=False, help='Debug mode resticts code to work on a handful of objects')
    parser.add_argument('-v','--verbose', action='store', type=int, default=0, help='Verbosity (defualt:0; currently values up to 2)')

    args = parser.parse_args()
    if (args.verbose > 0):
        print "Args: ",args

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
        dbSchema="%s." % (args.Schema)

    BandList=['u','g','r','i','z','Y']
    DetBand=['r','i','z']

########################################################
#
#   Setup a DB connection
#
    try:
        desdmfile = os.environ["des_services"]
    except KeyError:
        desdmfile = None
    dbh = despydb.desdbi.DesDbi(desdmfile,args.section)

    print "Using tag_constraint of ",args.proctag
#
    AttemptID=get_tile_attempt(args.tile,args.proctag,dbh,dbSchema,Timing=args.Timing,verbose=verbose)
    ObjList=get_ObjList(AttemptID,args.dbTable,dbh,dbSchema,Timing=args.Timing,debug=args.debug,verbose=verbose)
    ImgDict=get_ImgDict(AttemptID,dbh,dbSchema,Timing=args.Timing,verbose=verbose)

#
#   Some checks for debugging purposes.
#
#    print AttemptID
#    print len(ObjList)
#    for Band in ImgDict:
#        print Band,len(ImgDict[Band]['Y']),len(ImgDict[Band]['N'])
#
#    for Img in ImgDict['i']['N']:
#        print(" {:6d} {:2d} {:10.6f} {:10.6f} {:10.6f} {:10.6f} ".format(Img['expnum'],Img['ccdnum'],Img['racmin'],Img['racmax'],Img['deccmin'],Img['deccmax']))
#        print(" {:6d} {:2d} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(Img['expnum'],Img['ccdnum'],Img['racmin'],Img['racmax'],Img['decc1'],Img['decc2'],Img['decc3'],Img['decc4']))
#

    ObjList=find_object_images(ObjList,ImgDict,Timing=args.Timing,verbose=verbose) 

    exit(0)
