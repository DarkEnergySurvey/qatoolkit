#!/usr/bin/env python
"""
Calculates mean and stddev for a set of rows for each ccd/amp along the readout axis.
It determines if there is a strong gradient in the flat in the
shutter direction and it check for large deviations from a fitted function to the data.
It uses 4 tables in ricardoc schema of desoper.

FLATS_STATS_QA_Y2 : Contains mean and stddev for each box along readout axis
ESPOSURE_REGIONS_QA_Y2 : contains the definition of the box where the stats are calculated
GAIN_RDNOISE_QA_Y2 : quick and dirty (just for a rough idea) for gain and rdnoise
EXPOSURE_TEL_QA_Y2 : Several parameters from telescope written in the header of each exposure.


.. moduleauthor:: Ricardo Covarrubias <riccov@illinois.edu>

  Needs to setup the following before running:
  python
  pyfits
  despydb
  numpy

"""
# imports
import despydb.desdbi
import cx_Oracle
import argparse
import os
import sys
import numpy as np
import pyfits as pf
import scipy as scp
import os
import math
import re
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.colorbar as cbar
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.mplot3d.axes3d import Axes3D
import time


parser = argparse.ArgumentParser(description="Check if flats is appropiate for using it")
parser.add_argument("-n", "--nite", required=True, dest="night", help="Night to select the bias and flats automatically. This will override input flats and biases")
parser.add_argument("-f", "--filter", required=False, dest="filter", help="Filter to select flats")
parser.add_argument("-t", "--filetype", required=False, dest="type",help="Input file type. raw (multi Extension compressed fits file) or precal (crosstalked fits files)")
parser.add_argument("--section", "-s", required=False, dest="section", help="section of .desservices file with connection info (db-desoper, db-destest)", default="db-desoper")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")
parser.add_argument("--forceingest", required=False, action="store_true", help="Use this to run all ingestions again")


args = parser.parse_args()
night = args.night
filter = args.filter
type = args.type
section = args.section
verbose = args.verbose
forceIngest = args.forceingest

__version__ = "$Revision: 187 $"


def connectDB(dbase):
    """
    Connect to database, query and extract data from table, and close connection to database
    dbase can be: 
    db-destest
    db-desoper
    """
    try:
        desdmfile = os.environ["DES_SERVICES"]
    except KeyError:
        desdmfile  = None

    dbh = despydb.desdbi.DesDbi(desdmfile,dbase)
        
    return dbh


def query_to_cur(dbh, query, verbose=verbose):
    """Execute a query and return a cursor to a query
    :param dbh: Database connection handler
    :param query: query to execute
    :param debug: verbosity
    
    """
    if verbose : 
        print query
    cur = dbh.cursor()
    cur.execute(query)

    return cur


def insert_dictionary_2Db(dbh, query, dictionary, verbose=verbose):
    """Execute a query and return a cursor to a query
    :param dbh: Database connection handler
    :param query: query to execute
    :param debug: verbosity
    
    """
    #if verbose : 
    #    print query
    #cur = dbh.cursor()
    #cur.execute(query, dictionary)
    #dbh.commit()
    #print "Executing query :", query_flats
    #dbcursor.arraysize = 1000 # get 1000 at a time when fetching
    #dbcursor.execute(query_flats)


    try:      
        #dbh = connectDB('db-desoper')
        cur = dbh.cursor()        
        cur.execute(query,dictionary)
        #dbh.commit()
        if verbose:
            print query
    except cx_Oracle.IntegrityError as e:
    #except cx_Oracle.DatabaseError as e:
        error, = e.args
    #    if error.code == 955:
    #        print('Table already exists')
    #    elif error.code == 1031:
    #        print("Insufficient privileges")
        print(error.code)
        print(error.message)
        print(error.context)
        print "error while inserting into : ", e.args  
    #    # Raise the exception.
        #raise


        #print "error while inserting into Table: ", e 


    #return cur


 
def medsigclip(indata, clipsig=3.0, maxiter=10, converge_num=0.001):
   """
   Computes an iteratively sigma-clipped median on a
   data set. Clipping is done about median, but media and mean
   are returned.
  
   :History:
       * 21/10/1998 Written by RSH, RITSS
       * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
       * 24/11/2009 Converted to Python. PLL.
 
   Examples
   --------
   >>> median, mean, sigma = medsigclip(indata)
 
   Parameters
   ----------
   indata: array_like
       Input data.
 
   clipsig: float
       Number of sigma at which to clip.
 
   maxiter: int
       Ceiling on number of clipping iterations.
 
   converge_num: float
       If the proportion of rejected pixels is less than
       this fraction, the iterations stop.
 
   verbose: {0, 1}
       Print messages to screen?
 
   Returns
   -------
   median: float
       N-sigma clipped median.

   mean: float
       N-sigma clipped mean.
 
   sigma: float
       Standard deviation of remaining pixels.
 
   """
   # Flatten array
   skpix = indata.reshape( indata.size, )
 
   ct = indata.size
   iter = 0; c1 = 1.0 ; c2 = 0.0
 
   while (c1 >= c2) and (iter < maxiter):
       lastct = ct
       medval = np.median(skpix)
       sig = np.std(skpix)
       wsm = np.where( abs(skpix-medval) < clipsig*sig )
       ct = len(wsm[0])
       if ct > 0:
           skpix = skpix[wsm]
 
       c1 = abs(ct - lastct)
       c2 = converge_num * lastct
       iter += 1
   # End of while loop
 
   mean  = np.mean( skpix )
   sigma = np.std( skpix )
   median = np.median( skpix )
   
   #if verbose:
   #    prf = 'MEDSIGCLIP:'
   #    print '%s %.1f-sigma clipped mean' % (prf, clipsig)
   #    print '%s Mean computed in %i iterations' % (prf, iter)
   #    print '%s , Median = %.6f, Mean = %.6f, sigma = %.6f' % (prf, medval, mean, sigma)
 
   return median, mean, sigma


def robust_mean(x):
    y = x.flatten()
    if len(y) < 6:
        return -999.
    if len(np.unique(y))==1:
        meany=y[0]
    else:
        n = len(y)
        y.sort()
        ind_qt1 = round((n+1)/4.)
        ind_qt3 = round((n+1)*3/4.)
        IQR = y[ind_qt3]- y[ind_qt1]
        lowFense = y[ind_qt1] - 1.*IQR
        highFense = y[ind_qt3] + 1.*IQR
        if lowFense == highFense:
            meany=lowFense
        else:
            ok = (y>lowFense)*(y<highFense)
            yy=y[ok]
            meany=yy.mean(dtype='double')
    return meany

def robust_median(x):
    y = x.flatten()
    if len(y) < 6:
        return -999.
    if len(np.unique(y))==1:
        meany=y[0]
    else:
        n = len(y)
        y.sort()
        ind_qt1 = round((n+1)/4.)
        ind_qt3 = round((n+1)*3/4.)
        IQR = y[ind_qt3]- y[ind_qt1]
        lowFense = y[ind_qt1] - 1.*IQR
        highFense = y[ind_qt3] + 1.*IQR
        if lowFense == highFense:
            meany=lowFense
        else:
            ok = (y>lowFense)*(y<highFense)
            yy=y[ok]
            mediany=yy.mean(dtype='double')
    return mediany


def robust_std(x):
    y = x.flatten()
    if len(y) < 6:
        return -999
    if len(np.unique(y))==1:
        meany=y[0]
    else:
        n = len(y)
        y.sort()
        ind_qt1 = round((n+1)/4.)
        ind_qt3 = round((n+1)*3/4.)
        IQR = y[ind_qt3]- y[ind_qt1]
        lowFense = y[ind_qt1] - 1.5*IQR
        highFense = y[ind_qt3] + 1.5*IQR
        if lowFense == highFense:
            meany=lowFense
        else:
            ok = (y>lowFense)*(y<highFense)
            yy=y[ok]
            meany = yy.mean(dtype='double')
            stdy = yy.std(dtype='double')
    return stdy


def calculate_mean_ccd(ccd_flat, ccdnum, detsecA, detsecB, ccdsecA, ccdsecB,trimsecA = None, trimsecB = None):
    """
    Calculates the mean of each ccd for a flat exposure
    for a set of boxes along the readout axis (columns)
    
    :param ccd_flat: sci data from flat to calculate mean
    :param ccdnum: ccdnum header keyword
    :param detsecA: detector section A header keyword
    :param detsecB: detector section B header keyword
    :param ccdsecA: ccd section A header keyword
    :param ccdsecB: ccd section B header keyword
    :param trimsecA: trim section A header keyword
    :param trimsecB: trim section B header keyword
    """
    
    #box size to calculate mean and stddev
    box_size_row = 200
    box_size_col = 200

    #print "shape os ccd\n \n \n \n ",ccd_flat.shape
    #this dictionary contains the mean counts of 
    median_ccd_values = {}
    median_ccd_values_A = {}
    median_ccd_values_B ={}
    std_ccd_values = {}
    std_ccd_values_A = {}
    std_ccd_values_B ={}
    det_box_A = {}
    det_box_B = {}
    localSectionA = []
    localSectionB = []
    mean_A = []
    mean_B = []
    globalSectionA = []
    globalSectionB = []
    
    mean_counts = 0

    #print "trimsecA ", trimsecA
    #print "trimsecB ", trimsecB
    #print "ccdsecA ", ccdsecA
    #print "ccdsecB ", ccdsecB
    #print "detsecA ", detsecA
    #print "detsecB ", detsecB
    #ccdnum 1-31: ccdsecA [1025:2048,1:4096], ccdsecB [1:1024,1:4096]
    if int(ccdnum) <= 31:
        if trimsecA is not None and trimsecB is not None:                        
            trimsecA_coords = re.search('(\d+):(\d+),(\d+):(\d+)',trimsecA)    
            x1A = int(trimsecA_coords.group(1))
            x2A = int(trimsecA_coords.group(2))
            y1A = int(trimsecA_coords.group(3))
            y2A = int(trimsecA_coords.group(4))

                
            trimsecB_coords = re.search('(\d+):(\d+),(\d+):(\d+)',trimsecB)    
            x1B = int(trimsecB_coords.group(1))
            x2B = int(trimsecB_coords.group(2))
            y1B = int(trimsecB_coords.group(3))
            y2B = int(trimsecB_coords.group(4))

        else:
            ccdsecA_coords = re.search('(\d+):(\d+),(\d+):(\d+)',ccdsecA)    
            x1A = int(ccdsecA_coords.group(1))
            x2A = int(ccdsecA_coords.group(2))
    
            ccdsecB_coords = re.search('(\d+):(\d+),(\d+):(\d+)',ccdsecB)    
            x1B = int(ccdsecB_coords.group(1))
            x2B = int(ccdsecB_coords.group(2))

        detsecA_coords = re.search('(\d+):(\d+),(\d+):(\d+)',detsecA)    
        det_x1A = int(detsecA_coords.group(1))
        det_x2A = int(detsecA_coords.group(2))
        det_y1A = int(detsecA_coords.group(3))
        det_y2A = int(detsecA_coords.group(4))
    
        detsecB_coords = re.search('(\d+):(\d+),(\d+):(\d+)',detsecB)    
        det_x1B = int(detsecB_coords.group(1))
        det_x2B = int(detsecB_coords.group(2))
        det_y1B = int(detsecB_coords.group(3))
        det_y2B = int(detsecB_coords.group(4))
        
    else:
        if trimsecA is not None and trimsecB is not None:                        
            trimsecA_coords = re.search('(\d+):(\d+),(\d+):(\d+)',trimsecA)    
            x1A = int(trimsecA_coords.group(1))
            x2A = int(trimsecA_coords.group(2))
            y1A = int(trimsecA_coords.group(3))
            y2A = int(trimsecA_coords.group(4))
                
            trimsecB_coords = re.search('(\d+):(\d+),(\d+):(\d+)',trimsecB)    
            x1B = int(trimsecB_coords.group(1))
            x2B = int(trimsecB_coords.group(2))
            y1B = int(trimsecB_coords.group(3))
            y2B = int(trimsecB_coords.group(4))

        else:
            #ccdnum 32-62: ccdsecA [1:1024,1:4096], ccdsecB [1025:2048,1:4096]
            ccdsecB_coords = re.search('(\d+):(\d+),(\d+):(\d+)',ccdsecB)    
            x1B = int(ccdsecB_coords.group(1))
            x2B = int(ccdsecB_coords.group(2))
    
            ccdsecA_coords = re.search('(\d+):(\d+),(\d+):(\d+)',ccdsecA)    
            x1A = int(ccdsecA_coords.group(1))
            x2A = int(ccdsecA_coords.group(2))


        detsecA_coords = re.search('(\d+):(\d+),(\d+):(\d+)',detsecA)    
        det_x1A = int(detsecA_coords.group(1))
        det_x2A = int(detsecA_coords.group(2))
        det_y1A = int(detsecA_coords.group(3))
        det_y2A = int(detsecA_coords.group(4))
    
        detsecB_coords = re.search('(\d+):(\d+),(\d+):(\d+)',detsecB)    
        det_x1B = int(detsecB_coords.group(1))
        det_x2B = int(detsecB_coords.group(2))
        det_y1B = int(detsecB_coords.group(3))
        det_y2B = int(detsecB_coords.group(4))
    
    #print x1A,x2A,y1A,y2A
    #print x1B,x2B,y1B,y2B
    
    if int(ccdnum) <= 31:
        #print "CCDNUM:",ccdnum
        for i,row in enumerate(range(300,4096, box_size_row), start=1):
            if row == 300:
                the_box_coords_A = '[%s:%s,%s:%s]' % (det_x1A,det_x2A,det_y1A,row+box_size_row+det_y1A)
                the_box_coords_key_A = 'MEDIAN_GBOX' + str(i) + '_A'
                stdev_the_box_coords_key_A = 'STDEV_GBOX' + str(i) + '_A'
                globalSectionA.append(the_box_coords_A)
                #mean_counts_A = np.mean(np.asarray(ccd_flat[row:row+box_size_row,300:1024]))
                secA = '[%s:%s,%s:%s]' % (x1A,x2A-300,row,row+box_size_row+1)
                localSectionA.append(secA)
                #mean_counts_A = robust_mean(ccd_flat[row:row+box_size_row,x1A:x2A-300])
                #std_counts_A = robust_std(ccd_flat[row:row+box_size_row,x1A:x2A-300])
                median_counts_A, mean_counts_A, std_counts_A = medsigclip(ccd_flat[row:row+box_size_row,x1A:x2A-300], clipsig=3.0, maxiter=10, converge_num=0.001)

                #median_counts_A = robust_median(ccd_flat[row:row+box_size_row,x1A:x2A-300])
                
                the_box_coords_B = '[%s:%s,%s:%s]' % (det_x1B,det_x2B,det_y1B,row+box_size_row+det_y1B)
                the_box_coords_key_B = 'MEDIAN_GBOX' + str(i)+ '_B'
                stdev_the_box_coords_key_B = 'STDEV_GBOX' + str(i)+ '_B'
                globalSectionB.append(the_box_coords_B)
                #mean_counts_B = np.mean(np.asarray(ccd_flat[row:row+box_size_row,1030:1748]))
                secB = '[%s:%s,%s:%s]' % (x1B+300,x2B,row,row+box_size_row+1)
                localSectionB.append(secB)
                #mean_counts_B = robust_mean(ccd_flat[row:row+box_size_row,x1B+300:x2B])
                #std_counts_B = robust_std(ccd_flat[row:row+box_size_row,x1B+300:x2B])
                median_counts_B, mean_counts_B, std_counts_B = medsigclip(ccd_flat[row:row+box_size_row,x1B+300:x2B], clipsig=3.0, maxiter=10, converge_num=0.001)
                #median_counts_B = robust_median(ccd_flat[row:row+box_size_row,x1B+300:x2B])

                ratio = median_counts_A/median_counts_B
                new_B = ratio*median_counts_B
                    
                #print mean_counts_A, the_box_coords_A, ratio, new_B, the_box_coords_B, mean_counts_B
                #print ccdnum, the_box_coords, mean_counts
            elif row+box_size_row > 3800:
                the_box_coords_A = '[%s:%s,%s:%s]' % (det_x1A,det_x2A,row+det_y1A,4096+det_y1A)
                the_box_coords_key_A = 'MEDIAN_GBOX' + str(i) + '_A'
                stdev_the_box_coords_key_A = 'STDEV_GBOX' + str(i) + '_A'
                globalSectionA.append(the_box_coords_A)
                #mean_counts_A = np.mean(np.asarray(ccd_flat[3800:3900,300:1024]))
                secA = '[%s:%s,%s:%s]' % (x1A,x2A-300,3700,3800)
                localSectionA.append(secA)
                median_counts_A, mean_counts_A, std_counts_A = medsigclip(ccd_flat[3700:3800,x1A:x2A-300], clipsig=3.0, maxiter=10, converge_num=0.001)
                #mean_counts_A = robust_mean(ccd_flat[3700:3800,x1A:x2A-300])
                #std_counts_A = robust_std(ccd_flat[3700:3800,x1A:x2A-300])
                #median_counts_A = robust_median(ccd_flat[3700:3800,x1A:x2A-300])
                the_box_coords_B = '[%s:%s,%s:%s]' % (det_x1B,det_x2B,row+det_y1B,4096+det_y1B)
                the_box_coords_key_B = 'MEDIAN_GBOX' + str(i) + '_B'
                stdev_the_box_coords_key_B = 'STDEV_GBOX' + str(i) + '_B'
                globalSectionB.append(the_box_coords_B)
                #mean_counts_B = np.mean(np.asarray(ccd_flat[3800:3900,1030:1748]))
                secB = '[%s:%s,%s:%s]' % (x1B+300,x2B,3700,3800)
                localSectionB.append(secB)
                median_counts_B, mean_counts_B, std_counts_B = medsigclip(ccd_flat[3700:3800,x1B+300:x2B], clipsig=3.0, maxiter=10, converge_num=0.001)
                #mean_counts_B = robust_mean(ccd_flat[3700:3800,x1B+300:x2B])
                #std_counts_B = robust_std(ccd_flat[3700:3800,x1B+300:x2B])
                #median_counts_B = robust_median(ccd_flat[3700:3800,x1B+300:x2B])
                ratio = median_counts_A/median_counts_B
                new_B = ratio*median_counts_B

                    
                #print mean_counts_A, the_box_coords_A, ratio,new_B, the_box_coords_B,mean_counts_B
            else:
                the_box_coords_A = '[%s:%s,%s:%s]' % (det_x1A,det_x2A,row+det_y1A,row+box_size_row+det_y1A)
                the_box_coords_key_A = 'MEDIAN_GBOX' + str(i) + '_A'
                stdev_the_box_coords_key_A = 'STDEV_GBOX' + str(i) + '_A'
                globalSectionA.append(the_box_coords_A)
                #mean_counts_A = np.mean(np.asarray(ccd_flat[row:row+box_size_row,300:1024]))
                secA = '[%s:%s,%s:%s]' % (x1A,x2A-300,row,row+box_size_row+1)
                localSectionA.append(secA)
                median_counts_A, mean_counts_A, std_counts_A = medsigclip(ccd_flat[row:row+box_size_row,x1A:x2A-300], clipsig=3.0, maxiter=10, converge_num=0.001)
                #mean_counts_A = robust_mean(ccd_flat[row:row+box_size_row,x1A:x2A-300])
                #std_counts_A = robust_std(ccd_flat[row:row+box_size_row,x1A:x2A-300])
                #medclipdata_A, meanclipdata_A, sigmaclipdata_A = medsigclip(ccd_flat[row:row+box_size_row,x1A:x2A-300])
                #median_counts_A = robust_median(ccd_flat[row:row+box_size_row,x1A:x2A-300])
                the_box_coords_B = '[%s:%s,%s:%s]' % (det_x1B,det_x2B,row+det_y1B-1,row+box_size_row+det_y1B)
                the_box_coords_key_B = 'MEDIAN_GBOX' + str(i) + '_B'
                stdev_the_box_coords_key_B = 'STDEV_GBOX' + str(i) + '_B'
                globalSectionB.append(the_box_coords_B)
                #mean_counts_B = np.mean(np.asarray(ccd_flat[row:row+box_size_row,1030:1748]))
                secB = '[%s:%s,%s:%s]' % (x1B+300,x2B,row,row+box_size_row+1)
                localSectionB.append(secB)
                median_counts_B, mean_counts_B, std_counts_B = medsigclip(ccd_flat[row:row+box_size_row,x1B+300:x2B], clipsig=3.0, maxiter=10, converge_num=0.001)
                #mean_counts_B = robust_mean(ccd_flat[row:row+box_size_row,x1B+300:x2B])
                #std_counts_B = robust_std(ccd_flat[row:row+box_size_row,x1B+300:x2B])
                #medclipdata_B, meanclipdata_B, sigmaclipdata_B = medsigclip(ccd_flat[row:row+box_size_row,x1B+300:x2B])
                #median_counts_B = robust_median(ccd_flat[row:row+box_size_row,x1B+300:x2B])
                ratio = median_counts_A/median_counts_B
                new_B = ratio*median_counts_B

            
            #print "The Box coords A %s, meanCOUNTSA %s, SecA %s, trimsecA %s, detsecA %s ccdsecA %s" % (the_box_coords_A, mean_counts_A,secA, trimsecA, detsecA, ccdsecA)
            #print "The Box coords B %s, meanCOUNTSB %s, SecB %s, trimsecB %s, detsecB %s ccdsecB %s" % (the_box_coords_B, new_B,secB, trimsecB, detsecB, ccdsecB)
            #print "MEDIANSIGCLIP %.6f, MEAN_CLIP %.6f, ROBUST_MEAN %.6f, SIGMACLIP %.6f, ROBUST_STD %.6f" % (medclipdata_A, meanclipdata_A, mean_counts_A, sigmaclipdata_A, std_counts_A)
            
            median_ccd_values_A[the_box_coords_key_A] = median_counts_A
            median_ccd_values_B[the_box_coords_key_B] = median_counts_B
            median_ccd_values[the_box_coords_A] = median_counts_A
            median_ccd_values[the_box_coords_B] = new_B
            std_ccd_values_A[stdev_the_box_coords_key_A] = std_counts_A
            std_ccd_values_B[stdev_the_box_coords_key_B] = std_counts_B
            #det_box_A['DET_BOX_A'].append(the_box_coords_A,mean_counts_A)
            #det_box_B['DET_BOX_B'].append(the_box_coords_B,mean_counts_B)


    else: #ccdnum > 31
        #print "CCDNUM:", ccdnum
        for i,row in enumerate(range(300,4096, box_size_row), start=1):
            #print "ROW:", row, row+box_size_row        
            if row == 300:
                the_box_coords_A = '[%s:%s,%s:%s]' % (det_x1A,det_x2A,det_y1A,row+box_size_row+det_y1A-1)
                the_box_coords_key_A = 'MEDIAN_GBOX' + str(i) + '_A'
                stdev_the_box_coords_key_A = 'STDEV_GBOX' + str(i) + '_A'
                globalSectionA.append(the_box_coords_A)
                #mean_counts_A = np.mean(np.asarray(ccd_flat[row:row+box_size_row,300:1024]))
                secA = '[%s:%s,%s:%s]' % (x1A+300,x2A-1,row,row+box_size_row+1)
                localSectionA.append(secA)
                median_counts_A, mean_counts_A, std_counts_A = medsigclip(ccd_flat[row:row+box_size_row,x1A+300:x2A-1], clipsig=3.0, maxiter=10, converge_num=0.001)
                #mean_counts_A = robust_mean(ccd_flat[row:row+box_size_row,x1A+300:x2A-1])
                #std_counts_A = robust_std(ccd_flat[row:row+box_size_row,x1A+300:x2A-1])
                #median_counts_A = robust_median(ccd_flat[row:row+box_size_row,x1A+300:x2A-1])
                the_box_coords_B = '[%s:%s,%s:%s]' % (det_x1B,det_x2B,det_y1B,row+box_size_row+det_y1B-1)
                the_box_coords_key_B = 'MEDIAN_GBOX' + str(i) + '_B'
                stdev_the_box_coords_key_B = 'STDEV_GBOX' + str(i) + '_B'
                globalSectionB.append(the_box_coords_B)
                #mean_counts_B = np.mean(np.asarray(ccd_flat[row:row+box_size_row,1030:1748]))
                secB = '[%s:%s,%s:%s]' % (x1B+1,x2B-300,row,row+box_size_row+1)
                localSectionB.append(secB)
                median_counts_B, mean_counts_B, std_counts_B = medsigclip(ccd_flat[row:row+box_size_row,x1B+1:x2B-300], clipsig=3.0, maxiter=10, converge_num=0.001)
                #mean_counts_B = robust_mean(ccd_flat[row:row+box_size_row,x1B+1:x2B-300])
                #std_counts_B = robust_std(ccd_flat[row:row+box_size_row,x1B+1:x2B-300])
                #median_counts_B = robust_median(ccd_flat[row:row+box_size_row,x1B+1:x2B-300])
                ratio = median_counts_A/median_counts_B
                new_B = ratio*median_counts_B
                

            elif row+box_size_row > 3800:
                the_box_coords_A = '[%s:%s,%s:%s]' % (det_x1A,det_x2A,row+det_y1A,4096+det_y1A-1)
                the_box_coords_key_A = 'MEDIAN_GBOX' + str(i) + '_A'
                stdev_the_box_coords_key_A = 'STDEV_GBOX' + str(i) + '_A'
                globalSectionA.append(the_box_coords_A)
                #mean_counts_A = np.mean(np.asarray(ccd_flat[3800:3900,300:1024]))
                secA = '[%s:%s,%s:%s]' % (x1A+300,x2A-1,3700,3800)
                localSectionA.append(secA)
                median_counts_A, mean_counts_A, std_counts_A = medsigclip(ccd_flat[3700:3800,x1A+300:x2A-1], clipsig=3.0, maxiter=10, converge_num=0.001)
                #mean_counts_A = robust_mean(ccd_flat[3700:3800,x1A+300:x2A-1])
                #std_counts_A = robust_std(ccd_flat[3700:3800,x1A+300:x2A-1])
                #median_counts_A = robust_median(ccd_flat[3700:3800,x1A+300:x2A-1])
                the_box_coords_B = '[%s:%s,%s:%s]' % (det_x1B,det_x2B,row+det_y1B,4096+det_y1B-1)
                the_box_coords_key_B = 'MEDIAN_GBOX' + str(i) + '_B'
                stdev_the_box_coords_key_B = 'STDEV_GBOX' + str(i) + '_B'
                globalSectionB.append(the_box_coords_B)
                #mean_counts_B = np.mean(np.asarray(ccd_flat[3800:3900,1030:1748]))
                secB = '[%s:%s,%s:%s]' % (x1B+1,x2B-300,3700,3800)
                localSectionB.append(secB)
                median_counts_B, mean_counts_B, std_counts_B = medsigclip(ccd_flat[3700:3800,x1B+1:x2B-300], clipsig=3.0, maxiter=10, converge_num=0.001)
                #mean_counts_B = robust_mean(ccd_flat[3700:3800,x1B+1:x2B-300])
                #std_counts_B = robust_std(ccd_flat[3700:3800,x1B+1:x2B-300])
                #median_counts_B = robust_median(ccd_flat[3700:3800,x1B+1:x2B-300])
                ratio = median_counts_A/median_counts_B
                new_B = ratio*median_counts_B
            
            
            else:
                the_box_coords_A = '[%s:%s,%s:%s]' % (det_x1A,det_x2A,row+det_y1A,row+box_size_row+det_y1A-1)
                the_box_coords_key_A = 'MEDIAN_GBOX' + str(i) + '_A'
                stdev_the_box_coords_key_A = 'STDEV_GBOX' + str(i) + '_A'
                globalSectionA.append(the_box_coords_A)
                #mean_counts_A = np.mean(np.asarray(ccd_flat[row:row+box_size_row,300:1024]))
                secA = '[%s:%s,%s:%s]' % (x1A+300,x2A-1,row,row+box_size_row+1)
                localSectionA.append(secA)
                median_counts_A, mean_counts_A, std_counts_A = medsigclip(ccd_flat[row:row+box_size_row,x1A+300:x2A-1], clipsig=3.0, maxiter=10, converge_num=0.001)
                #mean_counts_A = robust_mean(ccd_flat[row:row+box_size_row,x1A+300:x2A-1])
                #std_counts_A = robust_std(ccd_flat[row:row+box_size_row,x1A+300:x2A-1])
                #median_counts_A = robust_median(ccd_flat[row:row+box_size_row,x1A+300:x2A-1])
                the_box_coords_B = '[%s:%s,%s:%s]' % (det_x1B,det_x2B,row+det_y1B,row+box_size_row+det_y1B-1)
                the_box_coords_key_B = 'MEDIAN_GBOX' + str(i) + '_B'
                stdev_the_box_coords_key_B = 'STDEV_GBOX' + str(i) + '_B'
                globalSectionB.append(the_box_coords_B)
                #mean_counts_B = np.mean(np.asarray(ccd_flat[row:row+box_size_row,1030:1748]))
                secB = '[%s:%s,%s:%s]' % (x1B+1,x2B-300,row,row+box_size_row+1)
                localSectionB.append(secB)
                median_counts_B, mean_counts_B, std_counts_B = medsigclip(ccd_flat[row:row+box_size_row,x1B+1:x2B-300], clipsig=3.0, maxiter=10, converge_num=0.001)
                #mean_counts_B = robust_mean(ccd_flat[row:row+box_size_row,x1B+1:x2B-300])
                #std_counts_B = robust_std(ccd_flat[row:row+box_size_row,x1B+1:x2B-300])
                #median_counts_B = robust_median(ccd_flat[row:row+box_size_row,x1B+1:x2B-300])
                ratio = median_counts_A/median_counts_B
                new_B = ratio*median_counts_B
                
            #print "The Box coords A %s, meanCOUNTSA %s, SecA %s, trimsecA %s, detsecA %s ccdsecA %s" % (the_box_coords_A, mean_counts_A, secA,trimsecA, detsecA, ccdsecA)
            #print "The Box coords B %s, meanCOUNTSB %s, SecB %s, trimsecB %s, detsecB %s ccdsecB %s" % (the_box_coords_B, new_B, secB, trimsecB, detsecB, ccdsecB)

            #local_box_A[ccdnum] = localSectionA
            median_ccd_values_A[the_box_coords_key_A] = median_counts_A
            median_ccd_values_B[the_box_coords_key_B] = median_counts_B
            median_ccd_values[the_box_coords_A] = median_counts_A
            median_ccd_values[the_box_coords_B] = new_B
            std_ccd_values_A[stdev_the_box_coords_key_A] = std_counts_A
            std_ccd_values_B[stdev_the_box_coords_key_B] = std_counts_B
            #det_box_A['DET_BOX_A'].append(the_box_coords_A,mean_counts_A)
            #det_box_B['DET_BOX_B'].append(the_box_coords_B, mean_counts_B)
            
            #print mean_counts_A, the_box_coords_A, ratio, new_B, the_box_coords_B, mean_counts_B
            #print ccdnum, the_box_coords, mean_counts

        #if row == 300 or row+200 > 3800:
        #    continue
        #else:
        #    print "The Box coords A %s, SecA %s, detsecA %s ccdsecA %s" % (the_box_coords_A, secA, detsecA, ccdsecA)
        #    print "The Box coords B %s, SecB %s, detsecB %s ccdsecB %s" % (the_box_coords_B, secB, detsecB, ccdsecB)
        #mean_ccd_values[the_box_coords_A] = mean_counts_A
        #mean_ccd_values[the_box_coords_B] = new_B


    #local_box_A[ccdnum] = localSectionA
    #local_box_B[ccdnum] = localSectionB
    #global_box_A[ccdnum] = global_section_A
    #global_box_B[ccdnum] = global_section_B
    
    #for i, values in mean_ccd_values_A.values():
    #    print i,values
        
    #sys.exit()
    return median_ccd_values, median_ccd_values_A, median_ccd_values_B, std_ccd_values_A, std_ccd_values_B, localSectionA, localSectionB, globalSectionA, globalSectionB

    
def get_detsec(fitsfile):
    """
    Get the DETSEC header keyword from image
    
    :param fitsfile: Image fits file 
    
    """
    
    if fitsfile.endswith('.fz'):
        hdu_fitsfile = pf.open(fitsfile)
        fitsfilename = fitsfile.rstrip('.fz')
        ccdnum = fitsfilename.split('_')[-1]
        ccdnum = ccdnum.split('.')[0]
        print "Opened file %s %s" % (fitsfile, ccdnum)
        imheader = hdu_fitsfile[1].header
        detsec = imheader['DETSEC']
        detsecA = imheader['DETSECA']
        detsecB = imheader['DETSECB']
        ccdsecA = imheader['CCDSECA']
        ccdsecB = imheader['CCDSECB']
        gaina = imheader['GAINA']
        gainb = imheader['GAINB']

    elif fitsfile.endswith('fits'):
        ccdnum = fitsfile.split('_')[-1]
        ccdnum = ccdnum.split('.')[0]
        print "Opened file %s %s" % (fitsfile, ccdnum)
        hdu_fitsfile = pf.open(fitsfile)
        imheader = hdu_fitsfile[0].header
        detsec = imheader['DETSEC'] 
        detsecA = imheader['DETSECA']
        detsecB = imheader['DETSECB']
        ccdsecA = imheader['CCDSECA']
        ccdsecB = imheader['CCDSECB']
        gaina = float(imheader['GAINA'])
        gainb = float(imheader['GAINB'])

    
    hdu_fitsfile.close()
            
    return ccdnum, detsec, detsecA,detsecB,ccdsecA,ccdsecB,hdu_fitsfile, gaina,gainb

def get_gain_from_file():
    gainA = {}
    gainB = {}
    #gain and rdnoise dictionaries
    with open('./gain_rdnoise.txt','r') as gains:
        all_lines = gains.readlines()
        for i,line in enumerate(all_lines):
            if i == 0:
                continue
            line = line.rstrip()
            ccd, gain_b, gain_a = line.split(',')
            if int(ccd) < 10:
                ccd = '0' + str(ccd)
            if ccd > int(31):                 
                gainA[ccd] = float(gain_b)
                gainB[ccd] = float(gain_a)
            else:
                gainA[ccd] = float(gain_a)
                gainB[ccd] = float(gain_b)

    return gainA,gainB

def meanForMultiExtensionFlat(hdulist):
    """
    Function to calculate the mean of all the ccds from a src multiextention fits file from DECam.
    fits file is read from DTS directory and is compressed and not Xtalked.

    :param hdulist: pyfits object contain the headers and science data from multiextension image

    """

    all_flats_stats = {}
    median_flats_stats_A = {}
    median_flats_stats_B = {}
    std_flats_A = {}
    std_flats_B = {}
    all_coords_indexes = {}
    all_det_box_A = {}
    all_det_box_B = {}
    local_box_A = {}
    local_box_B = {}
    global_box_A = {}
    global_box_B = {}


    #all ccd detector sections (detsec)
    allDetSections = {}
    allSections = {}
    allDetSecA = {}
    allDetSecB = {}
    headerKeywords = {}

    #get some telescope and istrument values from the exposure header    
    imheader = hdulist[0].header
    #print imheader
    try:
        headerKeywords['DOMEFLOR'] = float(imheader['DOMEFLOR'])
    except:
        print "Header Keyword DOMEFLOR Not present. Inserting NULL"
        headerKeywords['DOMEFLOR'] = None
    try:
        headerKeywords['DOMEAZ'] = float(imheader['DOMEAZ'])
    except:
        print "Header Keyword DOMEAZ Not present. Inserting NULL"
        headerKeywords['DOMEAZ'] = None
    try:
        headerKeywords['MSURTEMP'] = float(imheader['MSURTEMP'])
    except:
        print "Header Keyword MSURTEMP Not present. Inserting NULL"
        headerKeywords['MSURTEMP'] = None
    try:
        headerKeywords['MAIRTEMP'] = float(imheader['MAIRTEMP'])
    except:
        print "Header Keyword MAIRTEMP Not present. Inserting NULL"
        headerKeywords['MAIRTEMP'] = None
    try:
        headerKeywords['PMW_TEMP'] = float(imheader['PMW-TEMP'])
    except:
        print "Header Keyword PMW_TEMP Not present. Inserting NULL"
        headerKeywords['PMW_TEMP'] = None
    try:
        headerKeywords['UPTRTEMP'] = float(imheader['UPTRTEMP'])
    except:
        print "Header Keyword UPTRTEMP Not present. Inserting NULL"
        headerKeywords['UPTRTEMP'] = None
    try:
        headerKeywords['HA'] = str(imheader['HA'])
    except:
        print "Header Keyword HA Not present. Inserting NULL"
        headerKeywords['HA'] = None
    try:
        headerKeywords['UTE_TEMP'] = float(imheader['UTE-TEMP'])
    except:
        print "Header Keyword UTE_TEMP Not present. Inserting NULL"
        headerKeywords['UTE_TEMP'] = None
    try:
        headerKeywords['UTN_TEMP'] = float(imheader['UTN-TEMP'])
    except:
        print "Header Keyword UTN_TEMP Not present. Inserting NULL"
        headerKeywords['UTN_TEMP'] = None
    try:
        headerKeywords['WINDDIR'] = float(imheader['WINDDIR'])
    except:
        print "Header Keyword WINDDIR Not present. Inserting NULL"
        headerKeywords['WINDDIR'] = None
    try:
        headerKeywords['OUTTEMP'] = float(imheader['OUTTEMP'])
    except:
        print "Header Keyword OUTTEMP Not present. Inserting NULL"
        headerKeywords['OUTTEMP'] = None
    try:
        headerKeywords['DIMMSEE'] = float(imheader['DIMMSEE'])
    except:
        print "Header Keyword DIMSEE Not present. Inserting NULL"
        headerKeywords['DIMMSEE'] = None
    try:
        headerKeywords['PME_TEMP'] = float(imheader['PME-TEMP'])
    except:
        print "Header Keyword PME-TEMP Not present. Inserting NULL"
        headerKeywords['PME_TEMP'] = None
    try:
        headerKeywords['DOMEHIGH'] = float(imheader['DOMEHIGH'])
    except:
        print "Header Keyword DOMEHIGH Not present. Inserting NULL"
        headerKeywords['DOMEHIGH'] = None
    try:
        headerKeywords['ZD'] = float(imheader['ZD'])
    except:
        print "Header Keyword ZD Not present. Inserting NULL"
        headerKeywords['ZD'] = None
    try:
        headerKeywords['PMN_TEMP'] = float(imheader['PMN-TEMP'])
    except:
        print "Header Keyword PMN_TEMP Not present. Inserting NULL"
        headerKeywords['PMN_TEMP'] = None
    try:
        headerKeywords['DOMELOW'] = float(imheader['DOMELOW'])
    except:
        print "Header Keyword DOMELOW Not present. Inserting NULL"
        headerKeywords['DOMELOW'] = None
    try:
        headerKeywords['PRESSURE'] = float(imheader['PRESSURE'])
    except:
        print "Header Keyword PRESSURE Not present. Inserting NULL"
        headerKeywords['PRESSURE'] = None
    try:
        headerKeywords['AZ'] = float(imheader['AZ'])
    except:
        print "Header Keyword AZ Not present. Inserting NULL"
        headerKeywords['AZ'] = None
    try:
        headerKeywords['UTS_TEMP'] = float(imheader['UTS-TEMP'])
    except:
        print "Header Keyword UTS_TEMP Not present. Inserting NULL"
        headerKeywords['UTS_TEMP'] = None
    try:
        headerKeywords['AIRMASS'] = float(imheader['AIRMASS'])
    except:
        print "Header Keyword AIRMASS Not present. Inserting NULL"
        headerKeywords['AIRMASS'] = None
    try:
        headerKeywords['PMOSTEMP'] = float(imheader['PMOSTEMP'])
    except:
        print "Header Keyword PMOSTEMP Not present. Inserting NULL"
        headerKeywords['PMOSTEMP'] = None
    try:
        headerKeywords['UTW_TEMP'] = float(imheader['UTW-TEMP'])
    except:
        print "Header Keyword UTW_TEMP Not present. Inserting NULL"
        headerKeywords['UTW_TEMP'] = None
    try:
        headerKeywords['HUMIDITY'] = float(imheader['HUMIDITY'])
    except:
        print "Header Keyword HUMIDITY Not present. Inserting NULL"
        headerKeywords['HUMIDITY'] = None
    try:
        headerKeywords['WINDSPD']  = float(imheader['WINDSPD'])
    except:
        print "Header Keyword WINDSPD Not present. Inserting NULL"
        headerKeywords['WINDSPD'] = None
    try:
        headerKeywords['LWTRTEMP'] = float(imheader['LWTRTEMP'])
    except:
        print "Header Keyword LWTRTEMP Not present. Inserting NULL"
        headerKeywords['LWTRTEMP'] = None
    try:
        headerKeywords['PMS_TEMP'] = float(imheader['PMS-TEMP'])
    except:
        print "Header Keyword PMS_TEMP Not present. Inserting NULL"
        headerKeywords['PMS_TEMP'] = None
    
        
    for ccd in range(1,63):
        imheader = hdulist[ccd].header
        ccdnum = imheader['CCDNUM']
        num = int(ccdnum)
        #skip ccd 2 and 61 (bad ccds)
        if num == 61 or num == 2:
            continue
        #read science data from image
        sci_fitsfile = hdulist[ccd].data
        detsec = imheader['DETSEC']
        detsecA = imheader['DETSECA']
        detsecB = imheader['DETSECB']
        #if exposure is from src then:
        #TRIMSECA, TRIMSECB / Good section from amp A, B 
        trimsecA = imheader['TRIMSECA']
        trimsecB = imheader['TRIMSECB']
        ccdsecA = imheader['CCDSECA']
        ccdsecB = imheader['CCDSECB']

        header_gainA = imheader['GAINA']
        header_gainB = imheader['GAINB']
        
        if num < 10:
            ccdnum = str('0'+str(ccdnum))
        
        ccdnum = str(ccdnum)
        allDetSections[ccdnum] = detsec
        allDetSecA[ccdnum] = detsecA
        allDetSecB[ccdnum] = detsecB

        #dictionary with basic detector sections for each ccd read from image header
        allSections[ccdnum] = [detsec,detsecA,detsecB,trimsecA,trimsecB,ccdsecA,ccdsecB]
        
        #gainA, gainB = get_gain_from_file()


        #calcualte mean and std for each ccd/amp
        all_flats_stats[ccdnum], median_flats_stats_A[ccdnum], median_flats_stats_B[ccdnum], std_flats_A[ccdnum], std_flats_B[ccdnum],local_box_A[ccdnum], local_box_B[ccdnum], global_box_A[ccdnum], global_box_B[ccdnum]  = calculate_mean_ccd(sci_fitsfile, 
                                                                                                                                                                                                                                                 ccdnum, detsecA, detsecB, ccdsecA, ccdsecB,trimsecA = trimsecA, 
                                                                                                                                                                                                                                                 trimsecB=trimsecB)
        
        
        print "Processing ccdnum ", ccdnum
        
        #print median_flats_stats_A
        #print median_flats_stats_B
        #print allDetSections
        #print allSections            

    return all_flats_stats, allDetSections, allSections, median_flats_stats_A, median_flats_stats_B, std_flats_A, std_flats_B, headerKeywords, local_box_A, local_box_B, global_box_A, global_box_B

def plotExposureData(allDetSections, all_flats_stats, rootExposureName,std_flats_A, std_flats_B):
    """
    Plot the data calcualted for the whole exposure
    :param allDetSections: The detector coords for each ccd in the focal plane
    :param all_flats_stats: The mean counts for all the ccds in the exposure
    :param rootExposureName: The name of the exposure to use for saving the figures
    """

    #setting the figure
    fig = plt.figure()
    plt.rc("font", size=11)
    ax = fig.add_subplot(111)
    #ax = plt.subplot(111)
    ax.set_title(rootExposureName)
    ax.set_xlabel('Focal Plane Coords')
    ax.set_ylabel('Focal Plane Coords')

    
    for ccd, section in allDetSections.items():
        det_coords = re.search('(\d+):(\d+),(\d+):(\d+)',section)
        x1 = int(det_coords.group(1))
        x2 = int(det_coords.group(2))
        y1 = int(det_coords.group(3))
        y2 = int(det_coords.group(4))
        ax.text(y1+2048,x1+1024,ccd, ha='center', va='center',color='white')
        #rect = patches.Rectangle((y1,x1),4096,2048,color='red')
        #ax.add_patch(rect)


    #get the full range of colors for all values:
    allCountsExposureList = []
    for ccd, val in all_flats_stats.items():
        if ccd == '61' or ccd == '02':
            continue
        for individualCount in val.values():
            allCountsExposureList.append(individualCount)


    allCountsExposureArray = np.asarray(allCountsExposureList,dtype=float)

    mean = np.mean(allCountsExposureArray)
    std = np.std(allCountsExposureArray)    
    norm = colors.Normalize(mean - 2*std,mean+2*std)

    
    cax, _ = cbar.make_axes(ax)
    cb2 = cbar.ColorbarBase(cax, cmap=cm.jet,norm=norm) 

    for ccd, stats_boxes in all_flats_stats.items():
        #print ccd, stats_boxes.items()
        #print stats_boxes.values()
        
        if ccd == '61' or ccd == '02':
            continue
        tempX1 = []
        tempX2 = []
        tempY1 = []
        tempY2 = []
        dx = []
        dy = []
        
        for section, values in stats_boxes.items():
            det_coords = re.search('(\d+):(\d+),(\d+):(\d+)',section)
            detx1 = int(det_coords.group(1))
            tempX1.append(int(det_coords.group(1)))
            detx2 = int(det_coords.group(2))
            tempX2.append(int(det_coords.group(2)))
            dety1 = int(det_coords.group(3))
            tempY1.append(int(det_coords.group(3)))
            dety2 = int(det_coords.group(4))
            tempY2.append(int(det_coords.group(4)))
            dx.append(detx2 - detx1 -40)
            dy.append(dety2 - dety1 -40)
            
    
        myCounts = np.asarray(stats_boxes.values())
        theCounts = stats_boxes.values()
        minCounts = min(theCounts)
        maxCounts = max(theCounts)

        normal = colors.Normalize(mean - 2*std,mean+2*std)
        thecolors = cm.jet(normal(myCounts))
        
        for x,y,w,h,c in zip(tempX1,tempY1,dx,dy,thecolors):
            #print x1_ccd,y1_ccd,mycounts   
            rect = patches.Rectangle((y,x),h,w,color=c)
            ax.add_patch(rect)


    #ax = plt.gca()
    #ax.set_xlim(ax.get_xlim()[::-1])
    #ax.set_ylim(ax.get_ylim()[::-1])
    #ax.set_xlim(-100,29050)
    #ax.set_ylim(-100,29400)
    ax.set_xlim(29050,0)
    ax.set_ylim(29400,0)

    #ax2 = fig.add_subplot(111)

    allCounts_8_12 = []
    allCounts_13_18 = []
    allCounts_19_24 = []
    allCounts_25_31 = []
    allCounts_32_38 = []
    allCounts_39_44 = []
    allCounts_45_50 = []
            
    x_13_18 = []
    x_19_24 = []
    x_25_31 = []
    x_32_38 = []
    x_39_44 = []
    x_45_50 = []

    exposure_13_18 = [13,14,15,16,17,18]
    exposure_19_24 = [19,20,21,22,23,24]
    exposure_25_31 = [24,26,27,28,29,30,31]      
    exposure_32_38 = [32,33,34,35,36,37,38]
    exposure_39_44 = [39,40,41,42,43,44]
    exposure_45_50 = [45,46,47,48,49,50]
    
    for ccd, val in all_flats_stats.items():
        if int(ccd) in exposure_25_31:
            #print "ccd",ccd
            for coords, counts in val.items():
                det_coords = re.search('(\d+):(\d+),(\d+):(\d+)',coords)
                x_25_31.append(int(det_coords.group(3)))
                allCounts_25_31.append(counts)
        elif int(ccd) in exposure_19_24:
            print "ccd",ccd
            for coords, counts in val.items():
                det_coords = re.search('(\d+):(\d+),(\d+):(\d+)',coords)
                x_19_24.append(int(det_coords.group(3)))
                allCounts_19_24.append(counts)
        elif int(ccd) in exposure_13_18:
            print "ccd",ccd
            for coords, counts in val.items():
                det_coords = re.search('(\d+):(\d+),(\d+):(\d+)',coords)
                x_13_18.append(int(det_coords.group(3)))
                allCounts_13_18.append(counts)            
        elif int(ccd) in exposure_32_38:
            print "ccd",ccd
            for coords, counts in val.items():
                det_coords = re.search('(\d+):(\d+),(\d+):(\d+)',coords)
                x_32_38.append(int(det_coords.group(3)))
                allCounts_32_38.append(counts)
        elif int(ccd) in exposure_39_44:
            print "ccd",ccd
            for coords, counts in val.items():
                det_coords = re.search('(\d+):(\d+),(\d+):(\d+)',coords)
                x_39_44.append(int(det_coords.group(3)))
                allCounts_39_44.append(counts)
        elif int(ccd) in exposure_45_50:
            print "ccd",ccd
            for coords, counts in val.items():
                det_coords = re.search('(\d+):(\d+),(\d+):(\d+)',coords)
                x_45_50.append(int(det_coords.group(3)))
                allCounts_45_50.append(counts)

    fig2 = plt.figure()
    #plt.title('S and N CCDs')

    plt.rc("font", size=11)
    ax2 = fig2.add_subplot(211)
    SandN = 'S and N CCDs ' + rootExposureName
    ax2.set_title(SandN)
    #ax2 = plt.subplot(211)

    #print len(x_32_38), len(allCounts_32_38)
    ax2.plot(x_25_31,allCounts_25_31,'ro',label='ccd 25-31')
    ax2.plot(x_19_24,allCounts_19_24,'go', label ='ccd 19-24')    
    ax2.plot(x_13_18,allCounts_13_18,'ko',label='ccd 13-18')
    ax2.legend(loc='lower left',frameon=True,scatterpoints=1,ncol=2,prop={'size':8})
    allCountsArray = np.asarray(allCounts_25_31)
    meanAllCountsArray = np.mean(allCountsArray)
    stdAllCountsArray = np.std(allCountsArray)
    ax2.set_ylim(mean - 4*stdAllCountsArray, mean + 3.5*stdAllCountsArray)
    ax2.set_xlabel('Focal Plane Coords')
    ax2.set_ylabel('Focal Plane Coords')
 
    ax3 = fig2.add_subplot(212)   
    ax3.plot(x_32_38,allCounts_32_38,'ro',label='ccd 32-38')
    ax3.plot(x_39_44,allCounts_39_44,'go',label='ccd 39-44')
    ax3.plot(x_45_50,allCounts_45_50,'ko',label='ccd 45-50')
    ax3.legend(loc='lower left',frameon=True,scatterpoints=1,ncol=2,prop={'size':8})
    allCountsArray = np.asarray(allCounts_32_38)
    meanAllCountsArray = np.mean(allCountsArray)
    stdAllCountsArray = np.std(allCountsArray)
    ax3.set_ylim(mean - 4*stdAllCountsArray, mean + 3.5*stdAllCountsArray)
    ax3.set_xlabel('Focal Plane Coords')
    ax3.set_ylabel('Focal Plane Coords')

    
    figExposurePlane = rootExposureName + '_exposure.png'
    fig.savefig(figExposurePlane)
    figExposureLines = rootExposureName + '_lines.png'
    fig2.savefig(figExposureLines)

    #####################################################
    #Plot the STDEV of each ccd/amp for the exposure
    stdFig = plt.figure()
    plt.rc("font", size=11)
    ax4 = stdFig.add_subplot(111)
    ax4.set_title(rootExposurename)
    ax4.set_xlabel('Focal Plane Coords')
    ax4.set_ylabel('Focal Plane Coords')

    
    for ccd, section in allDetSections.items():
        det_coords = re.search('(\d+):(\d+),(\d+):(\d+)',section)
        x1 = int(det_coords.group(1))
        x2 = int(det_coords.group(2))
        y1 = int(det_coords.group(3))
        y2 = int(det_coords.group(4))
        ax4.text(y1+2048,x1+1024,ccd, ha='center', va='center',color='white')

    #get the full range of colors for all values:
    allStdExposureList = []
    for ccd, val in std_flats_A.items():
        if ccd == '61' or ccd == '02':
            continue
        for individualStd in val.values():
            allStdExposureList.append(individualStd)

    for ccd, val in std_flats_B.items():
        if ccd == '61' or ccd == '02':
            continue
        for individualStd in val.values():
            allStdExposureList.append(individualStd)


    allStdExposureArray = np.asarray(allStdExposureList,dtype=float)

    mean = np.mean(allStdExposureArray)
    std = np.std(allStdExposureArray)    
    norm = colors.Normalize(mean - 2*std,mean+2*std)
    
    cax, _ = cbar.make_axes(ax4)
    cb2 = cbar.ColorbarBase(cax, cmap=cm.jet,norm=norm) 

    for ccd, stats_boxes in std_flats_A.items():
        #print ccd, stats_boxes.items()
        #print stats_boxes.values()
        
        if ccd == '61' or ccd == '02':
            continue
        tempX1 = []
        tempX2 = []
        tempY1 = []
        tempY2 = []
        dx = []
        dy = []
        
        for section, values in stats_boxes.items():
            det_coords = re.search('(\d+):(\d+),(\d+):(\d+)',section)
            detx1 = int(det_coords.group(1))
            tempX1.append(int(det_coords.group(1)))
            detx2 = int(det_coords.group(2))
            tempX2.append(int(det_coords.group(2)))
            dety1 = int(det_coords.group(3))
            tempY1.append(int(det_coords.group(3)))
            dety2 = int(det_coords.group(4))
            tempY2.append(int(det_coords.group(4)))
            dx.append(detx2 - detx1 -40)
            dy.append(dety2 - dety1 -40)
            
    
        myCounts = np.asarray(stats_boxes.values())
        theCounts = stats_boxes.values()
        minCounts = min(theCounts)
        maxCounts = max(theCounts)

        normal = colors.Normalize(mean - 2*std,mean+2*std)
        thecolors = cm.jet(normal(myCounts))
        
        for x,y,w,h,c in zip(tempX1,tempY1,dx,dy,thecolors):
            #print x1_ccd,y1_ccd,mycounts   
            rect = patches.Rectangle((y,x),h,w,color=c)
            ax4.add_patch(rect)

    for ccd, stats_boxes in std_flats_B.items():
        #print ccd, stats_boxes.items()
        #print stats_boxes.values()
        
        if ccd == '61' or ccd == '02':
            continue
        tempX1 = []
        tempX2 = []
        tempY1 = []
        tempY2 = []
        dx = []
        dy = []
        
        for section, values in stats_boxes.items():
            det_coords = re.search('(\d+):(\d+),(\d+):(\d+)',section)
            detx1 = int(det_coords.group(1))
            tempX1.append(int(det_coords.group(1)))
            detx2 = int(det_coords.group(2))
            tempX2.append(int(det_coords.group(2)))
            dety1 = int(det_coords.group(3))
            tempY1.append(int(det_coords.group(3)))
            dety2 = int(det_coords.group(4))
            tempY2.append(int(det_coords.group(4)))
            dx.append(detx2 - detx1 -40)
            dy.append(dety2 - dety1 -40)
            
    
        myCounts = np.asarray(stats_boxes.values())
        theCounts = stats_boxes.values()
        minCounts = min(theCounts)
        maxCounts = max(theCounts)

        normal = colors.Normalize(mean - 2*std,mean+2*std)
        thecolors = cm.jet(normal(myCounts))
        
        for x,y,w,h,c in zip(tempX1,tempY1,dx,dy,thecolors):
            #print x1_ccd,y1_ccd,mycounts   
            rect = patches.Rectangle((y,x),h,w,color=c)
            ax4.add_patch(rect)


    #ax = plt.gca()
    #ax.set_xlim(ax.get_xlim()[::-1])
    #ax.set_ylim(ax.get_ylim()[::-1])
    #ax.set_xlim(-100,29050)
    #ax.set_ylim(-100,29400)
    ax4.set_xlim(29050,0)
    ax4.set_ylim(29400,0)

    figExposurePlaneStd = rootExposureName + '_exposure_stdev.png'
    stdFig.savefig(figExposurePlaneStd)

    
    print "Saved images for ", rootExposureName
    #plt.show()

    return 0


def getHeaderKeywords(hdulist):
    
    headerKeywords = {}

    #get some telescope and istrument values from the exposure header    
    imheader = hdulist[0].header
    #print imheader
    try:
        headerKeywords['DOMEFLOR'] = float(imheader['DOMEFLOR'])
    except:
        print "Header Keyword DOMEFLOR Not present. Inserting NULL"
        headerKeywords['DOMEFLOR'] = None
    try:
        headerKeywords['DOMEAZ'] = float(imheader['DOMEAZ'])
    except:
        print "Header Keyword DOMEAZ Not present. Inserting NULL"
        headerKeywords['DOMEAZ'] = None
    try:
        headerKeywords['MSURTEMP'] = float(imheader['MSURTEMP'])
    except:
        print "Header Keyword MSURTEMP Not present. Inserting NULL"
        headerKeywords['MSURTEMP'] = None
    try:
        headerKeywords['MAIRTEMP'] = float(imheader['MAIRTEMP'])
    except:
        print "Header Keyword MAIRTEMP Not present. Inserting NULL"
        headerKeywords['MAIRTEMP'] = None
    try:
        headerKeywords['PMW_TEMP'] = float(imheader['PMW-TEMP'])
    except:
        print "Header Keyword PMW_TEMP Not present. Inserting NULL"
        headerKeywords['PMW_TEMP'] = None
    try:
        headerKeywords['UPTRTEMP'] = float(imheader['UPTRTEMP'])
    except:
        print "Header Keyword UPTRTEMP Not present. Inserting NULL"
        headerKeywords['UPTRTEMP'] = None
    try:
        headerKeywords['HA'] = str(imheader['HA'])
    except:
        print "Header Keyword HA Not present. Inserting NULL"
        headerKeywords['HA'] = None
    try:
        headerKeywords['UTE_TEMP'] = float(imheader['UTE-TEMP'])
    except:
        print "Header Keyword UTE_TEMP Not present. Inserting NULL"
        headerKeywords['UTE_TEMP'] = None
    try:
        headerKeywords['UTN_TEMP'] = float(imheader['UTN-TEMP'])
    except:
        print "Header Keyword UTN_TEMP Not present. Inserting NULL"
        headerKeywords['UTN_TEMP'] = None
    try:
        headerKeywords['WINDDIR'] = float(imheader['WINDDIR'])
    except:
        print "Header Keyword WINDDIR Not present. Inserting NULL"
        headerKeywords['WINDDIR'] = None
    try:
        headerKeywords['OUTTEMP'] = float(imheader['OUTTEMP'])
    except:
        print "Header Keyword OUTTEMP Not present. Inserting NULL"
        headerKeywords['OUTTEMP'] = None
    try:
        headerKeywords['DIMMSEE'] = float(imheader['DIMMSEE'])
    except:
        print "Header Keyword DIMSEE Not present. Inserting NULL"
        headerKeywords['DIMMSEE'] = None
    try:
        headerKeywords['PME_TEMP'] = float(imheader['PME-TEMP'])
    except:
        print "Header Keyword PME-TEMP Not present. Inserting NULL"
        headerKeywords['PME_TEMP'] = None
    try:
        headerKeywords['DOMEHIGH'] = float(imheader['DOMEHIGH'])
    except:
        print "Header Keyword DOMEHIGH Not present. Inserting NULL"
        headerKeywords['DOMEHIGH'] = None
    try:
        headerKeywords['ZD'] = float(imheader['ZD'])
    except:
        print "Header Keyword ZD Not present. Inserting NULL"
        headerKeywords['ZD'] = None
    try:
        headerKeywords['PMN_TEMP'] = float(imheader['PMN-TEMP'])
    except:
        print "Header Keyword PMN_TEMP Not present. Inserting NULL"
        headerKeywords['PMN_TEMP'] = None
    try:
        headerKeywords['DOMELOW'] = float(imheader['DOMELOW'])
    except:
        print "Header Keyword DOMELOW Not present. Inserting NULL"
        headerKeywords['DOMELOW'] = None
    try:
        headerKeywords['PRESSURE'] = float(imheader['PRESSURE'])
    except:
        print "Header Keyword PRESSURE Not present. Inserting NULL"
        headerKeywords['PRESSURE'] = None
    try:
        headerKeywords['AZ'] = float(imheader['AZ'])
    except:
        print "Header Keyword AZ Not present. Inserting NULL"
        headerKeywords['AZ'] = None
    try:
        headerKeywords['UTS_TEMP'] = float(imheader['UTS-TEMP'])
    except:
        print "Header Keyword UTS_TEMP Not present. Inserting NULL"
        headerKeywords['UTS_TEMP'] = None
    try:
        headerKeywords['AIRMASS'] = float(imheader['AIRMASS'])
    except:
        print "Header Keyword AIRMASS Not present. Inserting NULL"
        headerKeywords['AIRMASS'] = None
    try:
        headerKeywords['PMOSTEMP'] = float(imheader['PMOSTEMP'])
    except:
        print "Header Keyword PMOSTEMP Not present. Inserting NULL"
        headerKeywords['PMOSTEMP'] = None
    try:
        headerKeywords['UTW_TEMP'] = float(imheader['UTW-TEMP'])
    except:
        print "Header Keyword UTW_TEMP Not present. Inserting NULL"
        headerKeywords['UTW_TEMP'] = None
    try:
        headerKeywords['HUMIDITY'] = float(imheader['HUMIDITY'])
    except:
        print "Header Keyword HUMIDITY Not present. Inserting NULL"
        headerKeywords['HUMIDITY'] = None
    try:
        headerKeywords['WINDSPD']  = float(imheader['WINDSPD'])
    except:
        print "Header Keyword WINDSPD Not present. Inserting NULL"
        headerKeywords['WINDSPD'] = None
    try:
        headerKeywords['LWTRTEMP'] = float(imheader['LWTRTEMP'])
    except:
        print "Header Keyword LWTRTEMP Not present. Inserting NULL"
        headerKeywords['LWTRTEMP'] = None
    try:
        headerKeywords['PMS_TEMP'] = float(imheader['PMS-TEMP'])
    except:
        print "Header Keyword PMS_TEMP Not present. Inserting NULL"
        headerKeywords['PMS_TEMP'] = None

    return headerKeywords


def insertTelescopeExposureInfo(dbh,headerKeywords,filename,filetype):
    """ Insert into DB only telescope info"""
    #prepare the insert statement
    
    sqlInsertIntoExposureQA = """insert into ricardoc.exposure_tel_qa_y2
    (FILENAME,FILETYPE,DOMEFLOR,MSURTEMP,MAIRTEMP, PMW_TEMP,UPTRTEMP,
    HA,DOMEAZ,UTE_TEMP,UTN_TEMP,WINDDIR,OUTTEMP,DIMMSEE,PME_TEMP,
    DOMEHIGH,ZD,PMN_TEMP,DOMELOW,PRESSURE,AZ,UTS_TEMP,AIRMASS,PMOSTEMP,
    UTW_TEMP,HUMIDITY,WINDSPD,LWTRTEMP,PMS_TEMP) 
    VALUES (:FILENAME,:FILETYPE,:DOMEFLOR,:MSURTEMP,:MAIRTEMP,:PMW_TEMP,:UPTRTEMP,
    :HA,:DOMEAZ,:UTE_TEMP,:UTN_TEMP,:WINDDIR,:OUTTEMP,:DIMMSEE,:PME_TEMP,
    :DOMEHIGH,:ZD,:PMN_TEMP,:DOMELOW,:PRESSURE,:AZ,:UTS_TEMP,:AIRMASS,:PMOSTEMP,
    :UTW_TEMP,:HUMIDITY,:WINDSPD,:LWTRTEMP,:PMS_TEMP)"""

    #merge statement example
    #
    #merge into ricardoc.exposure_qa exqa using (select id from exposure) eid on (eid.id = exqa.exposureid) when matched then update set exqa.MSURTEMP = 40.0;
    #
    
    #insert into headerKeywords dictionary the exposure.id value
    headerKeywords['FILENAME'] = filename
    headerKeywords['FILETYPE'] = filetype
    
    
    #print headerKeywords
    
    #print headerKeywords.keys()
    #print headerKeywords.values()
    
    #print sqlInsertIntoExposureQA
    insert_dictionary_2Db(dbh, sqlInsertIntoExposureQA, headerKeywords)
    
    dbh.commit()


def saveDataToDatabase(dbh, filename, filetype, band, night, median_flats_stats_A, median_flats_stats_B, std_flats_A, std_flats_B, allSections, headerKeywords, local_box_A, local_box_B, global_box_A, global_box_B):
    """
    Save all counts measured for each exposure/ccd to the database table ricardoc.flat_stats_qa
    :param dbh: Database connection
    :param expid: Exposure ID from exposure table
    :param all_flats_stats: Dictionary with detector sections with measurements for each ccd
    :param allDetSections: The detector sections on the focal plane
    :param headerKeyword: All headers keywords to insert in table
    """
    """    
    ricardoc.exposure_regions_qa_y2
    ricardoc.exposure_tel_qa_y2
    ricardoc.flats_stats_qa_y2
    

    """

    #prepare the insert statement
    sqlInsertIntoExposureQA = """insert into ricardoc.exposure_tel_qa_y2
    (FILENAME,FILETYPE,DOMEFLOR,MSURTEMP,MAIRTEMP, PMW_TEMP,UPTRTEMP,
    HA,DOMEAZ,UTE_TEMP,UTN_TEMP,WINDDIR,OUTTEMP,DIMMSEE,PME_TEMP,
    DOMEHIGH,ZD,PMN_TEMP,DOMELOW,PRESSURE,AZ,UTS_TEMP,AIRMASS,PMOSTEMP,
    UTW_TEMP,HUMIDITY,WINDSPD,LWTRTEMP,PMS_TEMP) 
    VALUES (:FILENAME,:FILETYPE,:DOMEFLOR,:MSURTEMP,:MAIRTEMP,:PMW_TEMP,:UPTRTEMP,
    :HA,:DOMEAZ,:UTE_TEMP,:UTN_TEMP,:WINDDIR,:OUTTEMP,:DIMMSEE,:PME_TEMP,
    :DOMEHIGH,:ZD,:PMN_TEMP,:DOMELOW,:PRESSURE,:AZ,:UTS_TEMP,:AIRMASS,:PMOSTEMP,
    :UTW_TEMP,:HUMIDITY,:WINDSPD,:LWTRTEMP,:PMS_TEMP)"""

    #merge statement example
    #
    #merge into ricardoc.exposure_qa exqa using (select id from exposure) eid on (eid.id = exqa.exposureid) when matched then update set exqa.MSURTEMP = 40.0;
    #
    
    #insert into headerKeywords dictionary the exposure.id value
    headerKeywords['FILENAME'] = filename
    headerKeywords['FILETYPE'] = filetype
    
    
    #print headerKeywords
    
    #print headerKeywords.keys()
    #print headerKeywords.values()
    
    #print sqlInsertIntoExposureQA
    insert_dictionary_2Db(dbh, sqlInsertIntoExposureQA, headerKeywords)
    
    detector_sections = ['DETSEC', 'DETSEC_A', 'DETSEC_B', 'TRIMSEC_A', 'TRIMSEC_B', 'CCDSEC_A', 'CCDSEC_B']

    local_exp_regions_keys_A = ['LBOX1_A','LBOX2_A','LBOX3_A','LBOX4_A','LBOX5_A','LBOX6_A','LBOX7_A','LBOX8_A','LBOX9_A','LBOX10_A','LBOX11_A','LBOX12_A',
                                'LBOX13_A','LBOX14_A','LBOX15_A','LBOX16_A','LBOX17_A','LBOX18_A','LBOX19_A']
    local_exp_regions_keys_B = ['LBOX1_B','LBOX2_B','LBOX3_B','LBOX4_B',
                                'LBOX5_B','LBOX6_B','LBOX7_B','LBOX8_B','LBOX9_B','LBOX10_B','LBOX11_B','LBOX12_B','LBOX13_B','LBOX14_B','LBOX15_B',
                                'LBOX16_B','LBOX17_B','LBOX18_B','LBOX19_B']
    global_exp_regions_keys_A = ['GBOX1_A','GBOX2_A','GBOX3_A','GBOX4_A','GBOX5_A','GBOX6_A','GBOX7_A','GBOX8_A','GBOX9_A','GBOX10_A','GBOX11_A','GBOX12_A',
                                'GBOX13_A','GBOX14_A','GBOX15_A','GBOX16_A','GBOX17_A','GBOX18_A','GBOX19_A']
    global_exp_regions_keys_B = ['GBOX1_B','GBOX2_B','GBOX3_B','GBOX4_B',
                                'GBOX5_B','GBOX6_B','GBOX7_B','GBOX8_B','GBOX9_B','GBOX10_B','GBOX11_B','GBOX12_B','GBOX13_B','GBOX14_B','GBOX15_B',
                                'GBOX16_B','GBOX17_B','GBOX18_B','GBOX19_B']
    
    insertDictionaryExposureRegionsQa = {}

    counter = 0
    #I'm going to try to create one dictionary for all
    for ccd in allSections.keys():
        insertDictionaryExposureRegionsQa['FILENAME'] = filename
        insertDictionaryExposureRegionsQa['FILETYPE'] = filetype
        insertDictionaryExposureRegionsQa['NITE'] = night
        insertDictionaryExposureRegionsQa['BAND'] = band
        insertDictionaryExposureRegionsQa['CCDNUM'] = ccd
        for key, val in zip(detector_sections, allSections[ccd]):
            insertDictionaryExposureRegionsQa[key] = val
        for key, val in zip(local_exp_regions_keys_A, local_box_A[ccd]):
            insertDictionaryExposureRegionsQa[key] = val
        for key, val in zip(local_exp_regions_keys_B, local_box_B[ccd]):
            insertDictionaryExposureRegionsQa[key] = val
        for key, val in zip(global_exp_regions_keys_A, global_box_A[ccd]):
            insertDictionaryExposureRegionsQa[key] = val
        for key, val in zip(global_exp_regions_keys_B, global_box_B[ccd]):
            insertDictionaryExposureRegionsQa[key] = val
        
        insertData = """ insert into ricardoc.exposure_regions_qa_y2 
                    (FILENAME, FILETYPE, CCDNUM, NITE, BAND, DETSEC, DETSEC_A, 
                    DETSEC_B, TRIMSEC_A, TRIMSEC_B, CCDSEC_A, CCDSEC_B,
                    LBOX1_A, LBOX2_A, LBOX3_A, LBOX4_A, LBOX5_A, LBOX6_A,
                    LBOX7_A, LBOX8_A, LBOX9_A, LBOX10_A, LBOX11_A, LBOX12_A,
                    LBOX13_A, LBOX14_A, LBOX15_A, LBOX16_A, LBOX17_A, LBOX18_A,
                    LBOX19_A, 
                    LBOX1_B, LBOX2_B, LBOX3_B, LBOX4_B, LBOX5_B, LBOX6_B,
                    LBOX7_B, LBOX8_B, LBOX9_B, LBOX10_B, LBOX11_B, LBOX12_B,
                    LBOX13_B, LBOX14_B, LBOX15_B, LBOX16_B, LBOX17_B, LBOX18_B,
                    LBOX19_B, 
                    GBOX1_A, GBOX2_A, GBOX3_A, GBOX4_A, GBOX5_A, GBOX6_A,
                    GBOX7_A, GBOX8_A, GBOX9_A, GBOX10_A, GBOX11_A, GBOX12_A,
                    GBOX13_A, GBOX14_A, GBOX15_A, GBOX16_A, GBOX17_A, GBOX18_A,
                    GBOX19_A, 
                    GBOX1_B, GBOX2_B, GBOX3_B, GBOX4_B, GBOX5_B, GBOX6_B,
                    GBOX7_B, GBOX8_B, GBOX9_B, GBOX10_B, GBOX11_B, GBOX12_B,
                    GBOX13_B, GBOX14_B, GBOX15_B, GBOX16_B, GBOX17_B, GBOX18_B,
                    GBOX19_B) VALUES 
                    (:FILENAME, :FILETYPE, :CCDNUM, :NITE, :BAND, :DETSEC, :DETSEC_A, 
                    :DETSEC_B, :TRIMSEC_A, :TRIMSEC_B, :CCDSEC_A, :CCDSEC_B,
                    :LBOX1_A, :LBOX2_A, :LBOX3_A, :LBOX4_A, :LBOX5_A, :LBOX6_A,
                    :LBOX7_A, :LBOX8_A, :LBOX9_A, :LBOX10_A, :LBOX11_A, :LBOX12_A,
                    :LBOX13_A, :LBOX14_A, :LBOX15_A, :LBOX16_A, :LBOX17_A, :LBOX18_A,
                    :LBOX19_A,
                    :LBOX1_B, :LBOX2_B, :LBOX3_B, :LBOX4_B, :LBOX5_B, :LBOX6_B,
                    :LBOX7_B, :LBOX8_B, :LBOX9_B, :LBOX10_B, :LBOX11_B, :LBOX12_B,
                    :LBOX13_B, :LBOX14_B, :LBOX15_B, :LBOX16_B, :LBOX17_B, :LBOX18_B,
                    :LBOX19_B, 
                    :GBOX1_A, :GBOX2_A, :GBOX3_A, :GBOX4_A, :GBOX5_A, :GBOX6_A,
                    :GBOX7_A, :GBOX8_A, :GBOX9_A, :GBOX10_A, :GBOX11_A, :GBOX12_A,
                    :GBOX13_A, :GBOX14_A, :GBOX15_A, :GBOX16_A, :GBOX17_A, :GBOX18_A,
                    :GBOX19_A,
                    :GBOX1_B, :GBOX2_B, :GBOX3_B, :GBOX4_B, :GBOX5_B, :GBOX6_B,
                    :GBOX7_B, :GBOX8_B, :GBOX9_B, :GBOX10_B, :GBOX11_B, :GBOX12_B,
                    :GBOX13_B, :GBOX14_B, :GBOX15_B, :GBOX16_B, :GBOX17_B, :GBOX18_B,
                    :GBOX19_B)"""                 
                        
        
        insert_dictionary_2Db(dbh, insertData, insertDictionaryExposureRegionsQa)

    #insert all calculated values for each bin/ccd    
    for key in median_flats_stats_A.keys():
        dictForKey = dict(median_flats_stats_A[key].items() + median_flats_stats_B[key].items() + std_flats_A[key].items() + std_flats_B[key].items())
        dictForKey['FILENAME'] = filename
        dictForKey['NITE'] = night
        dictForKey['BAND'] = band
        dictForKey['CCDNUM'] = key
        dictForKey['FILETYPE'] = filetype 
        
        firstInsert = """insert into ricardoc.flats_stats_qa_y2
                        (FILENAME, FILETYPE, CCDNUM, NITE, BAND, MEDIAN_GBOX1_A, MEDIAN_GBOX2_A, MEDIAN_GBOX3_A, MEDIAN_GBOX4_A,
                        MEDIAN_GBOX5_A, MEDIAN_GBOX6_A, MEDIAN_GBOX7_A, MEDIAN_GBOX8_A, MEDIAN_GBOX9_A, MEDIAN_GBOX10_A, MEDIAN_GBOX11_A,
                        MEDIAN_GBOX12_A, MEDIAN_GBOX13_A, MEDIAN_GBOX14_A, MEDIAN_GBOX15_A, MEDIAN_GBOX16_A, MEDIAN_GBOX17_A, MEDIAN_GBOX18_A,
                        MEDIAN_GBOX19_A,
                        MEDIAN_GBOX1_B, MEDIAN_GBOX2_B, MEDIAN_GBOX3_B, MEDIAN_GBOX4_B,
                        MEDIAN_GBOX5_B, MEDIAN_GBOX6_B, MEDIAN_GBOX7_B, MEDIAN_GBOX8_B, MEDIAN_GBOX9_B, MEDIAN_GBOX10_B, MEDIAN_GBOX11_B,
                        MEDIAN_GBOX12_B, MEDIAN_GBOX13_B, MEDIAN_GBOX14_B, MEDIAN_GBOX15_B, MEDIAN_GBOX16_B, MEDIAN_GBOX17_B, MEDIAN_GBOX18_B,
                        MEDIAN_GBOX19_B,
                        STDEV_GBOX1_A, STDEV_GBOX2_A, STDEV_GBOX3_A, STDEV_GBOX4_A,
                        STDEV_GBOX5_A, STDEV_GBOX6_A, STDEV_GBOX7_A, STDEV_GBOX8_A, STDEV_GBOX9_A, STDEV_GBOX10_A, STDEV_GBOX11_A,
                        STDEV_GBOX12_A, STDEV_GBOX13_A, STDEV_GBOX14_A, STDEV_GBOX15_A, STDEV_GBOX16_A, STDEV_GBOX17_A, STDEV_GBOX18_A,
                        STDEV_GBOX19_A,
                        STDEV_GBOX1_B, STDEV_GBOX2_B, STDEV_GBOX3_B, STDEV_GBOX4_B,
                        STDEV_GBOX5_B, STDEV_GBOX6_B, STDEV_GBOX7_B, STDEV_GBOX8_B, STDEV_GBOX9_B, STDEV_GBOX10_B, STDEV_GBOX11_B,
                        STDEV_GBOX12_B, STDEV_GBOX13_B, STDEV_GBOX14_B, STDEV_GBOX15_B, STDEV_GBOX16_B, STDEV_GBOX17_B, STDEV_GBOX18_B,
                        STDEV_GBOX19_B) VALUES (:FILENAME, :FILETYPE, :CCDNUM, :NITE, :BAND, :MEDIAN_GBOX1_A, :MEDIAN_GBOX2_A, 
                        :MEDIAN_GBOX3_A, :MEDIAN_GBOX4_A, :MEDIAN_GBOX5_A, :MEDIAN_GBOX6_A, :MEDIAN_GBOX7_A, :MEDIAN_GBOX8_A,
                        :MEDIAN_GBOX9_A,:MEDIAN_GBOX10_A, :MEDIAN_GBOX11_A, :MEDIAN_GBOX12_A, :MEDIAN_GBOX13_A, :MEDIAN_GBOX14_A,
                        :MEDIAN_GBOX15_A, :MEDIAN_GBOX16_A, :MEDIAN_GBOX17_A, :MEDIAN_GBOX18_A, :MEDIAN_GBOX19_A,
                        :MEDIAN_GBOX1_B, :MEDIAN_GBOX2_B, 
                        :MEDIAN_GBOX3_B, :MEDIAN_GBOX4_B, :MEDIAN_GBOX5_B, :MEDIAN_GBOX6_B, :MEDIAN_GBOX7_B, :MEDIAN_GBOX8_B,
                        :MEDIAN_GBOX9_B,:MEDIAN_GBOX10_B, :MEDIAN_GBOX11_B, :MEDIAN_GBOX12_B, :MEDIAN_GBOX13_B, :MEDIAN_GBOX14_B,
                        :MEDIAN_GBOX15_B, :MEDIAN_GBOX16_B, :MEDIAN_GBOX17_B, :MEDIAN_GBOX18_B, :MEDIAN_GBOX19_B,
                        :STDEV_GBOX1_A, :STDEV_GBOX2_A,
                        :STDEV_GBOX3_A, :STDEV_GBOX4_A, :STDEV_GBOX5_A, :STDEV_GBOX6_A, :STDEV_GBOX7_A, :STDEV_GBOX8_A,
                        :STDEV_GBOX9_A,:STDEV_GBOX10_A, :STDEV_GBOX11_A, :STDEV_GBOX12_A, :STDEV_GBOX13_A, :STDEV_GBOX14_A,
                        :STDEV_GBOX15_A, :STDEV_GBOX16_A, :STDEV_GBOX17_A, :STDEV_GBOX18_A, :STDEV_GBOX19_A,
                        :STDEV_GBOX1_B, :STDEV_GBOX2_B, 
                        :STDEV_GBOX3_B, :STDEV_GBOX4_B, :STDEV_GBOX5_B, :STDEV_GBOX6_B, :STDEV_GBOX7_B, :STDEV_GBOX8_B,
                        :STDEV_GBOX9_B,:STDEV_GBOX10_B, :STDEV_GBOX11_B, :STDEV_GBOX12_B, :STDEV_GBOX13_B, :STDEV_GBOX14_B,
                        :STDEV_GBOX15_B, :STDEV_GBOX16_B, :STDEV_GBOX17_B, :STDEV_GBOX18_B, :STDEV_GBOX19_B)"""
            
        insert_dictionary_2Db(dbh, firstInsert, dictForKey)
    
    dbh.commit()
    

           
    
def main():
    

    #start = time.time()
    #if necessary, I can use the gain from a file.
    gainA, gainB = get_gain_from_file()
    
    dbh = connectDB(section)
    #print "band", filter
    
    #src means raw data. It is compressed and nothing has done to it (no xtalk, etc)
    if type == 'raw':
        #root_dir = '/archive_data/Archive/DTS/src/'+night+'/src/' This is the root directory before refactored system
        root_dir = '/archive_data/desarchive/DTS/raw/'+night
        if filter is None:
            #Query for all bands but VR given as input
            print "no band given. Using all band present"
            
            #This is the new query I should use. There are not id anymore in the database
            #select band,filename,exptime,nite,expnum from exposure where obstype = 'dome flat' and band = 'g' and band not like '%Empty%' and nite = '20140807' order by filename;               
            queryitems = ["filename","band","filetype","exptime","nite"]
            querylist = ",".join(queryitems)
            query_flats = """select %s from exposure where 
                    obstype = 'dome flat'
                    and band not like '%%Empty%%'
                    and nite = '%s' order by filename""" % (querylist, night)
        else:
            #Query with specific band given as input
            queryitems = ["filename","band","filetype","exptime","nite"]
            querylist = ",".join(queryitems)
            query_flats = """select %s from exposure where 
                    obstype = 'dome flat' 
                    and band = '%s'
                    and band not like '%%Empty%%'
                    and nite = '%s' order by filename""" % (querylist, filter, night)


        cur_flats = query_to_cur(dbh, query_flats, verbose)

    #if it is raw, then data has been xtalked, bias corrected.
    elif type == 'precal':
        #type should be precal
        #(This still needs to be fixed. The query needs refinement)
        root_dir = '/archive_data/desarchive/'
        #Need to add nite-r{reqnum}/p{attnum}
        query_flats = """select a.path||'/'||a.filename,c.band,c.filetype,c.nite from file_archive_info a,calibration c 
                        where c.filename=a.filename 
                        and c.filetype='xtalked_dflat'
                        and c.band = '%s'
                        and c.nite = %s order by c.filename, c.ccdnum""" % (filter, night)
        
        #queryitems = ["band","nite"]
        #querylist = ",".join(queryitems)
        #query_flats = """select %s from exposure where 
        #        obstype = 'dome flat' and 
        #        band = '%s' and 
        #        project='DTS' 
        #        and nite = '%s' order by expnum""" % (querylist, filter, night)

        cur_flats = query_to_cur(dbh, query_flats, verbose)
        #print "Executing query :", query_flats
        #dbcursor.arraysize = 1000 # get 1000 at a time when fetching
        #dbcursor.execute(query_flats)
        
    else: #If type is none, then do it locally in the given directory below
        
        #dire = '/des005/red/20131125091730_20131121/raw/DECam_00257145/'
        root_dir = '/Users/ricardo/DESDM/DES-SNe/Manifest_database/image'
    #os.chdir(root_dir)
    #filenames = os.listdir(root_dir)
    #print filenames               

    #This dictionary contains the ccdnum as a key and a dictionary with all the mean values
    #for a set of boxes in each ccd
    median_flats_stats_A = {}
    median_flats_stats_B = {}
    std_flats_A = {}
    std_flats_B = {}
    all_flats_stats = {}
    all_coords_indexes = {}

    #all ccd detector sections (detsec)
    allDetSections = {}
    allDetSecA = {}
    allDetSecB = {}
    
 

    #loops throught the flats found for a given night
    for flats in cur_flats:
        #insert into ricardo.exposure_stats_qa_y2 the values of the exposure found
        filename, band, filetype, exptime, nite = flats
        
        print filename, band
        
        
        #Check if the exposure has already been ingested into the database table ricardoc.flats_stats_qa
        #if it has been ingested and processed for all ccds, then move to the next flat
        #checkExposure =  """select distinct(filename) from ricardoc.flats_stats_qa_y2 where FILENAME = '%s'""" % filename
        checkExposure =  """select count(distinct ccdnum) from ricardoc.flats_stats_qa_y2 where FILENAME = '%s'""" % filename
        #dbh = connectDB(section)
        exposureIngested = query_to_cur(dbh, checkExposure, verbose)
        
        #print exposureIngested
        #if all exposure were ingested, then ccdnum=60
        exposureWasIngested = False
        
        for val in exposureIngested:
            numCcdProc = val[0]
            print "Checking Exposure \n", filename
            if numCcdProc == 60:
                print "Exposure %s has %s ccd processed. \n" % (filename, numCcdProc)
                if forceIngest:
                    exposureWasIngested = False
                else:
                    exposureWasIngested = True
            else:
                print "Only %s ccd has been processed! reprocessing filename %s \n" % ( numCcdProc, filename)
        if verbose:
            print "Has Exposure Been ingested? \n", exposureWasIngested
        if exposureWasIngested:
            continue

        
                
        #SqlExposureInfo = 'select expnum, nite, band, exposurename, detsize, object, exptime, obstype, project from exposure where id=\'%s\' and project=\'%s\' ' % (expid, project)
        #if verbose:
        #    print "SqlExposureInfo", SqlExposureInfo
        #cur_exposure = query_to_cur(dbh, SqlExposureInfo, verbose)
        #for junk in cur_exposure:
        #    print junk
            
        #sys.exit()
        #loops trhough all the src exposures found
        #root_dir contains the path to the src file
        #for each src file, a rootExposureName can be acquaried
        #the order for each line in the loop is:
        #(494987031, 'g', 'DECam_00276299.fits', 30.0, '20140120', 00229034)
        #rootExposureName = str(exposurename.split('.')[0])
        #rootExposureName = rootExposureName+'_'+str(band)+'_'+str(nite)
        
        #if exposurename.endswith('.fz') or exposurename.endswith('fits'):
        if type=='raw':
            file = str(filename)+'.fz'
        print "Opening file ", file
        file = os.path.join(root_dir,file)
        hdu_fitsfile = pf.open(file, mode='readonly')
        print "openning filename ", file
        print "filename has headers ", len(hdu_fitsfile)
        if len(hdu_fitsfile) == 71:
            print "found raw file \n", filename
            
            #if forceIngest:
            #    if verbose:
            #        print "getting headers"
            #    headerkeywords = getHeaderKeywords(hdu_fitsfile)
            #    insertTelescopeExposureInfo(dbh,headerkeywords,filename,filetype)
            #else:
            
            all_flats_stats, allDetSections, allSections, median_flats_stats_A, median_flats_stats_B, std_flats_A, std_flats_B, headerKeywords,local_box_A, local_box_B, global_box_A, global_box_B = meanForMultiExtensionFlat(hdu_fitsfile)
            
            #plotExposureData(allDetSections, all_flats_stats, rootExposureName,std_flats_A, std_flats_B)
            hdu_fitsfile.close()
                           
            saveDataToDatabase(dbh, filename, filetype, band, nite, median_flats_stats_A, median_flats_stats_B, std_flats_A, std_flats_B, allSections, headerKeywords, local_box_A, local_box_B, global_box_A, global_box_B)
        
        else: #image is precal (already crosstalked)
            ccdnum, detsec, detsecA, detsecB, ccdsecA, ccdsecB, hdu_fitsfile = get_detsec(file)
            allDetSections[ccdnum] = detsec
            allDetSecA[ccdnum] = detsecA
            allDetSecB[ccdnum] = detsecB
    
            sci_fitsfile = hdu_fitsfile[1].data
            all_flats_stats[ccdnum], median_flats_stats_A[ccdnum], median_flats_stats_B[ccdnum] = calculate_mean_ccd(sci_fitsfile, ccdnum, detsecA, detsecB, ccdsecA, ccdsecB)
            hdu_fitsfile.close()
        
        #else: #found a non fits file
        #    print "found a file that doesn't look like it is a fits: ", exposurename
        #    continue
        #if verbose:
        #    print "CCDNUM %s , DETSEC %s " % (ccdnum, detsec)
  
    #plot the data when we have a raw exposures (already crosstalked)
    
    
    #if len(hdu_fitsfile) != 71:
    #    plotExposureData(allDetSections, all_flats_stats)
    #    saveDataToDatabase(all_flats_stats,allDetSections)
        #if len(hdulist) != 71:
        #    allDetSections[ccdnum] = detsec
        #    allDetSecA[ccdnum] = detsecA
        #    allDetSecB[ccdnum] = detsecB
    
        #    sci_fitsfile = hdu_fitsfile[1].data
        #    all_flats_stats[ccdnum] = calculate_mean_ccd(sci_fitsfile, gainA[ccdnum],gainB[ccdnum], header_gainA, header_gainB, ccdnum, detsecA, detsecB, ccdsecA, ccdsecB)
        
 

    #here is where all the plotting part goes

    #trying to do a 3d plot
    
    #print x1_ccd, y1_ccd, values.values(), len(x1_ccd), len(y1_ccd), len(values.values())
    #p = ax.plot_wireframe(x1_ccd, y1_ccd, values.values(), rstride=4,cstride=4)
    #ax.view_init(0, 45)
    #else:
    #    continue
    
    return 0

    
if  __name__ == '__main__':
    status = main()
    sys.exit(status)

    
    
