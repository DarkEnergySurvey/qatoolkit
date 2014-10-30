#!/usr/bin/env python

# imports
import coreutils.desdbi
import cx_Oracle
import argparse
import os
import sys
import numpy as np
import numpy.polynomial.polynomial as poly
import pyfits as pf
import scipy
import math
import re
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.colorbar as cbar
import matplotlib.colors as colors
import matplotlib.cm as cm
#from mpl_toolkits.mplot3d.axes3d import Axes3D
import itertools
import time
import datetime as dt
import matplotlib.dates as mdates


parser = argparse.ArgumentParser(description="Analysis of flats")
parser.add_argument("-d", "--dates", required=True, dest="dates", help="Range of dates or single date with format: first_nite-last_nite (20130801-20131120")
parser.add_argument("-b", "--band", required=True, dest="band", help="band to use. A list of band comma separated is accepted")
parser.add_argument("-c", "--ccds", required=True, dest="ccdrange", help="ccd to use. Range of ccds can be given (1-10)")
parser.add_argument("-a", "--amps", required=False, dest="amplifier", default="A,B", help="Amplifier to process (A or B). Default is to process both amps (A,B)")
parser.add_argument("-sp", "--saveplots", required=False, action="store_true", help="Flag to save plots. If given, then all plots will be saved")
parser.add_argument("-pa", "--plotall", required=False, action="store_true", help="Flag to plot  data for all Range of dates included")
parser.add_argument("--section", '-s', required=False, dest="section", default="db-desoper", help='section of .desservices file with connection info (db-desoper, db-destest)')
parser.add_argument("-t", "--flatstag", required=False, action="store_true", help="If given, only a query to the flats_stats_qa_y2 table will be done and update all exposures FLAT_STATUS column with corresponding tag (VETO, GOOD)")
parser.add_argument("-v", "--verbose", required=False,action="store_true",help="Use --verbose for verbose output")


args = parser.parse_args()
dateRange = args.dates
band = args.band
ccdRange = args.ccdrange
section = args.section
tosave = args.saveplots
plotall = args.plotall
theamp = args.amplifier
verbose = args.verbose
flatsTag = args.flatstag

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

    dbh = coreutils.desdbi.DesDbi(desdmfile,dbase)
        
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

def medsigclip(indata, clipsig=3.0, maxiter=10, converge_num=0.001, verbose=verbose):
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
   >>> median, mean, sigma = meanclip(indata)
 
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
   
   if verbose:
       prf = 'MEDSIGCLIP:'
       print '%s %.1f-sigma clipped mean' % (prf, clipsig)
       print '%s Mean computed in %i iterations' % (prf, iter)
       print '%s , Median = %.6f, Mean = %.6f, sigma = %.6f' % (prf, medval, mean, sigma)
 
   return median, mean, sigma

def saveMultiExpCcd(counter, axarr, x, medians, x_new, y_new, expnum, R2, chi2red, yerr, tosave=False, figurename = None, mytitle=None, flaglast = None):

    #print len(x), len(medians), len(yerr)
    start = time.time()
    #print "counter: ", counter
    # Four axes, returned as a 2-d array
    if counter == 0:
        axarr[0, 0].plot(x, medians, 'ro')
        axarr[0, 0].errorbar(x, medians, yerr=yerr, fmt = 'ro', ecolor='r', capsize=0)
        axarr[0, 0].plot(x_new, y_new, '-')
        axarr[0, 0].set_title('$%d  [R^2=%.2f,\, \chi^2_{red}=%.3f]$' % (expnum, R2, chi2red), fontsize=11, color=[0,0,0])
    elif counter == 1:
        axarr[0, 1].plot(x, medians, 'ro')
        axarr[0, 1].errorbar(x, medians, yerr=yerr, fmt = 'ro', ecolor='r', capsize=0)
        axarr[0, 1].plot(x_new, y_new, '-')
        #axarr[0, 1].set_title(expnum)
        axarr[0, 1].set_title('$%d  [R^2=%.2f,\, \chi^2_{red}=%.3f]$' % (expnum, R2, chi2red), fontsize=11, color=[0,0,0])
    elif counter == 2:
        axarr[1, 0].plot(x, medians, 'ro')
        axarr[1, 0].errorbar(x, medians, yerr=yerr, fmt = 'ro', ecolor='r', capsize=0)
        axarr[1, 0].plot(x_new, y_new, '-')
        #axarr[1, 0].set_title(expnum)
        axarr[1, 0].set_title('$%d  [R^2=%.2f,\, \chi^2_{red}=%.3f]$' % (expnum, R2, chi2red), fontsize=11, color=[0,0,0])
    else:
        axarr[1, 1].plot(x, medians, 'ro')
        axarr[1, 1].errorbar(x, medians, yerr=yerr, fmt = 'ro', ecolor='r', capsize=0)
        axarr[1, 1].plot(x_new, y_new, '-')
        #axarr[1, 1].set_title(expnum)
        axarr[1, 1].set_title('$%d  [R^2=%.2f,\, \chi^2_{red}=%.3f]$' % (expnum, R2, chi2red), fontsize=11, color=[0,0,0])
        # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
        plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
        plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
    
    
    if counter == 3 or flaglast == True:
        if tosave:
            plt.suptitle(mytitle)
            plt.savefig(figurename)
            plt.close()
        #else:
        #    plt.show()

    end = time.time()
    
    if figurename is not None:
        print "time to Create %s is %f" % (figurename, end-start) 

    #fig = plt.figure()
    #plt.rc("font", size=11)
    #ax1 = fig.add_subplot(222)
    #ax = plt.subplot(111)
    #ax.set_title(rootExposureName)
    #ax1.set_xlabel('exposureid')
    #ax1.set_ylabel('fit values 2 order polynomial')
    #ax2 = fig.a
 
def saveSingleExpCcd(ouptutFilename, expnum, x, medians, x_new, y_new, tosave=False):

        
    plt.figure()
    #figfit, ax5 = plt.subplots(1)    
    plt.plot(x_new, y_new, '-')
    plt.plot(x, medians, 'ro')
    plt.title(expnum) 
    if tosave:
        plt.savefig(ouptutFilename)
        plt.close()
    #else:
    #    plt.show()


def filenameToExpnum(filenames):
    """Extract the integer part of the filename, which is the expnum"""
    
    intexpnums = []
    for name in filenames:
        #print name
        digit = re.match('DECam_(\d+).fits',name)
        intDigit = int(digit.group(1))
        intexpnums.append(intDigit)
    
    return intexpnums

#def MyMultiPlot(expnum, nite, ccd, zero, first, second, third, fourth, tosave=False, figureName = None):
def MyMultiPlot(filenames, ccdnum, allfits, nite = None, tosave=False, figureName = None):

    start = time.time()
    
    f, (ax, ax1, ax2, ax3, ax4) = plt.subplots(5, sharex=True, sharey=False)

    #plt.title('20131227')

    #print resultFits
    zero =[]
    first = []
    second = []
    third = []
    fourth = []
    
    for fits in allfits:
        zero.append(fits[0])
        first.append(fits[1])
        second.append(fits[2])
        third.append(fits[3])
        fourth.append(fits[4])

    intexpnums = filenameToExpnum(filenames)
            
    minx = min(intexpnums)
    newexp = []
    for exp in intexpnums:
        val = exp - minx
        newexp.append(val)
    
    
    
    for x,z,fi,s,t,u in itertools.izip_longest(newexp, zero, first, second, third, fourth):
        ax.grid(), ax1.grid(), ax2.grid(), ax3.grid(), ax4.grid()
        ax.plot(x,z, 'ro')
        if nite is None:
            myTitle = 'range of nites'
        else:
            myTitle = 'Night = ' + nite + ', ccd = ' + ccdnum
        ax.set_title(myTitle)
        ax1.plot(x,fi, 'bs')
        ax2.plot(x,s, 'g^')
        ax3.plot(x,t, 'ro')
        ax4.plot(x,u, 'go')

        #ax.set_ylabel
    f.subplots_adjust(hspace=0, bottom = 0.1)
    xlabel = 'EXPNUM + ' + str(minx)
    ax4.set_xlabel(xlabel, fontsize=9)


    #plt.xticks(rotation=45)
    plt.rc("font", size=9)
    
    if tosave and figureName is not None:
        plt.savefig(figureName)
        plt.close()
    #else:
    #    plt.show()

    end = time.time()
    print "Time to create %s is %f" % (figureName, end-start)
    
def MySinglePlotRatios(allfits, ratios, nite, xlimits = None, title = None, tosave=False, figureName = None):


    start= time.time()
    
    #print resultFits
    zero =[]
    first = []
    second = []
    third = []
    fourth = []
    
    for fits in allfits:
        zero.append(fits[0])
        first.append(fits[1])
        second.append(fits[2])
        third.append(fits[3])
        fourth.append(fits[4])

    plt.figure(1)
    plt.clf()
    plt.rc("font", size=10)
    #plt.grid()
    #plt.xticks(rotation=45)
    plt.subplots_adjust(bottom=0.1)
    plt.xlabel('Counts at first Bin Mean y=[300:500] (fit)')
    plt.ylabel('Ratio y = 300:500 / 3800:4000 bins')
    if xlimits is not None:
        plt.xlim(xlimits[0], xlimits[1])
    if title is not None:
        theTitle = title
        plt.title(theTitle)
    else:
        theTitle = 'Ratios vs Mean Counts ' + nite
        plt.title(theTitle)

    for e,s in itertools.izip_longest(zero, ratios):
        plt.plot(e,s,'ro')
        plt.grid()
    if tosave and figureName is not None:
        plt.savefig(figureName)
        plt.close()
    #else:
    #    plt.show()
    
    end = time.time()
    
    print "Time to create %s is %f" % (figureName, end-start)
        
def MySinglePlotExptimes(exptimes, ratios, nite, xlimits = None, title = None, tosave=False, figureName = None):


    start =  time.time()
    
    plt.figure(1)
    plt.clf()
    plt.rc("font", size=10)
    plt.grid()
    #plt.xticks(rotation=45)
    plt.subplots_adjust(bottom=0.1)
    if xlimits is not None:
        plt.xlim(xlimits[0], xlimits[1])
    if title is not None:
        theTitle = title
        plt.title(theTitle)
    else:
        theTitle = 'Ratios vs Exposure Time ' + nite 
        plt.title(theTitle)

    plt.xlabel('Exptime (s)', fontsize=9)
    plt.ylabel('Ratio y = 300:500 / 3800:4000 bins')
    for e,s in itertools.izip_longest(exptimes, ratios):
        plt.plot(e,s,'ro')
    if tosave and figureName is not None:
        plt.savefig(figureName)
        plt.close()
    #else:
    #    plt.show()

    end = time.time()
    
    print "Time to create %s is %f" % (figureName, end-start)

        
def MySinglePlotSteps(filenames, allsteps, nite, tosave=False, figureName = None, title=None):
    

    start = time.time()
    
    #transform filename to expnums
    intexpnums = filenameToExpnum(filenames)
      
    #print expnums, allsteps
    minx = min(intexpnums)
    newexp = []
    for i,exp in enumerate(intexpnums):
        num = len(allsteps[i])
        #print "num = ", num
        val = exp - minx
        #print val
        x=[val+e for e in range(num)]
        #print x
        newexp.append(x)

    #print len(newexp), len(allsteps), len(newexp[0]), len(allsteps[0])
    plt.figure(1)
    plt.clf()
    #plt.xticks(rotation=45)
    plt.grid()
    plt.rc("font", size=9)
    plt.subplots_adjust(bottom=0.1)
    for s,c in itertools.izip_longest(newexp,allsteps):
        plt.plot(s,c,'bo')
    if title is not None:
        theTitle = title
    else:
        theTitle = 'Ratios of adjacent bin ' + nite
    plt.title(theTitle)
    xlabel = 'EXPNUM + ' + str(minx)
    plt.xlabel(xlabel, fontsize=9)
    plt.ylabel('Ratio of adjacent bins in Y direction', fontsize=9)

    if tosave and figureName is not None:
        plt.savefig(figureName)
        plt.close()
    #else:
    #    plt.show()

    end = time.time()
    
    print "Time to create %s is %f" % (figureName, end-start)

def plotExpnumRatios(filenames, allratios, nite, tosave=False, figureName = None, title=None):
    

    start = time.time()

    #transform filename to expnums    
    intexpnums = filenameToExpnum(filenames)
    
    minx = min(intexpnums)
    newexp = []
    for exp in intexpnums:
        val = exp - minx
        newexp.append(val)
      
    #print len(newexp), len(allsteps), len(newexp[0]), len(allsteps[0])
    plt.figure(1)
    plt.clf()
    #plt.xticks(rotation=45)
    plt.grid()
    plt.rc("font", size=9)
    plt.subplots_adjust(bottom=0.1)
    for s,c in itertools.izip_longest(newexp,allratios):
        plt.plot(s,c,'bo')
    if title is not None:
        theTitle = title
    else:
        theTitle = 'Ratios for each bin ' + nite
    plt.title(theTitle)
    xlabel = 'EXPNUM + ' + str(minx)
    plt.xlabel(xlabel, fontsize=9)
    plt.ylabel('Ratio y = 300:500 / 3800:4000 bin', fontsize=9)

    if tosave and figureName is not None:
        plt.savefig(figureName)
        plt.close()
    #else:
    #    plt.show()

    end = time.time()
    
    print "Time to create %s is %f" % (figureName, end-start)


def dontneed():
    #Fit data for each ccd
    for i in range(0,len(medians), 4):    
        f, axarr = plt.subplots(2, 2)
        for j in range(4):
            #print nites[i+j], expnums[i+j], i, j, i+j, nites[i]
            if i+j >= len(filenames):
                print "i+j= %d = len(filenames %d" % (i+j, len(filenames))
                continue
            elif nites[i+j] != nites[i]:
                print "nites", nites[i+j], nites[i]
                name_to_save = nites[i] + '_' + str(filenames[i]) + '_' + str(filenames[i+j]) + '_' + ccdnums[i+j] + '_' + band + '.png'
                title = str(nites[i]) + " ccd=" + ccdnums[i] + band
                plt.suptitle(title)
                plt.savefig(name_to_save)
                plt.close()
                break
            #    print "nite change %s %s" % (nites[i+j], nites[i])
            #    name_to_save = nites[i] + '_' + str(expnums[i]) + '_' + str(expnums[i+j]) + '_' + ccdnums[i+j] + '_' + band + '.png'
            #    title = str(nites[i]) + " ccd=" + ccdnums[i]
            #    flagsavelast = True 
            #    saveMultiExpCcd(j, axarr, x, medians[i+j], x_new, y_new, expnums[i+j], tosave=True, figurename = name_to_save, mytitle = title, flaglast = flagsavelast)
            else:
                #print expnums[j+i]
                
                xdata = len(medians[i+j]) + 1
                
                x = np.array(range(1,xdata))
                
                """
 
                #NOTE: The documentation on numpy recomends to not use np.polyfit, np.polyval, np.poly1d 
                fit = np.polyfit(x,medians[i+j],4)
                
                #this calculations are form 
                #http://nbviewer.ipython.org/github/duartexyz/BMC/blob/master/CurveFitting.ipynb
                yfit_new = np.polyval(fit, x)   # Evaluate polynomial at x
                R2 = np.corrcoef(x, medians[i+j])[0, 1]**2  #coefficient of determination between x and y
                resid = medians[i+j] - yfit_new #residuals
                chi2red = np.sum((resid/yerr[i+j])**2)/(medians[i+j].size - 2) #chisrq using stdev for each median bin
                print "chisqr = ", chi2red

                
                """
                #It is better to use the following:
                coefs = poly.polyfit(x, medians[i+j], 4)
                yfit_new = poly.polyval(x,coefs)
                R2 = np.corrcoef(x, medians[i+j])[0,1]**2 #coefficient of determination between x and y
                resid = medians[i+j] - yfit_new  #residuals
                chi2red = np.sum((resid/yerr[i+j])**2)/(medians[i+j].size - 2)
                #print "chisrq", chi2red
                
                #fitted = np.poly1d(fit) # Evaluate polynomial at x Using poly1d instead of polyval np function
                
                #calculate new x's and y's
                x_new = np.linspace(x[0], x[-1], 50)
                
                #using results from poly1d (both results for evaluation are from the resutlsfrom polyfit)
                #what it changed it poly1d or polyvall to evaluate data at certain points 
                #y_new = fitted(x_new)
                #y_polyval_new = np.polyval(fit,x_new)
                yfit_polyval_new = poly.polyval(x_new, coefs)

                #print expnums[i+j], R2, chi2red
                #if R2 == 1 and chi2red < 0.1:
                #    flatStatus[expnums[i+j]] = (ccdnums[i], 'VETO SHUTTER')
                #    print expnums[i+j],expnums[i], i,j,i+j,'VETO SHUTTER'
                #elif R2 < 1 and chi2red > 1:
                #    flatStatus[expnums[i+j]] = (ccdnums[i], 'VETO CHI2')
                #    print expnums[i+j],expnums[i],i,j,i+j,'VETO CHI2'
                #elif R2 < 1 and chi2red < 0.1:
                #    flatStatus[expnums[i+j]] = (ccdnums[i], 'GOOD')
                #    print expnums[i+j], expnums[i],i,j,i+j,'GOOD'

                
def insert_dictionary_2Db(dbh, query, dictionary, verbose=verbose):
    """Execute a query and return a cursor to a query
    :param dbh: Database connection handler
    :param query: query to execute
    :param debug: verbosity
    
    """
    
    try:      
        if verbose:
            print query
        #dbh = connectDB('db-desoper')
        cur = dbh.cursor()        
        cur.execute(query,dictionary)
        #dbh.commit()
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
    
    print "finished updateting..."

def update_flatstatus_2Db(dbh, query, verbose=verbose):
    """Execute a query and return a cursor to a query
    :param dbh: Database connection handler
    :param query: query to execute
    :param debug: verbosity
    
    """    

    try:      
        #dbh = connectDB('db-desoper')
        cur = dbh.cursor()
        cur.execute(query)
        #dbh.commit()
        if verbose:
            print query
    except cx_Oracle.IntegrityError as e:
        print "error while inserting into FLATS_STATS_QA_Y2 TABLE: ", e  

    if verbose:
        print "finished updating..."


def des_saturate(ccd, amp):


    saturateA = {}
    saturateB = {}
        
    saturateA['1'] =  38652.0 
    saturateA['2'] =  41945.0
    saturateA['3'] =  42449.0 
    saturateA['4'] =  42059.0 
    saturateA['5'] =  45237.0 
    saturateA['6'] =  50304.0 
    saturateA['7'] =  43219.0 
    saturateA['8'] =  50588.0
    saturateA['9'] =  46442.0 
    saturateA['10'] =  47817.0 
    saturateA['11'] =  44428.0 
    saturateA['12'] =  47000.0 
    saturateA['13'] =  43172.0 
    saturateA['14'] =  47577.0 
    saturateA['15'] =  46000.0 
    saturateA['16'] =  44000.0 
    saturateA['17'] =  48000.0 
    saturateA['18'] =  51319.0
    saturateA['19'] =  40534.0 
    saturateA['20'] =  43397.0
    saturateA['21'] =  47000.0 
    saturateA['22'] =  46380.0
    saturateA['23'] =  44823.0
    saturateA['24'] =  48000.0 
    saturateA['25'] =  40469.0 
    saturateA['26'] =  47435.0 
    saturateA['27'] =  33874.0 
    saturateA['28'] =  45246.0 
    saturateA['29'] =  46671.0
    saturateA['30'] =  47461.0
    saturateA['31'] =  37951.0
    saturateA['32'] =  50246.0
    saturateA['33'] =  44938.0 
    saturateA['34'] =  50010.0
    saturateA['35'] =  44144.0 
    saturateA['36'] =  45844.0 
    saturateA['37'] =  46251.0
    saturateA['38'] =  40983.0 
    saturateA['39'] =  44177.0
    saturateA['40'] =  51000.0

    saturateB['1'] =  37531.0 
    saturateB['2'] =  44078.0
    saturateB['3'] =  42345.0 
    saturateB['4'] =  42605.0
    saturateB['5'] =  41000.0
    saturateB['6'] =  43534.0 
    saturateB['7'] =  42163.0 
    saturateB['8'] =  48732.0 
    saturateB['9'] =  45402.0 
    saturateB['10'] =  49049.0
    saturateB['11'] =  46770.0 
    saturateB['12'] =  45537.0
    saturateB['13'] =  42871.0 
    saturateB['14'] =  44265.0 
    saturateB['15'] =  46370.0 
    saturateB['16'] =  44467.0 
    saturateB['17'] =  42982.0 
    saturateB['18'] =  49231.0 
    saturateB['19'] =  39960.0 
    saturateB['20'] =  42551.0 
    saturateB['21'] =  43939.0 
    saturateB['22'] =  46334.0 
    saturateB['23'] =  47770.0 
    saturateB['24'] =  50504.0 
    saturateB['25'] =  42151.0 
    saturateB['26'] =  44650.0 
    saturateB['27'] =  28164.0 
    saturateB['28'] =  42626.0 
    saturateB['29'] =  42720.0 
    saturateB['30'] =  44346.0 
    saturateB['31'] =  45307.0 
    saturateB['32'] =  45586.0 
    saturateB['33'] =  42834.0 
    saturateB['34'] =  49000.0 
    saturateB['35'] =  43448.0 
    saturateB['36'] =  43798.0 
    saturateB['37'] =  39364.0 
    saturateB['38'] =  42716.0 
    saturateB['39'] =  41854.0 
    saturateB['40'] =  48278.0 
    saturateA['41'] =  46390.0
    saturateB['41'] =  39973.0 
    saturateA['42'] =  45678.0 
    saturateB['42'] =  46060.0 
    saturateA['43'] =  39602.0 
    saturateB['43'] =  40252.0 
    saturateA['44'] =  46163.0
    saturateB['44'] =  43511.0 
    saturateA['45'] =  48694.0 
    saturateB['45'] =  42599.0 
    saturateA['46'] =  38696.0 
    saturateB['46'] =  37655.0 
    saturateA['47'] =  46559.0
    saturateB['47'] =  45650.0 
    saturateA['48'] =  41302.0 
    saturateB['48'] =  38771.0
    saturateA['49'] =  49133.0 
    saturateB['49'] =  42883.0 
    saturateA['50'] =  46642.0 
    saturateB['50'] =  40239.0 
    saturateA['51'] =  46415.0
    saturateB['51'] =  39981.0 
    saturateA['52'] =  40268.0 
    saturateB['52'] =  42913.0 
    saturateA['53'] =  40698.0 
    saturateB['53'] =  42397.0 
    saturateA['54'] =  50306.0 
    saturateB['54'] =  41791.0 
    saturateA['55'] =  46959.0 
    saturateB['55'] =  43185.0 
    saturateA['56'] =  49415.0
    saturateB['56'] =  42401.0 
    saturateA['57'] =  36677.0 
    saturateB['57'] =  33214.0 
    saturateA['58'] =  43135.0 
    saturateB['58'] =  41106.0 
    saturateA['59'] =  41544.0 
    saturateB['59'] =  42196.0 
    saturateA['60'] =  43295.0 
    saturateB['60'] =  41383.0
    saturateA['62'] =  47995.0
    saturateB['62'] =  43079.0 


    #print ccd
    if amp == 'A':
        value = saturateA[ccd]
        return value
    else:
        value = saturateB[ccd]
        return value


def calculateallvetoflats(dbh, band, dateRange, theccd, section, amp, tosave=False, plotall=False, verbose=False):

    """
    Calculate the statistics using the data in the database produced by check_flats.py
    Data is in table flats_stats_qa_y2
    :param band: Single band to process
    :param dateRange: Dante of ranges to process
    :param theccd: One single ccd to process
    """
    
    ############################################################
    #Cuts using fitting flat data:
    #
    # Anomalous flats:
    # Rsq = 0.99 and chi2red < 0.1  (this are flats that have a problem with the shutter)
    #
    # Rsq < 1 and chi2red > 1 (this are flats that have a problems somewhere in the fitting). They are not like a usual flat.
    #
    # Rsq < 0.94 and chi2red > 0.01 (and mean > 55,000 problem with saturation or exptime)
    #
    # Median of the flat is above 50,000 counts (or I could use the SATURATE keyword from the header as my limit)
    #
    # Ratio between consecutive bins. I call these steps.
    # 
    # Band     min      max (inclusive acceptable range)
    # Y        0.98    1.02
    # g,r,i,z  0.99    1.01
    #
    # any flat with a step outside these values will be vetoed
    #
    # Outer edge ratio (ratio between first and last bin along the long the readout axis)
    #
    # Band             min      max
    # g,r,i,z,Y        0.8      1.2
    # 
    # Exposure times. Typical Exposure times for each filter are:
    #
    # Band    exptime    counts (ADU)
    # g        30        16,000 - 18,000
    # r        10        20,000
    # i        22        18,000 - 20,000
    # z        10        21,000 - 24,000
    # Y        10        20,000 - 23,000
    #
    # How can we check if a flats is being taken for linearity tests?
    # When that happens, we look at the COMMENT from the header where it says Linearity.
    # 
    # If Mean counts of a ccd/Amp > 40,000, then Flag it as VETO SATURATE
    #
    ###############################################################
    #
    #
    # Tag bad flats as VETO (Didn't pass the cuts) or GOOD (pass the cuts), PTC, Linearity. A requirement is that
    # if one ccd is bad, then the whole exposure is declared bad.
    # 
    # The column for the Flat Status is FLATSTATUS 
    # Values: VETO, GOOD
    #
    #

    
    #tables to use:
    #ricardoc.exposure_regions_qa_y2 
    #ricardoc.exposure_tel_qa_y2
    #ricardoc.flats_stats_qa_y2
    #bands = band.split(',')
    #print ccdRange
    
    #print ccdRange.split('-')

    
    #try:
    #    ccdRange.split('-')
    #    firstCcd, lastCcd = ccdRange.split('-')
    #    ccdsToProcess = tuple(range(firstCcd,lastCcd+1))
    #    if 61 in ccdsToProcess:
    #        ccdsToProcess[61]
    #    print "ccds to process", ccdsToProcess
    #    flagCcdRange = True
    #except:
    #    ccdsToProcess = ccdRange
    #    flagCcdRange = False
    
    
    #print firstNight, lastNight
    
    #DOMEFLOR, UTE_TEMP, UTN_TEMP, PME_TEMP, DOMEHIGH, PMN_TEMP, DOMELOW, PRESSURE, UTS_TEMP,
    
    #domePosition = """select fs.MEDIAN_GBOX9_A, et.PRESSURE, fs.nite from flats_stats_qa fs, exposure_tel_qa et, exposure e
    #                where fs.EXPOSUREID = et.EXPOSUREID and fs.EXPOSUREID = e.id and e.id = fs.exposureid and fs.CCD = '10' and fs.nite >= '%s' and 
    #                fs.nite <= '%s' and fs.band = '%s' and e.exptime = 10 order by fs.nite""" % (firstNight, lastNight, band)
                    
    
    
    
    median_values_A = ['distinct(fs.FILENAME)','fs.nite', 'fs.ccdnum', 'e.exptime','fs.MEDIAN_GBOX1_A','fs.MEDIAN_GBOX2_A','fs.MEDIAN_GBOX3_A','fs.MEDIAN_GBOX4_A','fs.MEDIAN_GBOX5_A','fs.MEDIAN_GBOX6_A',
                     'fs.MEDIAN_GBOX7_A','fs.MEDIAN_GBOX8_A','fs.MEDIAN_GBOX9_A','fs.MEDIAN_GBOX10_A','fs.MEDIAN_GBOX11_A','fs.MEDIAN_GBOX12_A',
                     'fs.MEDIAN_GBOX13_A','fs.MEDIAN_GBOX14_A','fs.MEDIAN_GBOX15_A','fs.MEDIAN_GBOX16_A','fs.MEDIAN_GBOX17_A','fs.MEDIAN_GBOX18_A']

    stdev_values_A = ['distinct(fs.FILENAME)','fs.nite', 'fs.ccdnum', 'e.exptime','fs.STDEV_GBOX1_A','fs.STDEV_GBOX2_A','fs.STDEV_GBOX3_A','fs.STDEV_GBOX4_A','fs.STDEV_GBOX5_A','fs.STDEV_GBOX6_A',
                     'fs.STDEV_GBOX7_A','fs.STDEV_GBOX8_A','fs.STDEV_GBOX9_A','fs.STDEV_GBOX10_A','fs.STDEV_GBOX11_A','fs.STDEV_GBOX12_A',
                     'fs.STDEV_GBOX13_A','fs.STDEV_GBOX14_A','fs.STDEV_GBOX15_A','fs.STDEV_GBOX16_A','fs.STDEV_GBOX17_A','fs.STDEV_GBOX18_A']

    median_values = ['distinct(fs.FILENAME)','fs.nite', 'fs.ccdnum', 'e.exptime','fs.MEDIAN_GBOX1_','fs.MEDIAN_GBOX2_','fs.MEDIAN_GBOX3_','fs.MEDIAN_GBOX4_','fs.MEDIAN_GBOX5_','fs.MEDIAN_GBOX6_',
                     'fs.MEDIAN_GBOX7_','fs.MEDIAN_GBOX8_','fs.MEDIAN_GBOX9_','fs.MEDIAN_GBOX10_','fs.MEDIAN_GBOX11_','fs.MEDIAN_GBOX12_',
                     'fs.MEDIAN_GBOX13_','fs.MEDIAN_GBOX14_','fs.MEDIAN_GBOX15_','fs.MEDIAN_GBOX16_','fs.MEDIAN_GBOX17_','fs.MEDIAN_GBOX18_']

    stdev_values = ['distinct(fs.FILENAME)','fs.nite', 'fs.ccdnum', 'e.exptime','fs.STDEV_GBOX1_','fs.STDEV_GBOX2_','fs.STDEV_GBOX3_','fs.STDEV_GBOX4_','fs.STDEV_GBOX5_','fs.STDEV_GBOX6_',
                     'fs.STDEV_GBOX7_','fs.STDEV_GBOX8_','fs.STDEV_GBOX9_','fs.STDEV_GBOX10_','fs.STDEV_GBOX11_','fs.STDEV_GBOX12_',
                     'fs.STDEV_GBOX13_','fs.STDEV_GBOX14_','fs.STDEV_GBOX15_','fs.STDEV_GBOX16_','fs.STDEV_GBOX17_','fs.STDEV_GBOX18_']
    
    #Use to append amplifier name for the DB query
    median_values_amp = []
    for i,val in enumerate(median_values):
        if i <= 3: #Do not append the amp value to the 3 first elements of the list
            median_values_amp.append(val)

        else:
            boxName = val+amp #append the amp value on the rest elements
            median_values_amp.append(boxName)

    stdev_values_amp = []
    for i,val in enumerate(stdev_values):
        if i <= 3: #Do not append the amp value to the 4 first elements of the list
            stdev_values_amp.append(val)
        else:
            boxName = val+amp #append the amp value on the rest elements
            stdev_values_amp.append(boxName)

   
    #print "XXXXXX", median_values_amp
    
    
    query_median_values = ",".join(median_values_amp)
    query_stdev_values = ",".join(stdev_values_amp)
    
    #Assume the code is only processing one ccd
    flagCcdRange = False
 
    #Check if given a range of dates or a single date
    try: 
        #dateRange.split()
        firstNight, lastNight = dateRange.split('-')
        flagDateRange = True
        print "DateRange", dateRange
    except:
        nite = dateRange
        flagDateRange = False
        firstNight = nite
        lastNight = nite
        print "only One night provided", nite


    
    if flagDateRange:
        if flagCcdRange:
            if band == 'r':
                query_median = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum in %s 
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where band = '%s' and propid not like '%%2013A-9999%%' and exptime = 15 
                                    and nite >= '%s' and nite <= '%s') order by fs.filename"name""" % (query_median_values, theccd, band, firstNight, lastNight)

                query_stdev = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum in %s 
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where band = '%s' and propid not like '%%2013A-9999%%' and exptime = 15 
                                    and nite >= '%s' and nite <= '%s') order by fs.filename""" % (query_stdev_values, theccd, band, firstNight, lastNight)

            
            else:
                query_median = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s 
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where band = '%s' and propid not like '%%2013A-9999%%' and
                                    nite >= '%s' and nite <= '%s') order by fs.filename""" % (query_median_values, theccd, band, firstNight, lastNight)

                query_stdev = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s 
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where band = '%s' and propid not like '%%2013A-9999%%' and
                                    nite >= '%s' and nite <= '%s') order by fs.filename""" % (query_stdev_values, theccd, band, firstNight, lastNight)

        else:            
            if band == 'r':
                query_median = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s 
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where band = '%s' and propid not like '%%2013A-9999%%' and exptime = 15
                                    and nite >= '%s' and nite <= '%s') order by fs.filename""" % (query_median_values, theccd, band, firstNight, lastNight)

                query_stdev = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s 
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where band = '%s' and propid not like '%%2013A-9999%%' and exptime = 15 
                                    and nite >= '%s' and nite <= '%s') order by fs.filename""" % (query_stdev_values, theccd, band, firstNight, lastNight)

            else:
                query_median = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where band = '%s' and propid not like '%%2013A-9999%%' and
                                    nite >= '%s' and nite <= '%s') order by fs.filename""" % (query_median_values, theccd, band, firstNight, lastNight)

                query_stdev = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where band = '%s' and propid not like '%%2013A-9999%%' and
                                    nite >= '%s' and nite <= '%s') order by fs.filename""" % (query_stdev_values, theccd, band, firstNight, lastNight)
        
    else:
        if flagCcdRange:
            if band == 'r':
                query_median = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s 
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where nite = '%s' and band = '%s' and propid not like '%%2013A-9999%%' and 
                                    exptime = 15) order by fs.filename""" % (query_median_values_A, theccd,nite, band)

                query_stdev = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s 
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where nite = '%s' and band = '%s' and propid not like '%%2013A-9999%%' and 
                                    exptime = 15) order by fs.filename""" % (query_stdev_values_A, theccd,nite, band)

            else:
                query_median = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s 
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where nite = '%s' and band = '%s' and propid not like '%%2013A-9999%%') order by fs.filename""" % (query_median_values, theccd, nite, band)

                query_stdev = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s 
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where nite = '%s' and band = '%s' and propid not like '%%2013A-9999%%') order by fs.filename""" % (query_stdev_values, theccd, nite, band)

        
        else:     
            if band == 'r':
                query_median = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s 
                                    and fs.filename = e.filename and e.filename in (select filename from exposure 
                                    where nite = '%s' and band = '%s' and propid not like '%%2013A-9999%%' and 
                                    exptime = 15) order by fs.filename""" % (query_median_values, theccd, nite, band)

                query_stdev = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s 
                                    and fs.filename  = e.filename and e.filename in (select filename from exposure 
                                    where nite = '%s' and band = '%s' and propid not like '%%2013A-9999%%' and 
                                    exptime = 15) order by fs.filename""" % (query_stdev_values, theccd, nite, band)

            
            else:
                query_median = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s
                                    and fs.filename  = e.filename and e.filename in (select filename from exposure 
                                    where nite = '%s' and band = '%s' and propid not like '%%2013A-9999%%') order by fs.filename""" % (query_median_values, theccd, nite, band)
    
                query_stdev = """select %s from flats_stats_qa_y2 fs, exposure e where fs.ccdnum = %s
                                    and fs.filename  = e.filename and e.filename in (select filename from exposure 
                                    where nite = '%s' and band = '%s' and propid not like '%%2013A-9999%%') order by fs.filename""" % (query_stdev_values, theccd, nite, band)
    
    
    #print query_median_A
    #print query_stdev_A
    
    #print "query is:", query_median
    median_values = query_to_cur(dbh, query_median, verbose=False)
    stdev_values = query_to_cur(dbh, query_stdev, verbose=False)
    
    rowcounts = query_to_cur(dbh, query_median)
    
    for j, val in enumerate(rowcounts):
        j=j+1
        rows = j
        #print rows


    #print "rows",rows
        
    filenames = []
    nites = []
    ccdnums = []
    exptimes = []
    medians = []
    
    
    for i, val in enumerate(median_values,start=1):
        filenames.append(val[0])
        nites.append(val[1])
        ccdnums.append(str(val[2]))
        exptimes.append(val[3])
        medians.append(np.array(val[4:]))

    #stdev is Error for each median bin values
    yerr = [] #standard deviation calculated for each median bin
    for m, value in enumerate(stdev_values,start=1):
        #print len(value[4:])
        yerr.append(np.array(value[4:]))


    if len(filenames) != len(nites) and len(nites) != len(ccdnums) and len(ccdnums) != len(exptimes) and len(exptimes) != len(medians):
        sys.exit("missing data, can't do calculations...")
        
    
    
    #Dictionary with expnum and status
    #each dictionary created corresponds to single band, ccd and Amp and of course a date.
    flatStatus = {}
    insertDictionaryFitsDb = {}

    #how many distinct night have been provided to process     
    uniqueNights = set(nites)
    #uniqueNights = uniqueNights.sort()
    position = 0

    #print "unique nights:",uniqueNights 
       
    #Find the indexes of the exposures that correspond to each night
    for singleNite in uniqueNights:
        #Index position for all elements (filenames, nites, ccdnums,exptimes,medians) for a given unique night
        #It should be of a size = number of exposures in a given nite.
        thisNiteIndex = [i for i, x in enumerate(nites) if x == singleNite] 
        

        #for j in thisNiteIndex:
        #    print nites[j], expnums[j],ccdnums[j]
        
        #Fit data to each exposure/ccd and make a plot of each exposure with the fitted function.
        #Remember that the call to the current function is made per exposure/ccd/band, so this loop
        #is for a single ccd/band and all the exposures found for a given range of dates (or a single date)
        for i,idx in enumerate(thisNiteIndex):
            #print expnums[idx], nite
            
            #Mean of the ccd/amp using all the bins
            mean_ccd_amp = np.mean(medians[idx])
            
            xdata = len(medians[idx]) + 1   #medians is a np array of median values calculated for each exposure/ccd/band combination.
                                            #It contains the median values of each amp/ccd. 
            x = np.array(range(1,xdata))    #x array to calculate fit (it of size 18 in this case)
            #It is better to use the following:
            coefs, stats = poly.polyfit(x, medians[idx], 4, full=True)  #Fit fourth order polynomial
            yfitpredicted = poly.polyval(x,coefs)   #predicted values using fitting coefficients
            dev = medians[idx] - np.mean(medians[idx])  #deviation (measure of spread)
            SST = np.sum(dev**2)    #total variation to be accounted for
            resid = medians[idx] - yfitpredicted  #residuals
            SSE = np.sum(resid**2)  #Variation not accounted for (Total variation of the residuals)
            Rsq = 1 - SSE/SST   #Total percent of error (The closer to 1, the better the fit is to the data) 
            
            
            #stats contains a list with results in the following order: [residuals, rank, singular_values, rcond]
            #Sum of the squared residuals (SSR) of the least-squares fit (Want this small for a good fit)
            #the effective rank of the scaled Vandermonde matrix
            #it's singular values
            #the specified value of rcond (not given in my case)
            SSR = stats[0][0]
            rank = stats[1]
            singularValue = stats[2]
            
            #corrcoef = np.corrcoef(x, medians[idx])
            #print "corrcoef:", corrcoef
            
            #print len(resid), len(yerr[idx]), len(medians[idx]), medians[idx].size
            #R2 = np.corrcoef(x, medians[idx])[0,1]**2 #coefficient of determination between x and y
            chi2red = np.sum((resid/yerr[idx])**2)/(medians[idx].size -1 - 4)  #reduced chi squared. polynomial of order 4th, 4 degree of freedom

            #print "rsq= ", Rsq, "chi2red= ", chi2red

            x_new = np.linspace(x[0], x[-1], 50)
            
            #using results from poly1d (both results for evaluation are from the resutlsfrom polyfit)
            #what it changed it poly1d or polyvall to evaluate data at certain points 
            #y_new = fitted(x_new)
            #y_polyval_new = np.polyval(fit,x_new)
            yfit_polyval_new = poly.polyval(x_new, coefs)

            #print "expnum = %4d, R2 = %4.3f, chi2red = %4.4f" % (expnums[idx], R2, chi2red)
            #print SSR, rank, singularValue

            #ypred = polyval(coeff,x);   % predictions
            #dev = y - mean(y);          % deviations - measure of spread
            #SST = sum(dev.^2);          % total variation to be accounted for
            #resid = y - ypred;          % residuals - measure of mismatch
            #SSE = sum(resid.^2);        % variation NOT accounted for
            #Rsq = 1 - SSE/SST;          % percent of error explained


            
            if Rsq >= 0.98 and chi2red > 1:
                flatStatus[filenamess[idx]] = (ccdnums[idx], 'VETO CHI2')
                #flatStatus[expnums[idx]] = (ccdnums[idx], 'VETO')
                if verbose:
                    print filenames[idx],'VETO CHI2'
 
            if band == 'Y':
                if Rsq > 0.80 and chi2red < 0.2:
                    flatStatus[filenames[idx]] = (ccdnums[idx], 'GOOD')
                    if verbose:
                        print filenames[idx],'GOOD'
                else:
                    flatStatus[filenames[idx]] = (ccdnums[idx], 'VETO SATURATION')
                    if verbose:
                        print filenames[idx],'VETO SATURATION'                        
                
            if Rsq <= 0.95 and chi2red > 0.1 and band != 'Y':
                flatStatus[filenames[idx]] = (ccdnums[idx], 'VETO SATURATION')
                if verbose:
                    print filenames[idx],'VETO SATURATION'
                #flatStatus[expnums[idx]] = (ccdnums[idx], 'VETO')

            if mean_ccd_amp >= des_saturate(ccdnums[idx],amp):
                #print mean_ccd_amp, des_saturate(ccdnums[idx], amp), ccdnum[idx],amp
                flatStatus[filenames[idx]] = (ccdnums[idx], 'VETO SATURATION')
                if verbose:
                    print filenames[idx],'VETO SATURATION'
                #flatStatus[expnums[idx]] = (ccdnums[idx], 'VETO')
            if Rsq == 0 or chi2red == 0:
                flatStatus[filenames[idx]] = (ccdnums[idx], 'VETO SATURATION')
                if verbose:
                    print filenames[idx],'VETO SATURATION 3'                
            if Rsq > 0.99999 and chi2red > 0.009:
                flatStatus[filenames[idx]] = (ccdnums[idx], 'VETO SHUTTER')
                #flatStatus[expnums[idx]] = (ccdnums[idx], 'VETO')
            if Rsq < 0.99998 and chi2red < 0.05:
                #Rsq <= 0.9999 and chi2red <= 0.01:
                flatStatus[filenames[idx]] = (ccdnums[idx], 'GOOD')
                if verbose:
                    print filenames[idx],'GOOD'


            """    
            ricardoc.flats_stats_qa_y2
             FLAT_STATUS                             VARCHAR2(80)
             SSR_A                                  NUMBER(8,2)
             TOT_RSQ_B                             NUMBER(8,2)
             RED_CHI2_B                             NUMBER(8,2)
             BORDERS_RATIO_B                        NUMBER(8,2)
             STEPS_RATIO_B                             NUMBER(8,2)
             FIT_COEFF_ZERO_B                        NUMBER(8,2)
             FIT_COEFF_ONE_B                        NUMBER(8,2)
             FIT_COEFF_TWO_B                        NUMBER(8,2)
             FIT_COEFF_THREE_B                        NUMBER(8,2)
             FIT_COEFF_FOURTH_B                        NUMBER(8,2)
             MEAN_CCD_B                             NUMBER(8,2)
             SSR_B                                  NUMBER(8,2)
             TOT_RSQ_A                             NUMBER(8,2)
             RED_CHI2_A                             NUMBER(8,2)
             BORDERS_RATIO_A                        NUMBER(8,2)
             STEPS_RATIO_A                             NUMBER(8,2)
             FIT_COEFF_ZERO_A                        NUMBER(8,2)
             FIT_COEFF_ONE_A                        NUMBER(8,2)
             FIT_COEFF_TWO_A                        NUMBER(8,2)
             FIT_COEFF_THREE_A                        NUMBER(8,2)
             FIT_COEFF_FOURTH_A                        NUMBER(8,2)
             MEAN_CCD_A                             NUMBER(8,2)
             FLAT_STATUS_A                             VARCHAR2(80)
             FLAT_STATUS_B                             VARCHAR2(80)
            """
            
            #create all the strings for each dictionary. Most of them depends on the amp
            tot_rsq = 'TOT_RSQ_' + amp
            sumsqrred = 'SSR_' + amp
            red_chi2 = 'RED_CHI2_' + amp
            fitcoefzero = 'FIT_COEFF_ZERO_' + amp
            fitcoefone = 'FIT_COEFF_ONE_' + amp
            fitcoeftwo = 'FIT_COEFF_TWO_' + amp
            fitcoefthree = 'FIT_COEFF_THREE_' + amp
            fitcoeffourth = 'FIT_COEFF_FOURTH_' + amp
            meanccd = 'MEAN_CCD_' + amp
            flatstat = 'FLAT_STATUS_' + amp
            
            
            insertDictionaryFitsDb[tot_rsq] = Rsq
            insertDictionaryFitsDb[sumsqrred] = SSR
            insertDictionaryFitsDb[red_chi2] = chi2red 
            insertDictionaryFitsDb[fitcoefzero] = coefs[0]
            insertDictionaryFitsDb[fitcoefone] = coefs[1]
            insertDictionaryFitsDb[fitcoeftwo] = coefs[2]
            insertDictionaryFitsDb[fitcoefthree] = coefs[3]
            insertDictionaryFitsDb[fitcoeffourth] = coefs[4]
            insertDictionaryFitsDb[meanccd] = mean_ccd_amp
            insertDictionaryFitsDb[flatstat] = flatStatus[filenames[idx]][1]
            insertDictionaryFitsDb['FILENAME'] = filenames[idx]
            insertDictionaryFitsDb['CCDNUM'] = theccd
            insertDictionaryFitsDb['BAND'] = band
            insertDictionaryFitsDb['NITE'] = singleNite
            #insertDictionaryFitsDb['EXPNUM'] = expnums[idx]
            
            
            #insert into the database each dictionary value for each exposure/band/ccd value
            
            if amp == 'A':
                sqlInsertIntoDB = """ UPDATE ricardoc.flats_stats_qa_y2 SET 
                                    TOT_RSQ_A=:TOT_RSQ_A, SSR_A=:SSR_A, RED_CHI2_A=:RED_CHI2_A, FIT_COEFF_ZERO_A=:FIT_COEFF_ZERO_A,
                                    FIT_COEFF_ONE_A=:FIT_COEFF_ONE_A, FIT_COEFF_TWO_A=:FIT_COEFF_TWO_A, FIT_COEFF_THREE_A=:FIT_COEFF_THREE_A,
                                    FIT_COEFF_FOURTH_A=:FIT_COEFF_FOURTH_A, MEAN_CCD_A=:MEAN_CCD_A, FLAT_STATUS_A=:FLAT_STATUS_A
                                    WHERE FILENAME =:FILENAME  and CCDNUM=:CCDNUM and BAND=:BAND and NITE=:NITE""" 
            else:
                sqlInsertIntoDB = """ UPDATE ricardoc.flats_stats_qa_y2 SET 
                                    TOT_RSQ_B=:TOT_RSQ_B, SSR_B=:SSR_B, RED_CHI2_B=:RED_CHI2_B, FIT_COEFF_ZERO_B=:FIT_COEFF_ZERO_B,
                                    FIT_COEFF_ONE_B=:FIT_COEFF_ONE_B, FIT_COEFF_TWO_B=:FIT_COEFF_TWO_B, FIT_COEFF_THREE_B=:FIT_COEFF_THREE_B,
                                    FIT_COEFF_FOURTH_B=:FIT_COEFF_FOURTH_B, MEAN_CCD_B=:MEAN_CCD_B, FLAT_STATUS_B=:FLAT_STATUS_B
                                    WHERE FILENAME =:FILENAME  and CCDNUM=:CCDNUM and BAND=:BAND and NITE=:NITE""" 

            
            
            insert_dictionary_2Db(dbh, sqlInsertIntoDB, insertDictionaryFitsDb, verbose=verbose)
            
            #add all the other parameters for the fit and the mean to a dictionary to write them in the database
             
            #Making the plots
                
            axarr = plt.subplot(2, 2, position)
            #plt.clf()
            #plt.rc("font", size=9)

            
            if i%4 == 0:
                name_to_save_init = singleNite + '_' + str(filenames[idx])
            if i%4 == 3:
                name_to_save = name_to_save_init + '_' + str(filenames[idx]) + '_' + ccdnums[idx] + '_' + band + '_' + amp + '.png'
            elif i%4 > 0 and i%4 < 3:
                name_to_save = name_to_save_init + '_' + str(filenames[idx]) + '_' + ccdnums[idx] + '_' + band + '_' + amp + '.png'
                
            
            #print medians[idx], x
            #listMedians = medians[idx].tolist()
            #xvals = x.tolist()
            #print xvals, listMedians
            #x_new = x_new.tolist()
            #yfit_polyval_new = yfit_polyval_new.tolist()
            #yerr = yerr[idx].tolist()
            title = str(singleNite) + " ccd=" + ccdnums[i] + band
            
            axarr.plot(x, medians[idx], 'ro')
            axarr.plot(x_new, yfit_polyval_new, '-')
            axarr.errorbar(x, medians[idx], yerr[idx], fmt = 'ro', ecolor='r', capsize=0)
            axarr.set_title('$%s  [SSR=%4.2f, Rsq=%.4f,\, \chi^2_{red}=%.4f]$' % (filenames[idx], SSR, Rsq, chi2red), fontsize=9, color=[0,0,0])
            #print i, i%4, position, len(thisNiteIndex)
            if i%4 == 3 or i == len(thisNiteIndex)-1:
                
                if tosave:
                    if verbose:
                        print "saving image ", name_to_save
                    plt.suptitle(title)
                    plt.savefig(name_to_save)
                    plt.close()
                #else:
                #    plt.show()

                position = 0

            #elif i != 0 and i%4 == 0 or i == len(thisNiteIndex) - 1:
            #    plt.show()
            #    position = 0
            else:
                position = position + 1


    dbh.commit()
    
    allfits = []
    allratios = []
    xvals = []
    badIndex = [] #index for expnums/ccd which have a 0 at the edge of the ccd
    #Plot the fit to the data for each band and ccd
    for i in range(len(filenames)):
        edgeratioDic = {}
        xdata = len(medians[i]) + 1
        x = np.array(range(1,xdata))
        xvals.append(x)
        fit = poly.polyfit(x,medians[i],4)
        #print "fit ", fit
        allfits.append(fit)
        if medians[i][-1] != 0:
            edgeRatio = medians[i][0]/medians[i][-1]
            allratios.append(edgeRatio)
        else:
            badIndex.append(i)
            continue

        
        #Select good and bad flats
        if (edgeRatio > 0.8 and edgeRatio < 1.2):
            #print flatStatus[filenames[i]]
            if flatStatus[filenames[i]] == (ccdnums[i], 'VETO CHI2'):
                if verbose:
                    print "RATIO:",filenames[i], flatStatus[filenames[i]]
                flatStatus[filenames[i]] = (ccdnums[i], 'VETO CHI2 + GOOD RATIO')
                if verbose:
                    print "RATIO:",'VETO CHI2 + GOOD RATIO'
            if flatStatus[filenames[i]] == (ccdnums[i], 'VETO SATURATION'):
                if verbose:
                    print "RATIO:",filenames[i], flatStatus[filenames[i]]
                flatStatus[filenames[i]] = (ccdnums[i], 'VETO SATURATION + GOOD RATIO')
                if verbose:
                    print "RATIO:",'VETO SATURATION + GOOD RATIO'
            if flatStatus[filenames[i]] == (ccdnums[i], 'VETO SHUTTER'):
                if verbose:
                    print "RATIO:",filenames[i], flatStatus[filenames[i]]
                flatStatus[filenames[i]] = (ccdnums[i], 'VETO SHUTTER + GOOD RATIO')
                if verbose:
                    print "RATIO:",'VETO SHUTTER + GOOD RATIO'
            if flatStatus[filenames[i]] == (ccdnums[i], 'GOOD'):
                if verbose:
                    print "RATIO:",filenames[i], flatStatus[filenames[i]]
                flatStatus[filenames[i]] = (ccdnums[i], 'GOOD + GOOD RATIO')
                if verbose:
                    print "RATIO:",'GOOD + GOOD RATIO'
            #else:
            #    flatStatus[expnums[i]] = (ccdnums[i], 'GOOD RATIO')
            #if flatStatus[expnums[i]] == (ccdnums[i], 'VETO'):
            #    flatStatus[expnums[i]] = (ccdnums[i], 'VETO')
            #else:
            #    flatStatus[expnums[i]] = (ccdnums[i], 'GOOD')
        else: #VETO this exposure
            if flatStatus[filenames[i]] == (ccdnums[i], 'VETO CHI2'):
                if verbose:
                    print "RATIO:", filenames[i], flatStatus[filenames[i]]
                flatStatus[filenames[i]] = (ccdnums[i], 'VETO CHI2 + RATIO')
                if verbose:
                    print "RATIO:",'VETO CHI2 + RATIO'
            elif flatStatus[filenames[i]] == (ccdnums[i], 'VETO SATURATION'):
                if verbose:
                    print "RATIO:",filenames[i], flatStatus[filenames[i]]
                flatStatus[filenames[i]] = (ccdnums[i], 'VETO SATURATION + RATIO')
                if verbose:
                    print "RATIO:",'VETO SATURATION + RATIO'
            elif flatStatus[filenames[i]] == (ccdnums[i], 'VETO SHUTTER'):
                if verbose:
                    print "RATIO:",filenames[i], flatStatus[filenames[i]]
                flatStatus[filenames[i]] = (ccdnums[i], 'VETO SHUTTER + RATIO')
                if verbose:
                    print "RATIO:",'VETO SHUTTER + RATIO'
            else:
                if verbose:
                    print "RATIO:",filenames[i], flatStatus[filenames[i]]
                flatStatus[filenames[i]] = (ccdnums[i], 'VETO RATIO')
                if verbose:
                    print "RATIO:",'VETO RATIO'
            #flatStatus[expnums[i]] = (ccdnums[i], 'VETO')    

        borderratio = 'BORDERS_RATIO_' + amp
        flatstat = 'FLAT_STATUS_' + amp
        
        print "theccd", theccd
        edgeratioDic[borderratio] = edgeRatio
        edgeratioDic[flatstat] = flatStatus[filenames[i]][1]
        edgeratioDic['FILENAME'] = filenames[i]
        edgeratioDic['CCDNUM'] = theccd
        edgeratioDic['BAND'] = band
        edgeratioDic['NITE'] = singleNite
        

        #Insert edgratio into the database
        if amp == 'A':
            sqlInsertIntoDB = """ UPDATE ricardoc.flats_stats_qa_y2 SET
                                BORDERS_RATIO_A=:BORDERS_RATIO_A, FLAT_STATUS_A=:FLAT_STATUS_A
                                WHERE FILENAME=:FILENAME and CCDNUM=:CCDNUM and BAND=:BAND and NITE=:NITE"""
        else:
            sqlInsertIntoDB = """ UPDATE ricardoc.flats_stats_qa_y2 SET
                                BORDERS_RATIO_B=:BORDERS_RATIO_B, FLAT_STATUS_B=:FLAT_STATUS_B
                                WHERE FILENAME=:FILENAME and CCDNUM=:CCDNUM and BAND=:BAND and NITE=:NITE"""
                                
        insert_dictionary_2Db(dbh,sqlInsertIntoDB, edgeratioDic, verbose=verbose)


    
        #ingest into the database the values for each exposure/ccd/band the border_ratios values
    
    print "commiting borders"    
    dbh.commit()

    
    #Calculate the ratio for each bin along the Y axis of ccd

    stepsratioDic = {}
    allSteps = []
    myfilenames = []
    allfilenames = []
    for i, exp in enumerate(filenames):
        print i,exp
        listMedians = medians[i].tolist()
        steps = []
        #print listMedians, len(listMedians)
        for j, median in enumerate(listMedians):
            #print median
            if (j+1) < len(listMedians):
                if listMedians[j] == 0:
                    step1step2Ratio = 999
                    #print "XXXXXX Found value of Median of XXXXXXX = ", listMedians[j]
                else:
                    step1step2Ratio = listMedians[j+1]/listMedians[j]
                steps.append(step1step2Ratio)

        #Vetoing good and bad flats
        if band == 'Y': 
            #Steps are below or above the cut, so they are bad
            badSteps = [val for val in steps if (val <= 0.98 or val >= 1.02)]
        else:
            #steps are below or above the cut, so the are bad (cuts are slight different for Y than g,r,i,z
            badSteps = [val for val in steps if (val <= 0.99 or val >= 1.01)]
        print 'badsteps', badSteps
        
        #I might need to check if it has been veto by previous reason just to keep
        if len(badSteps) != 0: #If there is data with bad steps then VETO it
            
            #if flatStatus[exp] == (ccdnums[i], 'VETO CHI2'):
            #    flatStatus[expnums[i]] = (ccdnums[i], 'VETO CHI2 + STEPS')
            #elif flatStatus[exp] == (ccdnums[i], 'VETO SATURATION'):
            #    flatStatus[expnums[i]] = (ccdnums[i], 'VETO SATURATION + STEPS')
            #elif flatStatus[exp] == (ccdnums[i], 'VETO SHUTTER'):
            #    flatStatus[expnums[i]] = (ccdnums[i], 'VETO SHUTTER + STEPS')
            
            if flatStatus[exp] == (ccdnums[i], 'VETO RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO RATIO + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO SATURATION'):
                flatStatus[exp] = (ccdnums[i], 'VETO SATURATION + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO CHI2'):
                flatStatus[exp] = (ccdnums[i], 'VETO CHI2 + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO SHUTTER'):
                flatStatus[exp] = (ccdnums[i], 'VETO SHUTTER + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO CHI2 + RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO CHI2 + RATIO + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO SATURATION + RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO SATURATION + RATIO + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO SHUTTER + RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO SHUTTER + RATIO + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO CHI2 + GOOD RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO CHI2 + GOOD RATIO + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO SATURATION + GOOD RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO SATURATION + GOOD RATIO + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO SHUTTER + GOOD RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO SHUTTER + GOOD RATIO + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'GOOD'):
                flatStatus[exp] = (ccdnums[i], 'VETO GOOD + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'GOOD RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO GOOD RATIO + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'GOOD + GOOD RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO GOOD + GOOD RATIO + STEPS')
            else:
                flatStatus[exp] = (ccdnums[i], 'VETO STEPS')
            if verbose:
                print exp, flatStatus[exp], ccdnums[i]
        else:
            if flatStatus[exp] == (ccdnums[i], 'GOOD'):
                flatStatus[exp] = (ccdnums[i], 'GOOD + GOOD STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'GOOD RATIO'):
                flatStatus[exp] = (ccdnums[i], 'GOOD RATIO + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'GOOD + GOOD RATIO'):
                flatStatus[exp] = (ccdnums[i], 'GOOD + GOOD RATIO + STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO CHI2'):
                flatStatus[exp] = (ccdnums[i], 'VETO  CHI2 + GOOD STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO CHI2 + GOOD RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO  CHI2 + GOOD RATIO + GOOD STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO SATURATION'):
                flatStatus[exp] = (ccdnums[i], 'VETO SATURATION + GOOD STEPS')            
            elif flatStatus[exp] == (ccdnums[i], 'VETO SATURATION + GOOD RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO SATURATION + GOOD RATIO + GOOD STEPS')            
            elif flatStatus[exp] == (ccdnums[i], 'VETO SATURATION + CHI2'):
                flatStatus[exp] = (ccdnums[i], 'VETO SATURATION + CHI2 + GOOD STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO SHUTTER'):
                flatStatus[exp] = (ccdnums[i], 'VETO SHUTTER + GOOD STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO SHUTTER + GOOD RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO SHUTTER + GOOD RATIO + GOOD STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO RATIO + GOOD STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO CHI2 + RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO CHI2 + RATIO + GOOD STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO SATURATION + RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO SATURATION + RATIO + GOOD STEPS')
            elif flatStatus[exp] == (ccdnums[i], 'VETO SHUTTER + RATIO'):
                flatStatus[exp] = (ccdnums[i], 'VETO SHUTTER + RATIO + GOOD STEPS')
            else:
                flatStatus[exp] = (ccdnums[i], 'GOOD + GOOD RATIO + GOOD STEPS')
            if verbose:
                print exp, flatStatus[exp]
        
        print "flatsStatus[exp]", flatStatus[exp],exp    
        #print len(steps)
        allfilenames.append(exp)        
        allSteps.append(steps)
        
        #Insert into the database the list with steps for each expnum/amp
        flatstat = 'FLAT_STATUS_' + amp
        #stepratios = 'STEPS_RATIO_' + amp
        
        stepsratioDic['FILENAME'] = filenames[i]
        stepsratioDic['CCDNUM'] = theccd
        stepsratioDic['BAND'] = band        
        stepsratioDic[flatstat] = flatStatus[exp][1]
 
        
        #Insert edgratio into the database
        if amp == 'A':
            sqlInsertIntoDB = """ UPDATE ricardoc.flats_stats_qa_y2 SET
                                FLAT_STATUS_A=:FLAT_STATUS_A
                                WHERE FILENAME=:FILENAME and CCDNUM=:CCDNUM and BAND=:BAND"""
        else:
            sqlInsertIntoDB = """ UPDATE ricardoc.flats_stats_qa_y2 SET
                                FLAT_STATUS_B=:FLAT_STATUS_B
                                WHERE FILENAME=:FILENAME and CCDNUM=:CCDNUM and BAND=:BAND"""
        

        insert_dictionary_2Db(dbh,sqlInsertIntoDB, stepsratioDic, verbose=verbose)

        dbh.commit()
        
        #isVeto = flatStatus[exp][1]
        #if 'VETO' in isVeto:
        #     stepsratioDic['FLAT_STATUS'] = 'VETO'
        #else:
        #    stepsratioDic['FLAT_STATUS'] = 'GOOD'
    
        #print 'ccc'

        #dbh.commit()

    #dbh.close()
    #print flatStatus.items()

    
    uniqueNights = set(nites)
    #print "unique nights", uniqueNights, len(uniqueNights)
    #print uniqueNights
    #print flagDateRange
    
    if tosave:
        if flagDateRange and plotall is False:
            #create one file name for each date:
            foundFirst = False
            for val in uniqueNights:
                nite = val
                thisNiteIndexs = [i for i, x in enumerate(nites) if x == val] 
                fitDataFigurename = nite + '_' + str(ccdnums[0]) + '_' + band +  '_' + amp + '_fit.png'
                MyMultiPlot(filenames[thisNiteIndexs[0]:thisNiteIndexs[-1]], ccdnums[thisNiteIndexs[0]], allfits[thisNiteIndexs[0]:thisNiteIndexs[-1]], nite=nite, tosave=tosave, figureName = fitDataFigurename)    
                #print fitDataFigurename

                ratiosFigureName = nite +  '_' +str(ccdnums[0]) + '_' + band +  '_' + amp +'_ratios.png'
                stepsFigureName = nite + '_' +str(ccdnums[0]) + '_' + band +  '_' + amp +'_steps.png'
                exptimesRatiosFigurename = nite +  '_' +str(ccdnums[0]) + '_' + band +  '_' + amp +'_exptimes.png'
                expnumsRatiosFigurename = nite +  '_' +str(ccdnums[0]) + '_' + band +  '_' + amp +'_expnumsvsRatio.png'

                if band == 'r':
                    xLimRange = [5,25]
                else:
                    xLimRange = [min(exptimes) - 5, max(exptimes) + 5]
    
                MySinglePlotRatios(allfits[thisNiteIndexs[0]:thisNiteIndexs[-1]], allratios[thisNiteIndexs[0]:thisNiteIndexs[-1]], nite, tosave=tosave, figureName = ratiosFigureName)
                MySinglePlotExptimes(exptimes[thisNiteIndexs[0]:thisNiteIndexs[-1]], allratios[thisNiteIndexs[0]:thisNiteIndexs[-1]], nite, tosave=tosave, xlimits = xLimRange, figureName = exptimesRatiosFigurename)
                        
                MySinglePlotSteps(allfilenames[thisNiteIndexs[0]:thisNiteIndexs[-1]], allSteps[thisNiteIndexs[0]:thisNiteIndexs[-1]], nite, tosave=tosave, figureName = stepsFigureName)
                plotExpnumRatios(allfilenames[thisNiteIndexs[0]:thisNiteIndexs[-1]], allratios[thisNiteIndexs[0]:thisNiteIndexs[-1]], nite, tosave=tosave, figureName = ratiosFigureName)

        elif flagDateRange is True and plotall is True:
            if verbose:
                print "Plotting all data from the range of dates in one figure"
            nite = firstNight
            fitDataFigurename = firstNight + '-' + lastNight + '_' + str(ccdnums[0]) + '_' + band +  '_' + amp + '_fit.png'
            MyMultiPlot(filenames, ccdnums, allfits, nite = None, tosave=tosave, figureName = fitDataFigurename)    
            #print fitDataFigurename

            ratiosFigureName = firstNight + '-' + lastNight  +  '_' +str(ccdnums[0]) + '_' + band +  '_' + amp + '_ratios.png'
            stepsFigureName = firstNight + '-' + lastNight + '_' +str(ccdnums[0]) + '_' + band +  '_' + amp + '_steps.png'
            exptimesRatiosFigurename = firstNight + '-' + lastNight + '_' +str(ccdnums[0]) + '_' + band +  '_' + amp + '_exptimes.png'
            expnumsRatiosFigurename = firstNight + '-' + lastNight + '_' +str(ccdnums[0]) + '_' + band +  '_' + amp + '_expnumsvsRatio.png'

            if band == 'r':
                xLimRange = [5,25]
            else:
                xLimRange = [min(exptimes) - 5, max(exptimes) + 5]

            thedateRange = str(firstNight) + '-' + str(lastNight)
            MySinglePlotRatios(allfits, allratios, nite, title= thedateRange, tosave=tosave, figureName = ratiosFigureName)
            MySinglePlotExptimes(exptimes, allratios, nite, title= thedateRange, tosave=tosave, xlimits = xLimRange, figureName = exptimesRatiosFigurename)
            #MySinglePlotRatios(expnums, allratios, nite, title= thedateRange, tosave=tosave, figureName = ratiosFigureName)
                    
            MySinglePlotSteps(allfilenames, allSteps, nite, tosave=tosave, figureName = stepsFigureName)
            
            plotExpnumRatios(allfilenames, allratios, nite, title= thedateRange, tosave=tosave, figureName = expnumsRatiosFigurename)
                
        else:
            fitDataFigurename = nites[0] + '_' + str(ccdnums[0]) + '_' + band +  '_' + amp + '_fit.png'
            #print fitDataFigurename
            #print allfits
            MyMultiPlot(filenames, ccdnums[0], allfits, nite=nites[0], tosave=tosave, figureName = fitDataFigurename)
            
            ratiosFigureName = nites[0] +  '_' +str(ccdnums[0]) +  '_' + band +  '_' + amp + '_ratios.png'
            stepsFigureName = nites[0] + '_' +str(ccdnums[0]) + '_' + band +  '_' + amp + '_steps.png'
            exptimesRatiosFigurename = nites[0] +  '_' +str(ccdnums[0]) + '_' + band +  '_' + amp + '_exptimes.png'
            expnumsRatiosFigurename = nites[0] +  '_' +str(ccdnums[0]) + '_' + band +  '_' + amp + '_expnumsvsRatio.png'
            if band == 'r':
                xLimRange = [5,25]
            else:
                xLimRange = [min(exptimes) - 5, max(exptimes) + 5]
            
            MySinglePlotRatios(allfits, allratios, nites[0], tosave=tosave, figureName = ratiosFigureName)
            MySinglePlotExptimes(exptimes, allratios, nites[0], tosave=tosave, xlimits = xLimRange, figureName = exptimesRatiosFigurename)
            
            MySinglePlotSteps(allfilenames, allSteps, nites[0], tosave=tosave, figureName = stepsFigureName)
            
            plotExpnumRatios(allfilenames, allratios, nites[0], tosave=tosave, figureName = expnumsRatiosFigurename)




    print "End of Program:"

    print "TAG FLATS "
    toOrderDict = []
    allvetoexps = []
    summaryFilename = firstNight + '-' + lastNight + '_' +str(ccdnums[0]) + '_' + band +  '_' + amp + '_summary.txt'
    with open(summaryFilename, 'w') as f:

        for key, val in sorted(flatStatus.items()):
            toOrderDict.append(key)
            #create a list with all veto flats
            if 'VETO' in val[1]:
                allvetoexps.append(key)
                
            #print "expnum = %s and status = %s" % (key, val)
            
            f.write("expnum = %s and status = %s \n" % (key, val))
    
    #select the median values for each exposure
    
    print "Number of exposures %d and number of flags exposures %d " % (len(filenames), len(flatStatus.keys()))
    
    #Select all exposures with same exptime
    # Band    exptime    counts (ADU)
    # g        40        16,000 - 18,000
    # r        15        20,000
    # i        30        18,000 - 20,000
    # z        12        21,000 - 24,000
    # Y        12        20,000 - 23,000
 
    if band == 'g':
        expindx = [indx for indx, etime in enumerate(exptimes) if etime == 40]
        selallmedians = [medians[i] for i in expindx]
        selnites = [nites[i] for i in expindx]
        selexpnums = [filenames[i] for i in expindx]
    elif band == 'r':
        expindx = [indx for indx, etime in enumerate(exptimes) if etime == 15]
        selallmedians = [medians[i] for i in expindx]
        selnites = [nites[i] for i in expindx]
        selexpnums = [filenames[i] for i in expindx]
    elif band == 'i':
        expindx = [indx for indx, etime in enumerate(exptimes) if etime == 30]
        selallmedians = [medians[i] for i in expindx]
        selnites = [nites[i] for i in expindx]
        selexpnums = [filenames[i] for i in expindx]
    elif band == 'z' or band == 'Y':
        expindx = [indx for indx, etime in enumerate(exptimes) if etime == 12]
        selallmedians = [medians[i] for i in expindx]
        selnites = [nites[i] for i in expindx]
        selexpnums = [filenames[i] for i in expindx]
        
    
    intselexpnums = filenameToExpnum(selexpnums)
    
    allmedians = []
    allselexpnums = []
    allselnites = []
    #Calculate The median for all nights inputs
    for i, med in enumerate(selallmedians):
        for val in med:
            allselexpnums.append(intselexpnums[i])
            allselnites.append(selnites[i])
            allmedians.append(val)
    
    allmedians = np.asarray(allmedians,dtype=np.float32)
    allNightsMedians, allNightsMean, allNightsSigma = medsigclip(allmedians, clipsig=2.0, maxiter=10, converge_num=0.001, verbose=verbose)
        
    #print allNightsMedians, allNightsMean, allNightsSigma
    
    #Select the min expnum for plotting purposes
    minexpnum = min(allselexpnums)
    newexp = []
    for exp in allselexpnums:
        newexp.append(exp-minexpnum)
    
    #Transform dates to proper pyplot format
    xdates = [dt.datetime.strptime(str(int(date)),'%Y%m%d') for date in allselnites]
    
    if flagDateRange:
        firstNight, lastNight
        xrange = [dt.datetime.strptime(firstNight,'%Y%m%d'),dt.datetime.strptime(lastNight,'%Y%m%d')]
    else:
        xrange = [dt.datetime.strptime(str(int(20140821)),'%Y%m%d'),dt.datetime.strptime(str(int(20141023)),'%Y%m%d')]
        
    #print xdates
    #print dates

    #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d'))
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y%m/%d'))

    #find the difference between the 
    indxvetoexpnums = []
    mediansveto = []
    nitesveto = []
    for indx, exp in enumerate(allfilenames):
        if exp in allvetoexps:
            for med in medians[indx]:
                #indxvetoexpnums.append(indx)
                mediansveto.append(med)
                nitesveto.append(nites[indx])

    #print mediansveto
    #print nitesveto
    vetodates = [dt.datetime.strptime(str(int(date)),'%Y%m%d') for date in nitesveto]
    
    
    plt.figure(1)
    plt.clf()
    plt.grid()
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.15)

    plt.rc("font", size=9)
    #print len(selallmedians), len(allselexpnums)

    residuals = []
    for i,val in enumerate(allmedians):
        #print allselexpnums[i],val, allNightsMedians, val - allNightsMedians
        residuals.append((val - allNightsMedians)/allNightsMedians)
        #print len(residuals), (val-allNightsMedians), allselexpnums[i]
    medresids = allmedians - allNightsMedians
    
    residualsVeto = []
    for i, vetoval in enumerate(mediansveto):
        residualsVeto.append((vetoval - allNightsMedians)/allNightsMedians)
    
    
    #plt.plot(allmedians, medresids, 'ro')
    #plt.plot(xdates, residuals, 'bo')
    plt.plot_date(xdates, residuals, 'bo', xdate=True, ydate=False)
    plt.plot_date(vetodates, residualsVeto, 'ro', xdate=True, ydate=False)
    #loc = mdates.WeekdayLocator(byweekday=MO, tz=tz)
    #mdates.AutoDateLocator(minticks=10)
    #plt.gca().xaxis.set_major_locator(loc)
    plt.ylim([-0.2,0.2])
    plt.xlabel('Date')
    plt.ylabel('(Median Counts bin - Median_allNights_Y1) / Median_allNights_Y1', fontsize=9)
    figname = str(firstNight) + '-' + str(lastNight) + '_' + band + '_' + str(ccdnums[0]) +  '_' + amp + '_dateResidualFromMean.png'
    plt.savefig(figname)
    
    
def tagVetoFlats(dbh, band, verbose=verbose):
    

    
    update_query_1 = """update ricardoc.flats_stats_qa_y2 set 
                    flat_status='VETO' where filename in (
                    select filename from flats_stats_qa_y2 where band='%s' 
                    and flat_status_a like 'VETO%%' 
                    and flat_status_b like 'GOOD%%' 
                    )""" % (band)
                
    update_query_2 = """update ricardoc.flats_stats_qa_y2 set 
                    flat_status='VETO' where filename in (
                    select filename from flats_stats_qa_y2 where band='%s' 
                    and flat_status_a like 'GOOD%%' 
                    and flat_status_b like 'VETO%%' 
                    )""" % (band)

    update_query_3 = """update ricardoc.flats_stats_qa_y2 set 
                    flat_status='GOOD' where filename in (
                    select filename from flats_stats_qa_y2 where band='%s' 
                    and flat_status_a like 'GOOD%%' 
                    and flat_status_b like 'GOOD%%' 
                    )""" % (band)


    update_flatstatus_2Db(dbh, update_query_1, verbose=verbose)    
    update_flatstatus_2Db(dbh, update_query_2, verbose=verbose)    
    update_flatstatus_2Db(dbh, update_query_3, verbose=verbose)    

        #print flat_status
        
    #update_flatstatus_2Db(dbh, sqlInsertIntoDB, flat_status, verbose=verbose)
    
    dbh.commit()
        
        
    return 1
    
def main():
    
    
    
    #dateRange = args.dates
    #band = args.band
    #ccdRange = args.ccdrange
    #section = args.section
    #tosave = args.output
    #plotall = args.plotall
    #verbose = args.verbose
    if verbose:
        print "dateRange ", dateRange
        print "Band ", band
        print "DB section to use ", section
        print "ccd Amp to process ", theamp
        print "Single or range of ccds to process ", ccdRange
    if tosave:
        print "Saving all plots"
    if plotall:
        print "making a plot for all daterange"
        
    #If a ccdrange given it must contain firstccd-lastccd.
    #If a single ccd is given, then only process that
    #The code is not using ccd 2, 61 no matter what ccdrange and daterange is given.
    #If ccd = 2 needs to be processed for Y1 and SV data, do it by single ccd
    if len(ccdRange.split('-')) > 1:
        firstCcd, lastCcd = ccdRange.split('-')
        ccdsToProcess = list(range(int(firstCcd),int(lastCcd)+1))   #Creates a list with all ccds to process given by the user.
        if verbose:
            print "ccds to process", ccdsToProcess
        flagCcdRange = True
        if 61 in ccdsToProcess:
            ccdsToProcess.remove(61)
        if 2 in ccdsToProcess:
            ccdsToProcess.remove(2)
    #A single ccd is being processed
    else:
        singleccd = ccdRange
        flagCcdRange = False
        if verbose:
            print "Only One ccd provided", singleccd
        
    
    #If a list of comma separated bands is given, split and create a list
    if len(band.split(',')) > 1:
        bands = band.split(',')
        flagBands = True
        if verbose:
            print "bands", bands
    else:
        flagBands = False
        if verbose:
            print "only One band provided", band

    #If amp A and B are given, they must be comma separated. It creates a list with both amps.
    if len(theamp.split(',')) > 1:
        amps = theamp.split(',')
        flagAmp = True
        if verbose:
            print "Amps ", amps
    else:
        singleAmp = theamp
        flagAmp = False
        if verbose:
            print "only One Amp provided ", amp
    
    #make connection to database
    dbh = connectDB(section)


    #Check if only need to tag all flats (after the data has been ingested
    if flatsTag:
        if flagBands:
            for filter in bands:
                updated = tagVetoFlats(dbh, filter, verbose=verbose)
        else:
            updated = tagVetoFlats(dbh, theband, verbose=verbose)

        if updated:
            print "succesfully updated all flat tags"
            sys.exit(1)
        
    
    if flagBands: #A list of bands given comma separated
        for myband in bands:  #Loops trought the list of bands
            if flagCcdRange:  #Loops through the range of ccd given by user
                for ccd in ccdsToProcess:
                    if flagAmp:
                        for amp in amps:
                            if verbose:
                                print ccd,myband,amp
                            calculateallvetoflats(dbh, myband, dateRange, ccd, section, amp, tosave=tosave, plotall=plotall, verbose=verbose )
                            #update all expnum. If there is one ccd/amp in the exposure that is bad, then update to VETO the whole exposure
                            updated = tagVetoFlats(dbh, myband, verbose=verbose) 
                    else:
                        if verbose:
                            print ccd,myband,amp
                        calculateallvetoflats(dbh, myband, dateRange, ccd, section, singleAmp, tosave=tosave, plotall=plotall, verbose=verbose )
                        #update all expnum. If there is one ccd/amp in the exposure that is bad, then update to VETO the whole exposure
                        updated = tagVetoFlats(dbh, myband, verbose=verbose) 
            else:
                print myband
                if flagAmp:
                    for amp in amps:
                        if verbose:
                            print singleccd,myband,amp
                        calculateallvetoflats(dbh, myband, dateRange, singleccd, section, amp, tosave=tosave, plotall=plotall, verbose=verbose ) 
                        #update all expnum. If there is one ccd/amp in the exposure that is bad, then update to VETO the whole exposure
                        updated = tagVetoFlats(dbh, myband, verbose=verbose) 
                else:
                    if verbose:
                        print ccd,myband,amp
                    calculateallvetoflats(dbh, myband, dateRange, singleccd, section, singleAmp, tosave=tosave, plotall=plotall, verbose=verbose )
                    #update all expnum. If there is one ccd/amp in the exposure that is bad, then update to VETO the whole exposure
                    updated = tagVetoFlats(dbh, myband, verbose=verbose) 
    else: #One band given
        if flagCcdRange:
            for ccd in ccdsToProcess:
                if flagAmp:
                    for amp in amps:
                        print ccd,band,amp
                        calculateallvetoflats(dbh, band, dateRange, ccd, section, amp, tosave=tosave, plotall=plotall, verbose=verbose ) 
                        #update all expnum. If there is one ccd/amp in the exposure that is bad, then update to VETO the whole exposure
                        updated = tagVetoFlats(dbh, band, verbose=verbose) 
                else:
                    if verbose:
                        print ccd,band,amp
                    calculateallvetoflats(dbh, band, dateRange, ccd, section, singleAmp, tosave=tosave, plotall=plotall, verbose=verbose )
                    #update all expnum. If there is one ccd/amp in the exposure that is bad, then update to VETO the whole exposure
                    updated = tagVetoFlats(dbh, band, verbose=verbose) 
        else:
            print band
            if flagAmp:
                for amp in amps:
                    calculateallvetoflats(dbh, band, dateRange, singleccd, section, amp, tosave=tosave, plotall=plotall, verbose=verbose )
                    #update all expnum. If there is one ccd/amp in the exposure that is bad, then update to VETO the whole exposure
                    updated = tagVetoFlats(dbh, band, verbose=verbose) 
            else:
                calculateallvetoflats(dbh, band, dateRange, singleccd, section, singleAmp, tosave=tosave, plotall=plotall, verbose=verbose )
                #update all expnum. If there is one ccd/amp in the exposure that is bad, then update to VETO the whole exposure
                updated = tagVetoFlats(dbh, band, verbose=verbose) 



    dbh.close()
    
    #sys.exit()
    
    #plt.show()
    

    """
    #how to fit a polynomial using numpy polyfit
    # calculate polynomial
    z = np.polyfit(x, y, 3)
    f = np.poly1d(z)
    
    # calculate new x's and y's
    x_new = np.linspace(x[0], x[-1], 50)
    y_new = f(x_new)
    
    plt.plot(x,y,'o', x_new, y_new)
    plt.xlim([x[0]-1, x[-1] + 1 ])
    plt.show()
    

    
    domePositionCounts = query_to_cur(dbh, domePosition, verbose=verbose)
    counts = []
    domeAngle = []
    dates = []
    for val in domePositionCounts:
        singleCounts, singleDomeAngle, singleDate = val
        
        counts.append(singleCounts)
        domeAngle.append(singleDomeAngle)
        dates.append(singleDate)
        
    #plot all data. Use the median of the nite 20130806
    #I am using an arbitrary night as the zero point. This night is: 20130806
    dateIndex = []
    for i, date in enumerate(dates):
        if date == '20130806':
            dateIndex.append(i)
    
    referenceCounts = []
    referenceDomeAngle = []
    #calculate median for date
    for index in dateIndex:
        #print index
        #print dates[index], counts[index], domeAngle[index]
        referenceCounts.append(counts[index])
        referenceDomeAngle.append(domeAngle[index])
        
    #Calculate the Median and stddev for each nite
    medianReferenceCounts, MeanReferenceCounts, SigmaReferenceCounts = medsigclip(np.asarray(referenceCounts), clipsig=3.0, maxiter=10, converge_num=0.001)
    medianReferenceDomeAngle, MeanReferenceDomeAngle, SigmaReferenceDomeAngle = medsigclip(np.asarray(referenceDomeAngle), clipsig=3.0, maxiter=10, converge_num=0.001)
    
    print medianReferenceCounts, SigmaReferenceCounts, medianReferenceDomeAngle,SigmaReferenceDomeAngle
    
    #calculate the median for each night and plot relative to reference night
        #setting the figure
    fig = plt.figure()
    plt.rc("font", size=11)
    ax = fig.add_subplot(111)
    #ax = plt.subplot(111)
    #ax.set_title(rootExposureName)
    ax.set_xlabel('Dome angle')
    ax.set_ylabel('ratio of counts')
    
    for dateCounts, dateDomeAngle, dateDate in zip(counts, domeAngle, dates):
        if dateDomeAngle == -999:
            continue
        print dateCounts/medianReferenceCounts, dateDomeAngle, dateDate
        ratioCounts = dateCounts/medianReferenceCounts
        if dateDate == '20130806':
            ax.plot(dateDomeAngle, ratioCounts, 'bo')
        else:
            ax.plot(dateDomeAngle, ratioCounts, 'ro')
    
    ax.set_ylim(0.9,1.25)
    #ax.set_xlim(0,20)
    plt.show()
    """
    return 0

if  __name__ == '__main__':
    status = main()
    sys.exit(status)
    