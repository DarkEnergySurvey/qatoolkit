#! /usr/bin/env python3
# $Id: assess_SE.py 44620 2016-11-17 17:22:35Z rgruendl $
# $Rev: 44620 $:  # Revision of last commit.
# $LastChangedBy: rgruendl $:  # Author of last commit.

"""
Utilities internal to the single-epoch assessment script
"""
#from __future__ import print_function

import numpy
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')

########################
def QAplt_cloud( PlotFile, Data, astrom_good, verbose=False):
    """
    Function to produce a QA plot showing comparison between data and a reference
    catalog.  Plots show the diagnosed amount of atmospheric opacity (extinction)
    at the time of the observation.

    inputs:
        PlotFile:   filename where the output plot will be written (can include
                        relative path)
        Data:       Custom dictionary that holds data for the plots.  These include:
                        nmatch:  (integer) number of matches to the reference catalog

                        mag_des: list of rough measured DES magnitude
                        mag_cat: list of associated reference catalog data
                        mag_diff: the difference between the two
                        mag_des_reject: points in mag_des that were outlier rejects
                        mag_cat_reject: analogous points in mag_cat
                        mag_diff_reject: analogous points in mag_diff

                        mindes,maxdes: (float) range of data being plotted
                        med_magdiff:   (float) median magnitide offset/difference
                                       between objects and refernce catalog

                        band:          (string) Band of observations
                        magtype:       (string) Type of magnitude plotted
                        cat:           (string) Name of the reference catalog for plotting axis
                    For failure cases (number of stars < 100) the following are not needed:
                        mag_des_reject,mag_cat_reject, mag_diff_reject
                        mindes,maxdes,med_magdiff

        astrom_good: Flag that causes overplotting of text indicating an astrometric failure
        verbose:    Provide verbose output (curently there is none).

    OUTPUT:
        Nothing is returned to the calling program (a PNG file is written where directed)
    """

    plt.figure()
    plt.subplot(2,1,1)
    if (Data['nmatch']<100):
#       Case no fit
        if (Data['mag_des'].size > 0):
            plt.scatter(Data['mag_des'],Data['mag_cat'],marker='.',color='blue')
        plt.plot([10,18],[10,18],color='red',linewidth=1)
    else:
#       Case normal
        plt.scatter(Data['mag_des'],Data['mag_cat'],marker='.',color='blue')
        if (Data['mag_des_reject'].size > 0):
            plt.scatter(Data['mag_des_reject'],Data['mag_cat_reject'],marker='.',color='magenta')
        plt.plot([Data['mindes'],Data['maxdes']],[Data['mindes'],Data['maxdes']],color='red',linewidth=1)

    plt.xlabel('DES MAG_{:s}({:s})'.format(Data['magtype'],Data['band']))
    if (Data['band'] in ["u","g"]):
        plt.ylabel('{:s} g\''.format(Data['cat'].upper()) )
    elif (Data['band'] in ["r","VR"]):
        plt.ylabel('{:s} r\''.format(Data['cat'].upper()) )
    elif (Data['band'] in ["i","z","Y"]):
        plt.ylabel('{:s} i\''.format(Data['cat'].upper()) )
    else:
        plt.ylabel('{:s} [unknown]'.format(Data['cat'].upper()) )
    if ('mindes' in Data):
        xtxt=Data['mindes']+(0.025*(Data['maxdes']-Data['mindes']))
        ytxt=Data['maxdes']-(0.15*(Data['maxdes']-Data['mindes']))
    else:
        xtxt=10.0+(0.025*(18.-10.))
        ytxt=18.0-(0.15*(18.-10.))
    plt.text(xtxt,ytxt,'Number of {:s} matches: {:d}'.format(Data['cat'].upper(),Data['nmatch']))
    if (not(astrom_good)):
        if ('mindes' in Data):
            ytxt=Data['maxdes']-(0.075*(Data['maxdes']-Data['mindes']))
        else:
            ytxt=18.0-(0.075*(18.-10.))
        plt.text(xtxt,ytxt,'Probable failure to find an astrometric solution')
#
#   Second panel
#
    plt.subplot(2,1,2)
    if (Data['nmatch']<100):
#       Case no fit
        plt.plot([10,18],[0,0],color='red',linewidth=3)
    else:
#       Case normal
        plt.scatter(Data['mag_des'],Data['mag_diff'],marker='.',color='blue')
        if (Data['mag_des_reject'].size > 0):
            plt.scatter(Data['mag_des_reject'],Data['mag_diff_reject'],marker='.',color='magenta')
        plt.plot([Data['mindes'],Data['maxdes']],[Data['med_magdiff'],Data['med_magdiff']],color='red',linewidth=3)

    plt.xlabel('DES MAG_{:s}({:s})'.format(Data['magtype'],Data['band']))
    if (Data['band'] in ["u","g"]):
        plt.ylabel('DES - {:s} g\''.format(Data['cat'].upper()) )
    elif (Data['band'] in ["r","VR"]):
        plt.ylabel('DES - {:s} r\''.format(Data['cat'].upper()) )
    elif (Data['band'] in ["i","z","Y"]):
        plt.ylabel('DES - {:s} i\''.format(Data['cat'].upper()) )
    else:
        plt.ylabel('DES - {:s} [unknown]'.format(Data['cat'].upper()) )
    plt.savefig(PlotFile)

    return 0


########################
def QAplt_maghist(PlotFileName,Data,astrom_good,verbose=False):
    """
    Function to produce a QA plot showing luminosity function of objects in the
    exposure.

    inputs:
        PlotFile:   filename where the output plot will be written (can include
                        relative path)
        Data:       Custom dictionary that holds data for the plots.
                        nobj:        (integer) total number of objects
                                         contributing to the histogram
                        bins:        list of magnitude bins
                        mag_hist:    list of number of objects per bin
                        magerr_hist: list of median value of magerr for
                                        objects in each magnitude bin
                        band:        (string) Band of observations
                        magtype:       (string) Type of magnitude plotted
                    For failure cases (nobj < 100) no data from lists are plotted

        astrom_good: Flag that causes overplotting of text indicating an astrometric failure
        verbose:    Provide verbose output (curently there is none).

    OUTPUT:
        Nothing is returned to the calling program (a PNG file is written where directed)
    """


    plt.figure()
    plt.subplot(2,1,1)
#    plt.scatter(xplt,yplt1,marker='.',color='blue')
    plt.axis([8,26,0.5,10000])
    if (Data['nobj']>10):
        plt.semilogy(Data['bins'],Data['mag_hist'],marker='.',ls='None',color='blue')
    else:
        plt.text(9.,9000.,'Number of objects ({:d}) insufficient for histogram.'.format(Data['nobj']))
    if (not(astrom_good)):
        plt.text(9.,7800,'Probable failure to find an astrometric solution')
    plt.xlabel('DES {:s}({:s})'.format(Data['magtype'],Data['band']))
    plt.ylabel('# objects')
#
#   Second panel showing median MAGERR per bin
#
    plt.subplot(2,1,2)
    if (Data['nobj']>10):
        plt.scatter(Data['mag_raw'],Data['magerr_raw'],marker='.',color='gray')
        plt.scatter(Data['bins'],Data['magerr_hist'],marker='.',color='blue')
    plt.axis([8,26,-0.01,0.5])
    plt.xlabel('DES MAG_{:s}({:s})'.format(Data['magtype'],Data['band']))
    plt.ylabel('median(MAGERR_{:s}({:s})'.format(Data['magtype'],Data['band']))
    plt.savefig(PlotFileName)

    return 0
