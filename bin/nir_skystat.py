#! /usr/bin/env python3

# Importing necessary modules
#from __future__ import print_function
import argparse
import os
import despydb.desdbi
from pixcorrect.lightbulb_utils import medclip
from pixcorrect.corr_util import logger
import re
import time
import sys
import numpy as np
import fitsio
from scipy.optimize import curve_fit
from scipy.special import wofz

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

minorLocator = AutoMinorLocator()

###########################################
def get_data(filename,verbose=0):
    """Function to obtain header, image, mask data
    """
    
    ifits=fitsio.FITS(filename,'r') # Could be changed to 'rw' as needed
    ih=ifits['SCI'].read_header()
    isci=ifits['SCI'].read()
    imsk=ifits['MSK'].read()
#    iwgt=ifits['WGT'].read()
    ifits.close()
    if (verbose > 0):
        print("Successfully read {:s}".format(filename))

    return ih,isci,imsk


###########################################
def ingest_sky_entries(sDict,DBtable,DBorder,dbh,updateDB,verbose=0):
    """Ingest a set of sky entries"""

    t0=time.time()
    InsertCnt=0
    print("Preparing lists of list to ingest many zeropoints")
#
    new_data=[]
    for entry in sDict:
        AddRow=True
        new_row=[]
        for col in DBorder:
            if (col.lower() in sDict[entry]):
                 new_row.append(sDict[entry][col.lower()])
#                print(sDict[entry][col.lower()],type(sDict[entry][col.lower()]))
            elif (col.lower() == "filename"):
                new_row.append(entry)
            else:
                print("Warning: sDict[{:s}] entry does not have key {:s}.  Ingest will be skipped".format(entry,col))
                AddRow=False
        if (AddRow):
            new_data.append(new_row)

    t1=time.time()

    print("Successfully Formed list for Ingestion (Nrows={:d}, Timing: {:.2f})".format(len(new_data),(t1-t0)))
    print(DBorder)
    print(new_data[0])
#    print(new_data[1])
#    print(new_data[2])
#    for row in new_data:
#        print(row)

    if (updateDB):
        print("# Loading {:s} with {:d} entries".format(DBtable,len(new_data)))
        t1=time.time()
        dbh.insert_many(DBtable,DBorder,new_data)
        t2=time.time()
        print("Commit of {:d} rows, timing={:.2f} of {:.2f})".format(len(new_data),(t2-t1),(t2-t0)))
        dbh.commit()
    ningest=len(new_data)

    return ningest

###########################################
if __name__ == "__main__":


    t00=time.time()
    parser = argparse.ArgumentParser(description='Search for and mask pixels affected by light bulb')

    parser.add_argument('--img', action='store', type=str, required=True, help='Image filename (or list of images)')

    parser.add_argument('--trim',       action='store', type=int, default=200, help='Image border to trim before obtaining stats (default=200)')
    parser.add_argument('--patch_size', action='store', type=int, default=20,  help='Patch size for stats (default=20)')

    parser.add_argument('--pfw_attempt_id', action='store', type=int, default=None, help='PFW_ATTEMPT_ID (optional, used in DB ingest, default=NULL)')
    parser.add_argument('-D', '--DBTable', action='store', type=str, default='NIR_SKY_QA',  help='Tablename for QA ingest')
    parser.add_argument('--updateDB'    , action='store_true', default=False, help='Required flag for DB ingest (otherwise dry run)')
    parser.add_argument('-T','--Timing',  action='store_true', default=False, help='If set timing information accompanies output')
    parser.add_argument('--debug'       , action='store_true', default=False, help='Debug mode resticts code to work on a handful of objects')
    parser.add_argument('-v','--verbose', action='store', type=int, default=0, help='Verbosity (defualt:0; currently values up to 2)')
    parser.add_argument('-s','--section', action='store', type=str, default=None,   help='section of .desservices file with connection info')
    parser.add_argument('-S','--Schema',  action='store', type=str, default=None,   help='DB schema (do not include \'.\').')


    args = parser.parse_args()
    if (args.verbose > 0):
        print("Args: {:}".format(args))
    verbose = args.verbose

    if (args.Schema is None):
        dbSchema=""
    else:
        dbSchema="{:s}.".format(args.Schema)


#   Sigma clipping algorithm parameters
#
    maxiter=10
    converge_num=0.01
    clipsig=3.0

#
#   Blocking parameters
#
    trim = args.trim
    patch_size = args.patch_size

#
#   Setup DB connection
#
    try:
        desdmfile = os.environ["des_services"]
    except KeyError:
        desdmfile = None
    dbh = despydb.desdbi.DesDbi(desdmfile,args.section)
#
#   Obtain image/list to work through
#

    img_dict={}
    if (not(os.path.isfile(args.img))):
        print("Image/list {:s} not found!".format(args.img))
        print("Aborting!")
        exit(1)
    else:
        if ((args.img[-4:] == "fits")or
            (args.img[-7:] == "fits.fz")):
            tmp_fname=args.img.split('/')[-1]
            tmp_fname=re.sub('.fz','',tmp_fname)
            img_dict[tmp_fname]={'fpath':args.img}
        else:
            f_img=open(args.img,'r')
            for line in f_img:
                line=line.strip()
                columns=line.split(' ')
                if (columns[0] != "#"):
                    tmp_fname=columns[0].split('/')
                    img_fname=re.sub('.fz','',tmp_fname[-1])
                    img_dict[img_fname]={'fpath':columns[0]}
            f_img.close()

#
#   Main body:
#     - read each image/mask
#     - trim and then form statistic on individual patches
#     - accumulate into img_dict
#
    for img in img_dict:
        fname=img_dict[img]['fpath']
        if (not(os.path.isfile(fname))):
            print("Image {:s} not found!".format(fname))
            print("Aborting!")
            exit(1)
        else:
            # load files
            ih,isci,imsk=get_data(fname,verbose=args.verbose)

            # grab header args 
            img_dict[img]['pfw_attempt_id']=args.pfw_attempt_id
            img_dict[img]['expnum']=ih['EXPNUM']
            img_dict[img]['ccdnum']=ih['CCDNUM']
            img_dict[img]['band']=ih['BAND']
            img_dict[img]['nite']=ih['NITE']

            # trim the border
            isci = isci[trim:-1*trim,trim:-1*trim]
            imsk = imsk[trim:-1*trim,trim:-1*trim]

            # compute the patch means and RMSs
            means = []
            rms = []
            for i in range(int(isci.shape[0]/patch_size)-1):
                for j in range(int(isci.shape[1]/patch_size)-1):

                    # cut the patches
                    patch = isci[i*patch_size:(i+1)*patch_size, j*patch_size:(j+1)*patch_size]
                    pmask = imsk[i*patch_size:(i+1)*patch_size, j*patch_size:(j+1)*patch_size]

                    # test to see if there is >50% free pixels
                    # mask is non-zero where there are sources
                    if np.size(pmask[pmask==0]) >= 0.5*patch_size*patch_size:
                        means.append(np.mean(patch[pmask==0]))
                        rms.append(np.std(patch[pmask==0]))

            means = np.array(means)
            rms = np.array(rms)

            img_dict[img]['avg_means']=float(np.average(means))
            img_dict[img]['std_means']=float(np.std(means))
            img_dict[img]['avg_rms']=float(np.average(rms))
            img_dict[img]['std_rms']=float(np.std(rms))
            
            avgval,medval,stdval=medclip(means,clipsig=clipsig,maxiter=maxiter,converge_num=converge_num,verbose=args.verbose)

            p5s=avgval+(5.0*stdval)
            m5s=avgval-(5.0*stdval)
            wsm=np.where(np.logical_or(means>p5s,means<m5s))
            img_dict[img]['nsample']=int(means.size)
            img_dict[img]['nout_means']=int(means[wsm].size)
            img_dict[img]['cavg_means']=float(avgval)
            img_dict[img]['cmdn_means']=float(medval)
            img_dict[img]['cstd_means']=float(stdval)

            avgval,medval,stdval=medclip(rms,clipsig=clipsig,maxiter=maxiter,converge_num=converge_num,verbose=args.verbose)
            p5s=avgval+(5.0*stdval)
            m5s=avgval-(5.0*stdval)
            wsm=np.where(np.logical_or(rms>p5s,rms<m5s))
            img_dict[img]['nout_rms']=int(rms[wsm].size)
            img_dict[img]['cavg_rms']=float(avgval)
            img_dict[img]['cmdn_rms']=float(medval)
            img_dict[img]['cstd_rms']=float(stdval)

    DBorder=['FILENAME','EXPNUM','CCDNUM','BAND','NITE','PFW_ATTEMPT_ID','AVG_MEANS','STD_MEANS','AVG_RMS','STD_RMS',
             'NSAMPLE','NOUT_MEANS','CAVG_MEANS','CMDN_MEANS','CSTD_MEANS','NOUT_RMS','CAVG_RMS','CMDN_RMS','CSTD_RMS']
#
    ningest=ingest_sky_entries(img_dict,args.DBTable,DBorder,dbh,args.updateDB,verbose=verbose)


    exit(0)
        
    if (args.nite is not None):
#        try:
#            desdmfile = os.environ["des_services"]
#        except KeyError:
#            desdmfile = None
#        dbh = despydb.desdbi.DesDbi(desdmfile,args.section)
        query="""select 
            i.filename,
            i.expnum,
            i.ccdnum,
            oa.root,
            fai.path,
            fai.compression
        from proctag t, pfw_attempt_val av, image i, file_archive_info fai, ops_archive oa
        where t.tag='Y5N_FIRSTCUT'
            and t.pfw_attempt_id=av.pfw_attempt_id 
            and av.key='nite'
            and av.val='{nval:s}'
            and t.pfw_attempt_id=i.pfw_attempt_id
            and i.filetype='red_immask'
            and i.ccdnum=46
            and i.filename=fai.filename
            and fai.archive_name=oa.name
            """.format(nval=args.nite)

        iDict={}
        iDict=get_nite_img(iDict,query,dbh,dbSchema,Timing=args.Timing,verbose=verbose)

        imgList=[]
        for exp in iDict:
            if (iDict[exp]['compression'] is None):
                iDict[exp]['compression']=''
            imgList.append(iDict[exp]['root']+'/'+iDict[exp]['path']+'/'+iDict[exp]['filename']+iDict[exp]['compression'])

    print("Working from a list of {:d} images.".format(len(imgList)))

    LightBulbList=[46]
    LBD={46:{'xc':795,'yc':2620,'rad':500,'x1':595,'x2':995}}
    # form radial profile centered at 795,2620 (rad 350?) 
    # check columns 695,895 
#    DBorder=['EXPNUM','BULB','BULB_SIG','BULB_SB','RAD','RAD_OUTER','NUM_SAT','CMED','CSTD','Y_MEXTENT','Y_PEXTENT','X_WID','G_WID','L_WID','VG_WID','VL_WID','TIMING']
    DBorder=['EXPNUM','BULB','BULB_SIG','BULB_SB','RAD','RAD_OUTER','NUM_SAT','CMED','CSTD','Y_MEXTENT','Y_PEXTENT','X_WID',
             'G_AMP','G_WID','G_BKG','G_AMPERR','G_WIDERR','G_BKGERR','L_WID','VG_WID','VL_WID','TIMING','G_RAD',
             'Y1_AMP','Y1_WID','Y1_BKG','Y1_CEN','Y1_AMPERR','Y1_CENERR','X1_1','X1_2',
             'Y2_AMP','Y2_WID','Y2_BKG','Y2_CEN','Y2_AMPERR','Y2_CENERR','X2_1','X2_2']
#
    for img in imgList:
        print("Starting work on {:s}.".format(img))
        print("Time: {:.2f} ".format(time.time()-t00))


        ih,isci,imsk = get_data(img, verbose=verbose)
        print("Time: {:.2f} ".format(time.time()-t00))
        if (ih['CCDNUM'] in LBD):
            print("Checking for lightbulb on CCD {:d}".format(ih['CCDNUM']))
            iccd=ih['CCDNUM']

            if (('SKYVARA' in ih)and('SKYVARB' in ih)):
                print("SKYVAR A/B show expected noise level of: {:.3f} {:.3f} ".format(np.sqrt(ih['SKYVARA']),np.sqrt(ih['SKYVARB'])))
                print("   ")
            bulbDict=check_lightbulb(isci,LBD[iccd],verbose=verbose)
            bulbDict['expnum']=ih['EXPNUM']
            print("Time: {:.2f} ".format(time.time()-t00))

            if (bulbDict['isBulb']):
                bulbDict['bulb']='T'
            else:
                bulbDict['bulb']='F'
            print(bulbDict)

#            if (bulbDict['isBulb']):
            if(True):
                bulbDict=characterize_lightbulb(isci,imsk,LBD[iccd],bulbDict,qaplot=args.qa,verbose=verbose)

                bulbDict['bulb']='F'
                if (bulbDict['num_sat']>0):
                    if ((bulbDict['g_wid']>70.)and(bulbDict['g_amp']>100000)):
                        bulbDict['bulb']='S'
                if (bulbDict['bulb']=='F'):
                    if ((bulbDict['g_wid']>=70.)and(bulbDict['g_wid']<=140.)and
                        (bulbDict['g_widerr']>0.)and(bulbDict['g_widerr']<35.)):
                        if (bulbDict['g_amp']>=0.):
                            print(bulbDict['g_amp'])
                            if ((bulbDict['g_amperr']>0.)and(bulbDict['g_amperr']< np.sqrt(bulbDict['g_amp']))):
                                bulbDict['bulb']='T'

#                if (bulbDict['bulb']=='F'):
#                    if ((bulbDict['g_amp']>30.)and(bulbDict['g_wid']>70.65)and(bulbDict['g_wid']<140.)):
#                        bulbDict['bulb']='M'
                print(bulbDict)
                if (args.DBTable is not None):
                    bIngest={}
                    bIngest['1']=bulbDict
                    updateDB=True
                    ningest=ingest_bulb_entries(bIngest,args.DBTable,DBorder,dbh,updateDB,verbose=verbose)

        print("Finished work on {:s}".format(img))

    dbh.close()
 
    exit(0)
