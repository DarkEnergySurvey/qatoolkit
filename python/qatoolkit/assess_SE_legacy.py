#! /usr/bin/env python
# $Id: assess_SE.py 44620 2016-11-17 17:22:35Z rgruendl $
# $Rev: 44620 $:  # Revision of last commit.
# $LastChangedBy: rgruendl $:  # Author of last commit.

"""
Legacy functions from previous versions of the assessment scripts
"""
from __future__ import print_function

import numpy 

#######################################################################
def MkCatDict(CATcol,cols_retrieve):
    """
    For use with FITSIO.  Make a dictionary that relates the order 
    of the requested columns to the order they are retreived by FITSIO.

    INPUT: CATcol is the list of columns from the FITS table
                  e.g. CATcol=cat_fits[hdu].get_colnames() 
           cols_retrieve is the list of column names requested from the FITS table.
          e.g. CAT=cat_fits[hdu].read(columns=cols_retrieve)

    OUTPUT: catdict

            for row in CAT:
                row[catdict["COLUMN_NAME"]]
            where "COLUMN_NAME" was an element in cols_retrieve and a column in CAT
 
    """
#
#   Figure out which column in CAT corresponds to each requested column
#   Get all columns in table (form dictionary showing their order)
#   Form a reverse dictionary based on columns retrieved
#   Then sort to get the order of the columns retieved in CAT and 
#   form a dictionary/legend that contains their order
#
    full_catdict={}
    for index, item in enumerate(CATcol):
        full_catdict[item]=index
    rev_catdict={}
    catdict_list=[]
    for item in cols_retrieve:
        rev_catdict[full_catdict[item]]=item
        catdict_list.append(full_catdict[item])
    icnt=0
    catdict={}
    for item in sorted(list(set(catdict_list))):
        catdict[rev_catdict[item]]=icnt
        icnt=icnt+1
#
    return catdict
####################################################################
def binsrch_crosscorr_cat(ra,dec,catalog,match_rad):
    """
    Search through a catalog for match to a set of coordinates (only nearest match is returned).

    ra        = RA used for comparison
    dec       = Dec used for comparison
    catalog   = is a list containing records(dictionaries) that have an "ra" and "dec" entries 
                along with associated information.  NOTE: This routine assumes that the catalog
                is RA-ordered so that a binary search can be used to find the small range of 
                entries that need to be checked
    match_rad = minimum distance for an acceptable match (in seconds of arc).
    """

    pi=3.141592654
    halfpi=pi/2.0
    deg2rad=pi/180.0

#
#   Search radius value
#
#    match_rad2=4.*(match_rad/3600.)*(match_rad/3600.)
    match_rad2=4.*(match_rad/3600.)

    ra1_deg=ra
    dec1_deg=dec
    sinl2=numpy.sin(halfpi-(dec1_deg*deg2rad))
    cosl2=numpy.cos(halfpi-(dec1_deg*deg2rad))
    cosdec1=numpy.cos(dec1_deg*deg2rad)
    num_rec=len(catalog)

#
#   Work out maximum RA range that will need to be searched
#
#    print("DRA range calculation: {:.7f} {:.7f} ".format((match_rad2/cosdec1),(match_rad2/cosdec1*3600.0)))
    if (abs(dec1_deg)>89.0):
        ra_range1=0.0
        ra_range2=360.0
    else:
        ra_range1=ra1_deg-(match_rad2/cosdec1)
        if (ra_range1 > 360.):
            ra_range1=ra_range1-360.0
        if (ra_range1 < 0.0): 
            ra_range1=ra_range1+360.0
        ra_range2=ra1_deg+(match_rad2/cosdec1)
        if (ra_range2 > 360.):
            ra_range2=ra_range2-360.0
        if (ra_range2 < 0.0): 
            ra_range2=ra_range2+360.0
#
#   Create a list of RA ranges that will be searched
#   (two ranges are necessary for initial binary search)
#
    ra_range_list=[]
    if (ra_range1 > ra_range2):
#       print("Need two ranges for search around 0h-24h branch ",ra1_deg,dec1_deg)
        ra_range_list.append([0.0,ra_range2])
        ra_range_list.append([ra_range1,360.0])
#       print(ra_range_list)
    else:
#       print("Single range will suffice for ",ra1_deg,dec1_deg)
        ra_range_list.append([ra_range1,ra_range2])
#
#   Use each range to search for a counterpart.
#
    match_record=[]
    rad_dist=[]
    for ra_range in ra_range_list:
#
#   Initiate a binary search to find the appropriate entries in the catalog
#
        klo=0
        khi=num_rec-1
        while (khi-klo > 1):
            knew=(khi+klo)/2
#           print(klo,khi,knew)
            if (catalog[knew]["ra"] > ra_range[0]):
                khi=knew
            else:
                klo=knew
#
#       Should have initial spot to begin a serious comparison
#       Loop through catalog values until RA exceeds upper limit of range
#           or the end of the catalog has been reached
#
#       print(klo,khi,catalog[klo]["ra"])
        ncheck=0
        while ((klo < num_rec-1)and(catalog[klo]["ra"]<ra_range[1])):
            ncheck=ncheck+1
            record=catalog[klo]
            ra2_deg=record["ra"]
            dec2_deg=record["dec"]
            dra=abs(ra2_deg-ra1_deg)/cosdec1
            ddec=abs(dec2_deg-dec1_deg)
            if ((dra < match_rad2)and(ddec < match_rad2)):
#               print(ra2_deg,dec2_deg)
                ang1=halfpi-(dec2_deg*deg2rad)
                ang2=deg2rad*(ra2_deg-ra1_deg)
                if (ang2 < 0.0):
                    ang2=-ang2
                val=numpy.cos(ang1)*cosl2+numpy.sin(ang1)*sinl2*numpy.cos(ang2)
                if (val > 1.):
                    dist=0.0
                else:
                    val2=numpy.sqrt(1.0-val*val)
                    dist=3600.0*(halfpi-numpy.arctan2(val,val2))/deg2rad
                if (dist <= match_rad):
                    if (len(rad_dist) > 0):
                        if (dist < rad_dist[0]):
                            rad_dist=[]
                            match_record=[]
                            match_record.append(record)
                            rad_dist.append(dist)
                    else:
                        match_record.append(record)
                        rad_dist.append(dist)
            klo=klo+1
        if (ncheck > 1000):
            print("Warning possible brute force use detected (RA,DEC,ncheck,range): {:.7f},{.7f},{:d},{:s}".format(ra1_deg,dec1_deg,ncheck,map(str,ra_range)))
#       if (ncheck == 0):
#           print("Warning no checks made (RA,DEC,ncheck,range): ",ra1_deg,dec1_deg,ncheck,ra_range)

    return match_record
