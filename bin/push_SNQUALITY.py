#! /usr/bin/env python
"""
Create/write FIRSTCUT_EVAL table entries based on SNQUALITY table.

Syntax:
    push_SNQUALITY.py -e expnum -l list -n night -N -m -y [-D filename] [-u] [-C comment] [-a analyst] [-s section] [-v]

    Must have one or more of the following:
       -e (expnum(s)),
       -l file containing a list of eposures,
       -n night
       -N (flag to look for NEW entries)

    Optional outputs:
       -D generates a file with INSERT statements which can be manually inserted in the database.
       -u causes the database update to occur as part of the program execution.

Arguments:
     
"""

if __name__ == "__main__":

    import argparse
    import os
    import despydb.desdbi
    import stat
    import time
    import re
    import sys
    import datetime
    import numpy
#    import scipy
    
    svnid="$Id$"
    svnrev=svnid.split(" ")[2]
#    db_table="firstcut_eval"

    parser = argparse.ArgumentParser(description='Obtain new SN quality assessments and push them into the FIRSTCUT_EVAL table.')
    parser.add_argument('-a', '--analyst',   action='store', type=str, default='SNQUALITY', help='Provides override value for analyst (default is SNQUALITY)')
    parser.add_argument('-C', '--Comment',   action='store', type=str, default=None,        help='Optional comment/explanation')
    parser.add_argument('-D', '--DB_file',   action='store', type=str, default=None,        help='Optional output of DB update file')
    parser.add_argument('-e', '--expnum',    action='store', type=str, default=None,        help='Restrict search to an exposure number (or comma delimited list of exposure numbers)')
    parser.add_argument('-l', '--list',      action='store', type=str, default=None,        help='File containing a list of exposure numbers')
    parser.add_argument('-n', '--night',     action='store', type=str, default=None,        help='Restict search to a specific night')
    parser.add_argument('-N', '--NewEntry',  action='store_true',      default=False,       help='Use SNQUALITY and FIRSTCUT_EVAL to figure out new entries')
    parser.add_argument('-m', '--month',     action='store_true',      default=False,       help='Restricts search for NEW entries to those from the last month.')
    parser.add_argument('-y', '--year',      action='store_true',      default=False,       help='Restricts search for NEW entries to those from the last year.')
    parser.add_argument('-u', '--updateDB',  action='store_true',      default=False,       help='Flag for program to DIRECTLY update DB (firstcut_eval).')
    parser.add_argument('-v', '--verbose',   action='store_true',      default=False,       help='Print progress messages to stdout')
    parser.add_argument('-s', '--section',   action='store', type=str, default=None,        help='Section of .desservices file with connection info')
    parser.add_argument('-S', '--Schema',    action='store', type=str, default=None,        help='DB schema (do not include \'.\').')
    args = parser.parse_args()
    if (args.verbose):
        print "##### Initial arguments #####"
        print "Args: ",args

    argsOK=True
    if ((args.expnum is None)and(args.night is None)and(args.list is None)and(not(args.NewEntry))):
        print " "
        print " "
        print "ERROR: Must provide one (or more) of the following."
        print " 1) a night (-n), "
        print " 2) an expnum (-e), "
        print " 3) a file containg a list of exposure numbers (-l), or "
        print " 4) use (-N) to specify a search for new entries in SNQUALITY (can be further constrained by other flags)"
        print "Aborting!"
        print " "
        parser.print_help()
        exit(1)

#
#   Automatically get analyst name (from os.getlogin()) unless an override value has been specified
#
    if ((args.analyst is None)or(args.analyst == "None")):
        analyst=os.getlogin()
    else:
        analyst=args.analyst
#
#   Obtain Schema (if user specified).
#
    if (args.Schema is None):
        db_Schema=""
    else:
        db_Schema="%s." % (args.Schema)

#
#   If an optional DB_file (containing INSERT commands) is chosen the set flag...
#
    if (args.DB_file is None):
        use_DB_file=False
    else:
        use_DB_file=True

#
#   Assessment to quality conversion
#
    AssessVsQual={'True': 2, 'False': 1, 'Unknown': 0}
    AssessAsInt={2: 'True', 1: 'False', 0: 'Unknown'}
    SNQualAsInt={2: 'Good', 1: 'Bad', 0: 'Unknown'}

#
#   Setup database connection information based on services file.
#
    try:
        desdmfile = os.environ["des_services"]
    except KeyError:
        desdmfile = None
    dbh = despydb.desdbi.DesDbi(desdmfile,args.section)
    cur = dbh.cursor()
    curw = dbh.cursor()

    db_table="firstcut_eval"
#
#   populate the expnum_list with values
#

    prelim_expnum_list=[]
    check_expnum_for_change=[]

    if (not(args.expnum is None)):
        tmp_list=args.expnum.split(',')
        for tmp_entry in tmp_list:
            if (tmp_entry.strip() != ''):
                prelim_expnum_list.append(int(tmp_entry.strip()))

    if (not(args.list is None)):
        if os.path.isfile(args.list):
            f1=open(args.list,'r')
            for line in f1:
                line=line.strip()
                columns=line.split(',')
                if (columns[0] != "#"):
                    prelim_expnum_list.append(int(columns[0]))
            f1.close()
#
#   Preliminary search using night sepecified to find supernova exposures
#
    if (not(args.night is None)):
        queryitems = ["sq.expnum"]
        coldict={}
        for index, item in enumerate(queryitems):
            coldict[item]=index
        querylist = ",".join(queryitems)
        query = """select %s from %ssnquality sq where sq.nite=%s """ % ( querylist, db_Schema, args.night )
        if args.verbose:
            print query
        cur.arraysize = 1000 # get 1000 at a time when fetching
        cur.execute(query)

        for item in cur:
            prelim_expnum_list.append(int(item[coldict["sq.expnum"]]))

#
#   Preliminary search for entries in SNQUALITY with/and without associated entries in FIRSTCUT_EVAL
#
    if (args.NewEntry):
        sq_restrict_date=''
        ex_restrict_date=''
        if ((args.year)or(args.month)):
            cur_night=datetime.date.today().strftime('%Y%m%d')
            year=int(cur_night[0:4])
            month=int(cur_night[4:6])
            day=int(cur_night[6:8])
            if (args.year):
                year=year-1
            if (args.month):
                month=month-1
                if (month < 1):
                    year=year-1
                    month=12
            sq_restrict_date=" and sq.nite>%4d%02d%02d" % (year,month,day)
            ex_restrict_date=" and e.nite>%4d%02d%02d" % (year,month,day)
        if (args.night):
            sq_restrict_date=" and sq.nite=%s" % (args.night)
            ex_restrict_date=" and e.nite=%s" % (args.night)
            

        queryitems = ["sq.expnum","sq.snqual"]
        coldict={}
        for index, item in enumerate(queryitems):
            coldict[item]=index
        querylist = ",".join(queryitems)
        query = """select %s from %ssnquality sq where snqual>0 %s """ % ( querylist, db_Schema, sq_restrict_date )
        if args.verbose:
            print query
        cur.arraysize = 1000 # get 1000 at a time when fetching
        cur.execute(query)

        snqual_expnum=[]
        snqual_state={}
        for item in cur:
            snqual_expnum.append(int(item[coldict["sq.expnum"]]))
            snqual_state[int(item[coldict["sq.expnum"]])]=int(item[coldict["sq.snqual"]])

        queryitems = ["f.expnum","f.accepted","f.analyst_comment"] 
        coldict={}
        for index, item in enumerate(queryitems):
            coldict[item]=index
        querylist = ",".join(queryitems)
        query = """select %s from %s%s f, %sexposure e where f.program='supernova' and f.analyst='SNQUALITY' and f.exposurename=e.filename %s order by f.lastchanged_time """ % ( querylist, db_Schema, db_table, db_Schema, ex_restrict_date )
        if args.verbose:
            print query
        cur.arraysize = 1000 # get 1000 at a time when fetching
        cur.execute(query)

        firstcut_eval_expnum=[]
        firstcut_eval_state={}
        firstcut_eval_field={}
        for item in cur:
            firstcut_eval_expnum.append(int(item[coldict["f.expnum"]]))
            firstcut_eval_state[int(item[coldict["f.expnum"]])]=AssessVsQual[item[coldict["f.accepted"]]]
            firstcut_eval_field[int(item[coldict["f.expnum"]])]=item[coldict["f.analyst_comment"]]

        for expnum in snqual_expnum:
            if (expnum not in firstcut_eval_expnum):
                prelim_expnum_list.append(int(expnum))
            else:
                if (not(firstcut_eval_state[expnum] == snqual_state[expnum])):
                    print("Found previous SNQUALITY entry ({:s}) different from FIRSTCUT_EVAL({:s}) for expnum={:d} ({:s})".format(SNQualAsInt[snqual_state[expnum]],AssessAsInt[firstcut_eval_state[expnum]],expnum,firstcut_eval_field[expnum]))
                    prelim_expnum_list.append(int(expnum))
#                check_expnum_for_change.append(int(expnum))
#                else:
#                    print "Found previous was same for ",expnum

    uniq_prelim=sorted(list(set(prelim_expnum_list)))
    prelim_expnum_list=uniq_prelim
    print "Formed exposure list for processing: ",len(prelim_expnum_list)," preliminary exposures found."
#
#   So you have a preliminary list of exposures numbers to work on.
#   Now retrieve their entries from the SNQUALITY table
#
    snqual_info=[]

    queryitems = ["sq.id","sq.nite","sq.expnum","sq.sequence","sq.snqual","sq.task_id","l.unitname","l.reqnum","l.attnum","sq.camsym"] 
    coldict={}
    for index, item in enumerate(queryitems):
        coldict[item]=index
    querylist = ",".join(queryitems)

    expnum_count=len(prelim_expnum_list)
    expnum_done=0

    while (expnum_done < expnum_count):
        expnum_limit=expnum_done+200
        if (expnum_limit > expnum_count):
            expnum_limit=expnum_count
        print "Querying for SNQUALITY assessment (up though first ",expnum_limit," expnums)"
        expnum_list_for_query=",".join(str(expnum) for expnum in prelim_expnum_list[expnum_done:expnum_limit])
        expnum_constraint="and sq.expnum in (%s)" % expnum_list_for_query
    
        query = """select %s from %ssnquality sq, %slatest_snsubmit l where sq.snqual>0 %s and sq.task_id!=0 and sq.task_id=l.task_id and l.status=0 order by sq.expnum """ % ( querylist, db_Schema, db_Schema, expnum_constraint )
        if args.verbose:
            print query
        cur.arraysize = 1000 # get 1000 at a time when fetching
        cur.execute(query)

        for item in cur:
            tmp_dict={}
            tmp_dict["expnum"]=int(item[coldict["sq.expnum"]])
            tmp_dict["night"]=int(item[coldict["sq.nite"]])
            tmp_dict["sequence"]=item[coldict["sq.sequence"]]
            tmp_dict["field"]=item[coldict["sq.sequence"]][0:2]
            tmp_dict["band"]=item[coldict["sq.sequence"]][2:3]
            tmp_dict["unitname"]=item[coldict["l.unitname"]].strip()
            tmp_dict["reqnum"]=int(item[coldict["l.reqnum"]])
            tmp_dict["attnum"]=int(item[coldict["l.attnum"]])
            tmp_dict["camsym"]=item[coldict["sq.camsym"]]
            tmp_dict["qual"]=int(item[coldict["sq.snqual"]])
            skip_it=False
            if (not(tmp_dict["band"] in ["u","g","r","i","z"])):
                print("Encountered malformed band ({:s}) in SNQUALITY entry for exposure number {:d}.  Skipping record.".format(tmp_dict["band"],tmp_dict["expnum"] ))
                skip_it=True
            if (not(tmp_dict["field"] in ["C1","C2","C3","E1","E2","S1","S2","X1","X2","X3"])):
                print("Encountered malformed field name ({:s}) in SNQUALITY entry for exposure number {:d}.  Skipping record.".format(tmp_dict["field"],tmp_dict["expnum"] ))
                skip_it=True
            if (not(skip_it)):          
                snqual_info.append(tmp_dict)
        expnum_done=expnum_limit
 
    print "Found complete SNQUALITY entries for: ",len(snqual_info)," exposure numbers"
#    for snqual_rec in snqual_info:
#        print snqual_rec
#    exit(0)

#
#    Prepare for optional DB update file output if args.DBupdate is TRUE
#
    if (use_DB_file):
        fdbout=open(args.DB_file, 'w')
#    print snqual_info

    queryitems = ["e.filename","e.expnum","e.band","e.camsym","e.object"]
    coldict={}
    for index, item in enumerate(queryitems):
        coldict[item]=index
    querylist = ",".join(queryitems)

    num_insert=0
    for exp_record in snqual_info:
        dm_process="True"
        if (exp_record["qual"]==1):
            dm_accept="False"
        else:
            dm_accept="True"
        field_constraint="%"+exp_record["field"]+"%"

        print exp_record["night"],exp_record["field"],exp_record["band"],dm_accept
    
        query = """select %s from %sexposure e where e.band='%s' and e.nite=%d and e.object like '%s' and e.camsym='%s' order by expnum """ % ( querylist, db_Schema, exp_record["band"],exp_record["night"],field_constraint,exp_record["camsym"]) 
#
#       Have supressed typical query here so that verbosity can be used to see what the script would push.
#
#        if args.verbose:
#            print query
        cur.execute(query)

        for item in cur:
#
#           Write DB_file (or update database) if command line option(s) are present 
#
            if((args.updateDB)or(use_DB_file)):
                db_cname="INSERT INTO %s%s(EXPOSURENAME,EXPNUM,CAMSYM,PROCESSED,ACCEPTED,PROGRAM,ANALYST,LASTCHANGED_TIME" % ( db_Schema, db_table )
                db_value=""" VALUES ('%s', %d, '%s', '%s', '%s', '%s', '%s', sysdate """ % (item[coldict["e.filename"]],item[coldict["e.expnum"]],item[coldict["e.camsym"]],dm_process,dm_accept,'supernova',analyst)
                if (args.Comment is None):
                    db_cname=db_cname+','+'ANALYST_COMMENT'
                    db_value=db_value+','+"'"+exp_record["sequence"]+" "+str(exp_record["night"])+"'" 
                else:
                    db_cname=db_cname+','+'ANALYST_COMMENT'
                    db_value=db_value+','+"'"+exp_record["sequence"]+" "+str(exp_record["night"])+" "+args.Comment+"'" 
                db_cname=db_cname+")"
                db_value=db_value+")"
                insert_command=db_cname+db_value
                num_insert=num_insert+1
                if (args.verbose):
                    if (num_insert == 1):
                        print db_cname
                    print db_value
    
                if(args.updateDB):
                    curw.execute(insert_command)
                if (use_DB_file):
                    fdbout.write("{:s};\n".format(insert_command))

            
    if(args.updateDB):
        dbh.commit()
        print "DB update complete and committed"
    if (use_DB_file):
        fdbout.close()
