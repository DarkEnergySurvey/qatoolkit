import os,sys
import string,popen2,re
import math
from despyutils import tableio

"""
 Provides a colecction of tools to interact with ds9 from python using xpa.
 Author:
     Felipe Menanteau
"""

def prepare(nf, # Number of fields
            cmap     = 'Grey',
            scale    = 'linear',
            scalemode= '99.0',
            tile     = None,
            tilemode = 'column', # column,grid or row
            limits   = None,
            verb     = None):
            
    """ Starts up ds9 in XPA mode and several options """

    # Check env variables and set to proper values
    check_env(verb=verb);
    # Check if DS9 is already running
    check_ds9(verb=verb);
    # Start reseting
    os.system("xpaset -p ds9 scale "     + scale)
    os.system("xpaset -p ds9 scale mode" + scalemode)
    os.system("xpaset -p ds9 cmap "      + cmap)
    os.system("xpaset -p ds9 frame delete all")
    # Create the frames
    for i in range(nf):
        os.system("xpaset -p ds9 frame new")
        os.system("xpaset -p ds9 scale mode" + scalemode)
    # Tile up
    if tile:
        os.system ("xpaset -p ds9 tile ")
        os.system ("xpaset -p ds9 tile " + tilemode)

    return


def one_ellipse(shape,out,color='green',text='',options='',units='deg',system="image"):

    """ Draws just one ellipse of shape=(x,y,a,b,theta) on ds9
    -- redundant function that performs the task as ellipse() --- """

    if units == 'rad':
        theta  = shape[4]*180.0/math.pi
    else:
        theta  = shape[4]
    format = system +"; ELLIPSE(%8.2f,%8.2f,%7.2f,%7.2f,%6.1f) # color=%s text={%s} %s\n";
    out.write(format % (shape[0],shape[1],shape[2],shape[3],theta,color,text,options))
    return



def ellipse(shape,color='green',text='',options='',units='deg',out=None):

    """ Draw one ellipse of shape=(x,y,a,b,theta) on ds9 """

    if units == 'rad':
        theta  = shape[4]*180.0/math.pi
    else:
        theta  = shape[4]

    format = "image; ELLIPSE(%8.2f,%8.2f,%7.2f,%7.2f,%6.1f) # color=%s text={%s} %s\n";
    
    if out:
        reg = os.popen(out,'w')
        reg.write(format % (shape[0],shape[1],shape[2],shape[3],theta,color,text,options))
        reg.flush()
        reg.close()
    else:
        reg = os.popen('xpaset ds9 regions', 'w')
        reg.write(format % (shape[0],shape[1],shape[2],shape[3],theta,color,text,options))
        reg.flush()
        reg.close()
    return

def ellipses(shapes,color='green',text='',options='',units='deg',out=None):

    """ Draws ellipses on ds9 from a list of shapes with structure:
     shapes = [ (x1,y1,a1,b1,theta1),
                (x2,y2,a2,b2,theta2),
                ....
                (xn,yn,an,bn,thetan)]
    """

    format = "image; ELLIPSE(%8.2f,%8.2f,%7.2f,%7.2f,%6.1f) # color=%s text={%s} %s\n";
    
    check_env()
    check_ds9()
    
    if out:
        reg = os.popen(out,'w')
    else:
        reg = os.popen('xpaset ds9 regions', 'w')

    # For many shapes
    for shape in shapes:
        if units == 'rad':
            theta  = shape[4]*180.0/math.pi
        else:
            theta  = shape[4]
        reg.write(format % (shape[0],shape[1],shape[2],shape[3],theta,color,text,options))
            
        
    reg.flush()
    reg.close()
    return


def cat_ellipses(cat,IDs=None, color='green',label=None, text='',options='',units='deg',out=None):

    """
    Draw  all ellipses of shape=(x,y,a,b,theta) on ds9 from a SEx catalog dictionary,
    such as:
    (x,y,a,b,theta) = cat['SHAPE'][ID], where ID is ID of catalog
    """
    
    format = "image; ELLIPSE(%8.2f,%8.2f,%7.2f,%7.2f,%6.1f) # color=%s text={%s} %s\n";
    
    check_env()
    check_ds9()
    
    if out:
        reg = os.popen(out,'w')
    else:
        reg = os.popen('xpaset ds9 regions', 'w')
    
        
    if IDs == None:
        IDs = cat['SHAPE'].keys()

    for ID in IDs:
        
        shape = cat['SHAPE'][ID]
        
        if units == 'rad':
            theta  = shape[4]*180.0/math.pi
        else:
            theta  = shape[4]

        if label:
            text = ID

        reg.write(format % (shape[0],shape[1],shape[2],shape[3],theta,color,text,options))

    reg.flush()
    reg.close()
    return



def circles(x,y,radius=20,color='green',system='physical',text='',options=''):

    """ Draw circles at x,y on ds9 """


    format = "%s;circle(%8.2f,%8.2f,%s) # color=%s text={%s} %s\n"
    reg = os.popen('xpaset ds9 regions', 'w')

    for i in range(len(x)):
        reg.write(format % (system,x[i],y[i],radius,color,text,options))
    reg.flush()
    reg.close()
    return

def circles_fk5(x,y,radius=20,color='green',system='fk5',text='',options=''):

    """ Draw circles at x,y on ds9 """

    format = "%s;circle(%s,%s,%s\") # color=%s text={%s} %s\n"
    reg = os.popen('xpaset ds9 regions', 'w')

    for i in range(len(x)):
        reg.write(format % (system,x[i],y[i],radius,color,text,options))
    reg.flush()
    reg.close()
    return


def crosses(x,y,color='green',system='physical',text='',options=''):

    """ Draw an point cross at x,y on ds9 """
    if system == 'physical':
        format = "%s;point(%8.2f,%8.2f) # point=x color=%s text={%s} %s\n"
    else:
        format = "%s;point(%12.8f,%12.8f) # point=x color=%s text={%s} %s\n"
        
    reg = os.popen('xpaset ds9 regions', 'w')
    for i in range(len(x)):
        reg.write(format % (system,x[i],y[i],color,text,options))
    reg.flush()
    reg.close()
    return


def cross(x,y,color='green',system='physical',text='',options=''):

    """ Draw an point cross at x,y on ds9 """
    if system == 'physical':
        format = "%s;point(%8.2f,%8.2f) # point=x color=%s text={%s} %s\n"
    else:
        format = "%s;point(%12.8f,%12.8f) # point=x color=%s text={%s} %s\n"
        
    reg = os.popen('xpaset ds9 regions', 'w')
    reg.write(format % (system,x,y,color,text,options))
    reg.flush()
    reg.close()
    return



def circles_ids(x,y,ids,radius=20,color='green',system='physical',text='',options='',exclude=''):

    """ Draw circles at x,y on ds9 """

    format = "%s;%scircle(%f,%f,%s) # color=%s text={%s} %s\n"
    reg = os.popen('xpaset ds9 regions', 'w')

    for i in range(len(x)):
        reg.write(format % (system,exclude,x[i],y[i],radius,color,ids[i],options))
    reg.flush()
    reg.close()
    return


def circle(x,y,radius=20,color='green',system='physical',text='',options=''):

    """ Draws one circle at x,y on ds9 """

    format = "%s;circle(%f,%f,%s) # color=%s text={%s} %s\n"
    reg = os.popen('xpaset ds9 regions', 'w')
    reg.write(format % (system,x,y,radius,color,text,options))
    reg.flush()
    reg.close()
    return

def cross(x,y,color='green',system='physical',text='',options=''):

    """ Draw an point cross at x,y on ds9 """
    if system == 'physical':
        format = "%s;point(%8.2f,%8.2f) # point=x color=%s text={%s} %s\n"
    else:
        format = "%s;point(%12.8f,%12.8f) # point=x color=%s text={%s} %s\n"
        
    reg = os.popen('xpaset ds9 regions', 'w')
    reg.write(format % (system,x,y,color,text,options))
    reg.flush()
    reg.close()
    return


def frame(fr=1):
    ''' Go to frame fr'''
    os.system('xpaset -p ds9 frame ' + str(fr))
    return

def put(fitsfile,fr=None,zoom='yes',
        scale     = 'linear',
        scalemode = '99.0',
        cmap      = 'Gray',
        advance   = None,
        limits    = None):

    """ Display fits image into ds9 """

    #check_env()
    #check_ds9()

    if fr:
        job = popen2.Popen4("xpaset -p ds9 frame " + str(fr))
        job.wait()

        
    # Display the file
    job = popen2.Popen4("xpaset -p ds9 file "  + fitsfile)
    job.wait()

    if zoom:
        job_zoom = popen2.Popen4("xpaset -p ds9 zoom to fit " + str(fr))
        job_zoom.wait()

    if advance:
        os.system("xpaset -p ds9 frame next")

    if limits:
        os.system("xpaset -p ds9 scale limits %s %s" % (limits[0],limits[1]))

    job_scale = popen2.Popen4("xpaset -p ds9 scale "     + scale)
    job_scale.wait()

    job_mode = popen2.Popen4("xpaset -p ds9 scale mode " + scalemode)
    job_mode.wait()

    job_cmap = popen2.Popen4("xpaset -p ds9 cmap "      + cmap)
    job_cmap.wait()


    #else:
        #if verb:
        #    print "Scaling at",scalemode
        ##os.system("xpaset -p ds9 scale mode" + scalemode)
        #scalejob = popen2.Popen4("xpaset -p ds9 scale mode" + scalemode)
        #scalejob.wait()

    return

def close():
    """ Close the ds9 session and window """
    print "## Closing ds9 ";
    os.system ("xpaset -p ds9 exit")
    print "## Bye ...";
    return

def check_env(var="XPA_METHOD",value="local",verb=None):

    """ Check for the XPA_METHOD env variable to be set """

    if os.environ.has_key(var):
        if verb: print "# %s is already set to: %s" % (var, os.environ[var])
    else:
        os.environ[var] = value
        if verb: print "# %s is NOW set to: %s" % (var,value)

def check_ds9(verb=None):

    ''' Check if ds9 is running in -xpa mode, if not we start it and
    wait until is ready to listen'''
    
    ans = get_xpaaccess()
    
    if ans == "yes":
        if verb: print "# ds9 -xpa is running OK"
        answer = ans
        
    if ans == "no":
        if verb: print "# ds9 xpa is *NOT* running... trying to start it now"
        ds9job = popen2.Popen4("ds9 -xpa local")


        answer = get_xpaaccess()
        i = 1
        while answer != "yes":

            #print "waiting...", i, answer
            answer = get_xpaaccess()
            i = i + 1

        if ds9job.poll() == -1:
            if verb: print "# ds9 -xpa started OK"

    return answer


def get_xpaaccess():
    ''' To get the status of ds9 with xpa '''
    (stdout) = os.popen("xpaaccess ds9")
    ans = stdout.readline()
    ans = string.rstrip(ans)
    return ans


def inpath(program,verb=None):
    """ Checks if program is in the user's path """
    import os.path
    for path in os.environ['PATH'].split(':'):
        if os.path.exists( os.path.join(path,program) ):
            if verb: print "# program: %s found in: %s" % (program , os.path.join(path,program))
            return 1
    if verb: print "# program: %s NOT found in user's path " % program
    return 0


# Centroid selected regions
def centroid(fr=1):
    print >>sys.stderr,"Will try to centroid regions"
    os.system('xpaset -p ds9 frame %s' % fr)
    os.system('xpaset -p ds9 regions select all')
    os.system("xpaset -p ds9 regions centroid")
    os.system("xpaset -p ds9 regions select none")
    return

def read_and_parse_regions(counter=1):
    
    """ Read all the current drawn regions on ds9 and re-format the new x,y points
    as new regions in the ds9 format"""

    os.system("xpaset -p ds9 regions color blue")
    os.system("xpaset -p ds9 regions format ds9")
    (stdout,stdin) = popen2.popen2("xpaget ds9 regions")

    regexp_point = re.compile(r"(?P<type>physical;point)\("
                              r"(?P<x>[0-9]+),"
                              r"(?P<y>[0-9]+)\)\s*"
                              r"(?P<comments>#.+)?"
                              )
    x = []
    y = []
    i=0
    for line in stdout.readlines():
        
        # Check if point or polygon
        point = regexp_point.search(line)
        
        # Extract the information if a point was selected
        if point:
            type     = point.group('type')
            x.insert(i,point.group('x'))
            y.insert(i,point.group('y'))
            comments = point.group('comments')

            # Extract the comments and clean then up
            if (comments):
                comments = string.split(comments,"#")[1]
                comments = string.strip(comments)

            i = i + 1
            

    # After we've read all the lines and regions, we remove all the point
    # as regions are undeletable
    os.system("xpaset -p ds9 regions deleteall")

    # Format the coordinates and add comments
    comment = "color=blue text={region %d} delete = 0" % (counter)  
    regions = "physical;polygon(%s) # %s\n" % (pairs_into_list(x,y),comment)


    (stdout,stdin) = popen2.popen2("xpaset ds9 regions")
    #(io,w) = popen2.popen4("xpaset ds9 regions")
    stdin.write(regions)
    stdin.close()
    return

def split_into_pairs(list):

    """ To split a list into x,y pairs"""
 
    i=0
    n=0
    x = []
    y = []
    while (i<len(list)):
        x.insert(n,float(list[i]))
        y.insert(n,float(list[i+1]))
        n = n + 1
        i = i + 2
    return x,y

def pairs_into_list(x,y):

    """ Pairs into list"""

    i=0
    list = ""
    while (i<len(x)):
        if (i == 0):
            pair = "%s,%s"  % (x[i],y[i])
        else:
            pair = ",%s,%s" % (x[i],y[i])

        list = list + pair
        i = i + 1
    return list


