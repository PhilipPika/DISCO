# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/interpolate.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
from error import *

def interpolate(x,xlist,ylist,extrapol=0):
    '''
    Interpolate for given x and returns the interpolated value of y.
    extrapol == 1 then first or last given y value is return, when out of x-range.
    '''
    if (len(xlist) != len(ylist)):
           raise MyError ("Interpolation not possible because length of xlist differs from length ylist.")

    if (x < xlist[0]):
        if (extrapol != 0):
            return ylist[0]
        else:
            raise MyError ("Interpolation not possible for x: "+str(x)+" because list begins at: "+str(xlist[0]))
    elif (xlist[-1] < x):
        if (extrapol != 0):
            return ylist[-1]
        else:
            raise MyError ("Interpolation not possible for x: "+str(x)+" because list ends at: "+str(xlist[-1]))
    else:
        # Interpolation
        xbegin = xlist[0]
        xend   = xlist[-1]
        ibegin = 0
        iend = len(xlist)-1
        for item in range(1,len(xlist)):
            if (x >= xlist[item]):
                xbegin = xlist[item]
                ibegin = item
            else:
                break
        for item in range(len(xlist)-2,0,-1):            
            if (x < xlist[item]):
                xend = xlist[item]
                iend = item
            else:
                break
        # Now interpolation
        if (ibegin == iend):
            return ylist[ibegin]
        else:
            return ylist[ibegin] + (x - xlist[ibegin])*(ylist[iend] - ylist[ibegin])/float((xlist[iend] - xlist[ibegin]))
