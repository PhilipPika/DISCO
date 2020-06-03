# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/interpolate_list.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import own general modules
from error import *

def calculate(x,list,extrapol=0):
    '''
    Interpolate for given x and returns the interpolated value of all the values of list.
    List is a list of lists with the first entry the time.
    extrapol == 1 then first or last given y value is return, when out of x-range.
    '''
    time = [x]
    if (x < list[0][0]):
        if (extrapol != 0):
            time.extend(list[0][1:])
            return time
        else:
            raise MyError ("Interpolation not possible for x: "+str(x)+" because list begins at: "+str(list[0][0]))
    elif (list[-1][0] < x):
        if (extrapol != 0):
            time.extend(list[-1][1:])
            return time
        else:
            raise MyError ("Interpolation not possible for x: "+str(x)+" because list ends at: "+str(list[-1][0]))
    else:

        # Interpolation
        xbegin = list[0][0]
        xend   = list[-1][0]
        ibegin = 0
        iend = len(list)-1
        for item in range(1,len(list)):
            if (x >= list[item][0]):
                xbegin = list[item][0]
                ibegin = item
            else:
                break
        for item in range(len(list)-2,0,-1):            
            if (x < list[item][0]):
                xend = list[item][0]
                iend = item
            else:
                break
        # Now interpolation
        if (ibegin == iend):
            time.extend(list[ibegin][1:])
            return time
        else:
            fac = (x - list[ibegin][0])/float((list[iend][0] - list[ibegin][0]))
            output = [x]
            for i in range(1,len(list[0])):
                output.append(list[ibegin][i]+(list[iend][i] - list[ibegin][i])*fac)
            return output

def stepwise(x,list,extrapol=0):
    '''
    Stepwise for given x and returns the last given value before the given point.
    List is a list of lists with the first entry the time.
    '''
    time = [x]
    if (x < list[0][0]):
        if (extrapol != 0):
            time.extend(list[0][1:])
            return time
        else:
            raise MyError ("Stepwise not possible for x: "+str(x)+" because list begins at: "+str(list[0][0]))
    elif (list[-1][0] < x):
        if (extrapol != 0):
            time.extend(list[-1][1:])
            return time
        else:
            raise MyError ("Interpolation not possible for x: "+str(x)+" because list ends at: "+str(list[-1][0]))
    else:
        # Search for the last time which given and return this values.
        ibegin = 0
        for item in range(1,len(list)):
            if (x >= list[item][0]):
                ibegin = item
            else:
                break

        # Now stepwise
        time.extend(list[ibegin][1:])
        return time


def find_closest(x,list):
    '''
    Closest value for given x and returns the interpolated value of all the values of list.
    List is a list of lists with the first entry the time.
    '''
    time = [x]
    if (x < list[0][0]):
        time.extend(list[0][1:])
    elif (list[-1][0] < x):
        time.extend(list[-1][1:])
    else:
        iclosest = 0
        for item in range(1,len(list)):
            if (x <= list[item][0]):
                if (abs(x-list[item][0]) < abs(x-list[iclosest][0])):
                    iclosest = item
                break
            else:
                iclosest = item 
        time.extend(list[iclosest][1:])

    #print('Find closest to ' + str(x) + ' in ' + str(list))
    #print(str(time))
    return time
