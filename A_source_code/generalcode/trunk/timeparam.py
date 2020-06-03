# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/timeparam.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
'''
timeparam.py

Module with the definition of the timeparam class.
'''

import interpolate
from error import *
from iround import *

# python modules (for moving average)
import numpy as np
import copy

       
class TimeParam:
    """ Data class with years and a value.""" 
    def __init__(self,year=None,val=None):
        try:
            self.year = self._str2num(year)
        except (ValueError,TypeError):
            self.year = None
        try:
            self.val = float(val)
        except (ValueError,TypeError):
            self.val = None
    
    def copy(self):
        return TimeParam(year=self.year,val=self.val)

    def get_year(self):
        return self.year
    
    def set_year(self,value):
        self.year = value
        
    def get_val(self):
        return self.val

    def set_val(self,value):
        self.val = value
    
    def write(self):
        print(("("+str(self.year)+","+str(self.val)+")"))

    def _str2num(self, s):
        '''
        checks if a string is a float or integer. Returns a float or integer.
        @param s: string with a numerical value
        @type s: STRING
        '''
        try:
            return int(s)
        except (ValueError,TypeError):
            return float(s)

def get_val(list1,year,extrapol=1):
    '''
    Returns the value for a given year in a list1 of TimeParam objects
    Interpolation is used, when the year is not in the list1.
    extrapol == 1 then first or last given value is return, when out of year-range.    
    '''
    
    years = []
    vals  = []
    if (len(list1) > 0):
        prev_year = float(list1[0].get_year())-1.
    for item in range(len(list1)):
        if (list1[item].get_val() != None):
            years.append(float(list1[item].get_year()))
            if (years[-1] < prev_year):
                raise MyError("TIMEPARAM.get_val: list is not sorted.")
            prev_year = years[-1]
            vals.append(list1[item].get_val())

    return interpolate.interpolate(year,years,vals,extrapol=extrapol)
    
def moving_average(list1,num_of_year,trendlineyears=0):
    '''
    Returns the moving average of the timeparam list. Num_of_year means the number of years around a year.
    The number in the list are also changing.
    So list from 1961 upto 2010 with an moving average of 5 years, which means 5 years in front, and 5 years after and the current year.
    will start at 1966 upto 2005.
    If trendlineyears is non-zero, then the num_of_years which normally is lost (begin and end of series) by calculating the moving average, 
    is added again by an trendline over the number of years that are given in trendlineyears.
    '''
    # It is possible that the original data is changed, therefore we copy the data
    list2 = copy.deepcopy(list1)
    
    if (trendlineyears > 0):
        # We should add num_of_year points at the beginning and the end of the dataserie.
        # Adding in the beginning of the dataserie
        xlist = []
        ylist = []
        beginyear = iround(list2[0].get_year()) 
        for year in range(beginyear,beginyear+iround(trendlineyears)):
            xlist.append(year)
            ylist.append(get_val(list2,year))
        # Determine the regressionline for this dataset
        xlist = np.array(xlist)
        ylist = np.array(ylist)
        #print "START xlist: ",xlist
        #print "START ylist: ",ylist
        slope,intercept = np.polyfit(xlist,ylist,1)
        #print "slope,intercept: ",slope,intercept
        # Add values to the list
        for item in range(num_of_year):
            list2.insert(0,TimeParam(year = beginyear-item-1,val = (beginyear-item-1)*slope + intercept))
            #print "ADDED: ",beginyear-item-1,(beginyear-item-1)*slope + intercept

        # Adding to the end of the dataserie
        xlist = []
        ylist = []
        end_year = iround(list2[-1].get_year()) 
        for year in range(end_year-iround(trendlineyears)+1,end_year+1):
            xlist.append(year)
            ylist.append(get_val(list2,year))
        # Determine the regressionline for this dataset
        xlist = np.array(xlist)
        ylist = np.array(ylist)
        #print "START xlist: ",xlist
        #print "START ylist: ",ylist
        slope,intercept = np.polyfit(xlist,ylist,1)
        #print "slope,intercept: ",slope,intercept
        # Add values to the list
        for item in range(num_of_year):
            list2.append(TimeParam(year = end_year+item+1,val = (end_year+item+1)*slope + intercept))
            #print "ADDED: ",end_year+item+1,(end_year+item+1)*slope + intercept
        
            
    beginyear = iround(list2[0].get_year())
    end_year = iround(list2[-1].get_year())
    outlist = []
    total = 0.0
    ftot = 2.* num_of_year + 1.
    for year in range(beginyear,beginyear+2*num_of_year+1):
        try:
            total += get_val(list2,year)
        except TypeError:
            pass
    outlist.append(TimeParam(year = beginyear+num_of_year,val = total/ftot))
    
    for year in range(beginyear+num_of_year+1,end_year-num_of_year+1):
        val_start = get_val(list2,year-num_of_year-1)
        val_end = get_val(list2,year+num_of_year)
        total += val_end - val_start
        #print ("MOV: ",year,val_start,val_end,total)
        outlist.append(TimeParam(year = year,val = total/ftot)) 

    return outlist
    
if __name__ == "__main__":
    
    list1 = [TimeParam(year=2000,val=2000)]
    list1.append(TimeParam(year=2001,val=2001))
    list1.append(TimeParam(year=2002,val=2002))
    list1.append(TimeParam(year=2003,val=2003))
    list1.append(TimeParam(year=2004,val=2004))
    list1.append(TimeParam(year=2005,val=2005))
    list1.append(TimeParam(year=2006,val=2006))
    list1.append(TimeParam(year=2007,val=2007))
    list2 = moving_average(list1,2)
    for item in range(len(list2)):
        list2[item].write()
