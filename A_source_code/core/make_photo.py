# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/make_photo.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import general python modules
import os
import sys
import traceback
import math
import optparse
import cProfile as profile
import pickle
import multiprocessing
from scipy import *
import scipy.integrate as itg

# Import own general modules
import ascraster
import interpolate
import my_sys
from error import *

# Import local modules.
import general_func

# Convert daynumbers to floats per year.
fdays = [i/365.0 for i in range(365)]
 
def make_one_photosyn(arg):
    '''
    Makes the photosythese for one specie which needs it.
    Calculate this for each day and for just 1 year.
    Calculation is only done, when this file does not exist.
    '''
    params = arg[0]
    species = arg[1]
    ispec = arg[2]
    specname = species[ispec].get_name()
    inputfile = os.path.splitext(params.photosyn_file)[0] + "_" + specname + ".pkl"
    outputfile = os.path.join(params.outputdir,os.path.splitext(os.path.basename(\
                              params.photosyn_file))[0] + "_" + specname + ".pkl")

    if (os.path.isfile(inputfile)):
        # Check input file whether this file is belonging to this run by comparing W, kmax, J1 and eta.
        fp = open(inputfile, 'rb')
        # Read header and check this
        header = pickle.load(fp)
        fp.close()

        # First value is the alpha and second value is the kmax, J1 and eta.
        lerror = (header[0] != species[ispec].get_alpha())
        lerror = lerror or (header[1] != species[ispec].get_kmax())
        lerror = lerror or (header[2] != params.J1)
        lerror = lerror or (header[3] != params.eta)
        if (lerror):
            raise MyError("Photosyn_file: " + inputfile + " is not correct.",\
                          " Values in the file: " + " ".join(map(str,header)),\
                          " Values expected   : " + str(species[ispec].get_alpha()) + " " +\
                          str(species[ispec].get_kmax()) + " " + str(params.J1) + " " + str(params.eta))
        # Copy file to outputdirectory
        my_sys.my_copyfile(inputfile,outputfile)
        print("File " + inputfile + " copied.")
    else:
        print("Start calculating the photosynthese file for specie: ",species[ispec].get_name())
        # Calculate the photsyn for each lattitude and each day and store this on the output directory.
        # Determine the cellsize of the grid by reading a grid inputfile
        basin = ascraster.Asciigrid(ascii_file=params.basin,numtype=int)
        cellsize = basin.cellsize
        yllmin = basin.yllcorner
        yllmax = yllmin + cellsize*basin.nrows
        del basin

        # Open output file
        fp = open(outputfile, 'wb')
        # Make header of outputfile
        header = [species[ispec].get_alpha(),species[ispec].get_kmax(),params.J1,params.eta]
        pickle.dump(header,fp, -1)
        yll = yllmax - 0.5 * cellsize
        while yllmin < yll:
            # First value is the lattitude
            values = [yll]
            for day in range(1,366):
                val = photosynday(params.depth,day,species[ispec].get_kmax(),\
                      species[ispec].get_alpha(),yll,params.eta,params.J1)
                values.append(val)
            pickle.dump(values,fp, -1)       
            # Go to next cell (lattitude)
            yll -= cellsize

        # Close file
        fp.close()

def make_photosyn_lat(params,species,lat_list,filename):
    '''
    Calculate the photosyn for a list of lattitudes.
    Output is written to a output file.
    '''
    # Open output file
    fp = open(filename, 'wb')
    for lat in lat_list:
        # First value is the lattitude
        values = [lat]
        for day in range(1,366):
            val = photosynday(params.depth,day,species[ispec].get_kmax(),\
                  species[ispec].get_alpha(),yll,params.eta,params.J1)
            values.append(val)
        pickle.dump(values,fp, -1)       

    # Close file
    fp.close()

def make_photosyn(params,species):
    lmultiproc = True
    if (lmultiproc):

        list_arg = []
        for ispec in params.totalphotindex:
            list_arg.append([params,species,ispec])

        p = multiprocessing.Pool()
        p.map(make_one_photosyn,list_arg)

    else:
        for ispec in params.totalphotindex:
            make_one_photosyn([params,species,ispec])

def photosynday(depth,day,kmax,alpha,lat,eta,J1):    
    '''
    Photosynday returns the average depth and time integrated photosyntesis for a given species for 1 day
    ''' 
    def tzintegrateday(t,z,kmax,alpha,w,daylight,day,lat,eta,J1):
        '''
        Photosyntesis equation calculation From Billen et al Hydrobiologia 289: 119-137, 1994 
        '''
        J2 = lat/90.    #seasonal changes    
        
        return kmax*(1.-math.exp(-alpha/kmax*math.exp(-eta*z)*math.pi/2.*J1*(1.-J2*math.cos(w*day))*math.sin(math.pi*t/(daylight))))  

    #Parameters used in the daylight calculation, after "Ecological Modeling", volume 80 (1995) pp. 87-95
    #These "O" numers should not be changed, they are needed to correctly calculate the daylight for each latitude and day
    O18 = 0.39795
    O19 = 0.2163108
    O20 = 0.9671396
    O21 = 0.0086
    O23 = 0.8333
 
    
    # Converts Julian days to radians   
    w = 2.*math.pi/365.        
    
    #Calculates the daylight hours for each Julian day    
    daylight = 24.-24./math.pi*math.acos(max(min((math.sin(O23*math.pi/180.)+math.sin(lat*math.pi/181.)*math.sin(math.asin(O18 \
        *math.cos(O19+2.*math.atan(O20*math.tan(O21*(day-186.)))))))/(math.cos(lat*math.pi/180.*math.cos(math.asin(O18 \
        *math.cos(O19+2.*math.atan(O20*math.tan(O21*(day-186.)))))))),1.),-1.))        

    # No daylight, no photosynthese.
    if (daylight < 0.01):
        return 0.0     
    #print('daylight: ', day, daylight)

    #limits of the temporal photosynthesis integration (Here we assume integration over the entire daylight hours)
    hrmin = 0.0
    hrmax = daylight
    
    f = itg.dblquad(tzintegrateday,0,depth,lambda x: hrmin,lambda x: hrmax,args=(kmax,alpha,w,daylight,day,lat,eta,J1))

    
    return f[0]/depth*365.             #return the photosynthesis term in yr-1   

def photo_average(values,time_start,time_end):
    '''
    This function delivers the average of the daily photosynthese of the time period time_start - time_end.
    Time_start and time_end are given in years. The list with photosynthese is given for each day and for only one year. 
    '''
    rel_time1 = time_start - int(time_start)
    rel_time2 = time_end - int(time_end)
    #print("REL TIME: ", rel_time1,rel_time2)

    if (time_end - time_start >= 1.0):
        avg_year = sum(values)/float(len(values))
        factor = float(int(time_end - time_start))
    else:
        avg_year = 0.0
        factor = 0.0

    #print("AVG_YEAR,FACTOR,TIMEDIFF: ",avg_year,factor,  time_end - time_start)
    # This is part of a year.
    if (rel_time1 == rel_time2):
        return interpolate.interpolate(rel_time1,fdays,values)
    elif (rel_time1 < rel_time2):
        #print("START SINGLE: ",len(fdays),len(values))
        avg_period = general_func.avg_within_range_2dim(fdays,values,rel_time1,rel_time2)
        #print("SINGLE: ",avg_period)
        try:
            return ((rel_time2-rel_time1)*avg_period + factor*avg_year)/(rel_time2 - rel_time1 + factor)
        except ZeroDivisionError:
            return general_func.sum_within_range_2dim(fdays,values,rel_time1,rel_time2) * 0.5
    else:
        # It is an average of two periods
        #print("START DOUBLE: ")
        period1 = general_func.avg_within_range_2dim(fdays,values,0.0,rel_time2)
        period2 = general_func.avg_within_range_2dim(fdays,values,rel_time1,fdays[-1]+0.0001)
        #print("DOUBLE: ",period1,period2)
        return (period1*rel_time2 + period2*(1.0-rel_time1) + factor*avg_year)/(1.0-rel_time1 + rel_time2 + factor)

if __name__ == '__main__':

    # Set the general path for the own python modules
    import general_path

    # Import own general modules
    from iround import *

    # Get the command line arguments.
    parfile = sys.argv[1]

    # Get the general settings.   
    fp = open(parfile,"rb")
    params =  pickle.load(fp)
    species = pickle.load(fp)
    proc =    pickle.load(fp)
    fp.close()


    values = list(range(1,366))
    print(photo_average(values,0.99/365,1.01/365), "Value expected: ",float(values[1]))
    print(photo_average(values,7.01/365,7.01/365), "Value expected: ",sum(range(8,10))/float(len(list(range(8,10)))))
    print(photo_average(values,7.01/365,7.9/365), "Value expected: ", sum(range(8,10))/float(len(list(range(8,10)))))
    print(photo_average(values,7.01/365,8.01/365), "Value expected: ", sum(range(8,11))/float(len(list(range(8,11)))))
    # One year
    print("One year average: ",photo_average(values,7.01/365,372.01/365), "Value expected: ",sum(values)/float(len(values)))
    print(photo_average(values,0.0,7.01/365), "Value expected: ",sum(range(1,10))/float(len(list(range(1,10)))))
    print(photo_average(values,360.01/365,1.0), "Value expected: ", sum(values[360:])/float(len(values[360:])))
    QQ = ((7.01/365.)*5 + (1.0 - 360.01/365.)*362.5 )/(1.0 - (360.01/365.) + (7.01/365.))
    print(photo_average(values,360.01/365,372.01/365), "Value expected: ",QQ)
    if (len(params.phytoindex) > 0):
        make_photosyn(params,species)

        # Checking
        basin = ascraster.Asciigrid(ascii_file=params.basin,numtype=int)
        specname = species[params.phytoindex[0]].get_name()
        outputfile = os.path.join(params.outputdir,os.path.splitext(os.path.basename(params.photosyn_file))[0] + "_" + specname + ".pkl")
        fp = open(outputfile,"rb")
        photo = list(general_func.pickleLoader(fp))
        print(len(photo), params.photosyn_file)
        fp.close()
        lats= [89.75,87.25,-60.25,0.25,-89.75]
        for lat in lats:
            ind = iround((90.0 - lat+0.5*basin.cellsize)/basin.cellsize)
            print(lat,photo[ind][0])
    else:
        print("No species with photosyntheses.")



    # Test average
    lat = 60.25
    ind = iround((90.0 - lat+0.5*basin.cellsize)/basin.cellsize)
    print("Average: ",photo_average(photo[ind][1:],0.9,0.9999999999))
