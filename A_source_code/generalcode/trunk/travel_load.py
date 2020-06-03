# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/travel_load.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import ascraster
import aggregate

def calculate(mouth_dict,basin,load,traveltime,sumtime,avgtime,sumload):
    '''
    This function calculates the the sum and the average of the product of residence_time and load.
    This must be an index for the travel time of the load in a river basin.
    sumtime,avgtime and sumload are attributes which are in mouth_dict. They must be given as text.
    '''

    loadtime = ascraster.duplicategrid(traveltime)

    # Calculate the travel time for the N load
    loadtime.multiply(load,default_nodata_value = 0.0)
    # Check whether sumtime is a valid attribute of the mouth_dict dictionairy.
    for key in list(mouth_dict.keys()):
        try:
            sumtime_riv = getattr(mouth_dict[key],sumtime)
        except AttributeError:
            # Make attribute in the mouth_dict dictionary
            setattr(mouth_dict[key],sumtime, 0.0)

    # Aggregate traveltime over the river basin.
    aggregate.aggregate_grid(basin,loadtime,mouth_dict,item=sumtime)

    # Calculate the average traveltime for each river basin
    for key in list(mouth_dict.keys()):
        try:
            sumtime_riv = getattr(mouth_dict[key],sumtime)
            sumload_riv = getattr(mouth_dict[key],sumload)
            setattr(mouth_dict[key],avgtime, sumtime_riv/sumload_riv)
        except ZeroDivisionError:
            setattr(mouth_dict[key],avgtime, 0.0)


