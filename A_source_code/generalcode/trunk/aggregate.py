# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/aggregate.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

from error import *
from iround import *
import sys
                      
def aggregate_grid(grid_ids,grid_val,dict,item,method=1,weighing=None,min_value=None):
    '''
    grid_ids is ascraster object with id's for each cell
    grid_val is ascraster object with values for each cell
    If method == 1: Summation
    If method == 2: Average
    If method == 3: Weighing summation with raster "weighing"
    If method == 4: Weighing Average with raster "weighing"
    If min_value is given, then the aggregation is only performed on the cells which 
    have a value which is higher than min_value
    '''
    lweighing = (method == 3 or method == 4)
    
    # Check input
    if (lweighing):
        if (weighing == None):
            # The function call is not correct
            raise MyError("AGGREGATE_GRID: Error in function call. Method 3 or 4 needs a weighing grid.")   
    
    dict_out = {}
    if (lweighing): dict_weigh={}

    if (method == 2):
        dict_count = {}
 
    if (min_value == None):
        lmin_value = False
    else:
        lmin_value = True

    try:
        # Check length of the two lists:
        if (grid_ids.length != grid_val.length):
            raise MyError("AGGREGATE_GRID: Average is not possible on two grids of different length.",\
                           "Length of first list: "+str(grid_ids.length),\
                           "Length of second list: "+str(grid_val.length))

        if (lweighing):
            if (grid_ids.length != weighing.length):
                raise MyError("AGGREGATE_GRID: Length of weighing map is not correct.",\
                              "Length of value list: "+str(grid_ids.length),\
                              "Length of weighing list: "+str(weighing.length))
                           
        for i in range(grid_ids.length):
            ids = grid_ids.get_data(i)
            if ids == None:
                continue
            # Make from ids an integer.
            try:
                ids = iround(ids)
            except ValueError:
                pass
            val = grid_val.get_data(i)
            if (val == None):
                continue
            elif (lmin_value):
                if (val < min_value):
                    continue
            if (lweighing):
                weigh = weighing.get_data(i)
                if (weigh == None):
                    if (val == 0.0):
                        # In this case it does not matter.
                        weigh = 0.0
                    else:
                        print("Weighing map has nodata values on essential places.")
                        raise MyError("Weighing map is not correct. There are nodata values on essential places.")
                if (weigh < 0.0):
                    print("Weighing map has negative values on essential places.")
                    raise MyError("Weighing map is not correct. There are negative values on essential places.")

            key = str(ids)
            value = float(val)
            try:
                if (method == 1 or method == 2):
                    dict_out[key] += value
                if (method == 2): dict_count[key] += 1
                if (lweighing): 
                    dict_out[key] += value * weigh
                    if (method == 4): dict_weigh[key] += weigh
            except KeyError:
                # Make a new entry in dict_out
                if (method == 1 or method == 2):
                    dict_out[key] = value
                if (method == 2): dict_count[key] = 1
                if (lweighing): 
                    dict_out[key] = value * weigh
                    if (method == 4): dict_weigh[key] = weigh

            
        if (method == 2):
            # Calculate average:
            for key in list(dict_out.keys()):
                if (dict_count[key] > 0):
                    dict_out[key] /= dict_count[key]
        elif (method == 4):
            # Weighed average
            for key in list(dict_out.keys()):
                if (dict_weigh[key] > 0):
                    dict_out[key] /= dict_weigh[key]
            
        for key in list(dict.keys()):
            try:
                setattr(dict[key],str(item),dict_out[str(key)])
            except KeyError:
                pass

    except MyError as val:
        val.write()
        raise MyError()       
    except:
        raise MyError("AGGREGATE_GRID: Something goes wrong with error type: " + str(sys.exc_info()[0]),\
                      "ERROR in AGGREGATE_GRID and item " + str(item))


def aggre_list(list_values,method=1,weighing=None,lwarning=True):
    '''
    aggre_list aggregates the list of values to one value.
    If method == 1: Summation
    If method == 2: Average
    If method == 3: Weighing summation with "weighing" list
    If method == 4: Weighing Average with "weighing" list
    If method == 5: Choose dominant class (values must be integers (not checked), otherwise integer is returned!).
    list_values and weighing are lists.
    lwarning: when set to True a raise is given when list_values contains None values.
    When it is not possible to aggregate the list, then None is returned.
    '''
    lwarn = False
    # Set weighing on or off
    lweighing = (method == 3 or method == 4)

    # Check input
    if (lweighing):
        if (weighing == None):
            # The function call is not correct
            raise MyError("AGGRE_LIST: Error in function call. Method 3 or 4 needs a weighing list.")   
    
        if (len(list_values) != len(weighing)):
            raise MyError("AGGRE_LIST: Length of weighing list is not correct.",\
                          "Length of values list: "+str(len(list_values)),\
                          "Length of weighing list: "+str(len(weighing)))

    # Do the calculation.
    total = 0.0
    if (lweighing): 
        total_weighing = 0.0

    if (method == 2):
        total_count = 0
    if (method == 5):
        dict_key = {}

    for item in range(len(list_values)):
        value = list_values[item]
        if (value == None):
            # Nodata value in list_values.
            lwarn = True
            continue
        if (lweighing):
            weigh = weighing[item]
            if (weigh == None):
                if (value == 0.0):
                    # In this case it does not matter.
                    weigh = 0.0
                else:
                    print("Weighing map has nodata values on essential places.")
                    raise MyError("Weighing list is not correct. There are nodata values on essential places.")
            if (weigh < 0.0):
                print("Weighing map has negative values on essential places.")
                raise MyError("Weighing map is not correct. There are negative values on essential places.")

        if (method == 1 or method == 2):
            total += value
            if (method == 2): 
                total_count += 1
        if (lweighing): 
            total += value * weigh
            if (method == 4): 
                total_weighing += weigh
        if (method == 5):
            # Put value as integer in the dictionary
            try:
                dict_key[iround(value)] += 1
            except KeyError:
                dict_key[iround(value)] = 1

    # Calculate the average.
    if (method == 2):
        # Calculate average:
        if (total_count > 0):
            total /= float(total_count)
        else:
            total = None
    elif (method == 4):
        # Weighed average
        if (total_weighing > 0):
            total /= float(total_weighing)
        else:
            total = None
    elif (method == 5):
        max_count = None
        max_key = None
        for key in dict_key:
            if (max_count == None):
                max_count = dict_key[key]
                max_key = key
            else:
                if (dict_key[key] > max_count):
                    max_count = dict_key[key]
                    max_key = key
        if (max_key == None):
            total = None
        else:
            total = max_key
            
    if (lwarn and lwarning):
        raise MyError("AGGRE_LIST: Some cells of values list have nodata values.")

    return total


# Test script
if (__name__ == "__main__"):

    try:
        values = [1,2,3,4,5,6,7,8,9,10]
        print("Expected value 55.0, returned: ",aggre_list(values))
        print("Expected value 5.5, returned: ",aggre_list(values,method=2))
        values = [10,2,2,4,5,5,7,8,10,10]
        print("Expected value 10, returned: ",aggre_list(values,method=5))
        values   = [10,2,2,4 ,5,5,6,8,10,10]
        weighing = [2 ,2,2,10,2,2,5,9,8 ,8 ]
        print("Expected value 350.0, returned: ",aggre_list(values,method=3,weighing=weighing))
        print("Expected value 7.0, returned: ",aggre_list(values,method=4,weighing=weighing))

        # Check error messages
        try:
            # Check lenght of lists
            a = [10,11]
            aggre_list(values,method=3,weighing=a)
            print("This is wrong!!!!")
        except MyError as val:
            val.write()
        try:
            # Check without weighing
            aggre_list(values,method=3)
            print("This is wrong!!!!")
        except MyError as val:
            val.write()
        try:
            # None in weighing list
            b = [2 ,2,2,None,2,2,5,9,8 ,8 ]
            aggre_list(values,method=3,weighing=b)
            print("This is wrong!!!!")
        except MyError as val:
            val.write()
        try:
            # Negative weighing
            c = [2 ,2,2,-1.9,2,2,5,9,8 ,8 ]
            aggre_list(values,method=3,weighing=c)
            print("This is wrong!!!!")
        except MyError as val:
            val.write()
        try:
            # None in list_values
            d = [10,2,2,None ,5,5,6,8,10,10]
            aggre_list(d,method=1)
            print("This is wrong!!!!")
        except MyError as val:
            val.write()
        try:
            # None in list_values
            d = [10,2,2,None ,5,5,6,8,10,10]
            print("Expected value 58.0, returned: ",aggre_list(d,method=1,lwarning=False))
            print("This is wrong!!!!")
        except MyError as val:
            val.write()

    except MyError as val:
        val.write()
    except:
        import traceback
        traceback.print_exc()
