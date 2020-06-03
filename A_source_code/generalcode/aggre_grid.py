# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/aggre_grid.py $"
# ******************************************************
'''
Test script to test the functionality of ascraster.
'''
import os
import sys
__general = r"/home/beusena/globalnutrients/generalcode/trunk"
if os.path.exists(__general): 
    sys.path.insert(0, __general)
    print(__general + " is added to the python search path for modules.") 

import ascraster
import my_sys

isocodes_grid_file = sys.argv[1]
value_grid_file = sys.argv[2]
outputfile = sys.argv[3]

isogrid = ascraster.Asciigrid(ascii_file=isocodes_grid_file)
value_grid = ascraster.Asciigrid(ascii_file=value_grid_file)

dict_out = {}
for icell in range(isogrid.length):
    iso = isogrid.get_data(icell)
    if (iso != None):
        iso = int(iso+0.5)
        val = value_grid.get_data(icell,0.0)
        try:
            dict_out[iso] += val
        except KeyError:
            dict_out[iso] = val

fp = open(outputfile,"w")
for key in sorted(dict_out):
    fp.write(str(key) + ";" + str(dict_out[key]) + "\n")
fp.close()

