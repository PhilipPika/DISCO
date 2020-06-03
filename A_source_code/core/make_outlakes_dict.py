# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/make_outlakes_dict.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import Python modules
import os

# Import own general modules
from iround import *

# Import local modules.
import get_dynamic_gridinfo
import interpolate_list

def calculate(params,pointer1):
    
    water_root = params.water_inputdir

    # Read lake/reservoir id
    lakeid_grid = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.lakeid),params.lakeid_varname,temp_distrib=None,outputfile=None)
    
    # Read outlet cell of lake/reservoir id. Only outlet cell is filled
    outlakeid_grid = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.outlakeid),params.outlakeid_varname,temp_distrib=None,outputfile=None)

    #
    water_storage_grid = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.water_storage),params.water_storage_varname,temp_distrib=None,outputfile=None) #LV 09-05-2018

    # Read endorheic lake ids
    endo_lakeid_grid = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.endo_lakes),params.endo_lakes_varname,temp_distrib=None,outputfile=None)

    outlakes_dict = {}

    ndates = len(lakeid_grid[0])

    #print(lakeid_grid)

    # Loop over all timesteps:
    for itime in range(ndates):
        year = lakeid_grid[0][itime][0]
    
        # Loop over all cells:    
        for item in range(len(pointer1)):
            icell = pointer1[item].get_local_index()
            # Check whether cell is a lake/reservoir.
            lakeid = iround(interpolate_list.stepwise(year,lakeid_grid[icell],extrapol=1)[-1])
            if (lakeid > 0):
                if not(lakeid in outlakes_dict):
                    outlakes_dict[lakeid] = []
                outlakeid = interpolate_list.stepwise(year,outlakeid_grid[icell],extrapol=1)[-1]
                endo_lake = interpolate_list.stepwise(year,endo_lakeid_grid[icell],extrapol=1)[-1]
                wst = interpolate_list.stepwise(year,water_storage_grid[icell],extrapol=1)[-1]
                llakeout = (outlakeid  > 0.)
                #lendolake = (endo_lake > 0.)
                lendolake = ((endo_lake > 0.) and (wst > 0.)) #LV 09-05-2018
                if (llakeout or lendolake):
                    outlakes_dict[lakeid].append([year,icell])
    #print(outlakes_dict)
            
    return outlakes_dict

