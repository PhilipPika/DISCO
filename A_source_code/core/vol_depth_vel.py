# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/vol_depth_vel.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Based on volume_depth.py, but also returns volume and depth of water in the floodplain

# Import Python modules
import os
import math
import pickle

# Import own general modules
import ascraster
from error import *
from iround import *

# Import local modules.
import general_func
import get_dynamic_gridinfo
import interpolate_list


def calculate(params,mask,tmpdir,pointer1,outlakes_dict,jobid=None,loutputfile=False):
    '''
    This function calculates the volume and the depth of water in the main channel and floodplain in the grid cells
    '''
    print('Starting volume/depth calculations')
    water_root = params.water_inputdir
    
    # Read discharge for velocity calculations - LV 14-02-2018
    discharge = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.discharge),params.discharge_varname,temp_distrib=None,outputfile=None)

    # Read channel width grid for velocity calculations - LV 14-02-2018
    channel_width_grid = ascraster.Asciigrid(ascii_file=os.path.join(water_root,params.channel_width),mask=mask,numtype=float)
    
    # Read water storage of rivers and lakes/reservoirs
    water_storage = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.water_storage),params.water_storage_varname,temp_distrib=None,outputfile=None)
    
    # Read rivers and lakes/reservoirs areas as fraction of the cellarea
    fraction_waterbody_grid = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.fraction_water),params.fraction_water_varname,temp_distrib=None,outputfile=None)
    
    # Read water depth of flooding areas
    flooding_depth = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.flooding_depth),params.flooding_depth_varname,temp_distrib=None,outputfile=None)  
    
    # Read flooding areas as fraction of the cellarea
    flooding_fraction = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.flooding_fraction),params.flooding_fraction_varname,temp_distrib=None,outputfile=None)  
    
    # Read the cellarea - constant (contained in ASCII file)
    cellarea_grid = ascraster.Asciigrid(ascii_file=os.path.join(water_root,params.cellarea),mask=mask,numtype=float)
    
    # Read lake/reservoir id
    lakeid_grid = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.lakeid),params.lakeid_varname,temp_distrib=None,outputfile=None)
    
    # Read outlet cell of lake/reservoir id. Only outlet cell is filled
    outlakeid_grid = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.outlakeid),params.outlakeid_varname,temp_distrib=None,outputfile=None)
    
    # Read area of lakes/reservoirs
    waterbody_area_grid = get_dynamic_gridinfo.get_dynamic_gridinfo(params,pointer1,os.path.join(water_root,params.water_area),params.water_area_varname,temp_distrib=None,outputfile=None)
    
    # Read maximal depth of channels - constant (contained in ASCII file)
    channel_depth_grid = ascraster.Asciigrid(ascii_file=os.path.join(water_root,params.channel_depth),mask=mask,numtype=float)

    # Check if all read time series have the same length and if all variables are available at the same dates
    ndates = len(water_storage[0])
    #print('ndates  ' + str(ndates))

    #print('water_storage  ' + str(water_storage))
    #print('lakeid_grid  ' + str(lakeid_grid))
            
    if not(len(flooding_depth[0]) == ndates) or\
       not(len(flooding_fraction[0]) == ndates):
        raise MyError("Water storage, flooding fraction and flooding depth files should have the same number of time steps\n")
            
    for itime in range(ndates):
        if not(flooding_depth[0][itime][0] == water_storage[0][itime][0]) or\
           not(flooding_fraction[0][itime][0] == water_storage[0][itime][0]):
            raise MyError("Dates available in water storage, flooding fraction and flooding depth files are not all the same\n")
            
    # Make empty lists of volumes and depths for all cells
    vol_c = {}
    depth_c = {}
    vol_fp = {}
    depth_fp = {}
    vel = {} #LV 14-02-2018
    
    for item in range(len(pointer1)):
        icell = pointer1[item].get_local_index()
        vol_c[icell] = []
        depth_c[icell] = []
        vol_fp[icell] = []
        depth_fp[icell] = []
        vel[icell] = [] #LV 14-02-2018

    # Calculate volume and depth for each cell
    for item in range(len(pointer1)):
        #print 'icell  ' + str(icell)
        icell = pointer1[item].get_local_index()
        cellarea = cellarea_grid.get_data(icell,0.0)
        max_depth = channel_depth_grid.get_data(icell,-1)

        for itime in range(ndates):
            year = water_storage[0][itime][0]
            #print('year  ' + str(year))

            vol_waterbody = 0.
            vol_flooding_waterbody = 0.
            width = 0. #LV 14-02-2018

            # When there is a reservoir or a lake, apply volume and surface area values to outlet cell only
            #id = lakeid_grid[item][itime][1]
            id = iround((interpolate_list.find_closest(year, lakeid_grid[icell]))[-1])
            fraction_waterbody = (interpolate_list.calculate(year, fraction_waterbody_grid[icell], extrapol=1))[-1] #data not available monthly
            area_waterbody = cellarea * fraction_waterbody          
           
            if (id > 0):
                # It is a lake or reservoir
                # Check whether it is a outlet of a lake or reservoir
                #id_lake = outlakeid_grid[item][itime][1]
                #id_lake = (interpolate_list.find_closest(year, outlakeid_grid[icell]))[-1]
                #print(outlakes_dict)
                try: #LV 09-05-2018: there was a problem if it is an endorheic lake with no volume (no "outlet")
                    icell_out = interpolate_list.find_closest(year, outlakes_dict[id])[-1]
                    #LV 29-11-2017
                    # Total lake volume * swarea_cell / total lake swarea
                    vol_waterbody = water_storage[icell_out][itime][1]
                    tot_lake_area = (interpolate_list.calculate(year, waterbody_area_grid[icell_out], extrapol=1))[-1]
                    vol_waterbody *= area_waterbody / tot_lake_area 
                    width = math.pow(tot_lake_area,0.5) #LV 14-02-2018

                except:
                    vol_waterbody = 0.
                    width = 0.
                #print(icell, id, icell_out)

                '''
                #LV 29-11-2017
                # Total lake volume * swarea_cell / total lake swarea
                vol_waterbody = water_storage[icell_out][itime][1]
                tot_lake_area = (interpolate_list.calculate(year, waterbody_area_grid[icell_out], extrapol=1))[-1]
                vol_waterbody *= area_waterbody / tot_lake_area 
                width = math.pow(tot_lake_area,0.5) #LV 14-02-2018
                '''

            else:
                flood_frac = flooding_fraction[icell][itime][1]
                width = channel_width_grid.get_data(icell,0.0) #LV 14-02-2018
            
                if (flood_frac > 0.0):
                    if (max_depth > 0.0):
                        # Water body volume 
                        total_volume = water_storage[icell][itime][1]
                        volume_waterbody = max_depth * area_waterbody
                  
                        flood_volume = total_volume - volume_waterbody
                        if (flood_volume < 0.0):
                            # There is no flooding
                            flood_volume = 0.0
                            volume_waterbody = total_volume
                        
                        vol_waterbody = volume_waterbody
                        vol_flooding_waterbody = flood_volume
                    
                    else:
                        # We need another method for calculation the flooding volume and waterbody volume
                        # Now use flooding fraction instead of channel depth
                        total_volume = water_storage[icell][itime][1]
                        fl_depth = flooding_depth[icell][itime][1]
                        fl_frac = flooding_fraction[icell][itime][1]
                        flood_volume = fl_frac * fl_depth * cellarea
                 
                        volume_waterbody = total_volume - flood_volume
                    
                        if (volume_waterbody < 0.0):
                            # We assume there is no flooding
                            flood_volume = 0.0
                            volume_waterbody = total_volume
                        
                        vol_waterbody = volume_waterbody
                        vol_flooding_waterbody = flood_volume
                
                else:
                    vol_waterbody = water_storage[icell][itime][1]
                    vol_flooding_waterbody = 0.
    
            # Water depth of river in metres
            depth_waterbody = 0.0
            if (area_waterbody > 0.0):
                depth_waterbody = vol_waterbody / area_waterbody

            # Velocity - LV 14-02-2018
            velocity = 0.0
            if ((depth_waterbody*width) > 0.0):
                velocity = discharge[icell][itime][1] * 1.0e6 / (depth_waterbody * width)

            #print("velocity  " + str(velocity))

            # Change water volume from m3 into km3
            vol_waterbody *= 1.0e-9            
            vol_flooding_waterbody *= 1.0e-9
               
            vol_c[icell].append([year, vol_waterbody])
            depth_c[icell].append([year,depth_waterbody])
            vol_fp[icell].append([year, vol_flooding_waterbody])
            depth_fp[icell].append([year,flooding_depth[icell][itime][1]])

            vel[icell].append([year,velocity])

        #print("vel[icell]  " + str(vel[icell]))

    # If loutputfile, print in temporary pkl file
    if (loutputfile):
        fp_vol_c = open(os.path.join(tmpdir,"volume_c_"+str(jobid)+".pkl"),"wb")
        pickle.dump(vol_c,fp_vol_c,-1)
        fp_vol_c.close()
        fp_depth_c = open(os.path.join(tmpdir,"depth_c_"+str(jobid)+".pkl"),"wb")
        pickle.dump(depth_c,fp_depth_c,-1)
        fp_vol_fp = open(os.path.join(tmpdir,"volume_fp_"+str(jobid)+".pkl"),"wb")
        pickle.dump(vol_fp,fp_vol_fp,-1)
        fp_vol_fp.close()
        fp_depth_fp = open(os.path.join(tmpdir,"depth_fp_"+str(jobid)+".pkl"),"wb")
        pickle.dump(depth_fp,fp_depth_fp,-1)
        fp_depth_fp.close()
        #LV 14-02-2018
        fp_vel = open(os.path.join(tmpdir,"velocity_"+str(jobid)+".pkl"),"wb")
        pickle.dump(vel,fp_vel,-1)
        fp_vel.close()

    # Else return vol & depth objects
    else:
        return vol_c,depth_c,vol_fp,depth_fp,vel
