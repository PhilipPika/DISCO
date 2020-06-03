# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# General python modules
import math

# Local modules
import interpolate_list

year2second = 365*24*3600.0

def calculate_subgrid_hydro(lsteady,params,species,timeint,yearstart,yearend,runoff_in,\
                            discharge_in,volume_in,water_area_in,\
                            depth_in,width_in,volume_fp_in,vel_in,depth_fp_in,\
                            llake,llakeout,lendolake,args_strahler_cell=None):

    # Get strahler characteristics
    if not(args_strahler_cell == None):
        riverarea = args_strahler_cell.get_val("riverarea")
        riverlength = args_strahler_cell.get_val("riverlength")
        slope = args_strahler_cell.get_val("slope")
    else:
        riverarea = params.riverarea
        riverlength = params.riverlength
        slope = [params.slope]*params.norder

    # construct list of times with startyear, endyear and timesteps
    timeint = []
    if (lsteady): # only for steady state calculation
        # steady state calculation is performed on start time.
        timeint.append(params.starttime)
    else:
        time = yearstart
        while time < yearend:
            timeint.append(time)
            time += params.timestep
        timeint.append(yearend)
        
    runoff_cell = []
    discharge_cell = []
    water_area_cell = []
    volume_cell = []
    depth_cell = []
    dvoldt_cell = [] 
    vel_cell = []
    volume_fp_cell = []
    depth_fp_cell = []
    dvoldt_fp_cell = []

    # volume at timestep 0 (= startyear)
    vol_t = interpolate_list.calculate(timeint[0],volume_in,extrapol=1) #for calculation of dvoldt
    vol_prev = vol_t
    vol_next = vol_t

    # volume of floodplain at timestep 0 (=startyear)
    vol_fp_t = interpolate_list.calculate(timeint[0],volume_fp_in,extrapol=1) #for calculation of dvoldt
    vol_fp_prev = vol_fp_t
    vol_fp_next = vol_fp_t

    # iterate through list of times (from startyear to endyear, with timesteps params.timestep)
    for time in timeint:
        runoff_cell.append(interpolate_list.calculate(time,runoff_in,extrapol=1))
        discharge_cell.append(interpolate_list.calculate(time,discharge_in,extrapol=1))
        water_area_cell.append(interpolate_list.calculate(time,water_area_in,extrapol=1))
        volume_cell.append(interpolate_list.calculate(time,volume_in,extrapol=1)) 
        depth_cell.append(interpolate_list.calculate(time,depth_in,extrapol=1))
        vel_cell.append(interpolate_list.calculate(time,vel_in,extrapol=1)) #LV 14-02-2018
        volume_fp_cell.append(interpolate_list.calculate(time,volume_fp_in,extrapol=1))
        depth_fp_cell.append(interpolate_list.calculate(time,depth_fp_in,extrapol=1))

        # Calculation of dvoldt
        dvoldt = 0.0
        # no dvoldt at steady state
        if (lsteady):
            dvoldt_cell.append([time,0.0])
        else:
            vol_next = interpolate_list.calculate(time+params.timestep,volume_in,extrapol=1)
            dvoldt = (vol_next[-1] - vol_prev[-1]) / 2*params.timestep
            dvoldt_cell.append([time,dvoldt])
            vol_prev = vol_t
            vol_t = vol_next

        # Calculation of dvoldt in floodplain
        dvoldt_fp = 0.0
        # no dvoldt in floodplain at steady state
        if (lsteady):
            dvoldt_fp_cell.append([time,0.0])
        else:
            vol_fp_next = interpolate_list.calculate(time+params.timestep,volume_fp_in,extrapol=1)
            dvoldt_fp = (vol_fp_next[-1] - vol_fp_prev[-1]) / 2*params.timestep
            dvoldt_fp_cell.append([time,dvoldt_fp])
            vol_fp_prev = vol_fp_t
            vol_fp_t = vol_fp_next

    # Calculate the hydrological properties for the lower Strahler orders
    Q0 = []
    Q = []
    for item in range(len(timeint)):
        # Q zero order in m3/s,area in km2, runoff is in mm/yr, so conversion of 1000.0
        runoff_local  = max(0.0,1000*runoff_cell[item][1]/year2second)
        Q0.append(params.A0 * runoff_local)
        # Q for all streams in m3/s
        Q.append([])
        for iorder in range(params.norder):
            Q[-1].append(riverarea[iorder] * runoff_local)

    # Determine the flux at the midpoint of the stream section [m3/s].
    Qmid = []
    depth = []
    dvoldt = []
    width = []
    vel = [] 
    volume = []
    volume_fp = []
    depth_fp = []
    dvoldt_fp = []

    vol_t = [None] * params.norder
    vol_prev = [None] * params.norder
    vol_next = [None] * params.norder
    Qmid_next = [None] * params.norder

    for item in range(len(timeint)):
        depth.append([])
        dvoldt.append([])
        width.append([])
        vel.append([]) #LV 14-02-2018
        Qmid.append([])   #for one stream only
        volume.append([]) #for one stream only

        volume_fp.append([])
        depth_fp.append([])
        dvoldt_fp.append([])

        # If the cell is a lake, set everything to zero, else calculate hydrology in small streams
        # LV added 21-02-2017
        if (llake):
            for iorder in range(params.norder):
                depth[-1].append(0.)
                dvoldt[-1].append(0.)
                width[-1].append(0.)
                vel[-1].append(0.) #LV 14-02-2018
                Qmid[-1].append(0.)
                volume[-1].append(0.)
                
        else:
            if (item == 0): #Qmid at 1st time step
                Qmid[-1].append(Q0[item] + 0.5 * Q[item][0])
                for iorder in range(1,params.norder):
                    Qmid[-1].append(Q[item][iorder-1] + 0.5 * Q[item][iorder])
            else:
                for iorder in range(0,params.norder):
                    Qmid[-1].append(Qmid_next[iorder])

            # Calculate Qmid at next time step
            if (item+1 < len(timeint)):
                Qmid_next[0] = Q0[item+1] + 0.5 * Q[item+1][0]
                for iorder in range(1,params.norder):
                    Qmid_next[iorder] = (Q[item+1][iorder-1] + 0.5 * Q[item+1][iorder])
          

            for iorder in range(params.norder):
                # Stream width and depth for each order [m] and volume in m3.
                if (Qmid[item][iorder] < 0.0):
                    print("Qmid is negative: ",Qmid[item][iorder]," Runoff is: ", runoff_local)
                    print("***********************************************************")
                width[-1].append(params.width_factor * (Qmid[item][iorder]**params.width_exponent))
                depth[-1].append(params.depth_factor * (Qmid[item][iorder]**params.depth_exponent))
                if (item == 0): #volume at 1st time step
                    vol_t[iorder] = params.width_factor * (Qmid[item][iorder]**params.width_exponent)\
                                    * params.depth_factor * (Qmid[item][iorder]**params.depth_exponent)\
                                    * params.riverlength[iorder] * 1000.

                volume[-1].append(vol_t[iorder])

                if (params.lmanning): #LV 01-03-2018

                    vel[-1].append(params.manning_ks*(depth[-1][-1]**(2./3.))*(slope[iorder]**0.5))
                else:
                    wetarea = (params.depth_factor * (Qmid[item][iorder]**params.depth_exponent)\
                               * params.width_factor * (Qmid[item][iorder]**params.width_exponent))
                    if (wetarea > 0.):
                        vel[-1].append(Qmid[item][iorder]/ wetarea) #LV 14-02-2018
                        #print('Qmid[item][iorder]=', Qmid[item][iorder])
                        #print('wetarea=', wetarea)
                        #print('vel=', vel[-1][-1]*3600)
                    else:
                        vel[-1].append(0.)

                # Calculate dvoldt
                if (lsteady):
                    dvoldt[-1].append(0.0)
                else:
                    if (item+1 < len(timeint)):
                        # Calculate volume at next time step
                        vol_next[iorder] = params.width_factor * ((Qmid_next[iorder]/params.number_of_rivers[iorder])**params.width_exponent)\
                                           * params.depth_factor * ((Qmid_next[iorder]/params.number_of_rivers[iorder])**params.depth_exponent)\
                                           * params.riverlength[iorder] * 1000.
                        dvoldt_i = (vol_next[iorder] - vol_t[iorder]) / params.timestep
                    else:
                        dvoldt_i = (vol_t[iorder] - vol_prev[iorder]) / params.timestep
                    dvoldt[-1].append(dvoldt_i) 
                    vol_prev[iorder] = vol_t[iorder]
                    vol_t[iorder] = vol_next[iorder]

            # Conversion of Qmid from m3/s to km3/yr and volume from m3 to km3
            for iorder in range(params.norder):
                Qmid[-1][iorder] *= year2second * 1.0e-9
                volume[-1][iorder] *= 1.0e-9
                dvoldt[-1][iorder] *= 1.0e-9
                # from m/s to km/yr
                vel[-1][iorder] *= year2second*1.0e-3
                


        # floodplains are only available for main stream
        volume_fp[item] = volume_fp_cell[item][-1]
        depth_fp[item] = depth_fp_cell[item][-1]
        dvoldt_fp[item] = dvoldt_fp_cell[item][-1]
                
    for iorder in range(params.norder):
        # Replace properties of sixth order streams by information of PCRGLOBWB
        if (iorder == params.norder-1):
            for item in range(len(timeint)):
                if (llake): #LV 29-11-2017
                    width[item][iorder]  = math.pow(water_area_cell[item][-1],0.5)
                else:
                    width[item][iorder]  = width_in
                Qmid[item][iorder] = discharge_cell[item][-1]
                volume[item][iorder] = volume_cell[item][-1]
                depth[item][iorder]  = depth_cell[item][-1]
                dvoldt[item][iorder]  = dvoldt_cell[item][-1]
                vel[item][iorder]  = vel_cell[item][-1] 

    return Qmid, volume, depth, width, vel, dvoldt, volume_fp, depth_fp, dvoldt_fp
