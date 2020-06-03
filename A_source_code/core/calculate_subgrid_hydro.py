# ******************************************************
## Revision "$LastChangedDate: 2018-11-19 11:14:53 +0100 (ma, 19 nov 2018) $"
## Date "$LastChangedRevision: 7 $"
## Author "$LastChangedBy: laurianevilmin $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/calculate_subgrid_hydro.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Local modules
import interpolate_list

year2second = 365*24*3600.0

def calculate_subgrid_hydro(lsteady,params,species,timeint,runoff_in,\
                            discharge_in,volume_in,water_area_in,\
                            depth_in,volume_fp_in,vel_in,depth_fp_in,\
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
        
    runoff_cell = []
    discharge_cell = []
    water_area_cell = []
    volume_cell = []
    depth_cell = []
    vel_cell = [] #LV 14-02-2018
    volume_fp_cell = [] #LV 06-02-2017 
    depth_fp_cell = []  #LV 06-02-2017 

    for time in timeint:
            
        runoff_cell.append(interpolate_list.calculate(time,runoff_in,extrapol=1))
        discharge_cell.append(interpolate_list.calculate(time,discharge_in,extrapol=1))
        water_area_cell.append(interpolate_list.calculate(time,water_area_in,extrapol=1))
        volume_cell.append(interpolate_list.calculate(time,volume_in,extrapol=1))        
        depth_cell.append(interpolate_list.calculate(time,depth_in,extrapol=1))
        vel_cell.append(interpolate_list.calculate(time,vel_in,extrapol=1)) #LV 14-02-2018
        volume_fp_cell.append(interpolate_list.calculate(time,volume_fp_in,extrapol=1))
        depth_fp_cell.append(interpolate_list.calculate(time,depth_fp_in,extrapol=1))

    # Calculate the hydrological properties for the lower Strahler orders
    Q0 = []
    Q = []
    for item in range(len(timeint)):
        # Q zero order in m3/s,area in km2, runoff is in mm/yr, so conversion of 1000.0
        runoff_local  = max(0.0,1000*runoff_cell[item][1]/year2second)
        #print 'runoff_local ' + str(runoff_local)
        Q0.append(params.A0 * runoff_local)
        # Q for all streams in m3/s
        Q.append([])
        for iorder in range(params.norder):
            Q[-1].append(riverarea[iorder] * runoff_local)

    # Determine the flux at the midpoint of the stream section [m3/s].
    Qmid = []
    depth = []
    vel = [] #LV 14-02-2018
    volume = []
    volume_fp = [] #LV 06-02-2017 
    depth_fp = []  #LV 06-02-2017

    Qmid_next = [None] * params.norder

    for item in range(len(timeint)):
        depth.append([])
        vel.append([]) #LV 14-02-2018
        Qmid.append([])   #for one stream only
        volume.append([]) #for one stream only
        volume_fp.append([None]*params.norder)
        depth_fp.append([None]*params.norder)

        # If the cell is a lake, set everything to zero, else calculate hydrology in small streams
        # LV added 21-02-2017
        if (llake):
            for iorder in range(params.norder):
                depth[-1].append(0.)
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
                    Qmid_next[iorder] = Q[item+1][iorder-1] + 0.5 * Q[item+1][iorder]

            for iorder in range(params.norder):
                # Stream width and depth for each order [m] and volume in m3.
                if (Qmid[item][iorder] < 0.0):
                    print("Qmid is negative: ",Qmid[item][iorder]," Runoff is: ", runoff_local)
                    print("***********************************************************")

                depth[-1].append(params.depth_factor * (Qmid[item][iorder]**params.depth_exponent))
                if (params.lmanning): #LV 01-03-2018

                    vel[-1].append(params.manning_ks*(depth[-1][-1]**(2./3.))*(slope[iorder]**0.5))
                else:
                    wetarea = (params.depth_factor * (Qmid[item][iorder]**params.depth_exponent)\
                               * params.width_factor * (Qmid[item][iorder]**params.width_exponent))
                    if (wetarea > 0.):
                        vel[-1].append(Qmid[item][iorder]/ wetarea) #LV 14-02-2018
                    else:
                        vel[-1].append(0.)
                volume[-1].append(params.width_factor * (Qmid[item][iorder]**params.width_exponent)\
                                    * params.depth_factor * (Qmid[item][iorder]**params.depth_exponent)\
                                    * riverlength[iorder] * 1000.) #conversion from km to m                           
            
            # Conversion of Qmid from m3/s to km3/yr and volume from m3 to km3
            for iorder in range(params.norder):
                Qmid[-1][iorder] *= year2second * 1.0e-9
                vel[-1][iorder] *= year2second * 1.0e-3
                volume[-1][iorder] *= 1.0e-9
                
    for iorder in range(params.norder):
        # Replace properties of sixth order streams by information of PCRGLOBWB
        if (iorder == params.norder-1):
            for item in range(len(timeint)):

                Qmid[item][iorder] = discharge_cell[item][-1]
                volume[item][iorder] = volume_cell[item][-1]
                depth[item][iorder]  = depth_cell[item][-1]

                vel[item][iorder]  = vel_cell[item][-1] #LV 14-02-2018
                volume_fp[item][iorder] = volume_fp_cell[item][-1] #LV 06-02-2017
                depth_fp[item][iorder] = depth_fp_cell[item][-1]   #LV 06-02-2017                

    return Qmid, volume, depth, vel, volume_fp, depth_fp
