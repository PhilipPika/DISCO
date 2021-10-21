# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# General python modules
from copy import deepcopy
import math
import numpy as np
import operator
import os
import pickle
import scipy.integrate as itg
import scipy.optimize as opt
import time as time_module
import sys

# Import own general modules
import general_class
import my_sys

# Local modules
import calculate_subgrid_hydro
import general_func
import interpolate_list
import make_y0

# Import carbon version specific modules
#import global_radiation
import make_args
import mocsy
import reactions

# Conversion from year to seconds
year2second = 365*24*3600.0
second2year = 1.0/year2second
m3_s_to_km3_y = year2second * 1.0e-9

def add(x,y):
    return x+y

def positive(x):
    return max(x,0.0)


def calculate_cell(lsteady,lock,icell,params,species,proc,sources,tmpdir,next_cell,\
              yearstart,yearend,timeperiod,runoff_in,load_src,temperature_in,globrad_in,\
              discharge_in,volume_in,water_area_in,depth_in,width_in,slope,volume_fp_in,depth_fp_in,vel_in,low_veg_fr,high_veg_fr,\
              llake,llakeout,lendolake,args_strahler_cell=None):
    '''
    Calculates all processes in the subgrid in all stream orders.
    #JW comment: load_src[isrc][itime][1] is the multiple-source load for one cell, time is the 0 element in [[[time,value]]].
    '''
    ###print ("icell " + str(icell))

    # construct list of times with startyear, endyear and timesteps
    timeint = []
    if (lsteady):
        # Steady state calculation is perfomed on starttime
        timeint.append(params.starttime)
    else:
        # Make a list of all the time steps that we want to keep in years
        time = yearstart
        while time < yearend:
            timeint.append(time)
            time += params.timestep
        timeint.append(yearend)

    # Get strahler characteristics
    if not(args_strahler_cell == None):
        load_distribution = args_strahler_cell.get_val("load_distribution")
        flowpath_distribution = args_strahler_cell.get_val("flowpath_distribution")
        riverlength = args_strahler_cell.get_val("riverlength")
        number_of_rivers = args_strahler_cell.get_val("number_of_rivers")
    else:
        load_distribution = params.load_distribution
        flowpath_distribution = params.flowpath_distribution
        riverlength = params.riverlength
        number_of_rivers = params.number_of_rivers
    #print ("Number of rivers: " + number_of_rivers)


    # initialize load lists
    load_out = []
    load_in  = []

    # initialize file pointer
    fp_in = None

    # initialize arguments
    arguments = make_args.do()



    length_spec = len(species)

    # For each time point interpolate runoff, temperature and loads for all the sources (values in main channel)
    Qmid,volume,depth,width,vel,dvoldt,volume_fp,depth_fp,dvoldt_fp = \
      calculate_subgrid_hydro.calculate_subgrid_hydro(lsteady,params,species,\
                                                      timeint,yearstart,yearend,runoff_in,discharge_in,volume_in,water_area_in,\
                                                      depth_in,width_in,volume_fp_in,vel_in,depth_fp_in,llake,\
                                                      llakeout,lendolake,args_strahler_cell=args_strahler_cell)

    # PAY ATTENTION: the arguments related to floodplain, i.e. volumn_fp, depth_fp, dvoldt_fp are only for iorder=5, and are only one-level list: volumn_fp[i] instead of volumn_fp[i][j]
    #JW add the lines above 21-2-2021.
    # JW add the above 4 lines.
    #if icell in [XXXXXX]:
    #  print ("islake  " + str(llake))
    #  print ("Qmid  " + str(Qmid))
    #  print ("volume  " + str(volume))
    #  print ("depth  " + str(depth))
    #  print ("flow_velocity  " + str(vel))
    #  print ("dvoldt  " + str(dvoldt))


    # Get the previous load
    data_prev_period = []
    if (lsteady):
        # Get the start data from the species out of the class information
        for iorder in range(params.norder):
            data_prev_period.append(make_y0.make_y0(params,species))
            # JW modify and add the following 2 lines from DISC_C 24-2-2021, for floodplain 21-2-2021.
            if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False):
                data_prev_period[-1].extend(make_y0.make_y0(params,species))
        # data_prev_period is [[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],....,[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]]
    else:
        data_prev_period = []
        fp_in = open(os.path.join(tmpdir,str(icell)+"_"+str(timeperiod-1)+".pkl"),"rb")
        for line in general_func.pickleLoader(fp_in):
            # Remove the time
            data_prev_period.append(line[1:])
        fp_in.close()
        # data_prev_period is [[val1,val2,...,val8],....,[val51,val52,...,val58]]

    # Initialize budgets
    if (params.lbudget):
        xbud = general_func.initialize_budget(params, species, proc)
    # Output list of cell: [[time0,src0,  ..., srcn],[time1,src0,  ..., srcn],....,[time100,src0,  ..., srcn]]
    load_src_cell = []
    temperature_cell = []
    globrad_cell = []
    high_veg_fr_cell = []
    low_veg_fr_cell = []

    for time in timeint:
        load_src_cell.append(interpolate_list.calculate(time,load_src[0],extrapol=1))
        # JW comment: load_src is a list[isrc][time][1] for load value.
        # JW comment: load_src[0] is one-source (1st) load: [[time1,val1],...,[timen,valn]]
        # JW comment: interpolate(time,load_src[0]) is [time,val] for the 1st source.
        # JW comment: load_src_cell now is [[time01,val01]].
        for isrc in range(1,len(sources)):
            val = interpolate_list.calculate(time,load_src[isrc],extrapol=1) # [time,val] for the every source from the 2nd.
            load_src_cell[-1].append(val[1])
        # JW comment: load_src_cell now is [[time01,val01,val02,...,val0n]]. val is for all sources.
        # JW comment: load_src_cell finally is [[time01,val01,val02,...,val0n],...,[timex,valx1,valx2,...,valxn]].

        temperature_cell.append(interpolate_list.calculate(time,temperature_in,extrapol=1))
        globrad_cell.append(interpolate_list.calculate(time,globrad_in,extrapol=1))
        high_veg_fr_cell.append(interpolate_list.calculate(time,high_veg_fr,extrapol=1))
        low_veg_fr_cell.append(interpolate_list.calculate(time,low_veg_fr,extrapol=1))

        # Correct air temperature above freezing as a proxy for water temperature
        temperature_cell[-1][1] = max(0.0,temperature_cell[-1][1])

    # Read the load from other cells in Mmol/yr
    load_other_stream = []
    linput_other_stream = False
    if (os.path.isfile(os.path.join(tmpdir,"load_"+str(icell)+"_"+str(timeperiod)+".pkl"))):
        # Read the load
        fp_in = open(os.path.join(tmpdir,"load_"+str(icell)+"_"+str(timeperiod)+".pkl"),"rb")
        for line in general_func.pickleLoader(fp_in):
            load_other_stream.append(line)
        linput_other_stream = True
        fp_in.close()

    # Open output file for concentration and load at the end of the time period
    fp_out = open(os.path.join(tmpdir,str(icell)+"_"+str(timeperiod)+".pkl"),"wb")
    fp_conc = open(os.path.join(tmpdir,"conc_"+str(icell)+"_"+str(timeperiod)+".pkl"),"wb")

    # filepointer for the cell budget pkl
    fp_bud = None
    if (params.lbudget == 1):
        fp_bud = open(os.path.join(tmpdir,"budget_"+str(icell)+"_"+str(timeperiod)+".pkl"),"wb")

    # filepointer for the cell argument list pkl
    fp_args = None
    if (params.largs == 1):
        fp_args = open(os.path.join(tmpdir,"arguments_"+str(icell)+"_"+str(timeperiod)+".pkl"),"wb")

    # Solve the system for each Strahler order
    load_out = [] #load_out has is defined earlier together with load_in.

    #print ('flowpath_distribution ' + str(params.flowpath_distribution))

    for iorder in range(params.norder):
        # This is for floodplain, currently only used in carbon.
        if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False):
              floodplain_arguments = deepcopy(arguments)

        # Determine the load for this Strahler order
        load = []
        srcload = []
        upstrload = []
        # Take the local load from the cell (load_cell) and put the direct load
        # to the streams of this iorder in load array.
        # Load: [[time0,spec0,  ..., specn],[time1,spec0,  ..., specn],....,[time100,spec0,  ..., specn]]
        dt = 0. #LV 19-11-2018
        for item in range(len(timeint)): #item is the 0st, 1st,....,last time in the list of timeint.
            if (item > 0):
                dt = timeint[item] - timeint[item-1]  #JW add comment: dt ==0 means lsteady and it is for the starttime; dt>0 means not lsteady (for not starttime).

            # Start with time.
            # JW comment: load_src_cell is [[time01,val01,val02,...,val0n],...,[timex,valx1,valx2,...,valxn]]. Totally n sources.
            load.append([load_src_cell[item][0]]) #load now is [[time01]].
            srcload.append([load_src_cell[item][0]]) #srcload now is [[time01]].
            upstrload.append([load_src_cell[item][0]]) #upstrload now is [[time01]].
            # Initialize loads of all species to zero - modif LV 22-07-2016
            for j in range(length_spec):
                load[-1].append(0.0) #load now is [[time01,0.0]].
                srcload[-1].append(0.0) #srcload now is [[time01,0.0]].
                #The following 2 lines are for floodplains, currently used only in carbon.
                if iorder==params.norder-1 and params.lfloodplains==1:
                    load[-1].append(0.0)
                    srcload[-1].append(0.0)
                upstrload[-1].append(0.0) #upstrload now is [[time01,0.0]].
            # Calculate loads of species, depending on sources and order
            # If the cell is a lake, no small orders
            if (iorder <= params.norder-1) and (not llake): # LV added 21-02-2017
                for isrc in range(len(sources)):
                    for j in range(len(species)):
                        specname = species[j].get_name()
                        if (('fr_'+specname) in sources[isrc].get_attrib()):
                            if ('orders' in sources[isrc].get_attrib()):
                                try:
                                  minorder = int(sources[isrc].get_val('orders'))-1
                                except(ValueError):
                                  minorder = sources[isrc].get_val('orders')
                            try:
                              dummy_bool = (iorder >= (minorder))
                            except(TypeError):
                              dummy_bool = (minorder=='floodplain')
                            if dummy_bool:

                              # Apply fraction
                              load_spec = load_src_cell[item][isrc+1] * sources[isrc].get_val('fr_'+specname)

                              if load_spec is not np.ma.masked:
                                # If cell isn't a lake, apply order distribution
                                if (minorder == 0):
                                        load_spec *= params.load_distribution[iorder]
                                        load[item][j+1] += load_spec
                                if minorder!='floodplain':
                                        load_spec *= params.flowpath_distribution[minorder-1][iorder]
                                        load[item][j+1] += load_spec
                                if (iorder==params.norder-1) and (params.lfloodplains==1) and (minorder=='floodplain'):
                                        load[item][j+1+len(species)] += load_spec

                                if (params.lbudget) and (dt > 0.):
                                    general_func.add_budget_load(xbud, iorder, j, "load_src", load_spec * dt)
                                if (iorder==params.norder-1) and (params.lfloodplains==1) and (minorder=='floodplain'):
                                  general_func.add_budget_load(xbud, iorder+1, j, "load_src", load_spec * dt)

            elif (iorder == params.norder-1) and (llake): # LV added 21-02-2017 #JW comment, loads in lake cell in mainstream will only flow to mainstream.
                for isrc in range(len(sources)):
                    for j in range(len(species)):
                        specname = species[j].get_name()
                        if (('fr_'+specname) in sources[isrc].get_attrib()):
                            if ('orders' in sources[isrc].get_attrib()):
                                try:
                                  minorder = int(sources[isrc].get_val('orders'))-1
                                except(ValueError):
                                  minorder = sources[isrc].get_val('orders')
                            try:
                              dummy_bool = (iorder >= (minorder))
                            except(TypeError):
                              dummy_bool = (minorder=='floodplain')
                            if dummy_bool:
                                load_spec = load_src_cell[item][isrc+1] * sources[isrc].get_val('fr_'+specname)
                                if np.isnan(load_spec) or type(load_spec)!=np.float64:
                                  load_spec=0
                                if (iorder==params.norder-1) and (params.lfloodplains==1) and (minorder=='floodplain'):
                                    load[item][j+1+len(species)] += load_spec
                                    general_func.add_budget_load(xbud, iorder+1, j, "load_src", load_spec * dt)
                                else:
                                    load[item][j+1] += load_spec
                                    general_func.add_budget_load(xbud, iorder, j, "load_src", load_spec * dt)



	# Add all the output of the lower order Strahler to the list
        if not(llake): # LV added 21-02-2017
            for i in range(iorder):
                for item in range(len(timeint)):
                    for j in range(length_spec):
                        try:
                          outflow1 = flowpath_distribution[i][iorder] * load_out[i][item][j] * dt
                        except(FloatingPointError):
                          outflow1 = 0.
                        try:
                          load[item][j+1] += flowpath_distribution[i][iorder] * load_out[i][item][j]
                        except(FloatingPointError):
                          load[item][j+1] += 0.
                        if (params.lbudget) and (dt > 0.):
                          general_func.add_budget_load(xbud, iorder, j, "load_hw", outflow1)

        # Add all the output of the other streams (only sixth order)
        #JW comment: load_other_stream is the total load from other streams (not specified into orders, only changing with time), and it will flow to the mainstream (6th order) of running this cell.
        if (linput_other_stream and (iorder == params.norder-1)):
            for item in range(len(timeint)):
                for j in range(length_spec):
                    load[item][j+1] += load_other_stream[item][j]
                    if (params.lbudget) and (dt > 0.):
                        try:
                          outflow2 = load_other_stream[item][j] * dt
                        except(FloatingPointError):
                          outflow2 = 0.
                        general_func.add_budget_load(xbud, iorder, j, "load_up", outflow2)


        # Take amounts for iorder out of data_prev_period to proceed.
        y0 = []
        for i in range(len(data_prev_period[iorder])):

          try:
            y0.append(data_prev_period[iorder][i]/number_of_rivers[iorder])
          except(FloatingPointError):
            y0.append(0.)
        #y0 = np.divide(data_prev_period[iorder], number_of_rivers[iorder]).tolist()
        if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False) and (not len(y0)==2*len(species)):
            for i in range(len(data_prev_period[iorder+1])):
                try:
                    y0.extend([data_prev_period[iorder+1][i]/number_of_rivers[iorder]])
                except(FloatingPointError):
                    y0.extend([0.])
                except:
                    print('i=', i)
                    print('data_prev_period[iorder+1]=',data_prev_period[iorder+1])
                #y0.extend(np.divide(data_prev_period[iorder+1], number_of_rivers[iorder]).tolist())

        # Add load to load_out for this order.
        load_out.append([]) #JW comment: load_out is currently [[]] for iorder0 (1st element in the list) or [...,[[val1t0,...,val8t0],...,[val1t0,...,val8tx]],[]].

        if (lsteady):
            istart = 0
        else:
            # Fill first element of load_out. Is not used.
            load_out[-1] = [len(y0) * [0.0]] #JW comment: when not steady, load_out[-1] starts from [[0.,0.,0.,0.,0.,0.,0.,0.]] for each iorder, load_out is [[[v1...v8],...,[vx1...vx8]],...,[[0.,0.,0.,0.,0.,0.,0.,0. for iorder time0]]].
            istart = 1


        outdict = None # dictionnary with outputs from odeint

        for item in range(istart,len(timeint)):

            # Make time range for this run
            timerange=[timeint[item-1],timeint[item]] #JW comment: timerange = [starttime, starttime] (lsteady,i.e. istart=0); = [last_timepoint,timepoint] (not steady, i.e., istart=1).
            dt = timeint[item] - timeint[item-1] #JW comment: dt == 0 (lsteady,i.e.istart=0); dt>0 (not lsteady,i.e.istart=1).

            # LV 28-11-2017 - Determine whether it's an endorheic water body
            # Other lake/reservoir cells are calculated as rivers
            if (iorder == params.norder-1):
                if (llake):
                    # This cell is a lake/reservoir
                    if (lendolake):
                        # Retention is set to 1.0. So nothing is in the system. load_out for lake is zero.
                        load_out[-1].append(len(load[item][1:])*[0.0])
                        #JW comment: For endolake cell when iorder5, load_out[-1] is currently made [[0.,0.,0.,0.,0.,0.,0.,0.]] (lsteady), [[0.,0.,0.,0.,0.,0.,0.,0.],...,[0.,0.,0.,0.,0.,0.,0.,0. for time_item]] (not steady).
                        #JW comment: load_out is currently [[[v1...v8] for iorder0],...,[[0.,0.,0.,0.,0.,0.,0.,0. for iorder5]]] (lsteady), [[[0,0,0,0,0,0,0,0],[v1...v8],...,[vx1,...,vx8] for iorder0],...,[[0,0,0,0,0,0,0,0 for time0],...,[0,0,0,0,0,0,0,0 for time_item] for iorder5]].
                        for item1 in range(len(Y[-1])):
                            Y[-1][item1] = 0.0
                        # Set y0 for the next small timestep
                        y0 = deepcopy(Y[-1])

                        # Save total load (no load out)
                        if (params.lbudget) and (dt > 0.):
                            for ispec in range(len(species)):
                                general_func.add_budget_load(xbud, iorder, ispec, "loadIN", load[item][ispec+1] * dt)
                        continue

            #if (depth[item][iorder]==0.0) or volume[item][iorder]<params.minimal_watervolume or Qmid[item][iorder]<1e-6: #JW comment this line 24-2-2021. PP main difference is depth.
            if (depth[item][iorder]<0.001) or (volume[item][iorder]< params.minimal_watervolume) or (Qmid[item][iorder]<1e-6): # depth smaller than 0.001 m means dry.
            # Note that 1) compare two floats make computer systems don't understand, e.g. float_a == float_b can raise problems sometimes. We take 0.001 instead of 0.0.
                #print("CHECK: ",icell,timeint[item],volume[item][iorder],depth[item][iorder],Qmid[item][iorder])
                # There is no water volume, so skip this calculation.
                load_out[-1].append(list(map(positive,load[item][1:])))
                Y = [[]]
                for item1 in range(len(y0)):
                    Y[-1].append(0.0)
                Y = np.array(Y)
                # Set y0 for the next small timestep
                y0 = deepcopy(Y[-1])

                # Save total load (in and out)
                if (params.lbudget) and (dt > 0.):
                    for ispec in range(len(species)):
                        general_func.add_budget_load(xbud, iorder, ispec, "loadIN", load[item][ispec+1] * dt) #JW comment: has itime element in the lowest level, so ispec+1 for specie ispece; without highest level as iorder.
                        # general_func.add_budget_load(xbud, iorder, ispec, "load_src", srcload[item][ispec+1] * dt) #JW comment: same as above line.
                        # general_func.add_budget_load(xbud, iorder, ispec, "load_up", upstrload[item][ispec+1] * dt) #JW comment: same as above line.
                        # general_func.add_budget_load(xbud, iorder, ispec, "load_hw", hwload[item][ispec+1] * dt) #JW comment: same as above line.
                        general_func.add_budget_load(xbud, iorder, ispec, "loadOUT", load_out[-1][item][ispec] * dt) #JW comment: has no itime element in the lowest level, so ispec for specie ispece; with highest level as iorder, so this line means use the mainstream order.
                continue

            arguments.set_val("temperature",temperature_cell[item][-1]+params.tempcorrection)
            arguments.set_val("windspeed", params.windspeed)
            arguments.set_val("globrad", globrad_cell[item][-1])
            arguments.set_val("dvoldt",dvoldt[item][iorder])
            arguments.set_val("vol",volume[item][iorder])
            arguments.set_val("width",width[item][iorder])
            arguments.set_val("slope",slope*1e-5)
            arguments.set_val("discharge",Qmid[item][iorder])
            arguments.set_val("residence_time", arguments.get_val("vol")/Qmid[item][iorder])
            arguments.set_val("area", arguments.get_val("vol")/(depth[item][iorder]*1e-3))
            arguments.set_val("depth",depth[item][iorder])
            arguments.set_val("length", (arguments.get_val("area")*1e6)/width[item][iorder])

            arguments.set_val("icell", icell)
            arguments.set_val("next_cell", next_cell)
            arguments.set_val("flow_velocity", vel[item][iorder])

            if params.lsensitivity==1:
              arguments.set_val("temperature",(temperature_cell[item][-1]+params.tempcorrection)+params.temperature_add)
              arguments.set_val("windspeed", params.windspeed*params.windspeed_factor)
              arguments.set_val("globrad", globrad_cell[item][-1]*params.global_radiation_factor)
              arguments.set_val("slope",slope*1e-5*params.slope_factor)
              Qmid[item][iorder] = Qmid[item][iorder]*params.discharge_factor
              arguments.set_val("discharge",Qmid[item][iorder])
              arguments.set_val("residence_time", arguments.get_val("vol")/Qmid[item][iorder])
              arguments.set_val("flow_velocity", Qmid[item][iorder]/(arguments.get_val("width")*arguments.get_val("depth")*1e-6))

            if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False):
              floodplain_arguments = deepcopy(arguments)
              floodplain_arguments.set_val("dvoldt",dvoldt_fp[item])
              floodplain_arguments.set_val("vol",volume_fp[item])
              floodplain_arguments.set_val("depth",depth_fp[item])
              if (floodplain_arguments.get_val("vol") >0.) and (floodplain_arguments.get_val("depth")>0.):
                floodplain_arguments.set_val("area",floodplain_arguments.get_val("vol")/(floodplain_arguments.get_val("depth")*1e-3))
                floodplain_arguments.set_val("width", (floodplain_arguments.get_val('area')*1e6) / arguments.get_val('length'))

                # flow velocity in floodplain is assumed to be a fraction (between 0 and 1) of the flow velocity of the main stream
                floodplain_arguments.set_val("flow_velocity",vel[item][iorder]*params.fp_vel_fraction)
                floodplain_arguments.set_val("discharge",(floodplain_arguments.get_val("depth")*floodplain_arguments.get_val("width")*1e-6)*floodplain_arguments.get_val('flow_velocity'))
                floodplain_arguments.set_val("windspeed", params.windspeed*(params.fp_wind_reduc_factor*high_veg_fr_cell[item][-1]+(1-high_veg_fr_cell[item][-1])))
              else:
                  floodplain_arguments.set_val("area", 0)
                  floodplain_arguments.set_val("width", 0)
                  floodplain_arguments.set_val("length", 0)
                  floodplain_arguments.set_val("discharge",0)
                  floodplain_arguments.set_val("flow_velocity",0)
              if params.lsensitivity==1:
                  floodplain_arguments.set_val("flow_velocity",arguments.get_val("flow_velocity")*params.fp_vel_fraction*params.fp_vel_factor)
                  floodplain_arguments.set_val("discharge",(floodplain_arguments.get_val("depth")*floodplain_arguments.get_val("width")*1e-6)*floodplain_arguments.get_val('flow_velocity'))
                  high_veg_fraction = params.high_veg_fr_factor*high_veg_fr_cell[item][-1]
                  floodplain_arguments.set_val("windspeed", params.windspeed*(params.fp_wind_reduc_factor*params.windspeed_reduction_factor*high_veg_fraction+(1-high_veg_fraction)))

            setattr(params, 'debug', False)

            x_args = list() # wj201709
            for arg in sorted(arguments.get_attrib()): # wj201709
                x_args.append(arguments.get_val(arg)) # wj201709

            #### start steady state calc
            if timeint[0] == params.starttime and lsteady:


                load_reach = [] #LV 01-01-2017: loads only for one river reach (not total length of one order)





                for ispec in range(len(species)):
                    load_reach.append(load[item][ispec+1]/number_of_rivers[iorder])
                # for floodplains, currently used only in carbon. #if name == 'carbon':
                if params.lfloodplains==1 and iorder==params.norder-1:
                    for ispec in range(len(species)):
                        load_reach.append(load[item][ispec+1+len(species)]/number_of_rivers[iorder])
                timerange=[timeint[0],timeint[0]+params.outputtime]
                dt = timerange[-1] - timerange[0]
                for i_y in range(len(species)):
                    y0[i_y] = species[i_y].get_val('amount')*arguments.get_val('vol')
                    if ((params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False)): #carbon
                      y0[i_y+len(species)] = species[i_y].get_val('amount')*arguments.get_val('vol')





                    if params.lfloodplains==0 or iorder<params.norder-1 or llake==True:
                        Y1,outdict = itg.odeint(reactions.dy, y0,\
                                                timerange,\
                                                args=(params,species,proc,load_reach,\
                                                      Qmid[item][iorder],arguments), \
                                             full_output=True,printmessg=True,mxstep=5000)
                        Y = []
                        Y.append(list(map(positive,Y1[-1])))
                        # Make numpy array
                        Y = np.array(Y)



            #### Non-steady state case
            else:
                load_reach = [] #LV 01-01-2017: loads only for one river reach (not total length of one order)
                for ispec in range(len(species)):
                    load_reach.append(load[item][ispec+1]/number_of_rivers[iorder])
                # for floodplains, currently used only in carbon. #if name == 'carbon':
                if params.lfloodplains==1 and iorder==params.norder-1:
                    for ispec in range(len(species)):
                        load_reach.append(load[item][ispec+1+len(species)]/number_of_rivers[iorder])

                if params.lfloodplains==0 or iorder<(params.norder-1) or llake==True:
                    try:
                        Yy = itg.solve_ivp(reactions.dy2, timerange, y0, method='BDF',args=(params,species,proc,load_reach, Qmid[item][iorder],arguments))
                        Y = (Yy.y).T # means Yy.status == 0:  i.e. integration succeeds.
                        if Yy.status != 0:  # If the integration fails (usually when the step is alrealy smaller than the different of y0 values between steps), take the Y values of timeini[0] (last year or last month).
                            Y = []
                            Y.append(list(map(positive,y0)))
                            Y = np.array(Y)

                    except FloatingPointError: # If the integration fails (usually when the step is alrealy smaller than the different of y0 values between steps), take the Y values of timeini[0] (last year or last month).
                        Y = []
                        Y.append(list(map(positive,y0)))
                        Y = np.array(Y)

                    #  print('JW DEBUG',icell, Yy.message,(Yy.t[-1]-Yy.t[-2]),y0,Y[-1],[species0.get_name() for species0 in species])

                elif (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False):
                    print ("Solving 2: icell_"+str(icell) + ': iorder_' + str(iorder)) #JW 30-1-2021.
                    Y,outdict = itg.odeint(reactions.dy, y0,\
                                           timerange,\
                                           args=(params,species,proc,load_reach,\
                                                 Qmid[item][iorder],arguments,floodplain_arguments), full_output=True,printmessg=False,mxstep=8760, rtol=1e-3, atol=1e-3)

                else:
                    #print ("Solving 3: icell_"+str(icell) + ': iorder_' + str(iorder)) #JW 30-1-2021.
                    Y,outdict = itg.odeint(reactions.dy, y0,\
                                           timerange,\
                                           args=(params,species,proc,load_reach,\
                                                Qmid[item][iorder],arguments), full_output=1,printmessg=True,mxstep=8760, rtol=1e-6, atol=1e-6) #full_output=True, rtol=1e-3, atol=1e-3
            #JW comment: Y is finially a list of time series: Y[itime][ispec]
            #### Ending solver part



            setattr(params, 'debug', False)
            # Convert to load
            # Convert amount Y(Mmol) using Qmid(runoff, km3/yr) and volume (km3) into load (Mmol/yr)
            # Y[-1][0:len(species)] is used to account for the transporting in streams only  and exclude the floodplain transport in the extended Y (in the highest order only)
            load_out[-1].append([])
            for i in range(len(Y[-1][:len(species)])):
              try:
                load_out[-1][-1].append((Qmid[item][iorder]/volume[item][iorder])*Y[-1][:len(species)][i]*params.number_of_rivers[iorder])
                #load_out[-1][-1].append((Qmid[item][iorder]/arguments.get_val('vol')*Y[-1][:len(species)][i]*params.number_of_rivers[iorder])
              except(FloatingPointError):
                  load_out[-1][-1].append(0.)
            # Set outflow to zero for benthic species - LV 17-07-2017
            for ispec in range(len(species)):
                if (species[ispec].get_name().endswith("_benth")):
                    load_out[-1][item][ispec] = 0.

            #print ("order " + str(iorder+1))
            #print ("load_out = " + str(load_out))

            for i in range(len(load_out[-1][-1])):
                if load_out[-1][-1][i]<0.:
                    load_out[-1][-1][i] = 0

            # Set y0 for the next small timestep
            y0 = deepcopy(Y[-1])

            # Add budget fluxes for current internal timestep
            if (params.lbudget) and (not lsteady): #lbudget==1 and not lsteady (dt > 0.)
                for ispec in range(len(species)):
                    try:
                        load_dt = load[item][ispec+1] * dt
                    except(FloatingPointError):
                        load_dt = 0.
                    general_func.add_budget_load(xbud, iorder, ispec, "loadIN", load_dt)
                    try:
                        load_out_dt = load_out[-1][item][ispec] * dt
                    except(FloatingPointError):
                      load_out_dt = 0.
                    general_func.add_budget_load(xbud, iorder, ispec, "loadOUT", load_out_dt)
                    #general_func.add_budget_load(xbud, iorder, ispec, "dvoldt", arguments.dvoldt * Y[-1][ispec]/volume[item][iorder] * dt)
                    # to add: src_load for floodplains
                    if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False) and volume_fp[item]>0:
                      general_func.add_budget_load(xbud, iorder, ispec, "loadIN", load_reach[ispec+len(species)])
                      #if (Y[-1][ispec+len(species)]>0) and (llake==False) and ('benth' not in species[ispec].get_val('name')):
                      #  try:
                      #    general_func.add_budget_load(xbud, iorder+1, ispec, "dvoldt", floodplain_arguments.dvoldt * Y[-1][ispec+len(species)]/volume_fp[item] * dt)
                      #  except(FloatingPointError):
                      #    general_func.add_budget_load(xbud, iorder+1, ispec, "dvoldt", 0.)
                      #else:
                      #  general_func.add_budget_load(xbud, iorder+1, ispec, "dvoldt", floodplain_arguments.dvoldt * 0 * dt)

                proc_rates = reactions.procfunc(Y[-1][0:len(species)],params,species,proc,\
                                                Qmid[item][iorder],arguments) #LV 02-08-2017
                general_func.add_budget_procs(species, xbud, iorder, number_of_rivers[iorder], dt, proc_rates)

                # for floodplains, currently used only in carbon. #if name == 'carbon':
                if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False):
                    proc_rates_fp = reactions.procfunc(Y[-1][-len(species):],params,species,proc,\
                                                0,floodplain_arguments) #LV 02-08-2017
                    general_func.add_budget_procs(species, xbud, iorder+1, 1, dt, proc_rates_fp)
                else:
                    general_func.add_budget_procs(species, xbud, iorder, 1, dt, len(proc_rates)*[0])

        # Save last state situation
        # Y is a numpy array, so first make it a list.
        #x = Y[-1].tolist() # This line is the same to the line below.
        x = Y[-1][0:len(species)].tolist() #Y is array([[a,b,...,h]]) Y[-1] is array([a,b,...,h]), Y[-1][0:len(species)] is the same as Y[-1].
        xtot = []
        for i in range(len(Y[-1][0:len(species)])):
            try:
                xtot.append(number_of_rivers[iorder]*Y[-1][i])
            except(FloatingPointError):
                xtot.append(0.)

        if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False):
            x_floodplains = Y[-1][-len(species):].tolist()

        # Calculate concentration for all species for last order.
        # Change Mmol/km3 to mg/l is 10^9 mg/10^12 l * molarmass = 0.001 * molarmass as conversion
        # !!! For benthic species, xconc is not the concentration but amount in the benthic layer per volume of water column
        xconc=[]
        for ispec in range(len(x)):
            if (volume[-1][iorder] > params.minimal_watervolume):
                try:
                    xconc.append(1e-3 * x[ispec] * species[ispec].get_molarmass()/volume[item][iorder]) # item means time, but the loop above is finished earlier.
                except(FloatingPointError):
                    xconc.append(0.)
            else:
                xconc.append(0)
        if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False):
            xconc_floodplains=[]
            for ispec in range(len(x_floodplains)):
                if volume_fp[-1]>0:
                    try:
                        xconc_floodplains.append(0.001 * x_floodplains[ispec] * species[ispec].get_molarmass()/volume_fp[-1])
                    except(FloatingPointError):
                        xconc_floodplains.append(0.)
                else:
                    xconc_floodplains.append(0)

        # Add time to the list
        x.insert(0,timerange[-1])
        xtot.insert(0,timerange[-1])
        # Write to file
        pickle.dump(xtot,fp_out,-1) # -1 means the highest binary protocol.
        # for floodplains, currently used only in carbon. #if name == 'carbon':
        if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False):
            x_floodplains.insert(0,timerange[-1])
            pickle.dump(x_floodplains,fp_out,-1)
        elif (params.lfloodplains==1) and (iorder==params.norder-1):
            x_floodplains = [0]*len(species)
            x_floodplains.insert(0,timerange[-1])
            pickle.dump(x_floodplains,fp_out,-1)

        # Add time to the list
        xconc.insert(0,timerange[-1])
        pickle.dump(xconc,fp_conc,-1)
        if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False):
            xconc_floodplains.insert(0,timerange[-1])
            pickle.dump(xconc_floodplains,fp_conc,-1)
        elif (params.lfloodplains==1) and (iorder==params.norder-1):
            xconc_floodplains = [0]*len(species)
            xconc_floodplains.insert(0,timerange[-1])
            pickle.dump(xconc_floodplains,fp_conc,-1)


        # Save budget
        if (params.lbudget == 1):
            xbud[iorder].insert(0,timerange[-1]) # add the last timing at the beginning of the budget row.
            pickle.dump(xbud[iorder],fp_bud,-1)
            # for floodplains, currently used only in carbon. #if name == 'carbon':
            if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False): #floodplain (mainstream) & not lake
                xbud[iorder+1].insert(0,timerange[-1])
                pickle.dump(xbud[iorder+1],fp_bud,-1)
            elif (params.lfloodplains==1) and (iorder==params.norder-1): #floodplain (mainstream) & lake
                xbud[iorder+1] = [0]*len(general_func.get_all_header_budget_names(species, proc, line = []))
                xbud[iorder+1].insert(0,timerange[-1])
                pickle.dump(xbud[iorder+1],fp_bud,-1)

        # Save arguments
        if (params.largs == 1):
            x_args=[]
            for arg in sorted(arguments.get_attrib()):
                x_args.append(arguments.get_val(arg))
            # Add time to the list
            x_args.insert(0, timerange[-1]) # add the last timing at the beginning of the argument calculation row.
            # Write to file
            pickle.dump(x_args, fp_args,-1)
            # for floodplains, currently used only in carbon. #if name == 'carbon':
            if (params.lfloodplains==1) and (iorder==params.norder-1) and (llake==False):
                x_floodplain_args=[]
                for flpl_arg in sorted(floodplain_arguments.get_attrib()):
                    x_floodplain_args.append(floodplain_arguments.get_val(flpl_arg))
                x_floodplain_args.insert(0, timerange[-1])
                pickle.dump(x_floodplain_args, fp_args,-1)
            elif (params.lfloodplains==1) and (iorder==params.norder-1):
                x_floodplain_args=[]
                for arg in sorted(arguments.get_attrib()):
                    x_floodplain_args.append(0)
                x_floodplain_args.insert(0, timerange[-1])
                pickle.dump(x_floodplain_args, fp_args,-1)

    fp_out.close()
    fp_conc.close()
    if (params.lbudget == 1):
        fp_bud.close()
    if (params.largs == 1):
        fp_args.close()

    # Remove the previous timeperiod file of this cell. It is not needed anymore.
    filename = os.path.join(tmpdir,str(icell)+"_"+str(timeperiod-1)+".pkl")
    if (os.path.isfile(filename)):
        my_sys.my_removefile(filename)
    filename = os.path.join(tmpdir,"load_"+str(icell)+"_"+str(timeperiod-1)+".pkl")
    if (os.path.isfile(filename)):
        my_sys.my_removefile(filename)
    filename = os.path.join(tmpdir,"conc_"+str(icell)+"_"+str(timeperiod-1)+".pkl")
    if (os.path.isfile(filename)):
        my_sys.my_removefile(filename)
    if (params.lbudget == 1):
        filename = os.path.join(tmpdir,"budget_"+str(icell)+"_"+str(timeperiod-1)+".pkl")
        if (os.path.isfile(filename)):
            my_sys.my_removefile(filename)
    if (params.largs == 1): #wj201709
        filename = os.path.join(tmpdir,"arguments_"+str(icell)+"_"+str(timeperiod-1)+".pkl") #wj201709
        if (os.path.isfile(filename)): #wj201709
            my_sys.my_removefile(filename) #wj201709

    # Write load output of this cell to the output file.
    # Acquire the lock.
    if lock != None:
        lock.acquire()

    # Open output file for the load of the next cell.
    filename = os.path.join(tmpdir,"load_"+str(next_cell)+"_"+str(timeperiod)+".pkl")

    if (os.path.isfile(filename)):
        #print (" Load is added for cell from ", icell,"to ",next_cell)
        # Read the file and add this load to it.
        fp_out = open(os.path.join(tmpdir,"load_"+str(next_cell)+"_"+str(timeperiod)+".pkl"),"rb")
        load_other = []
        for line in general_func.pickleLoader(fp_out):

            load_other.append(line)
        fp_out.close()
        # Write it again.
        fp_out = open(os.path.join(tmpdir,"load_"+str(next_cell)+"_"+str(timeperiod)+".pkl"),"wb")
        for item in range(len(timeint)):
            added_list = list(map(add,load_other[item],load_out[-1][item]))
            pickle.dump(added_list,fp_out, -1)
        fp_out.close()
    else:
        #print(" Load is stored for cell from ", icell,"to ",next_cell)
        fp_out = open(os.path.join(tmpdir,"load_"+str(next_cell)+"_"+str(timeperiod)+".pkl"),"wb")

        for item in range(len(timeint)):
            pickle.dump(load_out[-1][item],fp_out, -1)
        fp_out.close()

    # Release lock
    if lock != None:
        lock.release()
