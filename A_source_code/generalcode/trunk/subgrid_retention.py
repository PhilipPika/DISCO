# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/subgrid_retention.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
import os
import math

import aggregate
import ascraster

# Note that the following import takes another module for N or P!
import accuflux_retention

# Conversion from year to seconds
year2second = 365*24*3600.0
second2year = 1.0/year2second

def define_subgrid_streamorder(params):
    '''
    This function defines the constant factors for the subgrid_streamorder (cell independent).
    All output parameters are set in the params object.
    The order streams (one to six) are stored in a list.
    The stream orders 1 to params.norder are defined, where the last order is the first stream that is modelled by PCRGLOBWB.
    '''
    # Calculate the length of the streams (km)
    params.riverlength = [params.L1]
    for iorder in range(1,params.norder):
        params.riverlength.append(params.riverlength[-1]*params.Rl)

    # Calculate the area of the streams (km2)
    params.riverarea = [params.A1]
    for iorder in range(1,params.norder):
        params.riverarea.append(params.riverarea[-1]*params.Ra)

    # Calculate the number of the streams, first create the list and after that fill the list.
    params.number_of_rivers = []
    for iorder in range(params.norder):
        params.number_of_rivers.append(1.)

    for iorder in range(params.norder-2,-1,-1):
        params.number_of_rivers[iorder] = params.number_of_rivers[iorder+1]*params.Rb

    # Calculate the distribution of N or P load for each stream order as porpotion river*length of the total riverlength
    total_length = []
    for iorder in range(params.norder):
        total_length.append(params.riverlength[iorder] * params.number_of_rivers[iorder])
    params.load_distribution = [] 
    for iorder in range(params.norder):
        params.load_distribution.append(total_length[iorder]/sum(total_length))

    # Calculate the flow path probabilities from stream order i to order j as proportionof total length in order j relative to 
    # the total length of all stream orders greater than i. First create a matrix and after that fill the matrix (only the upper part).
    params.flowpath_distribution = []
    for iorder in range(params.norder):
        params.flowpath_distribution.append(params.norder*[0.0])
    for iorder in range(params.norder):
        sumlength = 0.0
        for iorder2 in range(iorder+1,params.norder):
            sumlength += total_length[iorder2]
        for iorder2 in range(iorder+1,params.norder):
            #print iorder,iorder2,total_length[iorder2],sumlength,total_length[iorder2]/sumlength
            params.flowpath_distribution[iorder][iorder2] = total_length[iorder2]/sumlength

    params.A0 = params.A1/float(params.Ra)

    #print "vf = 35.0,pnet=500",calculate_localretention(params,35.,500)
    #print "vf = 45.0,pnet=500",calculate_localretention(params,45.,500)
    #print "vf = 35.0,pnet=1500",calculate_localretention(params,35.,1500)
    #print "vf = 35.0,pnet=100",calculate_localretention(params,35.,100)
    #print "LENGTH",params.riverlength
    #print "AREA",params.riverarea
    #print "NUMBER",params.number_of_rivers
    #print "LOAD DISTRIBUTION",params.load_distribution
    #for iorder in range(params.norder):
    #    print "FLOWPATH DISTRIBUTION ",iorder," ",params.flowpath_distribution[iorder]

def calculate_localretention(params,runoff,vf,load_org,temperature):
    '''
    Calculate the total retention in the subgrid in the stream orders 1 upto 5.
    Runoff in mm/yr. Return value is the total retention of this subgrid.
    vf is dependent on temperature and for N also on concentration.
    '''
    # When there is no runoff in this area, there is no transport.
    if (runoff < params.minimal_soil_infiltration):
        return 1.0

    # Check whether there is a load.
    if (load_org < 0.001):
        # No load and so no retention.
        return 0.0

    # Q zero order in m3/s,area in km2, so conversion of 1000.0
    Q0 = params.A0 *1000*runoff/year2second
    # Q for all streams
    Q = []
    for iorder in range(params.norder):
        # runoff in mm/yr, area in km2, so conversion of 1000.0
        Q.append(params.riverarea[iorder] * 1000. * runoff/year2second)

    # Determine the flux at the midpoint of the stream section [m3/s].
    Qmid = []
    Qmid.append(Q0 + 0.5 * Q[0])
    for iorder in range(1,params.norder):
        Qmid.append(Q[iorder-1] + 0.5 * Q[iorder])

    # Stream width for each order [m]
    width = []
    for iorder in range(params.norder):
        width.append(params.width_factor * (Qmid[iorder]**params.width_exponent))

    # Calculate the hydraulic load [m/yr] for each order
    Hl = []
    for iorder in range(params.norder):
        # width in m, riverlength in km, so conversion of 1000.0
        Hl.append(year2second * Qmid[iorder]/(1000.0 * width[iorder]*params.riverlength[iorder]))

    # Calculate the retention in all the stream orders. Qmid must be transformed into km3/yr.
    m3_s_to_km3_y = year2second * 1.0e-9
    load_out = []
    load_in  = []
    retention = []
    for iorder in range(params.norder):
        load = params.load_distribution[iorder] * load_org
        for i in range(len(load_out)):
            load += params.flowpath_distribution[i][iorder] * load_out[i]
        load_in.append(load)
        # Calculate the retention
        ret = accuflux_retention.calculate_retention(params,vf,Hl[iorder],1,\
                        temperature,m3_s_to_km3_y * Qmid[iorder],load)
        retention.append(ret)
        load_out.append(load * (1.0 - retention[iorder]))

    #print "Q",Q
    #print "Qmid",Qmid
    #print "Width",width
    #print "Hl",Hl
    #print "retention",retention
    #for iorder in range(params.norder):
    #    print "load_in,load_out ",iorder," ",load_in[iorder],load_out[iorder]
    #print (load_org - load_in[-1])/load_org 
   
    # The final retention of the norder streams is the original load minus the in load of the last stream (which is modelled in PCRGLOBWB) 
    # devided by the original load 
    return (load_org - load_in[-1])/load_org

def calculate_subgrid_retention(params,mask,mouth_dict,basin,lake_icell_dict,load,\
                                water_body,pnet,vf):
    '''
    The subgrid retention is based on Wollheim (2006, GBC). This method is based on the stream order.
    In each cell, which is not a lake or reservoir, the subgrid retention is calculated.
    Only the local situation is used. For the load only the diffuse source from this cell are used. 
    The runoff is used (so only local water). Also the concentration effect or temperature effect on retention is not used here.
    '''
    # Make the subgrid stream order properties.
    # Setup the stream order river length,area and distribution and flowpath. Everything is stored in params.
    define_subgrid_streamorder(params)

    # Set the relation between concentration and net uptake velocity
    try:
        accuflux_retention.set_conc_relation_retention(params)
    except AttributeError:
        # This function is only for N and not for P
        pass

    # Make output rasters
    load_grid       = ascraster.duplicategrid(load)
    retention      = ascraster.duplicategrid(load)

    # Make nodata on all cells.
    retention.add_values(retention.length*[retention.nodata_value])
    retentionload  = ascraster.duplicategrid(retention)

    # Read water temperature
    temperature = ascraster.Asciigrid(ascii_file=params.water_temperature,mask=mask,numtype=float)

    # Loop over all cells:    
    for icell in range(load_grid.length):
        # Get mass flux at the cell.
        load_cell = load_grid.get_data(icell)
        if (load_cell == None):
            # No data value
            continue
        if load_cell < 0.0:
            print("Negative mass flow of " + str(load_cell) + " at cell: " + str(icell))           
            # Go to next cell.
            continue

        # Retention in cell
        water_body_cell = water_body.get_data(icell,0)
        if (water_body_cell > 0):
            runoff_cell = pnet.get_data(icell,0.0)
            # Check whether this cell is an outflow of a lake/reservoir
            try:            
                loutflow = lake_icell_dict[icell].outcell
            except KeyError:
                # No lake or reservoir
                loutflow = -1
            try:
                if (loutflow == -1):
                    # Normal cell, no lake or reservoir.
                    # Default temperature of 20C because then temperature has no effect.
                    temperature_cell = temperature.get_data(icell,20.0)
                    ret_cell = calculate_localretention(params,runoff_cell,vf,load_cell,temperature_cell)

                else:
                    #TODO: local retention when it is a lake or reservoir. For example as normal retention * fr_land.
                    # Lake or reservoir.
                    ret_cell = 0.0

            except OverflowError:
                raise OverflowError
        else:
            ret_cell = 0.0

        # Fill retention grid with retention value and retention load
        retention.set_data(icell,ret_cell)
        retentionload.set_data(icell,load_cell*ret_cell)

    # Write to output files
    retention.write_ascii_file(os.path.join(params.outputdir,"retention_subgrid.asc"))
    retentionload.write_ascii_file(os.path.join(params.outputdir,"removed_subgridload.asc"))

    # Write in mouth_dict database
    aggregate.aggregate_grid(basin,retentionload,mouth_dict,item="retentionload_subgrid")
    aggregate.aggregate_grid(basin,load_grid,mouth_dict,item="subgrid_load")
    # Put the retention_subgrid of the river in the mouth database
    for key in list(mouth_dict.keys()):
        try:
            mouth_dict[key].retention_subgrid = mouth_dict[key].retentionload_subgrid/mouth_dict[key].subgrid_load
        except ZeroDivisionError:
            mouth_dict[key].retention_subgrid = 0.0

    return retentionload,retention
