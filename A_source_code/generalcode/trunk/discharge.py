# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/discharge.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import os
import accuflux
import ascraster
                 
def calculate_discharge(params,mask,pnet):
    '''
    Calculate the accumulated water discharge 
    '''
    ldd = ascraster.Asciigrid(ascii_file=params.ldd,mask=mask,numtype=int)
    if params.ldebug: print(params.ldd  + " has been read.")  
    
    landarea = ascraster.Asciigrid(ascii_file=params.landarea,mask=mask,numtype=float)
    if params.ldebug: print(params.landarea  + " has been read.")    

    water = ascraster.duplicategrid(landarea)    
    # Calculate accumulated water flux of each grid cell to the mouth of the river basin. Make discharge in km3
    # Convert pnet [mm/yr] to km3 by multiplying with landarea and 10e-6 (conversion from mm to km)
    
    for icell in range(landarea.length):
        wt = pnet.get_data(icell,0.0) * landarea.get_data(icell,0.0) * 1.0e-6
        water.set_data(icell,wt)
    
    discharge = accuflux.accuflux(ldd,water,negative_flux_possible=1)
    
    # Write discharge to output file:
    discharge.write_ascii_file(os.path.join(params.outputdir,"discharge.asc")) 
    
    return discharge    
