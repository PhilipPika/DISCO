# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/create_outfiles.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# General python modules
import os

# Import own general modules
import general_func

def create_outfiles(params,species,proc,year,lrestart=False):
    '''
    Create the headers of the output pkl files 
    '''
    # In order to avoid too long names in the case of monthly simulations...
    strtime = general_func.stryear(year) 

    if (params.allorders == 1):
        general_func.make_header(os.path.join(params.outputdir,strtime+"_allorders.pkl"),species,year)
        general_func.make_header(os.path.join(params.outputdir,"conc_"+strtime+"_allorders.pkl"),species,year)
    else:
        general_func.make_header(os.path.join(params.outputdir,strtime+".pkl"),species,year)
        general_func.make_header(os.path.join(params.outputdir,"conc_"+strtime+".pkl"),species,year)
        
    general_func.make_header(os.path.join(params.outputdir,"start"+strtime+".pkl"),species,year)

    if (params.lbudget == 1):
        general_func.make_header_budget(os.path.join(params.outputdir,"budget_"+strtime+".pkl"),species,proc,year)
    
    if (params.largs == 1):
        general_func.make_header_arguments(os.path.join(params.outputdir,"arguments_"+strtime+".pkl"),year)

    if (lrestart):
        general_func.make_header(os.path.join(params.outputdir,"start"+strtime+".pkl"),species,year)
