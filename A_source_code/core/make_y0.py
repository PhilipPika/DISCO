# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/make_y0.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import local modules.
import general_func

def make_y0(params,species):

    y0 = general_func.get_amount(species)

    # Append spool and rpool 
    #for i in params.phytoindex:
        #y0.append(species[i].get_spool())
        #y0.append(species[i].get_rpool())
    
    #for i in params.pocindex:
        #y0.append(species[i].get_dissolved())  
  
    return y0
