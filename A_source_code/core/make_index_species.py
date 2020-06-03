# ******************************************************
## Revision "$LastChangedDate: 2018-07-08 18:08:17 +0200 (zo, 08 jul 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/make_index_species.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import local modules.
import general_func

def make_index_species(params,species,proc): 
   # Make an indexnumber in the species list with the name of the specie.
    # For example:
    #    params.io2   = general_func.find_spec(species,"O2")
    for i in range(len(species)):
        name = species[i].get_name().lower()
        setattr(params,"i"+name,i)

    # Create dictionary for all the available processes in proc.
    params.proc_dict = general_func.find_func_dict(proc)
    
    # Make an indexnumber in the proc list with the name of the process.
    params.ipsed   = general_func.get_val_dict(params.proc_dict,"sedimentation")
    params.ipero   = general_func.get_val_dict(params.proc_dict,"erosion")
    params.ipdiss  = general_func.get_val_dict(params.proc_dict,"pop_diss")
    params.ippup   = general_func.get_val_dict(params.proc_dict,"phot_p_uptake")
    params.ipsorp  = general_func.get_val_dict(params.proc_dict,"sorption")
    params.ipsed_pop   = general_func.get_val_dict(params.proc_dict,"sedim_pop")
    params.ipsed_pip   = general_func.get_val_dict(params.proc_dict,"sedim_pip")
    params.ipsed_ondet   = general_func.get_val_dict(params.proc_dict,"sedim_ondet")
    params.ipsed_onpp   = general_func.get_val_dict(params.proc_dict,"sedim_onpp")
    params.ipsed_opdet   = general_func.get_val_dict(params.proc_dict,"sedim_opdet")
    params.ipsed_oppp   = general_func.get_val_dict(params.proc_dict,"sedim_oppp")
    params.ipero_pop   = general_func.get_val_dict(params.proc_dict,"erosion_pop")
    params.ipero_pip   = general_func.get_val_dict(params.proc_dict,"erosion_pip")
    params.ipero_ondet   = general_func.get_val_dict(params.proc_dict,"erosion_ondet")
    params.ipero_opdet   = general_func.get_val_dict(params.proc_dict,"erosion_opdet")
    params.ipbdiss  = general_func.get_val_dict(params.proc_dict,"benth_pop_diss")
    params.ipbsorp  = general_func.get_val_dict(params.proc_dict,"benth_sorption")
    params.ipbsorp2  = general_func.get_val_dict(params.proc_dict,"benth_sorption2")
    params.ipbsorp3  = general_func.get_val_dict(params.proc_dict,"benth_sorption3")
    params.ipbdesorp3  = general_func.get_val_dict(params.proc_dict,"benth_desorption3")
    params.ipnpp   = general_func.get_val_dict(params.proc_dict,"phot_n_pp")
    params.ipppp   = general_func.get_val_dict(params.proc_dict,"phot_p_pp")
    params.ipnresp   = general_func.get_val_dict(params.proc_dict,"phot_n_resp")
    params.ippresp   = general_func.get_val_dict(params.proc_dict,"phot_p_resp")
    params.ipnmort   = general_func.get_val_dict(params.proc_dict,"phot_n_mort")
    params.ippmort   = general_func.get_val_dict(params.proc_dict,"phot_p_mort")
    params.ipnh4up   = general_func.get_val_dict(params.proc_dict,"phot_nh4_uptake")
    params.ipno3up   = general_func.get_val_dict(params.proc_dict,"phot_no3_uptake")
    params.ipdissox   = general_func.get_val_dict(params.proc_dict,"pon_diss_oxic")
    params.ipdissanox   = general_func.get_val_dict(params.proc_dict,"pon_diss_anoxic")
    params.ipnit   = general_func.get_val_dict(params.proc_dict,"nitrification")
    params.ipbnit   = general_func.get_val_dict(params.proc_dict,"benth_nitrification")
    params.ipbondiss   = general_func.get_val_dict(params.proc_dict,"benth_ondiss")
    params.ipbopdiss   = general_func.get_val_dict(params.proc_dict,"benth_opdiss")
    params.ipbammonr   = general_func.get_val_dict(params.proc_dict,"benth_ammonr")
    params.ipbpo4r   = general_func.get_val_dict(params.proc_dict,"benth_po4r")
    params.ipbresp   = general_func.get_val_dict(params.proc_dict,"benth_resp")
    params.ipbdenit   = general_func.get_val_dict(params.proc_dict,"benth_denit")
    params.ipcomp_pip   = general_func.get_val_dict(params.proc_dict,"compaction_pip")
    params.ipcomp_tss   = general_func.get_val_dict(params.proc_dict,"compaction_tss")

    # Make a list with the indexes in the species list.
    params.totalphotindex = []

    for i in range(len(species)):  
        #Create the index for all species which require photosynthesis
        if (species[i].get_type() == "SPEC" and species[i].get_kmax() > 0. and species[i].get_alpha() > 0.):
            params.totalphotindex.append(i)

