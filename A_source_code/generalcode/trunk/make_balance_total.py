# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/make_balance_total.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import ascraster

    
def make_balance_total(params,N_bal_grass,N_bal_arable,N_bal_natural,\
                       grassarea,croparea,natarea,mouth_dict,basin,check=1):
    '''
    Adds all positive balances to one total balance.
    Balances are changed when balance term is smaller than params.epsilon_load.
    When area does not match with the balances, than the balances are made zero.
    Input argument check determines whether the negative mass is counted for. 
    '''
    errorlevel = 1.0e-15
       
    if params.ldebug: print("Start with calculation of make_balance_total.")

    balance_total = ascraster.duplicategrid(N_bal_grass)
    if (check):
        N_bal_grass,N_bal_arable,N_bal_natural,removed_balance_arable,removed_balance_grs,removed_balance_nat = \
              clean_balance(params,N_bal_grass,N_bal_arable,N_bal_natural,\
                            grassarea,croparea,natarea,mouth_dict,basin)
    else:
        removed_balance_arable = ascraster.duplicategrid(N_bal_grass)
        removed_balance_arable.add_values(N_bal_grass.length * [0.0])
        removed_balance_grs = ascraster.duplicategrid(N_bal_grass)
        removed_balance_grs.add_values(N_bal_grass.length * [0.0])
        removed_balance_nat = ascraster.duplicategrid(N_bal_grass)
        removed_balance_nat.add_values(N_bal_grass.length * [0.0])
    
    for icell in range(N_bal_grass.length):
        N_bal_grs_cell = N_bal_grass.get_data(icell,0.0)
        N_bal_arable_cell = N_bal_arable.get_data(icell,0.0)
        N_bal_nat_cell = N_bal_natural.get_data(icell,0.0)
        total_load = N_bal_grs_cell +  N_bal_arable_cell + N_bal_nat_cell      
        balance_total.set_data(icell,total_load)

    return balance_total,removed_balance_arable,removed_balance_grs,removed_balance_nat


def clean_balance(params,N_bal_grass,N_bal_arable,N_bal_natural,\
                  grassarea,croparea,natarea,mouth_dict,basin):
    '''
    Adds all positive balances to one total balance.
    Balances are changed when balance term is smaller than params.epsilon_load.
    When area does not match with the balances, than the balances are made zero.
    Input argument check determines whether the negative mass is counted for. 
    '''
    errorlevel = 10
       
    if params.ldebug: print("Start with calculation of clean_balance.")

    removed_balance_arable = ascraster.duplicategrid(N_bal_grass)
    removed_balance_arable.add_values(N_bal_grass.length * [0.0])
    removed_balance_grs = ascraster.duplicategrid(N_bal_grass)
    removed_balance_grs.add_values(N_bal_grass.length * [0.0])
    removed_balance_nat = ascraster.duplicategrid(N_bal_grass)
    removed_balance_nat.add_values(N_bal_grass.length * [0.0])
    
    for icell in range(N_bal_grass.length):
        grass_cell = grassarea.get_data(icell,0.0)
        N_bal_cell = N_bal_grass.get_data(icell,0.0)
        if (grass_cell > 0.0):
            if (N_bal_cell < params.epsilon_load):
                # Make balance really zero
                removed_balance_grs.set_data(icell,N_bal_cell)
                N_bal_grass.set_data(icell,0.0)
        elif (N_bal_cell > 0.0):
            if (N_bal_cell > errorlevel):
                print("Grass area and grasbalance do not belong to each other. Gras area = 0.0 and balance_grass = " + str(N_bal_cell))
            N_bal_grass.set_data(icell,0.0)
            removed_balance_grs.set_data(icell,N_bal_cell)
        elif (N_bal_cell < 0.0):
            if (N_bal_cell < -errorlevel):
                print("Grass area and grasbalance do not belong to each other. Gras area = 0.0 and balance_grass = " + str(N_bal_cell))
            N_bal_grass.set_data(icell,0.0)
            removed_balance_grs.set_data(icell,N_bal_cell)
        
        crop_cell = croparea.get_data(icell,0.0)
        N_bal_cell = N_bal_arable.get_data(icell,0.0)        
        if (crop_cell > 0.0):
            if (N_bal_cell < params.epsilon_load):
                # Make balance really zero
                N_bal_arable.set_data(icell,0.0)
                removed_balance_arable.set_data(icell,N_bal_cell)
        elif (N_bal_cell > 0.0):
            if (N_bal_cell > errorlevel):
                print("Crop area and cropbalance do not belong to each other. Crop area = 0.0 and balance_crop = " + str(N_bal_cell))
            N_bal_arable.set_data(icell,0.0)
            removed_balance_arable.set_data(icell,N_bal_cell)

        elif (N_bal_cell < 0.0):
            if (N_bal_cell < -errorlevel):
                print("Crop area and cropbalance do not belong to each other. Crop area = 0.0 and balance_crop = " + str(N_bal_cell))
            N_bal_arable.set_data(icell,0.0)
            removed_balance_arable.set_data(icell,N_bal_cell)

        nat_cell = natarea.get_data(icell,0.0)
        N_bal_cell = N_bal_natural.get_data(icell,0.0)
        if (nat_cell > 0.0):
            if (N_bal_cell < params.epsilon_load):
                # Make balance really zero
                N_bal_natural.set_data(icell,0.0)
                removed_balance_nat.set_data(icell,N_bal_cell)
        elif (N_bal_cell > 0.0):
            if (N_bal_cell > errorlevel):
                print("Natural area and natural balance do not belong to each other. Natural area = 0.0 and balance_nat = " + str(N_bal_cell))
            N_bal_natural.set_data(icell,0.0)
            removed_balance_nat.set_data(icell,N_bal_cell)
        elif (N_bal_cell < 0.0):
            if (N_bal_cell < -errorlevel):
                print("Natural area and natural balance do not belong to each other. Natural area = 0.0 and balance_nat = " + str(N_bal_cell))
            N_bal_natural.set_data(icell,0.0)
            removed_balance_nat.set_data(icell,N_bal_cell)

    return N_bal_grass,N_bal_arable,N_bal_natural,removed_balance_arable,removed_balance_grs,removed_balance_nat
    
