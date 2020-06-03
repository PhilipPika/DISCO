# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
'''
Calculation of biogeochemical process rates
All species are represented with specie class
'''

import math
import mocsy_module
import mocsy
import numpy as np
import operator
import os
import sys
import time
from contextlib import contextmanager


# Import local modules.
import general_func

# Generalcode module
from error import *

class suppress_stdout_stderr(object):
    def __enter__(self):
        self.outnull_file = open(os.devnull, 'w')
        self.errnull_file = open(os.devnull, 'w')

        self.old_stdout_fileno_undup    = sys.stdout.fileno()
        self.old_stderr_fileno_undup    = sys.stderr.fileno()

        self.old_stdout_fileno = os.dup ( sys.stdout.fileno() )
        self.old_stderr_fileno = os.dup ( sys.stderr.fileno() )

        self.old_stdout = sys.stdout
        self.old_stderr = sys.stderr

        os.dup2 ( self.outnull_file.fileno(), self.old_stdout_fileno_undup )
        os.dup2 ( self.errnull_file.fileno(), self.old_stderr_fileno_undup )

        sys.stdout = self.outnull_file        
        sys.stderr = self.errnull_file
        return self

    def __exit__(self, *_):        
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr

        os.dup2 ( self.old_stdout_fileno, self.old_stdout_fileno_undup )
        os.dup2 ( self.old_stderr_fileno, self.old_stderr_fileno_undup )

        os.close ( self.old_stdout_fileno )
        os.close ( self.old_stderr_fileno )

        self.outnull_file.close()

##############################################################################
######################################

def specie_dy(proc,speciesname):
    '''
    Makes a list of the factors for each species.
    '''
    # Set all factors to zero.
    vals = len(proc)*[0.0]
    
    # Find all the processess which influence this species.
    for item in range(len(proc)):
        try:
            vals[item] = proc[item].get_val("dy_"+speciesname)
            
        except MyError:
            pass
    return vals

##############################################################################
######################################
##### GENERAL REACTION FUNCTIONS #####    
    
def gen_prim_prod(sp,amount,lim):
  if amount>0:
    np.seterr(all='raise')
    try:  
      return lim*sp.get_val('max_growth_r')*365*amount
    except(FloatingPointError):
      return 0.
  else:
    return 0.
   
def gen_respiration(params,sp,amount,temperature, vol):
  if amount>0 and vol>0:
    np.seterr(all='raise')
    try:
      fT = general_func.fT(sp.get_val("Topt"), sp.get_val("sigma"), temperature-params.tempcorrection)
      return fT*sp.get_val('max_resp_r')*365*amount
    except(FloatingPointError):
      return 0.
  else:
    return 0

def gen_mineralization(params, sp, amount, temperature):
  if amount>0:
    np.seterr(all='raise')
    try:
      if sp.get_val('q10')==None:
        return sp.get_val('max_oxi_r')*365*amount
      else:
        return amount*sp.get_val('min_oxi_r')*math.pow(sp.get_val('q10'), (temperature-params.tempcorrection-15)/10)*365
    except(FloatingPointError):
      return 0.
  else:
    return 0

def gen_mortality(params,sp,amount,temperature,vol, area):
  if amount > 0 and vol>0:
    np.seterr(all='raise')
    fT = general_func.fT(sp.get_val("Topt"), sp.get_val("sigma"), temperature-params.tempcorrection)
    mort = fT* sp.get_val('max_mort_r')
    try:
        phyconc = (amount/vol)*sp.get_val('molarmass')*1e-3 #mg/l
        factor = general_func.MM(phyconc, sp.get_val("paras_thres")*sp.get_val("spec_chl_ratio"))
    except(OverflowError):
        factor = 1
    except(FloatingPointError):
         factor = 0
    try:
      mort += mort + factor*mort*sp.get_val("vf")
    except(FloatingPointError):
      return 0.
        
    try:
        return mort*365*amount
    except(FloatingPointError):
        return 0.
  else:
    return 0


def gen_doc_excretion(params,sp,amount,temperature, vol):
  if amount>0 and vol>0.:
    np.seterr(all='raise')
    try:	
      fT = general_func.fT(sp.get_val("Topt"), sp.get_val("sigma"), temperature-params.tempcorrection)
      return fT*sp.get_val('max_excr_r')*365*amount
    except(FloatingPointError):
      return 0.
  else:
    return 0.

######################################
######################################
##############################################################################

##############################################################################
######################################
##### SPECIE SPECIFIC FUNCTIONS ######

# primary production of algae limited by C, light and temperature
def prim_prod_C_PHYTO(spec,params,species,vol,depth,temperature, glob_rad):
    # DIC effect 
    try:
      conc_DIC = spec[params.idic]/vol #mmol/m3 = um/L
      try:
        factor = 1/(1+math.exp(-10*conc_DIC))
      except(OverflowError):
        if conc_DIC>0:
          factor = 1
        else:
          factor = 0
      DIC_lim = general_func.MM(conc_DIC, species[params.iphyto].get_kDIC())*factor
    except(FloatingPointError):
      DIC_lim = 0


    # Turbidity
    eta_min = 0.8
    eta_max = 20
    eta_dum = 0.8
    for s in species:
        ispec = getattr(params, "i"+s.get_name().lower())
        if spec[ispec]>0:
          try:
            eta_dum += s.get_val('eta')*((spec[ispec])/(vol))*s.get_val('molarmass') #ug/l
          except(FloatingPointError):
            if spec[ispec]>1:
              eta_dum += 20
            else:
              eta_dum += 0
    eta_dum2 = max(eta_min, eta_dum) 
    eta = min(eta_max, eta_dum2)

    # Light attenuation in the water column (2 layers)
    I_lim = 0
    for d in np.linspace(0, depth, 2):
        try:
          I_lim += general_func.MM(glob_rad*math.exp(-d*eta), species[params.iphyto].get_kI())*0.5
        except(FloatingPointError):
          I_lim += 0
     
    # Temperature effect
    fT = general_func.fT(species[params.iphyto].get_val("Topt"), species[params.iphyto].get_val("sigma"), temperature-params.tempcorrection)

    # limitation of primary production
    try:
      lim = I_lim*DIC_lim*fT
    except(FloatingPointError):
      lim = 0

    try:
      if vol>0.:
        if spec[params.iphyto]>0:
          vol_density = spec[params.iphyto]/vol*species[params.iphyto].get_val('molarmass') #ug/l
          amount_lim = 1-general_func.MM(vol_density*1e-3, species[params.iphyto].get_val("paras_thres")*species[params.iphyto].get_val('spec_chl_ratio'))
          if amount_lim<0:
            amount_lim = 0
          elif amount_lim>1:
            amount_lim = 1
        else:
          amount_lim=1
      else:
        return 0.
    except(OverflowError):
      amount_lim = 0
    except(FloatingPointError):
      amount_lim = 1

    return gen_prim_prod(species[params.iphyto],spec[params.iphyto],lim)*amount_lim


# primary production of periphyton limited by C, light and temperature
def prim_prod_C_PERIPHYTON(spec,params,species,vol,depth,temperature, glob_rad, area):
    # DIC effect 
    try:
      conc_DIC = spec[params.idic]/vol #mmol/m3 = um/L
      try:
        factor = 1/(1+math.exp(-10*conc_DIC))
      except(OverflowError):
        if conc_DIC>0:
          factor = 1
        else:
          factor = 0

      DIC_lim = general_func.MM(conc_DIC, species[params.iperiphyton_benth].get_kDIC())*factor
    except(FloatingPointError):
      DIC_lim = 0

    # Turbidity
    eta_min = 0.8
    eta_max = 20
    eta_dum = 0.8
    for s in species:
        ispec = getattr(params, "i"+s.get_name().lower())
        if spec[ispec]>0:
          try:
            eta_dum += s.get_val('eta')*((spec[ispec])/(vol))*s.get_val('molarmass') #ug/l
          except(FloatingPointError):
            if spec[ispec]>1:
              eta_dum += 20
            else:
              eta_dum += 0
    eta_dum2 = max(eta_min, eta_dum) 
    eta = min(eta_max, eta_dum2)

    # Light attenuation at water bottom
    try:
      light = glob_rad*math.exp(-depth*eta)
      I_lim = general_func.MM(light, species[params.iperiphyton_benth].get_kI())
    except(FloatingPointError):
        I_lim = 0

    # temperature effect
    fT = general_func.fT(species[params.iperiphyton_benth].get_val("Topt"), species[params.iperiphyton_benth].get_val("sigma"), temperature-params.tempcorrection)
    
    # limitation of primary production
    try:
      lim = I_lim*DIC_lim*fT
    except(FloatingPointError):
      lim = 0
    
    try:
      if vol>0.:
        
        if spec[params.iperiphyton_benth]>0: 
          peri_vol_density = spec[params.iperiphyton_benth]/vol*species[params.iperiphyton_benth].get_val('molarmass') #ug/l
          amount_lim = 1-general_func.MM(peri_vol_density*1e-3, species[params.iperiphyton_benth].get_val("paras_thres")*species[params.iperiphyton_benth].get_val('spec_chl_ratio'))
          if amount_lim<0:
            amount_lim = 0
          elif amount_lim>1:
            amount_lim = 1
        else:
          amount_lim=1
      else:
        return 0.
    except(OverflowError):
      amount_lim = 0
    except(FloatingPointError):
      amount_lim = 1
    try:     
      return gen_prim_prod(species[params.iperiphyton_benth],spec[params.iperiphyton_benth],lim)*amount_lim
    except(FloatingPointError):
      return 0

# exchange of DIC with atmosphere    
def atmospheric_exchange_DIC(spec,params,species,temperature,windspeed,Q, flow_velocity, vol,width,length):
    # [mol/m3]
    try:    
      conc_DIC = max(spec[params.idic]*1e-3/vol, 0)
    except(FloatingPointError):
      conc_DIC = 0
    try:
      conc_ALK = max(spec[params.ialk]*1e-3/vol, 0)
    except(FloatingPointError):
      conc_ALK = 0

    if conc_DIC<999 and conc_ALK<999 and conc_DIC>0:
      carb_dat = np.dstack(mocsy.mvars(temp=temperature-params.tempcorrection, sal=0, alk=conc_ALK, dic=conc_DIC, sil=0, phos=0, patm=1., depth=1, lat=0, optcon='mol/m3', optt='Tinsitu', optp='db', optb="u74", optk1k2='l', optkf="dg", optgas='Pinsitu'))[0][0]
      # [micro atmosphere]    
      pCO2 = operator.itemgetter(1)(carb_dat)

    elif conc_DIC<=0:
      pCO2 = 0
    else:
      with suppress_stdout_stderr(): # suppressing printed output of mocsy module when concentrations of either ALK or DIC are more than 1M 
        carb_dat = np.dstack(mocsy.mvars(temp=temperature-params.tempcorrection, sal=0, alk=conc_ALK, dic=conc_DIC, sil=0, phos=0, patm=1., depth=1, lat=0, optcon='mol/m3', optt='Tinsitu', optp='db', optb="u74", optk1k2='l', optkf="dg", optgas='Pinsitu'))[0][0]
      # [micro atmosphere]    
      pCO2 = operator.itemgetter(1)(carb_dat)
      if pCO2>1e6:
        pCO2=1e6
   
    # [mol/m3]
    CO2_conc = pCO2*math.pow(10,-6)*species[params.idic].get_val("kH")*1000
    
    # [cm s-1]
    v = (flow_velocity*1e5 / (365.*86400.))
          
    # [km]
    length *= math.pow(10,-3)
    # [km]
    width *= math.pow(10,-3)

    if width< 100*math.pow(10,-3) and v>0.0:
        # [cm h-1]
        k600 = (species[params.idic].get_val("exch_a1") + species[params.idic].get_val("exch_b1")*v)
    else:
        # [cm h-1] 
        k600 = (species[params.idic].get_val("exch_a2") + species[params.idic].get_val("exch_b2")*windspeed)
      
    sc = 1911.1-118.11*(temperature-params.tempcorrection)+3.4527*(temperature-params.tempcorrection)**2-0.04132*(temperature-params.tempcorrection)**3
    try:
      k = k600/((600/sc)**(-0.5))
    except(FloatingPointError):
      k = 1e-9

    # [cm d-1]
    k*= 24.
    # [cm y-1]
    k*=365
    # [m y-1]
    k*=0.01


              
    #[mol/m3]
    dCO2 = CO2_conc-params.CO2_eq 

    if params.lsensitivity==1:
      #[mol/m3]
      dCO2 = CO2_conc-params.CO2_eq*params.co2_eq_factor
                 
    # [km y-1]
    k*=math.pow(10,-3)
    #[mol/km3]
    dCO2*=math.pow(10,9)
    #[Mmol/km3]
    dCO2*=math.pow(10,-6)
       
    #exchange in Mmol/year
    ex = dCO2*k*length*width

    return ex
   

def sedimentation(species, spec, ispec, proc, depth, vol, area):
    '''
    Returns sedimentation rate in Mmol.yr-1, tons.yr-1 or kg.yr-1
    '''   
    # Scheffer 2004, Ecology of shallow lakes
    try:
      return species[ispec].get_val("vsed") * spec[ispec] / (depth*1e-3) #(km/yr * tons / km = tons/yr)
    except(FloatingPointError):
      return 0.

def erosion(spec, ispec_benth, params, species, proc, width, depth, bedarea, vel, slope, area, vol, Q):
    '''
    Returns erosion rate in Mmol.yr-1 (tons.yr-1 for total sediments)
    '''

    np.seterr(all='raise')
    spec_benth_amount = spec[ispec_benth]
    if bedarea>0. and vol>0. and vel>0. and spec_benth_amount>0.:
      try:
        bsed = (spec[params.ipim_benth]+spec[params.idetritushighcn_benth]*24+spec[params.idetrituslowcn_benth]*24) / bedarea #sediment stock per bed surface area
      except(AttributeError):
        try:
          bsed = (spec[params.ipim_benth]+spec[params.idetritushighcn_benth]*24) / bedarea #sediment stock per bed surface area          
        except(AttributeError):
          bsed = spec[params.ipim_benth] / bedarea #sediment stock per bed surface area
    else:
      return 0.

    try: 
      factor = 1/(1+math.exp(-1*(bsed-species[params.itss_benth].get_val("k_sed"))))
    except(OverflowError):
      if bsed > 0:
        factor = 1
      else:
        factor = 0
      #print('overflow bsed =', bsed)

    ero_tss = factor*species[params.itss_benth].get_val("kero1")/(1e-3*6.66) * slope * vel * bedarea

    try:
      ero = ero_tss * spec_benth_amount/(spec[params.ipim_benth]+spec[params.idetritushighcn_benth]*24+spec[params.idetrituslowcn_benth]*24)
    except(AttributeError):
      try:
        ero = ero_tss * spec_benth_amount/(spec[params.ipim_benth]+spec[params.idetritushighcn_benth]*24)
      except(AttributeError):
        try:
          ero = ero_tss * spec_benth_amount/(spec[params.ipim_benth])
        except(FloatingPointError):
          return 0.
      except(FloatingPointError):
          return 0.
    except(FloatingPointError):
          return 0.

    return ero


def burial(spec, ispec, params, species, depth, vol):

    bur = 0.
    if spec[ispec]>0 and vol>0.:
        bedarea = vol * 1.e3/depth # km2

        try:
          bsed = (spec[params.ipim_benth]+spec[params.idetritushighcn_benth]*24+spec[params.idetrituslowcn_benth]*24) / bedarea #sediment stock per bed surface area
        except(AttributeError):
            try:
              bsed = (spec[params.ipim_benth]+spec[params.idetritushighcn_benth]*24) / bedarea #sediment stock per bed surface area
            except(AttributeError):
              bsed = (spec[params.ipim_benth]) / bedarea #sediment stock per bed surface area
        try:
            factor = 1/(1+math.exp(-1*(math.log(bsed)-math.log(species[params.itss_benth].get_val("bsedlim")))))
        except(OverflowError):
            if bsed > 0:
              factor = 1
            else:
              factor = 0
            #print('overflow bsed =', bsed)
        except(ValueError):
            if bsed > 0:
              factor = 1
            else:
              factor = 0
            #print('value bsed =', bsed)
        bur = species[params.itss_benth].get_val("compmax")
        bur *= factor
        try:
          bur *= spec[ispec]
        except(FloatingPointError):
          bur = 0.
    return bur


      
##### SPECIE SPECIFIC FUNCTIONS ######
######################################

#####
def procfunc(spec,params,species,proc,Q,arguments):
    '''
    Creates a list of processes
    '''
    out = []
    vol = arguments.get_val("vol")
    area = arguments.get_val("area")
    globrad_cell = arguments.get_val("globrad")
    temperature = arguments.get_val("temperature")
    windspeed = arguments.get_val("windspeed")
    vel = arguments.get_val("flow_velocity")
    width = arguments.get_val("width")
    depth = arguments.get_val("depth")
    length = arguments.get_val("length")
    slope =  arguments.get_val("slope")

    if vol>0.:
      for iproc in range(len(proc)):
        name = proc[iproc].get_val("name")
        if (name == "prim_prod_C_PHYTO"):
            out.append(prim_prod_C_PHYTO(spec,params,species,vol,depth,temperature, globrad_cell))
        elif (name == "prim_prod_C_PERIPHYTON"):
            out.append(prim_prod_C_PERIPHYTON(spec,params,species,vol,depth,temperature, globrad_cell, area))   

        elif (name == "respiration_PHYTO"):
            out.append(gen_respiration(params, species[params.iphyto],spec[params.iphyto],temperature, vol))
        elif (name == "respiration_PERIPHYTON"):
            out.append(gen_respiration(params,species[params.iperiphyton_benth],spec[params.iperiphyton_benth],temperature, vol))            

        elif (name == "mortality_PHYTO"):
            out.append(gen_mortality(params, species[params.iphyto],spec[params.iphyto],temperature,vol,area))
        elif (name == "mortality_PERIPHYTON"):
            out.append(gen_mortality(params,species[params.iperiphyton_benth],spec[params.iperiphyton_benth],temperature,vol,area)) 

        elif (name == "doc_excretion_PHYTO"):
            out.append(gen_doc_excretion(params,species[params.iphyto],spec[params.iphyto],temperature, vol))
        elif (name == "doc_excretion_PERIPHYTON"):
            out.append(gen_doc_excretion(params,species[params.iperiphyton_benth],spec[params.iperiphyton_benth],temperature, vol))  

        elif (name == "oxidation_DETRITUSlowCN"):
            out.append(gen_mineralization(params, species[params.idetrituslowcn], spec[params.idetrituslowcn], temperature))
        elif (name == "oxidation_DETRITUSlowCN_benth"):
            out.append(gen_mineralization(params, species[params.idetrituslowcn_benth], spec[params.idetrituslowcn_benth], temperature))
        elif (name == "oxidation_DETRITUShighCN"):
            out.append(gen_mineralization(params, species[params.idetritushighcn], spec[params.idetritushighcn], temperature))
        elif (name == "oxidation_DETRITUShighCN_benth"):
            out.append(gen_mineralization(params, species[params.idetritushighcn_benth], spec[params.idetritushighcn_benth], temperature))
        elif (name == "oxidation_DOC"):
            out.append(gen_mineralization(params, species[params.idoc], spec[params.idoc], temperature))
        elif (name == "sedimentation_TSS"):
            out.append(sedimentation(species, spec, params.itss, proc, depth, vol, area))
        elif (name == "sedimentation_PIM"):
            out.append(sedimentation(species, spec, params.ipim, proc, depth, vol, area))
        elif (name == "sedimentation_DETRITUShighCN"):
            out.append(sedimentation(species,spec,params.idetritushighcn,proc,depth, vol, area))
        elif (name == "sedimentation_DETRITUSlowCN"):
            out.append(sedimentation(species,spec,params.idetrituslowcn,proc,depth, vol, area))

        elif (name == "erosion_TSS"):
            out.append(erosion(spec, params.itss_benth, params, species, proc, width, depth, area, vel, slope, area, vol, Q))
        elif (name == "erosion_PIM"):
            out.append(erosion(spec, params.ipim_benth, params, species, proc, width, depth, area, vel, slope, area, vol, Q))
        elif (name == "erosion_DETRITUShighCN"):
            out.append(erosion(spec,params.idetritushighcn_benth,params,species,proc,width,depth, area, vel, slope, area, vol, Q))
        elif (name == "erosion_DETRITUSlowCN"):
            out.append(erosion(spec,params.idetrituslowcn_benth,params,species,proc,width,depth,area, vel, slope, area, vol, Q))

        elif (name == "burial_PIM"):
            out.append(burial(spec, params.ipim_benth, params, species, depth, vol))
        elif (name == "burial_DETRITUShighCN"):
            out.append(burial(spec, params.idetritushighcn_benth, params, species, depth, vol))
        elif (name == "burial_DETRITUSlowCN"):
            out.append(burial(spec, params.idetrituslowcn_benth, params, species, depth, vol))

        elif (name == "atmospheric_exchange_DIC"):
            dic_exch = atmospheric_exchange_DIC(spec,params,species,temperature,windspeed,Q,vel,vol,width,length)
            out.append(dic_exch)
        else:
            print(name)
            MyError.write("Function %s has not been defined." % name)
    else:
      return [0.]*len(proc)
     
    return out


# Total change rate due to biogeochemical processes 
def change_spec(spec, params, species, name, funcs, proc):
    '''
    Calculates the total change rate due to biogeochemical processes 
    for one species
    '''

    # Get the stoichiometry for the processes.
    vals = specie_dy(proc,name)

    # Calculate the product of the two lists and sum them.   
    total = 0.0

    for item in range(len(vals)):
        try:
          total += vals[item] * funcs[item]
        except(TypeError):
          total='error'

    return total 


def dy_list(spec,params,species,proc,load,Q,arguments,fp_arguments):
    '''
    This function returns a list of changing_rate function for each specie.
    '''
  
    # if these conditions meet it is a stream without floodplains
    if (params.lfloodplains == 0) or (fp_arguments==None) or len(spec)==len(species):
      out = []

      dvoldt = arguments.get_val("dvoldt")
      vol = arguments.get_val("vol")
 
      # Calculation of the change rates due to biogeochemical processing
      funcs = procfunc(spec,params,species,proc,Q,arguments)         
    
      # Calculation of total change rates of species in the water column
      for item in range(len(species)):
        try:      
          conc = max(0, spec[item] / vol)
        except(FloatingPointError):
          conc = 0.

        try:
          outflow = Q*conc
        except(FloatingPointError):
          outflow = 0.

        if not (species[item].get_name().endswith("_benth")):
            out.append(change_spec(spec, params, species, species[item].get_name(), funcs, proc) +\
                       load[item] - outflow)
        else:
              out.append(change_spec(spec, params, species, species[item].get_name(), funcs, proc))
      return out

    # if these conditions meet it is a stream with a floodplain
    elif (params.lfloodplains == 1):
      dvoldt = arguments.get_val("dvoldt")
      vol = arguments.get_val("vol")

      fp_dvoldt = fp_arguments.get_val("dvoldt")
      fp_vol = fp_arguments.get_val("vol")
      fp_depth = fp_arguments.get_val("depth")

      # Calculation of the change rates due to biogeochemical processing in mainstream
      funcs = procfunc(spec,params,species,proc,Q,arguments)


      if fp_vol>0.0:
        out = [0.]*len(spec)
        # Calculation of the change rates due to biogeochemical processing in floodplain
        fp_funcs = procfunc(spec[-len(species):],params,species,proc,0,fp_arguments)
        for item in range(len(species)):
          try:
            conc = max(0, spec[item] / vol)
          except(FloatingPointError):
            conc = 0.
          try:
            fp_conc = max(0, spec[item+len(species)] / fp_vol)
          except(FloatingPointError):
            fp_conc = 0.
    
          # all dissolved and suspended contents are calculated within this condition
          if not (species[item].get_name().endswith("_benth")):
            if fp_dvoldt>0. and fp_depth>0.:
              try:
                to_fp_dvoldt_amount = conc*fp_dvoldt
              except(FloatingPointError):
                to_fp_dvoldt_amount = 0.
              try:
                through_fp_volIN_amount = conc*fp_arguments.get_val('discharge')
                through_fp_volOUT_amount = fp_conc*fp_arguments.get_val('discharge')
                if through_fp_volIN_amount<0:
                  through_fp_volIN_amount = 0
                if through_fp_volOUT_amount<0:
                  through_fp_volOUT_amount= 0 
              except(FloatingPointError):
                through_fp_volIN_amount = 0.
                through_fp_volOUT_amount = 0.

              try:
                outflow = Q*conc
              except(FloatingPointError):
                outflow = 0.

              # in stream dY
              out[getattr(params, 'i'+species[item].get_name().lower())] = change_spec(spec[:len(species)], params, species, species[item].get_name(), funcs, proc)+load[item]\
                                                                           +through_fp_volOUT_amount-through_fp_volIN_amount\
                                                                           -outflow
              # floodplain dY
              out[getattr(params, 'i'+species[item].get_name().lower())+len(species)] = change_spec(spec[-len(species):], params, species, species[item].get_name(), fp_funcs, proc)+load[item+len(species)]\
                                                                           +through_fp_volIN_amount-through_fp_volOUT_amount

            else:
              try:
                from_fp_dvoldt_amount = fp_conc*fp_dvoldt
              except(FloatingPointError):
                from_fp_dvoldt_amount = 0.
              try:
                through_fp_volIN_amount = conc*fp_arguments.get_val('discharge')
                through_fp_volOUT_amount = fp_conc*fp_arguments.get_val('discharge')
                if through_fp_volIN_amount<0:
                  through_fp_volIN_amount = 0
                if through_fp_volOUT_amount<0:
                  through_fp_volOUT_amount= 0 
              except(FloatingPointError):
                through_fp_volIN_amount = 0.
                through_fp_volOUT_amount = 0.

              try:
                outflow = Q*conc
              except(FloatingPointError):
                outflow = 0.

              # in stream dY
              out[getattr(params, 'i'+species[item].get_name().lower())] = change_spec(spec[:len(species)], params, species, species[item].get_name(), funcs, proc) + load[item]\
                                                                           +through_fp_volOUT_amount-through_fp_volIN_amount\
                                                                           -outflow

              # floodplain dY
              out[getattr(params, 'i'+species[item].get_name().lower())+len(species)] = change_spec(spec[-len(species):], params, species, species[item].get_name(), fp_funcs, proc)+load[item+len(species)]\
                                                                                        +through_fp_volIN_amount-through_fp_volOUT_amount\

          # all benthic (surface attached) contents are calculated within this condition
          else:

            out[getattr(params, 'i'+species[item].get_name().lower())] = change_spec(spec[:len(species)], params, species, species[item].get_name(), funcs, proc)
            out[getattr(params, 'i'+species[item].get_name().lower())+len(species)] = change_spec(spec[-len(species):], params, species, species[item].get_name(), fp_funcs, proc)

        return out
      else:
        out = [0.]*len(spec)
        # Calculation of total change rates of species in the water column
        for item in range(len(species)):
          try:      
            conc = max(0., spec[item] / vol)
          except(FloatingPointError):
            conc = 0.

          try:
                outflow = Q*conc
          except(FloatingPointError):
                outflow = 0.
          if not (species[item].get_name().endswith("_benth")):
              out[getattr(params, 'i'+species[item].get_name().lower())] = change_spec(spec[:len(species)], params, species, species[item].get_name(), funcs, proc) + load[item] - outflow
          else:
              out[getattr(params, 'i'+species[item].get_name().lower())] = change_spec(spec[:len(species)], params, species, species[item].get_name(), funcs, proc)
       
        return out


def dy_steady(y,params,species,proc,load,Q,arguments,fp_arguments=None):
    '''
    This function is used in the solver, depending on y (list of species).
    It returns a list of changing_rate function for each specie.
    Process change rates are expressed in amounts (not in concentrations)
    '''
    return dy_list(y,params,species,proc,load,Q,arguments,fp_arguments)

    
def dy(y,t,params,species,proc,load,Q,arguments,fp_arguments=None):
    '''
    This function is used in the solver, depending on y (list of species) and t (time).
    It returns a list of changing_rate function for each species.
    Process change rates are expressed in amounts (not in concentrations)
    '''
    return dy_list(y,params,species,proc,load,Q,arguments, fp_arguments)


