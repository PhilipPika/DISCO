# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import general_class

'''
Class with general species specifications
'''
class Spec(general_class.General):
    def __init__(self, name=None, type="SPEC", unit=None, amount=0.0, lake_amount=0.0,\
                     molarmass=None, conv2dm=None, k_I=0.0, k_DIC=0.0, k_DIN=0.0, k_DIP=0.0, k_LOC=0.0, k_LOCbenth=0.0, k_space=0.0, k_DSi=0.0, k_sed=0.0, k_mort=0.0,\
                     Topt=0.0, sigma=0.0, benthic_uptake=0.0, \
                     Tref_growth=0.0, Q10_growth=0.0, Tref_resp=0.0, Q10_resp=0.0, Tref_mort=0.0, Q10_mort=0.0,\
                     max_growth_r=0.0, max_grazing_r=0.0, max_resp_r=0.0, max_excr_r=0.0, max_oxi_r=0.0, min_oxi_r=0.0, max_mort_r=0.0, min_mort_r=0.0, vf=0.0, paras_thres=0.0, spec_chl_ratio=0.0,\
                     sed_r=0.0, diss_r=0.0, bur_r=0.0, dens=0.0, poro=0.0,compmax=0.0,bsedlim=0.0,\
                     CNratio=0.0, CPratio=0.0, kmax=0.0, alpha=0.0, a_min=0, b_min=0, q10=None,\
                     basedlim=0.0, vsed=0.0, ldistrib=1, exch_rate=0.0, exch_a1=0.0, exch_b1=0.0, exch_a2=0.0, exch_b2=0.0, kH=0.0, eta=0.0, kero1=0.0, kdiss=0.0, \
                     atmospheric_conc=0, a=0, b=0, phylim=0.0, budget=[],
                     color=None):
        general_class.General.__init__(self)             
        
        if lake_amount==0.0:
            lake_amount=amount
        
        self.add_item('name', name)
        self.add_item('unit', unit)
        self.add_item('type',type)
        self.add_item('amount',amount)
        self.add_item('lake_amount', lake_amount)
        self.add_item('molarmass',molarmass)
        self.add_item('conv2dm',conv2dm)
        self.add_item('k_I',k_I)
        self.add_item('k_DIC', k_DIC)
        self.add_item('k_DIN',k_DIN)
        self.add_item('k_DIP',k_DIP)
        self.add_item('k_LOC', k_LOC)
        self.add_item('k_LOCbenth', k_LOCbenth)
        self.add_item('k_space', k_space)
        self.add_item('k_DSi', k_DSi)
        self.add_item('k_sed', k_sed)
        self.add_item('k_mort', k_mort)
        self.add_item('Topt',Topt)
        self.add_item('sigma',sigma)
        self.add_item('Tref_growth',Tref_growth)
        self.add_item('Q10_growth',Q10_growth)
        self.add_item('Tref_resp',Tref_resp)
        self.add_item('Q10_resp',Q10_resp)
        self.add_item('Tref_mort',Tref_mort)
        self.add_item('Q10_mort',Q10_mort)        
        self.add_item('max_growth_r',max_growth_r)
        self.add_item('max_grazing_r',max_grazing_r)
        self.add_item('max_resp_r',max_resp_r)
        self.add_item('max_excr_r',max_excr_r)
        self.add_item('max_oxi_r',max_oxi_r)
        self.add_item('min_oxi_r',min_oxi_r)
        self.add_item('max_mort_r',max_mort_r)
        self.add_item('min_mort_r', min_mort_r)
        self.add_item('paras_thres', paras_thres)
        self.add_item('vf', vf)
        self.add_item('spec_chl_ratio', spec_chl_ratio)
        self.add_item('benthic_uptake', benthic_uptake)
        self.add_item('CNratio',CNratio)
        self.add_item('CPratio',CPratio)
        self.add_item('exch_a1',exch_a1)
        self.add_item('exch_b1',exch_b1)
        self.add_item('exch_a2',exch_a2)
        self.add_item('exch_b2',exch_b2)  
        self.add_item('sed_r', sed_r)
        self.add_item('diss_r', diss_r)
        self.add_item('bur_r', bur_r)   
        self.add_item('dens', dens)   
        self.add_item('poro', poro)
        self.add_item('compmax', compmax)
        self.add_item('bsedlim', bsedlim)
        self.add_item('kH',kH)
        self.add_item('eta',eta)
        self.add_item('kmax',kmax)
        self.add_item('alpha',alpha)
        self.add_item('a_min', a_min)
        self.add_item('b_min', b_min)
        self.add_item('q10', q10)
        self.add_item('vsed',vsed)
        self.add_item('ldistrib',ldistrib)
        self.add_item('exch_rate',exch_rate)
        self.add_item('kero1',kero1)
        self.add_item('kdiss', kdiss)
        self.add_item('basedlim', basedlim)
        self.add_item('atmospheric_conc',atmospheric_conc) #atmospheric_concentration[Mmol/km3]
        self.add_item('phylim', phylim)
        self.add_item('a', a)
        self.add_item('b', b)
        self.add_item('budget',budget)
        self.add_item('color', color)

    def check_type(self):
        if (self.name         != None):    
            self.name          = str(self.name)
        if (self.amount       != None): 
            self.amount        = float(self.amount)
        if (self.molarmass    != None): 
            self.molarmass     = float(self.molarmass)
        if (self.conv2dm      != None): 
            self.conv2dm       = float(self.conv2dm)
        if (self.k_I          != None):
            self.k_I           = float(self.k_I)
        if (self.kmax         != None): 
            self.kmax          = float(self.kmax)
        if (self.alpha        != None): 
            self.alpha         = float(self.alpha)
        if (self.vsed         != None):
            self.vsed          = float(self.vsed)
        if (self.ldistrib     != None):
            self.ldistrib      = int(self.ldistrib)

    def get_name(self):
        return str(self.name)

    def get_type(self):
        return self.type

    def get_amount(self):
        return float(self.amount)

    def get_molarmass(self):
        return float(self.molarmass)

    def get_conv2dm(self):
        return float(self.conv2dm)

    def get_kmax(self):
        return float(self.kmax)

    def get_alpha(self):
        return float(self.alpha)
        
    def get_kI(self):
        return self.k_I
    
    def get_kDIC(self):
        return self.k_DIC

    def get_kDIN(self):
        return self.k_DIN

    def get_kDIP(self):
        return self.k_DIP
        
    def get_exch_a(self):
        return self.exch_a
        
    def get_exch_b(self):
        return self.exch_b        

    def get_kH(self):
        return self.kH

    def get_vsed(self):
        return float(self.vsed)
    
    def get_ldistrib(self):
        return int(self.ldistrib)

    def get_exch_rate(self):
        return float(self.exch_rate) 

    def set_exch_rate(self, val):
        self.exch_rate = val

    def add_exch_rate(self, val):
        self.exch_rate += val

    def reset_budget(self, nproc):
        #changes due to reactions processes + load + outflow + dvoldt + exchanges
        #self.budget = [0.0] * (nproc + 4)
        #changes due to load + outflow + dvoldt + exchanges
        self.budget = [0.0] * 4

    def set_budget(self, val, i):
        self.budget[i] = val

    def add_budget(self, val, i):
        self.budget[i] += val

    def get_budget(self, i):
        return float(self.budget[i])
