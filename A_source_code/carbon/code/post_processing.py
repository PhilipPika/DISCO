# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

try:
    # import general Anaconda3 Python modules
    import copy
    from netCDF4 import Dataset, MFDataset
    import numpy as np
    import operator
    import os
    import pickle
    import random
    import shutil
    import sys

    # import own general modules
    import ascraster
    import cmd_options_dgnm
    import define_subgrid_streamorder
    import directory
    import general_func
    import get_dynamic_gridinfo
    import interpolate_list
    import make_index_species
    import make_mask
    import manip
    import mocsy
    import output_conversion
    import pointer_class
    import read_parameter
except:    
    # set the general path for the own python modules
    import general_path

    # import general Anaconda3 Python modules
    from netCDF4 import Dataset, MFDataset
    import numpy as np
    import operator
    import os
    import pickle
    import random
    import shutil
    import sys

    # import own general modules
    import ascraster
    import cmd_options_dgnm
    import define_subgrid_streamorder
    import directory
    import general_func
    import get_dynamic_gridinfo
    import interpolate_list
    import make_index_species
    import manip
    import make_mask
    import mocsy
    import output_conversion
    import pointer_class
    import read_parameter

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

def calc_carb_dat(conc_DIC, conc_ALK, T):
    if (conc_DIC<1/100000) or (conc_ALK<1/100000):
        conc_DIC = 1e-6  
    with suppress_stdout_stderr():        
      carb_dat = np.dstack(mocsy.mvars(temp=T, sal=0, alk=conc_ALK/12, dic=conc_DIC/12, sil=0, phos=0, patm=1., depth=1, lat=0, optcon='mol/m3', optt='Tinsitu', optp='db', optb="u74", optk1k2='l', optkf="dg", optgas='Pinsitu'))[0][0]
    pH, pCO2 = operator.itemgetter(0,1)(carb_dat)
    if pH>14:
        pH = -999999999
    if pCO2>1e6:
        pCO2 = -999999999
    return pCO2, pH

def post_processing(params,dformat):
  print('POST PROCESSING DATA IN '+params.outputdir+' AS '+dformat)  
  np.seterr(all='ignore')
  folder = os.path.join(params.outputdir, '..', "STATES", "subgrid")
  if dformat=='NETCDF':
      alklist = sorted(directory.get_files_with_str(folder, 'conc*ALK*'))
      diclist = sorted(directory.get_files_with_str(folder, 'conc*DIC*'))
      Tlist = sorted(directory.get_files_with_str(os.path.join(params.outputdir, '..', "STREAM_ENV_CONDITIONS", "subgrid"), '*temperature*'))

      icell_list = ascraster.create_mask(params.file_mask, params.maskid, params.mask_bool_operator,numtype=int)
      dum_asc = ascraster.Asciigrid(ascii_file=params.file_mask)
      mask_2d = make_mask.do(params.file_mask, params.maskid, dum_asc, logical=params.mask_bool_operator, mask_type='np_grid')

      for order in range(len(alklist)):
          alk = Dataset(alklist[order], "a")
          dic = Dataset(diclist[order], "a")
          Tdat = Dataset(Tlist[order], "a")
          time_list = alk.variables['time'][:]
                  
          pCO2 = Dataset(os.path.join(folder,'pCO2_order'+str(order+1)+'.nc'), "w")
          pH = Dataset(os.path.join(folder,'pH_order'+str(order+1)+'.nc'), "w")
          pCO2_grid = np.ma.array(np.zeros((dum_asc.nrows,dum_asc.ncols)), mask=mask_2d).unshare_mask()
          pH_grid = np.ma.array(np.zeros((dum_asc.nrows,dum_asc.ncols)), mask=mask_2d).unshare_mask()
          
          output_conversion.init_ncdata(folder, pCO2, 'pCO2', dum_asc)
          output_conversion.init_ncdata(folder, pH, 'pH', dum_asc)      
 
          alk_matrix = alk['conc_ALK']
          dic_matrix = dic['conc_DIC']
          T_matrix = Tdat['temperature']
 
          for itime in range(len(time_list)):
              for icell in icell_list:
                  ilat, ilon = manip.calc_row_col_from_index(icell, dum_asc.ncols)
                  conc_alk = alk_matrix[itime,ilat,ilon]
                  conc_dic = dic_matrix[itime,ilat,ilon]
                  T = T_matrix[itime,ilat,ilon]-273
                  pCO2dat, pHdat = calc_carb_dat(conc_dic, conc_alk, T)
                  pH_grid[ilat, ilon] = pHdat
                  pCO2_grid[ilat, ilon] = pCO2dat
              manip.add_grid_time(pH, 'pH', pH_grid, time_list[itime])
              manip.add_grid_time(pCO2, 'pCO2', pCO2_grid, time_list[itime])
            
          alk.close()        
          dic.close()
          pCO2.close()
          pH.close()

      for fn in directory.get_files_with_str(folder, "*order6*"):
          shutil.copyfile(fn, os.path.join(folder, '..', os.path.basename(fn).replace("_order6", ""))) 

if __name__ == "__main__":
    # Set the general path for the own python modules
    import general_path 
    
    # Parse command-line arguments and set parameters for script
    try:
      param = cmd_options_dgnm.Input_dgnm(sys.argv)
      params = param.options
    except SystemExit:
      raise MyError("Error has occured in the reading of the commandline options.")      
    
    dformat = 'NETCDF'
    #post_processing(params, dformat)
    post_processing_states(params, dformat)
    #post_processing_budgets(params, dformat)	
    
			
