# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

from netCDF4 import Dataset
import numpy as np

import general_path
import accuflux
import ascraster
import get_surrounding_cells
import make_np_grid


def do(mask_asc_fn, mask_id, dum_asc, logical = "EQ", mask_type='np_grid'):

    dum_mask = ascraster.create_mask(mask_asc_fn, mask_id, logical = logical, numtype=int)
    mask=[]
    if mask_type=="rowcol":
      for i in dum_mask:
        mask.append(dum_asc.get_row_col_from_index(i))
    elif mask_type=="index":
      for i in dum_mask:
        mask.append(i)
    elif mask_type=="latlon":
      for i in dum_mask:
        mask.append(dum_asc.get_coord_from_index(i))  
    elif mask_type=="np_grid":
      mask = np.zeros((dum_asc.nrows, dum_asc.ncols), dtype=bool)
      mask[:,:] = True
      for i in dum_mask:
        row, col = dum_asc.get_row_col_from_index(i)
        mask[row,col]=False
    return mask
