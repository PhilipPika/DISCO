# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import numpy

def do(ascii_object):
  lons = numpy.arange(ascii_object.xllcorner, ascii_object.xllcorner+(ascii_object.cellsize*ascii_object.ncols), ascii_object.cellsize)
  lats = numpy.arange(ascii_object.yllcorner, ascii_object.yllcorner+(ascii_object.cellsize*ascii_object.nrows), ascii_object.cellsize)
  grid = numpy.zeros((len(lats),len(lons)))
  for ilat in range(grid.shape[0]):
    for ilon in range(grid.shape[1]):
      grid[ilat, ilon] = ascii_object.get_data(ascii_object.get_index_from_row_col(ilat, ilon))
  return grid
