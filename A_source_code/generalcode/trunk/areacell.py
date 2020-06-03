# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/areacell.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import math


#radius of the earth in km
R_earth = 6356.751 

def areacell(cellsize,lat_midpoint):
  '''
  http://mathforum.org/library/drmath/view/63767.html
  A = 2*pi*R^2 |sin(lat1)-sin(lat2)| |lon1-lon2|/360
    = (pi/180)R^2 |sin(lat1)-sin(lat2)| |lon1-lon2|
  R is from Gerard.
  Return value is area in km2.  
  '''
  lat1 = math.pi * (lat_midpoint + 0.5 * cellsize)/180.
  lat2 = math.pi * (lat_midpoint - 0.5 * cellsize)/180. 
  return (math.pi/180.) *  R_earth * R_earth * abs(math.sin(lat1) - math.sin(lat2)) * cellsize
