LIBRARY TO CALCULATE INORGANIC CARBON CHEMISTRY

## DESCRIPTION FOR LINUX ##

ROUTINES DOWNLOADED FROM:
http://ocmip5.ipsl.jussieu.fr/mocsy/index.html

IMPORTANT NOTE:
############################################
before compiling please change in vars.f90
lines that contain:
    !f2py logical optional, intent(in) :: verbose = .true.
to: 
    !f2py logical optional, intent(in) :: verbose = 1
############################################