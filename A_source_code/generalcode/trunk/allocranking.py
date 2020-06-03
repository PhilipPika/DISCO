# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/allocranking.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
#
# This program's roots can be found in program INDUSTRY.FOR,
# part of WARIBAS, Water Assessment on a RIver BAsin Scale, 
# Olivier Klepper and Gerard van Drecht, 1997, RIVM, Netherlands
# Gerard van Drecht, June, 1998
# September 2004:
# Major changes were made for the purpose of a more generic approach:
# no units and no fixed mapextensions, no internal weighting function,
# instead the user should provide a weighting value map.

# February 2006:
# no more ascii grid map output (too big)
# March 2006: allocation by ranking is implemented,
#
# ranking is done by numerical sorting of the weigting factor,
# thus ranking can be seen as a kind of extreme weighting.
# September 2012:
# Program converted from fortran into Python.
#
# USES:
# SORTIND(N,ARRIN,INDX) (Numerical Recipes, Press et al.)
#
#==============================================================================================

from operator import itemgetter

def allocranking(sq,sw,wReg,qmaxReg):

    #===================================================================================
    #INPUT (All input is changed during this function)
    # sq               Regional sum of values of the variable that must be allocated
    # sw               Regional sum of values of weighing factor
    # wReg             Weighting factor grid map
    # qmaxReg          Maximum value of allocation variable per grid cell

    #OUTPUT
    # qReg             Values of the variable that must be allocated 
    #===================================================================================

    qReg = len(wReg) * [0.0]
    if (sq > 0.0 and sw > 0.0 ):
        # fill list with index and weighting factor of the grid map.
        indx_old = []
        for icell in range(len(wReg)):
           indx_old.append([wReg[icell],icell])

        pool = sq

        # sort cells in ascending order, so the last one has the highest weight
        indx = sorted(indx_old, key=itemgetter(0))
        for k in range(len(wReg)-1,-1,-1):
            if (pool > 0.0):
                i = indx[k][1]
                qReg[i] = min(pool,qmaxReg[i])
                pool = pool - qReg[i]
                if (pool < 0.0):
                    pool = 0.0
    return qReg
