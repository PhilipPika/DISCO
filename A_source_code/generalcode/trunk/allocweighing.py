# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/allocweighing.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
#=============================================================================
#
# This program's roots can be found in program INDUSTRY.FOR,
# part of WARIBAS, Water Assessment on a RIver BAsin Scale, 
# Olivier Klepper and Gerard van Drecht, 1997, RIVM, Netherlands
# Gerard van Drecht, June, 1998
# September 2004:
# Major changes were made for the purpose of a more generic approach:
# no units and no fixed mapextensions, no internal weighting function,
# instead the user should provide a weighting value map.
# September 2012:
# Program converted from fortran into Python.
#
# PURPOSE:
# Allocate per country available data, for example industrial water use,
# in a map. Use a so called carrier or weigthing factor map, for a
# distributed weighting factor. For example a population numbers map.
#
# METHOD:
# The weighting factor map is used because the weighting factor W is
# expected to have simular spatial distribution properties as variable Q,
# which we are trying to allocate in a grid.
# A weighting factor may be derived in several ways from multiple maps.
# Therefore this program only requires the resulting weighting factor map
# as input.
# The weighting factor ranges from 0 (no chance to find any Q),
# to an arbitrairy maximum value (very likely to find Q),
# so W and Q must be positively correlated.
#
# Allocation by weighting is done by the following formula:
#
# Q(i) = [W(i)/SW(j)] . SQ(j)  <= Qmax(i)                (1)
#
# where:
# Q(i) = allocation variable in cell i [unit]
# W(i) = weighting variable in cell i [any unit]
# SW(j)= sum of weights W(i) of region j [any unit]
# SQ(j)= sum of all values of allocated variable Q(i) of region j [unit]
# Qmax(i) = maximum value of Q(i) in cell i [unit]
#
# Qmax also is a spatially explicit variable. 
# In the process of allocating SQ(j) to the grid cells(i) of region j,
# temporally Q(i) may exceed Qmax(i). If that is the case we have a
# residual Qres in grid cell i. The sum of all residual Q's of region j,
# SQres(j) must be allocated in the next round in cells, where
# we still have any allocation space left, in other words Qmax-Q > 0.
# The allocation process is ended after 10 rounds or when the following
# breakoff criterion is fullfilled:
# SQres(j) < 0.001 SQ(j)
#
#==============================================================================================
def allocweighing(sq,sw,wReg,qmaxReg):

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
       ermax = 0.001 * sq
       done = len(wReg) * [0]
       for it in range(10):
           # No check on sw (divding by zero) because it is not possible here.
           scale = sq /sw
           sw = 0.0

           # get row and column number
           for icell in range(len(wReg)):
               if (wReg[icell] > 0.0 and done[icell] != 1):
                   qalloc = scale * wReg[icell]
                   if (qmaxReg[icell] >= 0.):
                       # deficit is room for more q in this cell
                       # print reg,icell,qmaxReg[icell], qReg[icell]
                       deficit= max(0.0,qmaxReg[icell] - qReg[icell])
                       if (qalloc>=deficit):
                           # cell is/wordt volgeboekt en gewicht op nul gezet
                           qalloc = deficit
                           done[icell] = 1
                           wReg[icell] = 0.0
                   qReg[icell] = qReg[icell] + qalloc
                   sw = sw + wReg[icell]
                   sq = sq - qalloc
           if (sq <= ermax or sw <= 0.0):
               # Allocation is finished
               break
    
    return qReg
