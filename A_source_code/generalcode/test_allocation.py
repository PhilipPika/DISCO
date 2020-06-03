# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/test_allocation.py $"
# ******************************************************
'''
Test script to test the functionality of allocation functions.
'''
import os
import sys
__general = os.path.join(os.getcwd(), 'trunk')
if os.path.exists(__general): 
    sys.path.insert(0, __general)
    print(__general + " is added to the python search path for modules.") 

import allocranking
import allocweighing


#def allocranking(sq,sw,wReg,qmaxReg):

    #===================================================================================
    #INPUT (All input is changed during this function)
    # sq               Regional sum of values of the variable that must be allocated
    # sw               Regional sum of values of weighing factor
    # wReg             Weighting factor grid map
    # qmaxReg          Maximum value of allocation variable per grid cell
# Test 1
sq = 100.
wReg = [1.,1.,2.,1.5]
sw = sum(wReg)
qmaxReg = [10.,10.,75.,50.]
qReg =  allocranking.allocranking(sq,sw,wReg,qmaxReg)
qReg_exp = [0.0,0.0,75.0,25.0]
ltest1 = 1
if (qReg != qReg_exp):
    ltest1 = 0
    print("Test 1 is not a succes. Values found: " + str(qReg) + " and " + str(qReg_exp))
if (ltest1 == 1):
    print("Test1 passed.")
else:
    print("Allocation ranking method")
    print(qReg)
    print(qReg_exp)

# Test 2
sq = 100.
wReg = [1.,1.,2.,1.]
sw = sum(wReg)
qmaxReg = [100.,20.,10.,50.]
qReg =  allocranking.allocranking(sq,sw,wReg,qmaxReg)
qReg_exp = [20.0,20.0,10.0,50.0]
ltest1 = 1
if (qReg != qReg_exp):
    ltest1 = 0
    print("Test 2 is not a succes. Values found: " + str(qReg) + " and " + str(qReg_exp))
if (ltest1 == 1):
    print("Test2 passed.")
else:
    print("Allocation ranking method 2")
    print(qReg)
    print(qReg_exp)

# Test 3
sq = 100.
wReg = [1.,1.,2.,1.]
sw = sum(wReg)
qmaxReg = [1.,0.,0.,0.]
qReg =  allocranking.allocranking(sq,sw,wReg,qmaxReg)
qReg_exp = [1.0,0.0,0.0,0.0]
ltest1 = 1
if (qReg != qReg_exp):
    ltest1 = 0
    print("Test 3 is not a succes. Values found: " + str(qReg) + " and " + str(qReg_exp))
if (ltest1 == 1):
    print("Test3 passed.")
else:
    print("Allocation ranking method 3")
    print(qReg)
    print(qReg_exp)

# Testing allocweighing
print("Start allocation weighing method testing.") 
# Test 1
sq = 100.
wReg = [1.,1.,2.,6]
sw = sum(wReg)
qmaxReg = [100.,100.,100.,100.]
qReg =  allocweighing.allocweighing(sq,sw,wReg,qmaxReg)
qReg_exp = [10.0,10.0,20.0,60.0]
ltest1 = 1
if (qReg != qReg_exp):
    ltest1 = 0
    print("Test 1 is not a succes. Values found: " + str(qReg) + " and " + str(qReg_exp))
if (ltest1 == 1):
    print("Test1 passed.")
else:
    print("Allocation weighing method")
    print(qReg)
    print(qReg_exp)

# Test 2
sq = 100.
wReg = [1.,1.,2.,6]
sw = sum(wReg)
qmaxReg = [100.,100.,100.,50.]
qReg =  allocweighing.allocweighing(sq,sw,wReg,qmaxReg)
qReg_exp = [12.5,12.5,25.0,50.0]
ltest1 = 1
if (qReg != qReg_exp):
    ltest1 = 0
    print("Test 2 is not a succes. Values found: " + str(qReg) + " and " + str(qReg_exp))
if (ltest1 == 1):
    print("Test2 passed.")
else:
    print("Allocation weighing method 2")
    print(qReg)
    print(qReg_exp)

# Test 3
sq = 80.
wReg = [0,0,2.,6]
sw = sum(wReg)
qmaxReg = [100.,100.,100.,50.]
qReg =  allocweighing.allocweighing(sq,sw,wReg,qmaxReg)
qReg_exp = [0,0,30.0,50.0]
ltest1 = 1
if (qReg != qReg_exp):
    ltest1 = 0
    print("Test 3 is not a succes. Values found: " + str(qReg) + " and " + str(qReg_exp))
if (ltest1 == 1):
    print("Test3 passed.")
else:
    print("Allocation weighing method 3")
    print(qReg)
    print(qReg_exp)



sys.exit(12)
# First three tests are on a grid with no nodata
# Test 1
# Multiply grid1 with a scalar of type integer
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
factor = 1
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 1 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test1 passed.")
else:
    print("Multiplying with factor 1")
    print(grid1_old.values)
    print(grid1.values)

# Test 2
# Multiply grid1 with a scalar of type float
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
factor = 1.0
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 2 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test2 passed.")
else:
    print("Multiplying with factor 1.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 3
# Multiply grid1 with a grid of one
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
factor = 1.0
grid1_old = ascraster.duplicategrid(grid1)
gridone = ascraster.duplicategrid(grid1)
# Make all grid entries one.
for i in range(gridone.length):
    gridone.set_data(i,1.0)

grid1.multiply(gridone)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 3 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test3 passed.")
else:
    print("Multiplying with another grid with 1.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# First three tests are on a grid with no nodata. Now with grid with no nodata but in the header the nodata_value specified.
# Test 4
# Multiply grid1 with a scalar of type integer
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
factor = 1
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 4 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test4 passed.")
else:
    print("Multiplying with factor 1")
    print(grid1_old.values)
    print(grid1.values)

# Test 5
# Multiply grid1 with a scalar of type float
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
factor = 1.0
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 5 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test5 passed.")
else:
    print("Multiplying with factor 1.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 6
# Multiply grid1 with a grid of one
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
factor = 1.0
grid1_old = ascraster.duplicategrid(grid1)
gridone = ascraster.duplicategrid(grid1)
# Make all grid entries one.
for i in range(gridone.length):
    gridone.set_data(i,1.0)

grid1.multiply(gridone)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 6 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test6 passed.")
else:
    print("Multiplying with another grid with 1.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# First six tests are on a grid with no nodata. Now with grid with nodata.
# Test 7
# Multiply grid1 with a scalar of type integer
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
factor = 1
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 7 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test7 passed.")
else:
    print("Multiplying with factor 1")
    print(grid1_old.values)
    print(grid1.values)

# Test 8
# Multiply grid1 with a scalar of type float
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
factor = 1.0
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 8 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test8 passed.")
else:
    print("Multiplying with factor 1.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 9
# Multiply grid1 with a grid of one
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
factor = 1.0
grid1_old = ascraster.duplicategrid(grid1)
gridone = ascraster.duplicategrid(grid1)
# Make all grid entries one.
for i in range(gridone.length):
    gridone.set_data(i,1.0)

grid1.multiply(gridone)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 9 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test9 passed.")
else:
    print("Multiplying with another grid with 1.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 10
# Multiply grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(gridone)
# Grid1 must have the multiplication of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1*val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 10 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test10 passed.")
else:
    print("Multiplying with another grid with 1.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)    

# Test 11
# Multiply grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(gridone)
# Grid1 must have the multiplication of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1*val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 11 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test11 passed.")
else:
    print("Multiplying with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values) 

# Test 12
# Multiply grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(gridone)
# Grid1 must have the multiplication of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1*val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 12 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test12 passed.")
else:
    print("Multiplying with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 13
# Multiply grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(gridone)
# Grid1 must have the multiplication of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1*val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 13 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test13 passed.")
else:
    print("Multiplying with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 14
# Multiply grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(gridone)
# Grid1 must have the multiplication of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1*val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 14 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test14 passed.")
else:
    print("Multiplying with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 15
# Multiply grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.multiply(gridone)
# Grid1 must have the multiplication of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1*val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 15 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test15 passed.")
else:
    print("Multiplying with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)    

# Now tests for the summation of two objects (grid or scalar). 
# First three tests are on a grid with no nodata
# Test 1
# Add grid1 with a scalar of type integer
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
factor = 0
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 1 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test1 passed.")
else:
    print("Sum with factor 0")
    print(grid1_old.values)
    print(grid1.values)

# Test 2
# Add grid1 with a scalar of type float
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
factor = 0.0
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 2 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test2 passed.")
else:
    print("Sum with factor 0.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 3
# Add grid1 with a grid of one
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
factor = 0.0
grid1_old = ascraster.duplicategrid(grid1)
gridone = ascraster.duplicategrid(grid1)
# Make all grid entries one.
for i in range(gridone.length):
    gridone.set_data(i,0.0)

grid1.add(gridone)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 3 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test3 passed.")
else:
    print("Adding with another grid with 0.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# First three tests are on a grid with no nodata. Now with grid with no nodata but in the header the nodata_value specified.
# Test 4
# Multiply grid1 with a scalar of type integer
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
factor = 0
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 4 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test4 passed.")
else:
    print("Sum with factor 0")
    print(grid1_old.values)
    print(grid1.values)

# Test 5
# Multiply grid1 with a scalar of type float
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
factor = 0.0
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 5 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test5 passed.")
else:
    print("Sum with factor 0.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 6
# Sum grid1 with a grid of zeros
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
factor = 0.0
grid1_old = ascraster.duplicategrid(grid1)
gridone = ascraster.duplicategrid(grid1)
# Make all grid entries one.
for i in range(gridone.length):
    gridone.set_data(i,0.0)

grid1.add(gridone)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 6 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test6 passed.")
else:
    print("Multiplying with another grid with 1.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# First six tests are on a grid with no nodata. Now with grid with nodata.
# Test 7
# Sum grid1 with a scalar of type integer
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
factor = 0
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 7 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test7 passed.")
else:
    print("Sum with factor 0")
    print(grid1_old.values)
    print(grid1.values)

# Test 8
# Sum grid1 with a scalar of type float
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
factor = 0.0
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 8 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test8 passed.")
else:
    print("Summing with factor 0.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 9
# Sum grid1 with a grid of zeros
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
factor = 0.0
grid1_old = ascraster.duplicategrid(grid1)
gridone = ascraster.duplicategrid(grid1)
# Make all grid entries one.
for i in range(gridone.length):
    gridone.set_data(i,0.0)

grid1.add(gridone)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 9 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test9 passed.")
else:
    print("Sum with another grid with 0.0. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 10
# Multiply grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(gridone)
# Grid1 must have the sum of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1+val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 10 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test10 passed.")
else:
    print("Sum with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)    

# Test 11
# Multiply grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(gridone)
# Grid1 must have the sum of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1+val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 11 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test11 passed.")
else:
    print("Sum with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values) 

# Test 12
# Multiply grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(gridone)
# Grid1 must have the sum of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1+val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 12 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test12 passed.")
else:
    print("Summing with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 13
# Sum grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(gridone)
# Grid1 must have the multiplication of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1+val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 13 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test13 passed.")
else:
    print("Sum with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 14
# Multiply grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(gridone)
# Grid1 must have the multiplication of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1+val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 14 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test14 passed.")
else:
    print("Sum with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)

# Test 15
# Sum grid1 with nodata with a grid with no nodata
# Read ascii grid
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(gridone)
# Grid1 must have the sum of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        grid1_old.set_data(i,val1+val2)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 15 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test15 passed.")
else:
    print("Sum with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values)   

# Check range checkers
# Test 1
# Sum grid1 with nodata with a grid with no nodata
# Read ascii grid
xmin=4
xmax=7
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid2.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(gridone,minimum= xmin, maximum=xmax)
# Grid1 must have the sum of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        val3 = val1+val2
        if (val3 < xmin): val3 = xmin
        if (val3 > xmax): val3 = xmax
        grid1_old.set_data(i,val3)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 1 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test1 passed.")
else:
    print("Range checker. Sum with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values) 

# Test 2
# Sum grid1 with nodata with a grid with nodata
# Read ascii grid
xmin=4
xmax=7
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid4.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(gridone,minimum= xmin, maximum=xmax)
# Grid1 must have the sum of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        val3 = val1+val2
        if (val3 < xmin): val3 = xmin
        if (val3 > xmax): val3 = xmax
        grid1_old.set_data(i,val3)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 2 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test2 passed.")
else:
    print("Range checker. Sum with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values) 

# Test 3
# Sum grid1 with nodata with a grid with nodata
# Read ascii grid
xmin=4
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid4.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(gridone,minimum= xmin)
# Grid1 must have the sum of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        val3 = val1+val2
        if (val3 < xmin): val3 = xmin
        grid1_old.set_data(i,val3)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 3 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test3 passed.")
else:
    print("Range checker. Sum with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values) 

# Test 4
# Sum grid1 with nodata with a grid with nodata
# Read ascii grid
xmax=7
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
gridone = ascraster.Asciigrid(ascii_file='testgrid4.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.add(gridone, maximum=xmax)
# Grid1 must have the sum of grid1 and gridone
# Calculation is done on a different way
for i in range(grid1_old.length):
    val1 = grid1_old.get_data(i)
    val2 = gridone.get_data(i)
    if (val1 == None or val2 == None):
        grid1_old.set_data(i,grid1_old.nodata_value)
    else:
        val3 = val1+val2
        if (val3 > xmax): val3 = xmax
        grid1_old.set_data(i,val3)

ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 4 is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Test4 passed.")
else:
    print("Range checker. Sum with another grid. Changing int grid into float")
    print(grid1_old.values)
    print(grid1.values) 

# Test 1 with divide
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
factor = 1
grid1_old = ascraster.duplicategrid(grid1)
grid1.divide(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1_old.length):
    if (grid1.values[i] != grid1_old.values[i]):
        ltest1 = 0
        print('Test 1 divide is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " and " + str(grid1_old.values[i]))
if (ltest1 == 1):
    print("Divide Test1 passed.")
else:
    print("Divide with factor 1")
    print(grid1_old.values)
    print(grid1.values)

# Test 2 with divide
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
factor = ascraster.duplicategrid(grid1)
grid1.divide(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1.length):
    if (grid1.values[i] != 1.0):
        ltest1 = 0
        print('Test 2 divide is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be 1.0 ")
if (ltest1 == 1):
    print("Divide Test2 passed.")
else:
    print("Divide with itself")
    print(grid1.values)

# Test 3 with divide
grid1 = ascraster.Asciigrid(ascii_file='testgrid3.asc',numtype=int)
factor = ascraster.duplicategrid(grid1)
grid1.divide(factor)
# Grid1 and grid1_old must be the same
ltest1 = 1
for i in range(grid1.length):
    val = grid1.get_data(i)
    if (val != None):
        if (grid1.values[i] != 1.0):
            ltest1 = 0
            print('Test 3 divide is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be 1.0 ")
if (ltest1 == 1):
    print("Divide Test3 passed.")
else:
    print("Divide with itself")
    print(grid1.values)

# Test 4 with divide
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
grid2 = ascraster.Asciigrid(ascii_file='testgrid5.asc',numtype=int)
grid1.divide(grid2,default_nodata_value=-12)

ltest1 = 1
if (grid1.nodata_value != -12):
    print('Test 4 divide is not a succes. Setting of nodata goes wrong.')
    ltest1 = 0
for i in range(grid1.length):
    val = grid1.get_data(i)
    if (val != None):
        if (grid1.values[i] != 1.0):
            ltest1 = 0
            print('Test 4 divide is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be 1.0 ")
if (ltest1 == 1):
    print("Divide Test4 passed.")
else:
    print("Divide with itself")
    print(grid1.values)

# Test 5 with divide
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
grid2 = ascraster.Asciigrid(ascii_file='testgrid6.asc',numtype=int)
grid1.divide(grid2,default_nodata_value=-12)

ltest1 = 1
if (grid1.nodata_value != -12):
    print('Test 5 divide is not a succes. Setting of nodata goes wrong.')
    ltest1 = 0
for i in range(grid1.length):
    val = grid1.get_data(i)
    if (val != None):
        if (grid1.values[i] != 1.0):
            ltest1 = 0
            print('Test 5 divide is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be 1.0 ")
if (ltest1 == 1):
    print("Divide Test5 passed.")
else:
    print("Divide with itself")
    print(grid1.values)
