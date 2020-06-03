# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/test_ascraster.py $"
# ******************************************************
'''
Test script to test the functionality of ascraster.
'''
import os
import sys
__general = os.path.join(os.getcwd(), 'trunk')
if os.path.exists(__general): 
    sys.path.insert(0, __general)
    print(__general + " is added to the python search path for modules.") 

import ascraster
import my_sys

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

# Test 1 with power
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
grid2 = ascraster.Asciigrid(ascii_file='testgrid5.asc',numtype=int)
grid1.power(grid2,default_nodata_value=-12)

ltest1 = 1
if (grid1.nodata_value != -12):
    print('Test 6 power is not a succes. Setting of nodata goes wrong.')
    ltest1 = 0
for i in range(grid2.length):
    val = grid2.get_data(i)
    if (val == None):
        if (grid1.values[i] != -12.0):
            ltest1 = 0
            print('Test 6 power is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be -12.0 ")
if (ltest1 == 1):
    print("Power Test1 passed.")
else:
    print("Power with nodata values")
    print(grid1.values)

# Test 2 with power
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
# Testgrid 7 is zero grid without nodata
grid2 = ascraster.Asciigrid(ascii_file='testgrid7.asc',numtype=int)
grid1.power(grid2,default_nodata_value=-12)

ltest1 = 1
if (grid1.nodata_value != None):
    print('Test 2 power is not a succes. Setting of nodata goes wrong.')
    ltest1 = 0
for i in range(grid1.length):
    val = grid1.get_data(i)
    if (val != None):
        if (grid1.values[i] != 1.0):
            ltest1 = 0
            print('Test 2 power is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be 1.0 ")
if (ltest1 == 1):
    print("Power Test 2 passed.")
else:
    print("Power of zero, so outcome must be 1.")
    print(grid1.values)

# Test 3 with power
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
# Testgrid 7 is zero grid with nodata
grid2 = ascraster.Asciigrid(ascii_file='testgrid8.asc',numtype=int)
grid1.power(grid2,default_nodata_value=-12)

ltest1 = 1
if (grid1.nodata_value != -12):
    print('Test 3 power is not a succes. Setting of nodata goes wrong.')
    ltest1 = 0
for i in range(grid1.length):
    val = grid1.get_data(i)
    if (val != None):
        if (grid1.values[i] != 1.0):
            ltest1 = 0
            print('Test 3 power is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be 1.0 ")
if (ltest1 == 1):
    print("Power Test 3 passed.")
else:
    print("Power of zero, so outcome must be 1.")
    print(grid1.values)

# Test 4 with power
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
# Testgrid 7 is zero grid with nodata
grid2 = ascraster.Asciigrid(ascii_file='testgrid9.asc',numtype=int)
grid1_old = ascraster.duplicategrid(grid1)
grid1.power(grid2,default_nodata_value=-12)
import math
ltest1 = 1
if (grid1.nodata_value != -12):
    print('Test 4 power is not a succes. Setting of nodata goes wrong.')
    ltest1 = 0
for i in range(grid1.length):
    val = grid1.get_data(i)
    if (val != None):
        val_old = grid1_old.get_data(i)
        pow = grid2.get_data(i)
        if (grid1.values[i] != math.pow(val_old,pow)):
            ltest1 = 0
            print('Test 4 power is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be ",math.pow(val_old,pow))
if (ltest1 == 1):
    print("Power Test 4 passed.")
else:
    print("Power of zero, so outcome must be 1.")
    print(grid1.values)

# Test 5 with power
grid1 = ascraster.Asciigrid(ascii_file='testgrid9.asc',numtype=int)
grid2 = ascraster.Asciigrid(ascii_file='testgrid10.asc',numtype=int)
grid1.power(grid2,default_nodata_value=-12)

ltest1 = 1
if (grid1.nodata_value != -12):
    print('Test 5 power is not a succes. Setting of nodata goes wrong.')
    ltest1 = 0
for i in range(grid1.length):
    val = grid1.get_data(i)
    if (val != None):
        if (grid1.values[i] != -1.0):
            ltest1 = 0
            print('Test 5 power is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be -1.0 ")
if (ltest1 == 1):
    print("Power Test 5 passed.")
else:
    print("Power of zero, so outcome must be -1.")
    print(grid1.values)
    
# Test 6 with power
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
grid1.power(0.0,default_nodata_value=-12)
ltest1 = 1
if (grid1.nodata_value != None):
    print('Test 6 power is not a succes. Setting of nodata goes wrong.')
    ltest1 = 0
for i in range(grid1.length):
    val = grid1.get_data(i)
    if (val != None):
        if (grid1.values[i] != 1.0):
            ltest1 = 0
            print('Test 6 power is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be 1.0 ")
if (ltest1 == 1):
    print("Power Test 6 passed.")
else:
    print("Power of zero, so outcome must be 1.")
    print(grid1.values)    
    
# Test 7 with power
grid1 = ascraster.Asciigrid(ascii_file='testgrid5.asc',numtype=int)
nodat = grid1.nodata_value
grid1.power(0.0,default_nodata_value=-12)
ltest1 = 1
if (grid1.nodata_value != nodat):
    print('Test 7 power is not a succes. Setting of nodata goes wrong.')
    ltest1 = 0
for i in range(grid1.length):
    val = grid1.get_data(i)
    if (val != None):
        if (grid1.values[i] != 1.0):
            ltest1 = 0
            print('Test 7 power is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be 1.0 ")
if (ltest1 == 1):
    print("Power Test 7 passed.")
else:
    print("Power of zero, so outcome must be 1.")
    print(grid1.values)      
    
# Test 8 with power
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
grid1.add_values(grid1.length*[-1])
grid1.power(1.5,default_nodata_value=-12)
ltest1 = 1
if (grid1.nodata_value != -12):
    print('Test 8 power is not a succes. Setting of nodata goes wrong.')
    ltest1 = 0
for i in range(grid1.length):
    val = grid1.get_data(i)
    if (val != None):
        ltest1 = 0
        print('Test 8 power is not a succes for item: ' + str(i) + ". Values found: " + str(grid1.values[i]) + " Must be nodata ")
if (ltest1 == 1):
    print("Power Test 8 passed.")
else:
    print("Power of negative value to the power .., so outcome must be nodata")
    print(grid1.values)

# Test 1 write ascii file, all cells have data
grid1 = ascraster.Asciigrid(ascii_file='testgrid1.asc',numtype=int)
grid1.write_ascii_file("testout1.asc")
msg = os.system("diff testgrid1.asc testout1.asc")
if (msg != 0):
    print('Test 1 write ascii file is not a succes. Differences above.')
else:
    print('Test 1 write ascii file passed.')
    my_sys.my_removefile("testout1.asc")

# Test 2 write ascii file, grid with nodata values
grid1 = ascraster.Asciigrid(ascii_file='testgrid8.asc',numtype=int)
grid1.write_ascii_file("testout2.asc")
msg = os.system("diff testgrid8.asc testout2.asc")
if (msg != 0):
    print('Test 2 write ascii file is not a succes. Differences above.')
else:
    print('Test 2 write ascii file passed.')
    my_sys.my_removefile("testout2.asc")

# Test 3 write ascii file, grid with nodata values and mask
mask = ascraster.create_mask('testgrid8.asc',0, logical = "GE",numtype=int)
grid1 = ascraster.Asciigrid(ascii_file='testgrid8.asc',mask=mask,numtype=int)
grid1.write_ascii_file("testout3.asc")
msg = os.system("diff testgrid8.asc testout3.asc")
if (msg != 0):
    print('Test 3 write ascii file is not a succes. Differences above.')
else:
    print('Test 3 write ascii file passed.')
    my_sys.my_removefile("testout3.asc")

# Test 4 write ascii file, grid with nodata values and mask with all cells have a value
mask = ascraster.create_mask('testgrid8.asc',0, logical = "GE",numtype=int)
#Add two cells with nodata to the mask
mask.append(16)
mask.append(17)
grid1 = ascraster.Asciigrid(ascii_file='testgrid8.asc',mask=mask,numtype=int)
grid1.write_ascii_file("testout4.asc")
msg = os.system("diff testgrid8.asc testout4.asc")
if (msg != 0):
    print('Test 4 write ascii file is not a succes. Differences above.')
else:
    print('Test 4 write ascii file passed.')
    my_sys.my_removefile("testout4.asc")

# Test 5 write ascii file, grid with nodata values and mask with all cells have nodata value
mask = ascraster.create_mask('testgrid8.asc',0, logical = "GE",numtype=int)
#Add two cells with nodata to the mask
mask.append(16)
mask.append(17)
grid1 = ascraster.Asciigrid(ascii_file='testgrid8.asc',mask=mask,numtype=int)
grid1.add_values(grid1.length * [grid1.nodata_value])
grid1.write_ascii_file("testout5.asc")
msg = os.system("diff testgrid11.asc testout5.asc")
if (msg != 0):
    print('Test 5 write ascii file is not a succes. Differences above.')
else:
    print('Test 5 write ascii file passed.')
    my_sys.my_removefile("testout5.asc")


# Test 6 write ascii file, grid with nodata values with all cells have nodata value
grid1 = ascraster.Asciigrid(ascii_file='testgrid8.asc',numtype=int)
grid1.add_values(grid1.length * [grid1.nodata_value])
grid1.nodata_value = None
grid1.write_ascii_file("testout6.asc",output_nodata_value = -2)
msg = os.system("diff testgrid11.asc testout6.asc")
if (msg != 0):
    print('Test 6 write ascii file is not a succes. Differences above.')
else:
    print('Test 6 write ascii file passed.')
    my_sys.my_removefile("testout6.asc")

# Test 7 write ascii file, grid with nodata values and mask with all cells have nodata value
mask = ascraster.create_mask('testgrid8.asc',0, logical = "GE",numtype=int)
grid1 = ascraster.Asciigrid(ascii_file='testgrid8.asc',mask=mask,numtype=int)
grid1.add_values(grid1.length * [None])
grid1.nodata_value = None
grid1.write_ascii_file_slow("testout7.asc",output_nodata_value = -2)
msg = os.system("diff testgrid11.asc testout7.asc")
if (msg != 0):
    print('Test 7 write ascii file is not a succes. Differences above.')
else:
    print('Test 7 write ascii file passed.')
    my_sys.my_removefile("testout7.asc")

