# ******************************************************
## Revision "$LastChangedDate: 2018-12-04 10:26:58 +0100 (di, 04 dec 2018) $"
## Date "$LastChangedRevision: 14 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/mapcomparison.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# P(A): Fraction of Agreement.
# P(E): Expected fraction of Agreement subject to observed distribution.
# P(max): Maximal fraction of Agreement subject to observed distribution.
# Assume two maps A and B with N classes.
# p_ij = 0.1 means that 10% of all the cells of map A with class i have class j in map B.
# p_iT = fraction of all the cells of map A with class i.
# p_Ti = fraction of all the cells of map B with class i. 
# P(A) = SUM(over Nclasses)_(p_ii)
# P(E) = SUM(over Nclasses)_(p_iT * P_Ti)
# P(max) = SUM(over Nclasses)_(min(p_iT , P_Ti))

# Import general python modules
import math


# These methods are only used on integer maps!
from error import *
from iround import *
import ascraster

def write_diffmap(mapA,mapB,outputfilename):
    '''
    Create the cell-by-cell comparison of all cells and write this to ans ascii raster file.
    '''
    if (mapA.length != mapB.length):
        raise MyError("COMPARE_RASTERS: To compare two rasters, length should be equal.")

    # Duplicate output map.
    outgrid = ascraster.duplicategrid(mapA)

    # Set all data on nodata.
    if (not outgrid.nodata):
        outgrid.nodata_value = -9999.
    outgrid.add_values(outgrid.length *[outgrid.nodata_value])

    # Create the diffences map.
    for icell in range(mapA.length):
        valA = mapA.get_data(icell)
        valB = mapB.get_data(icell)
        if (valA == None):
            # Nodata is excluded.
            continue
        if (valB == None):
            # Nodata is excluded.
            continue

        valA = iround(valA)
        valB = iround(valB)

        # Set this cell.
        if (valA == valB):
            outgrid.set_data(icell, 1)
        else:
            outgrid.set_data(icell, 0)

    # Write to output file
    outgrid.write_ascii_file(outputfilename)

def calculate_chi2_matrix(mapA,mapB,matrixfile=None,weighmap=None):
    '''
    Calculate the cell-by-cell chi square matrix of all cells. The cells contain the number of cells with this combination.
    arg mapA: input map with classes for comparison.
    arg mapB: input map with classes for comparison.
    arg matrixfile: name of outputfile with the calculated matrix.
    arg weighmap: input map with weighing for each cell.
    Returned: matrix[i][j]: Chi-squared of map A with class i and class j in map B.
    Returned: ncount: total number of cells (or total sum of weighing) which are included in this analyses.
    '''
    if (mapA.length != mapB.length):
        raise MyError("CALCULATE_CHI2_MATRIX: To compare two rasters, length should be equal.")

    if (weighmap != None):
        if (mapA.length != weighmap.length):
            raise MyError("CALCULATE_CHI2_MATRIX: Length of weighing map should be equal to the input maps.")

    # Fill the matrix with the frequency values.
    ncount = 0
    count_dict = {}
    for icell in range(mapA.length):
        valA = mapA.get_data(icell)
        valB = mapB.get_data(icell)
        if (valA == None):
            # Nodata is excluded.
            continue
        if (valB == None):
            # Nodata is excluded.
            continue
        if (weighmap != None):
            valweigh = weighmap.get_data(icell)
            if (valweigh == None):
                # Nodata is excluded, also in weighing map.
                print("CALCULATE_CHI2_MATRIX: Nodata found in weighing map, which can influence the comparison.")
                continue
            else:
                # Check whether the weighing is positive everywhere.
                if (valweigh < 0.):
                    raise MyError("CALCULATE_CHI2_MATRIX: Negaive values found in weighing map, which influence the comparison.")
        else:
            valweigh = 1

        valA = iround(valA)
        valB = iround(valB)
        ncount += valweigh
        # Add this cell to the matrix.
        #matrix[valA-min_num_rows][valB-min_num_cols] += 1
        try:
           count_dict[str(valA)+"_" + str(valB)] += valweigh
        except KeyError:
           count_dict[str(valA)+"_" + str(valB)] = valweigh
    #print("MATRIX: ", matrix)
    #print("count_dict: ", count_dict)
    
    # Get the dimensions of the matrix by examing the dictionary.
    row_dict = {}
    col_dict = {}
    Nrows = 0
    Ncols = 0
    for key in count_dict:
        irow,icol = key.split("_")
        try:
            qq = row_dict[iround(irow)]
        except KeyError:
            row_dict[iround(irow)] = Nrows
            Nrows += 1            

        try:
            qq = col_dict[iround(icol)]
        except KeyError:
            col_dict[iround(icol)] = Ncols
            Ncols += 1  
    
    # Sort the keys in the matrix, so the one to one are on the diagonal.
    rowlist = sorted(row_dict)
    collist = sorted(col_dict)
    # Reassign the row and columnnumbers
    count = 0
    for item in rowlist:
        row_dict[item] = count
        count += 1
    count = 0
    for item in collist:
        col_dict[item] = count
        count += 1
        
    #print("COL_DICT: ", col_dict)
    #print("ROW_DICT: ", row_dict)
    # Make matrix
    # Create an matrix for all these values.
    matrix = []
    for irow in range(Nrows):
        matrix.append(Ncols*[0.0])
    #print "Nrows,Ncols",Nrows,Ncols
    for key in count_dict:
        irow,icol = key.split("_")
        irow_mat = row_dict[iround(irow)]
        icol_mat = col_dict[iround(icol)]
        #print key, irow_mat,icol_mat
        matrix[irow_mat][icol_mat] = count_dict[key]

    # Calculate the sum of the columns and the rows.
    P_iT = Nrows * [0.0]
    P_Ti = Ncols * [0.0]
    for irow in range(Nrows):
        for icol in range(Ncols):
            P_iT[irow] += matrix[irow][icol]
            P_Ti[icol] += matrix[irow][icol]

    # Write matrix to output file when asked for:
    if (matrixfile != None):
        sep = ";"
        fp = open(matrixfile,"w")
        # Write the column headers
        fp.write(sep)
        for icol in range(Ncols):
            fp.write(str(collist[icol]) + sep)
        fp.write(str("Total") + "\n")
        
        # Write each row to the output file
        for irow in range(Nrows):
            fp.write(str(rowlist[irow]) + sep)
            for icol in range(Ncols):
                fp.write(str(matrix[irow][icol]) + sep)
            fp.write(str(P_iT[irow]) + "\n")
        
        # Write the sum of the columns
        fp.write("Total" + sep)
        for icol in range(Ncols):
            fp.write(str(P_Ti[icol]) + sep)
        fp.write(str(ncount) + "\n")        
        


    # Based on this matrix the expected distribution for p_i,j (when they are not dependent on each other).
    # Ep_i,j = p_.,j * p_i,. /ncount with p_.,j is the sum of column j. and p_i,. the sum of row i.

    # The Chi-squared method is defined as:
    # chi2[i][j] = (matrix[i][j] - Ep_i,j)^2/Ep_i,j

    # Calculate the Chi-square matrix.
    for irow in range(Nrows):
        for icol in range(Ncols):
            Pexpect = P_iT[irow] * P_Ti[icol]/float(ncount)
            diff = matrix[irow][icol] - Pexpect
            try:
                matrix[irow][icol] = (diff * diff)/Pexpect
            except ZeroDivisionError:
                 matrix[irow][icol] = 0.0
            #if (irow == icol):
            #    print irow, ncount,P_iT[irow] ,P_Ti[icol],diff,Pexpect,matrix[irow][icol]
    #print("MATRIX: ", matrix)
    return matrix,ncount

def calculate_chi2_matrix_old(mapA,mapB,matrixfile=None):
    '''
    Calculate the cell-by-cell chi square matrix of all cells. The cells contain the number of cells with this combination.
    Returned: matrix[i][j]: Chi-squared of map A with class i and class j in map B.
    Returned: ncount: total number of cells which are included in this analyses.
    '''
    if (mapA.length != mapB.length):
        raise MyError("COMPARE_RASTERS: To compare two rasters, length should be equal.")

    # Dimension of matrix.
    #min_num_rows = iround(mapA.minimum())
    #min_num_cols = iround(mapB.minimum())
    #Nrows = iround(mapA.maximum()) - min_num_rows + 1
    #Ncols = iround(mapB.maximum()) - min_num_cols + 1

    # Create an matrix for all these values.
    #matrix = []
    #for irow in range(Nrows):
    #    matrix.append(Ncols*[0.0])

    # Fill the matrix with the frequency values.
    ncount = 0
    count_dict = {}
    for icell in range(mapA.length):
        valA = mapA.get_data(icell)
        valB = mapB.get_data(icell)
        if (valA == None):
            # Nodata is excluded.
            continue
        if (valB == None):
            # Nodata is excluded.
            continue

        valA = iround(valA)
        valB = iround(valB)
        ncount += 1
        # Add this cell to the matrix.
        #matrix[valA-min_num_rows][valB-min_num_cols] += 1
        try:
           count_dict[str(valA)+"_" + str(valB)] += 1
        except KeyError:
           count_dict[str(valA)+"_" + str(valB)] = 1
    #print("MATRIX: ", matrix)
    #print("count_dict: ", count_dict)
    
    # Get the dimensions of the matrix by examing the dictionary.
    row_dict = {}
    col_dict = {}
    Nrows = 0
    Ncols = 0
    for key in count_dict:
        irow,icol = key.split("_")
        try:
            qq = row_dict[iround(irow)]
        except KeyError:
            row_dict[iround(irow)] = Nrows
            Nrows += 1            

        try:
            qq = col_dict[iround(icol)]
        except KeyError:
            col_dict[iround(icol)] = Ncols
            Ncols += 1  
    
    # Sort the keys in the matrix, so the one to one are on the diagonal.
    rowlist = sorted(row_dict)
    collist = sorted(col_dict)
    # Reassign the row and columnnumbers
    count = 0
    for item in rowlist:
        row_dict[item] = count
        count += 1
    count = 0
    for item in collist:
        col_dict[item] = count
        count += 1
        
    #print("COL_DICT: ", col_dict)
    #print("ROW_DICT: ", row_dict)
    # Make matrix
    # Create an matrix for all these values.
    matrix = []
    for irow in range(Nrows):
        matrix.append(Ncols*[0.0])
    #print "Nrows,Ncols",Nrows,Ncols
    for key in count_dict:
        irow,icol = key.split("_")
        irow_mat = row_dict[iround(irow)]
        icol_mat = col_dict[iround(icol)]
        #print key, irow_mat,icol_mat
        matrix[irow_mat][icol_mat] = count_dict[key]

    # Calculate the sum of the columns and the rows.
    P_iT = Nrows * [0.0]
    P_Ti = Ncols * [0.0]
    for irow in range(Nrows):
        for icol in range(Ncols):
            P_iT[irow] += matrix[irow][icol]
            P_Ti[icol] += matrix[irow][icol]

    # Write matrix to output file when asked for:
    if (matrixfile != None):
        sep = ";"
        fp = open(matrixfile,"w")
        # Write the column headers
        fp.write(sep)
        for icol in range(Ncols):
            fp.write(str(collist[icol]) + sep)
        fp.write(str("Total") + "\n")
        
        # Write each row to the output file
        for irow in range(Nrows):
            fp.write(str(rowlist[irow]) + sep)
            for icol in range(Ncols):
                fp.write(str(matrix[irow][icol]) + sep)
            fp.write(str(P_iT[irow]) + "\n")
        
        # Write the sum of the columns
        fp.write("Total" + sep)
        for icol in range(Ncols):
            fp.write(str(P_Ti[icol]) + sep)
        fp.write(str(ncount) + "\n")        
        


    # Based on this matrix the expected distribution for p_i,j (when they are not dependent on each other).
    # Ep_i,j = p_.,j * p_i,. /ncount with p_.,j is the sum of column j. and p_i,. the sum of row i.

    # The Chi-squared method is defined as:
    # chi2[i][j] = (matrix[i][j] - Ep_i,j)^2/Ep_i,j

    # Calculate the Chi-square matrix.
    for irow in range(Nrows):
        for icol in range(Ncols):
            Pexpect = P_iT[irow] * P_Ti[icol]/float(ncount)
            diff = matrix[irow][icol] - Pexpect
            try:
                matrix[irow][icol] = (diff * diff)/Pexpect
            except ZeroDivisionError:
                 matrix[irow][icol] = 0.0
            #if (irow == icol):
            #    print irow, ncount,P_iT[irow] ,P_Ti[icol],diff,Pexpect,matrix[irow][icol]
    #print("MATRIX: ", matrix)
    return matrix,ncount

def cramer_v(matrix,ncount):
    '''
    Calculate the Cramer V statistics for these matrix.
    Cramer's V2 = chi2 / n * (min(r , k) -1) =
    # r: number of rows
    # k: number of columns
    # n: ncount, or number of cells that are included in this analyses.
    # chi2: sum of all chi_2 in the matrix.
    # V2: square of Cramer V.
    '''
    
    chi2 = 0.0
    for irow in range(len(matrix)):
        for icol in range(len(matrix[0])):
            chi2 += matrix[irow][icol]
    try:
        V2 = chi2/(ncount * (min(len(matrix),len(matrix[0]))-1))
    except ZeroDivisionError:
        raise MyError("Cramer V method is not possible for dimension smaller than 2.")
    return math.sqrt(V2)

def calculate_cramer_v(mapA,mapB,matrixfile=None,weighmap=None):
    '''
    Return the level of agreement between the two maps (cramer_V).
    '''
    # Calculate the Chi-square matrix of these two maps.
    matrix,ncount = calculate_chi2_matrix(mapA,mapB,matrixfile=matrixfile,weighmap=weighmap)
    #print("Chi-square matrix: ",matrix)

    # Calculate the Cramer V value
    cramer_value = cramer_v(matrix,ncount)

    return cramer_value
    

def calculate_matrix(mapA,mapB,weighmap=None):
    '''
    Calculate the cell-by-cell matrix of all cells. The cells contain the number of cells with this combination.
    Returned: matrix[i][j]: fraction of map A with class i that have class j in map B.
    Returned: the minimal value of all classes. 
    '''
    if (mapA.length != mapB.length):
        raise MyError("CALCULATE_MATRIX: To compare two rasters, length should be equal.")

    if (weighmap != None):
        if (mapA.length != weighmap.length):
            raise MyError("CALCULATE_MATRIX: Length of weighing map should be equal to the input maps.")

    # Find minimum and maximum values in both maps.
    max_num = iround(max(mapA.maximum(),mapB.maximum()))
    min_num = iround(min(mapA.minimum(),mapB.minimum()))

    # Dimension of matrix.
    N = max_num - min_num + 1

    # Create an matrix for all these values.
    matrix = []
    for irow in range(N):
        matrix.append(N*[0.0])

    # Fill the matrix
    ncount = 0
    for icell in range(mapA.length):
        valA = mapA.get_data(icell)
        valB = mapB.get_data(icell)
        if (valA == None):
            # Nodata is excluded.
            continue
        if (valB == None):
            # Nodata is excluded.
            continue
        if (weighmap != None):
            valweigh = weighmap.get_data(icell)
            if (valweigh == None):
                # Nodata is excluded, also in weighing map.
                print("CALCULATE_MATRIX: Nodata found in weighing map, which can influence the comparison.")
                continue
            else:
                # Check whether the weighing is positive everywhere.
                if (valweigh < 0.):
                    raise MyError("CALCULATE_MATRIX: Negaive values found in weighing map, which influence the comparison.")
        else:
            valweigh = 1

        valA = iround(valA)
        valB = iround(valB)
        ncount += valweigh
        # Add this cell to the matrix.
        matrix[valA-min_num][valB-min_num] += valweigh

    # Make a fraction of the matrix.
    for irow in range(N):
        for icol in range(N):
            matrix[irow][icol] /= float(ncount)

    return matrix,min_num

def fr_agreement(matrix):
    '''
    Calculate the fraction of agreement P(A).
    '''
    P_A = 0.0
    for irow in range(len(matrix)):
        P_A += matrix[irow][irow]
    return P_A

def fr_expected_agreement(matrix):
    '''
    Calculate the expected fraction of agreement P(E) and the maximal fraction of Agreement subject to observed distribution.
    '''

    # Calculate the sum of the columns and the rows.
    P_iT = len(matrix) * [0.0]
    P_Ti = len(matrix) * [0.0]
    for irow in range(len(matrix)):
        for icol in range(len(matrix)):
            P_iT[irow] += matrix[irow][icol]
            P_Ti[icol] += matrix[irow][icol]
    # print("P_iT", P_iT)
    # print("P_Ti", P_Ti)
    P_E = 0.0
    P_max = 0.0
    for irow in range(len(matrix)):
        P_E += P_iT[irow]*P_Ti[irow]
        # print("P_E_irow ", irow,P_iT[irow]*P_Ti[irow])
        P_max += min(P_iT[irow],P_Ti[irow])
    # print("P_E", P_E)
    # print("P_max", P_max)
    return P_E,P_max

def kappa_stat(P_A,P_E):
    '''
    Calculate the kappa statistics (level of agreement).
    Note: kappa_stat is equal to multiplication of kappa_location and kappa_histo.
    '''
    try:
        return (P_A - P_E)/(1.0 - P_E)
    except:
        return 0.0

def kappa_location(P_A,P_E,P_max):
    '''
    Calculate the kappa-location statistics (level of agreement).
    Kappa-location is a measure for the similarity of spatial allocation of categories
    for the two maps compared.
    '''
    try:
        return (P_A - P_E)/(P_max - P_E)
    except:
        return 0.0


def kappa_histogram(P_E,P_max):
    '''
    Calculate the kappa-histogram statistics (level of agreement).
    Kappa_histogram is a measure of the quantitive similarity of two maps. 
    '''
    try:
        return (P_max - P_E)/(1.0 - P_E)
    except:
        return 0.0

def calculate_statistics(matrix):
    '''
    Calculate the statistics based on the matrix provided.
    '''

    # Calculate the fraction of agreement.
    P_A = fr_agreement(matrix)

    # Calculate P_E and P_max.
    P_E,P_max = fr_expected_agreement(matrix)

    # Calculate kappa.
    kappa = kappa_stat(P_A,P_E)

    # Calculate kappa_location.
    kappa_loc = kappa_location(P_A,P_E,P_max)

    # Calculate kappa_histo (of histograms).
    kappa_histo = kappa_histogram(P_E,P_max)

    return kappa,kappa_loc,kappa_histo

def compare_rasters(mapA,mapB,weighmap=None):
    '''
    Return the level of agreement between the two maps (kappa,kapp_location and kappa_histo).
    Also weighing can play a role here.
    '''
    # Calculate the cell-by-cell difference.
    matrix,min_num = calculate_matrix(mapA,mapB,weighmap=weighmap)
    # print("MATRIX: ",matrix)

    # Calculate the statistics.
    kappa,kappa_loc,kappa_histo = calculate_statistics(matrix)

    return kappa,kappa_loc,kappa_histo

def compare_rasters_per_class(mapA,mapB,weighmap=None):
    '''
    Return the level of agreement between the two maps (kappa,kapp_location and kappa_histo).
    '''
    # Calculate the cell-by-cell difference.
    matrix,min_num = calculate_matrix(mapA,mapB,weighmap=weighmap)
    # print("MATRIX: ",matrix)

    # Create output dictionary
    output_dict = {}

    # Calculate the statistics for each class
    for iclass in range(len(matrix)):
        # print("ICLASS: ", iclass + min_num)
        # Create a new matrix with the class and not the class (binair).
        # Element zero represent the class choosen, element one is the rest.
        matrix_class = []
        for irow in range(2):
            matrix_class.append(2*[0.0])

        # Fill the matrix
        for irow in range(len(matrix)):
            for icol in range(len(matrix)):
                if (irow != iclass and icol != iclass):
                    matrix_class[1][1] += matrix[irow][icol]
                elif (irow == iclass and icol != iclass):
                    matrix_class[0][1] += matrix[irow][icol]
                elif (irow != iclass and icol == iclass):
                    matrix_class[1][0] += matrix[irow][icol]
                else:
                    matrix_class[0][0] += matrix[irow][icol]                

        # Calculate the statistics.
        kappa,kappa_loc,kappa_histo = calculate_statistics(matrix_class)

        output_dict[min_num + iclass] = [kappa,kappa_loc,kappa_histo]

    return output_dict


# Source code to test this functions.
           
if (__name__ == "__main__"):

    import ascraster

    mapA = ascraster.Asciigrid(ncols = 4,nrows = 4, xllcorner = 0.0, yllcorner = 0.0, \
                               cellsize = 1.0, nodata_value = -1, numtype = int)
    mapB = ascraster.duplicategrid(mapA)
    mapA.values = [-1,-1,-1,-1, 1, 1, 1, 1,\
                    2, 2, 2, 2, 3, 3, 3, 3]
    mapB.values = [-1,-1,-1,-1, 1, 1, 1, 1,\
                    3, 2, 2, 2, 3, 3, 3, 3]
    print("Raster are equal.")
    print(compare_rasters(mapA,mapA))
    print("Cramer: ",calculate_cramer_v(mapA,mapA))
    print("Raster are equal. PER CLASS")
    print(compare_rasters_per_class(mapA,mapA))
    print("Raster are equal except one cell (3 versus 2).")
    print(compare_rasters(mapA,mapB))
    print(compare_rasters_per_class(mapA,mapB))
    print("Cramer: ",calculate_cramer_v(mapA,mapA))
    mapA.values = [-1,-1,-1,-1, 1, 1, 1, 1,\
                    2, 3, 2, 2, 3, 3, 3, 3]
    mapB.values = [-1,-1,-1,-1, 1, 1, 1, 1,\
                    3, 2, 2, 2, 3, 3, 3, 3]
    print("Raster are equal except two cells (3 versus 2).")
    print(compare_rasters(mapA,mapB))
    print(compare_rasters_per_class(mapA,mapB))
    print("Cramer: ",calculate_cramer_v(mapA,mapA))
    mapA.values = [ 2, 3, 2, 2, 3, 3, 3, 3,\
                   -1,-1,-1,-1, 1, 1, 1, 1]
    mapB.values = [-1,-1,-1,-1, 1, 1, 1, 1,\
                    3, 2, 2, 2, 3, 3, 3, 3]
    print("Raster cells are all not equal, but the same distribution.")
    print(compare_rasters(mapA,mapB))
    print(compare_rasters_per_class(mapA,mapB))
    print("Cramer: ",calculate_cramer_v(mapA,mapB))
    mapA.values = [ 2, 2, 2, 2, 3, 3, 3, 3,\
                   -1,-1,-1,-1, 1, 1, 1, 1]
    mapB.values = [ 4, 4, 4, 4, 5, 5, 5, 5,\
                    6, 6, 6, 6, 7, 7, 7, 7]
    print("Raster cells are all not equal, and the different distribution.")
    print(compare_rasters(mapA,mapB))
    print(compare_rasters_per_class(mapA,mapB))
    print("Cramer: ",calculate_cramer_v(mapA,mapB))
    mapA.values = [ 2, 2, 2, 2, 2, 2, 2, 2,\
                    2, 2, 2, 2, 2, 2, 2, 2]
    mapB.values = [ 2, 2, 2, 2, 2, 2, 2, 2,\
                    1, 1, 2, 2, 2, 2, 2, 2]
    print("Raster cells are partly equal, and the different distribution, only two values.")
    matrix,min_num = calculate_matrix(mapA,mapB)
    print("MATRIX: ", matrix)
    print(compare_rasters(mapA,mapB))
    try:
        print("Cramer: ",calculate_cramer_v(mapA,mapB))
    except MyError:
        print("This is correct. No problem. Size of matrix is one.")

    mapA = ascraster.Asciigrid(ncols = 2,nrows = 2, xllcorner = 0.0, yllcorner = 0.0, \
                               cellsize = 1.0, nodata_value = -1, numtype = int)
    mapB = ascraster.duplicategrid(mapA)
    mapA.values = [1, 1, 1, 1]
    mapB.values = [1, 1, 1, 1]
    print("SAME_4_1: ",compare_rasters(mapA,mapB))
    mapA.values = [1, 2, 1, 1]
    mapB.values = [1, 2, 1, 1]
    print("SAME_4_12: ",compare_rasters(mapA,mapB))
    mapA.values = [1, 2, 3, 1]
    mapB.values = [1, 2, 2, 1]
    print("SAME_4_123: ",compare_rasters(mapA,mapB))
    mapA.values = [1, 2, 3, 4]
    mapB.values = [1, 2, 3, 4]
    print("SAME_4_1234: ",compare_rasters(mapA,mapB))
    mapA.values = [1, 2, 3, 4]
    mapB.values = [1, 2, 3, 3]
    print("SAME_4_1234_3: ",compare_rasters(mapA,mapB))
    mapA = ascraster.Asciigrid(ncols = 6,nrows = 6, xllcorner = 0.0, yllcorner = 0.0, \
                               cellsize = 1.0, nodata_value = -1, numtype = int)
    mapB = ascraster.duplicategrid(mapA)
    mapA.values = [-1,-1, 1,-1,-1,-1, \
                   -1, 2, 2, 2,-1,-1, \
                    3, 3, 3, 3, 3, 3, \
                    4, 4, 4, 4, 4, 4, \
                   -1,-1, 5, 5, 5,-1, \
                   -1,-1,-1, 6,-1,-1]
    mapB.values = [-1,-1, 1,-1,-1,-1, \
                   -1, 2, 2, 2,-1,-1, \
                    3, 3, 3, 3, 3, 3, \
                    4, 4, 5, 5, 4, 4, \
                   -1,-1, 5, 5, 5,-1, \
                   -1,-1,-1, 6,-1,-1]
    print("comparison 6by6: ",compare_rasters(mapA,mapB))
    matrix,min_num = calculate_matrix(mapA,mapB)
    print("MATRIX: ", matrix)
    mapA.write_ascii_file("mapA_6_6.asc")
    mapB.write_ascii_file("mapB_6_6.asc")
    write_diffmap(mapA,mapB,"map_diffAB.asc")
    print("comparison 6by6 per class: ",compare_rasters_per_class(mapA,mapB))

    # Check whether method is dependent on nodata. Results should be the same!
    mapC = ascraster.duplicategrid(mapA)
    mapD = ascraster.duplicategrid(mapB)
    mapE = ascraster.duplicategrid(mapA)
    weighmap = ascraster.duplicategrid(mapA)
    weighmap.values = 36 *[1.]

    # MapC and MapD, change one nodata cell into value.
    mapC.set_data(0,1)
    mapD.set_data(6,1)
    # MapE, change nodata value.
    mapE.nodata_value = 10
    print("comparison 6by6 A_B: ",compare_rasters(mapA,mapB))
    print("comparison 6by6 A_B equal weighing: ",compare_rasters(mapA,mapB,weighmap=weighmap))
    for i in range(36):
        weighmap.set_data(i,i+1.)
    print("comparison 6by6 A_B non-equal weighing: ",compare_rasters(mapA,mapB,weighmap=weighmap))
    print("comparison 6by6 C_B: ",compare_rasters(mapC,mapB))
    print("comparison 6by6 A_D: ",compare_rasters(mapA,mapD))
    print("comparison 6by6 C_D: ",compare_rasters(mapC,mapD))
    print("comparison 6by6 E_B: ",compare_rasters(mapE,mapB))
    print("comparison 6by6 E_D: ",compare_rasters(mapE,mapD))
    print("Cramer: ",calculate_cramer_v(mapA,mapB))

    # Calculation of the Chi-square matrix.
    # Example from the https://en.wikipedia.org/wiki/Chi-squared_test
    # I just fill the maps with the values from the example.
    mapA = ascraster.Asciigrid(ncols = 65,nrows = 10, xllcorner = 0.0, yllcorner = 0.0, \
                               cellsize = 1.0, nodata_value = -1, numtype = int)
    mapB = ascraster.duplicategrid(mapA)

    mapA.values = 90*[5]   + 30*[6]  + 30*[7] + \
                  60*[5]   + 50*[6]  + 40*[7] + \
                  104*[5]  + 51*[6]  + 45*[7] + \
                  95*[5]   + 20*[6]  + 35*[7]
    mapB.values = 150*[1]  + 150*[2] + 200*[3] + 150*[5]
    matrix,ncount = calculate_chi2_matrix(mapA,mapB)
    print("Chi-square matrix: ",matrix)
    print("Cramer: ",cramer_v(matrix,ncount))

    # Calculation of the Chi-square matrix.
    # Example from the http://www.let.leidenuniv.nl/history/RES/stat/html/les9.html (Politics and voting behaviour)
    # I just fill the maps with the values from the example.
    mapA = ascraster.Asciigrid(ncols = 77,nrows = 5, xllcorner = 0.0, yllcorner = 0.0, \
                               cellsize = 1.0, nodata_value = -1, numtype = int)
    mapB = ascraster.duplicategrid(mapA)

    mapA.values =  94*[3]  + 56*[4]  + \
                  100*[3]  + 135*[4]
    mapB.values = 150*[1]  + 235*[3]
    matrix,ncount = calculate_chi2_matrix(mapA,mapB)
    print("Chi-square matrix: ",matrix)
    print("Cramer: ",cramer_v(matrix,ncount))
    print("Cramer: ",calculate_cramer_v(mapA,mapB))
    weighmap = ascraster.duplicategrid(mapA)
    weighmap.values = 385*[3.]
    print("Cramer equal weighing: ",calculate_cramer_v(mapB,mapA,weighmap=weighmap))
    weighmap.values = 150*[2.5] + 235*[1.4]   
    print("Cramer non-equal weighing: ",calculate_cramer_v(mapB,mapA,weighmap=weighmap))
    print("Cramer equal matrix: ",calculate_cramer_v(mapB,mapB))
    # To produce an warning on the weighing map (nodata in weighing map where the system does not expect it.)
    weighmap.values = 150*[2.5] + 234*[1.4] + [-1.]
    print("There should be a warning on the occurence of nodata in the weighmap.")  
    print("Cramer non-equal weighing: ",calculate_cramer_v(mapB,mapA,weighmap=weighmap))
    # To produce an error on the weighing map (negative value in the weighing map)
    weighmap.values = 150*[2.5] + 234*[1.4] + [-2.]
    try:  
        val = calculate_cramer_v(mapB,mapA,weighmap=weighmap)
        print("Cramer non-equal weighing: ",val)
    except: 
        print("Pass on test of error message.")
