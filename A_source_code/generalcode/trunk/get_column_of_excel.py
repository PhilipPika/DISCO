## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/get_column_of_excel.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
## ---------------------------------------------------------------------------
## -
## -    Name:    get_column_of_excel.py
## -    Author:  Arthur Beusen PBL/IMP
## -    Date:    June 4 2009
## -    Purpose: Gets a column out of an excel file and puts this in a list.
## -             Headers are places on the first line of the sheet.
## -
## ---------------------------------------------------------------------------

import sys
import os
from error import *
import xlrd
# from readexcel import *

def get_column_of_excel(filename,worksheetname,colname,ldebug=0):
    '''
    Puts one column of a worksheet of a workbook  into a list.
       Arguments: filename: Name of excel file.
                  worksheetname: is the name of the worksheet.
                  colname: name of the column which values are put in the output list.
       Return value: List with values from the column colname of the worksheet worksheetname.
    '''

    # Make default output list
    result_list = []

    # Check whether excel file exist:
    if not os.path.isfile(filename):
        raise MyError("Excel file " + filename + " does not exist.")
        
    if (ldebug): print("Start reading excel file " + filename)

    # Open workbook:
    xls = xlrd.open_workbook(filename=filename,on_demand=True) 

    # Check whether worksheet worksheetname exist:
    sheet_names=xls.sheet_names()
    for item in range(0,xls.nsheets):
        if (str(sheet_names[item]) == worksheetname):
            break
    else:
        print("Sheetnames available in the excel file:")
        print(xls.sheet_names())
        raise MyError("Excel file " + filename + " has no worksheet "+worksheetname)
      
    # Load excel sheet in memory
    xlsTuple = xls.sheet_by_name(worksheetname)    

    # Number of rows available in the sheet.
    rows = xlsTuple.nrows
    if (ldebug): print("Number of rows: ", rows)
    if (rows < 2):
        raise MyError("Number of rows is too small in worksheet "+worksheetname)

    # Number of columns available in the sheet.
    cols = xlsTuple.ncols
    if (ldebug): print("Number of columns: ", cols)
    if (cols < 2): 
        raise MyError("Number of columns is too small in worksheet "+worksheetname)

    # Find column with the right header.
    col_id = -1
    COLNAME = str(colname).upper()
    header = xlsTuple.row_values(0)
    for column in range(0, len(header)):
        if (str(header[column]).upper() == COLNAME):
            col_id = column
        
    if (col_id == -1):
        print("Columnnames available in the worksheet:")
        for column in range(0, len(header)):
            print(header[column])     
        raise MyError("Excel file: " + filename + " Worksheet "+worksheetname+" has no column "+colname)        

    # Put data into a list
    val = xlsTuple.col_values(col_id)
    for row in range(1,rows):
        result_list.append(val[row])

    if (ldebug): print("End reading excel file " + filename + " with worksheet " + worksheetname + " for column " + colname)

    # Unload sheet from memory
    xls.unload_sheet(worksheetname)
    del xls
    del xlsTuple
        
    return result_list

