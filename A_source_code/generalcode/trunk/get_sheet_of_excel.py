# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/get_sheet_of_excel.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

## ---------------------------------------------------------------------------
## -
## -    Name:    get_sheet_of_excel.py
## -    Author:  Arthur Beusen PBL/IMP
## -    Date:    June 4 2009
## -    Purpose: Gets a sheet out of an excel file and puts it in a list.
## -
## ---------------------------------------------------------------------------


import sys
import os
from error import *
import xlrd

def get_sheet_of_excel(filename,worksheetname,ldebug=0):
    """Puts a worksheet of a workbook  into a list.
       Arguments: filename: Name of excel file.
                  worksheetname: is the name of the worksheet.
       Return value: List with of rows with a list of columns. """

    try:
        result_list = []

        # Check whether excel file exist:
        if not os.path.isfile(filename):
            raise MyError1("Excel file " + filename + " does not exist.")

        if (ldebug): print("Start reading excel file " + filename + " with worksheet " + worksheetname)

        # Open workbook:
        try:
            xls = xlrd.open_workbook(filename=filename,on_demand=True)
        except Exception as e:
            print(e)
        # Check whether worksheet worksheetname exist:
        sheet_names=xls.sheet_names()
        for item in range(0,xls.nsheets):
            if (str(sheet_names[item]) == worksheetname):
                break
        else:
            print("Sheetnames available in the excel file:")
            print(xls.sheet_names())
            raise MyError1("Excel file " + filename + " has no worksheet "+worksheetname)

        # Load excel sheet in memory
        xlsTuple = xls.sheet_by_name(worksheetname)

        # Number of rows available in the sheet.
        rows = xlsTuple.nrows
        if (ldebug): print("Number of rows: ", rows)
        if (rows < 2):
           raise MyError1("Number of rows is too small in worksheet "+worksheetname)


        # Number of columns available in the sheet.
        cols = xlsTuple.ncols
        if (ldebug): print("Number of columns: ", cols)
        if (cols < 2):
           raise MyError1("Number of columns is too small in worksheet "+worksheetname)


        # Put data into a list
        for row in range(0,rows):
            val = xlsTuple.row_values(row)
            result_list.append(val)

        if (ldebug): print("End reading excel file " + filename + " with worksheet " + worksheetname)

        # Unload sheet from memory
        xls.unload_sheet(worksheetname)
        del xls
        del xlsTuple

        return result_list

    except MyError0:
        raise MyError0()
    except MyError1 as val:
        val.write()
        raise MyError0()
    except:
        print("***** ERROR ******")
        print("An error has occurred in get_sheet_of_excel.py")
        print(sys.exc_info()[0])
        raise MyError0()
