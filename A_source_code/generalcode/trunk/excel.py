# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/excel.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

#! /usr/bin/env python
#encoding:UTF-8
'''
Created on 23 sep 2010

@author: warrinka
'''
import xlwt
import xlrd


class Worksheet:
#    def __init__(self, workbook, name="Sheet1"):
#        self.workbook = workbook
#        self.sheet = self.workbook.add_sheet(name, cell_overwrite_ok=True)
#        self.defaultstyle = Style()
    def __init__(self, workbook, name=None):
        self.workbook = workbook
        if (name == None):
            self.sheet = {}
        else:
            self.sheet = self.workbook.add_sheet(name, cell_overwrite_ok=True)
        self.defaultstyle = Style()
    
    def get_cell_value(self):
        pass
    
    def get_cell_style(self):
        pass
    
    def change_col_width(self, col_index, width):
        self.sheet.col(col_index).width = width
    
    def write(self, location, value, style=None):
#        label=""
#        style=Style.default_style
        if not style:
            style = self.defaultstyle.style
        r, c = location
        try:
            self.sheet.write(r, c, value, style)
        except UnicodeDecodeError as msg:
            print(msg, value)
    
    def place_picture(self, location, filename):
        r, c = location
        self.sheet.insert_bitmap(filename, r, c) #, x, y, scale_x, scale_y)


class Style:
    # font name, type, size, bold
    def __init__(self, fname='Arial', fbold=False, fitalic=False, fsize=10, borders=None, wrap=False,\
                 num_format_str='general'):
        self.style = xlwt.XFStyle()
        font = xlwt.Font()
        font.name = fname
        font.bold = fbold
        font.italic = fitalic
        font.height = self._getheight(fsize)
        self.style.font = font
        self.style.num_format_str = num_format_str
        if borders == 'bottom': self.style.borders.bottom = 2
        elif borders == 'top': self.style.borders.top = 2
        self.style.alignment.wrap = wrap
        
    def _getheight(self, size):
        heigth = size*20
        return heigth

class Workbook:
    def __init__(self, save_path, worksheet = None):
        self.save_path = save_path
        self.workbook = xlwt.Workbook('latin1')
        self.sheets = {}
        if (worksheet != None):
            self.add_sheet(worksheet)
    
    def add_sheet(self, name):
        #TODO: foutmelding bij al bestaande sheet
        sheet = Worksheet(self.workbook, name)
        self.sheets[name] = sheet
        
    def save_workbook(self):
        self.workbook.save(self.save_path)


def get_range(filepath, worksheetname, datarange):
    # TODO: open file, sheet en get data
    workbook = xlrd.open_workbook(filepath)
    worksheet = workbook.sheet_by_name(worksheetname)
    if datarange[0] == datarange[1]:
        data = worksheet.cell_value(datarange[0][1], datarange[0][0])
    elif datarange[0][0] == datarange[1][0]:
        data = worksheet.col_values(datarange[0][0], datarange[0][1], datarange[1][1]+1)
    elif datarange[0][1] == datarange[1][1]:
        data = worksheet.row_values(datarange[0][1], datarange[0][0], datarange[1][0]+1)
    return data

def coord2tuple(coord):
    col, row = '', ''
    alpha = 'abcdefghijklmnopqrstuvwxyz'.upper()
    pairs = [''.join((x,y)) for x in alpha for y in [''] + [z for z in alpha]]
    pairs = sorted(pairs, key=len)
    coord = coord.upper()
    for c in coord:
        if c in alpha:
            col += c
        else:
            row += c
    return (pairs.index(col), int(row)-1)


if __name__ == '__main__':
    pass
