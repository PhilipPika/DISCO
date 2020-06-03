# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/general_class.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import my_sys
from error import *

class General:
    """
    Empty data class.
    This class can be used as a template to create a data class with all kinds of attributes.
    For example: It can be used to read a matrix of data in a file, where the column headers are used as attributes
    with the function read_general_file.
    """
    def __init__(self):
        # List of attributes
        self.names = []
    
    def get_attrib(self):
        return self.names

    def copy(self):
        out = General()
        for name in self.names:
            out.names.append(name)
            setattr(out,name,getattr(self,name))
        return out

    def merge(self,C):       
        for name in C.names:
            self.add_item(name,C.get_val(name))

    def copy_into(self,C):
        '''
        Copy all attributes of this class to an existing class, when the attribute name is common in both classes.
        '''

        for name in self.names:
            try:
                qq = getattr(C,name)
                setattr(C,name,getattr(self,name))
            except AttributeError:
                # Attribute does not exist
                pass
        return C

    def add_item(self,attrib,val):
        # Replace all spaces in the attrib into underscores
        attrib = attrib.replace(' ', '_').replace('\"', '')
        # Add attrib to the list of attributes
        self.names.append(attrib)
        # Add attrib to class a set a value
        setattr(self,attrib,val)

    def set_val(self,attrib,val):
        # Replace all spaces in the attrib into underscores
        attrib1 = attrib.replace(' ', '_').replace('\"', '')
        try:
            # Check whether attribute is new
            qq = getattr(self,attrib1)
            setattr(self,attrib1,val)
        except AttributeError:
            raise MyError("New attribute is added to the General class. This is not allowed.",\
                          "Use function add_item.",\
                          "Attribute that is new: " + str(attrib1))        

    def get_val(self,attrib,val=None):
        '''
        Get value for a given attribute. When the attribute is NoData None is returned.
        When val is not None, val will be returned in stead of None (in case of NoData) 
        @param attrib: String with attribute name
        @param val: Default value in stead of None in case of Nodata
        
        '''
        # Replace all spaces in the attrib into underscores
        attrib1 = attrib.replace(' ', '_').replace('\"', '')
        try:
            # Check whether attribute is element of this class
            returnval = getattr(self,attrib1)
            if (returnval != None):
                return returnval
            else:
                if (val == None):
                    return returnval
                else:
                    return val
 
        except AttributeError:
            # Attribute is unkown
            raise MyError("General class does not contain attribute " + str(attrib1))

    def write(self,fp,sep=";",lheader=False,NoneValue=""):
        '''
        Write General class to output file.
        @fp: file pointer to a opened file
        '''       
        # Write header
        if (lheader):
            fp.write(str(self.names[0])+sep)
            for item in range(1,len(self.names)-1):
                fp.write(str(self.names[item]) + sep)
            fp.write(str(self.names[-1]) + "\n")

        # Write values to output file. None values are written as empty strings.
        fp.write(str(self.get_val(self.names[0],val=NoneValue))+sep)
        for item in range(1,len(self.names)-1):
            fp.write(str(self.get_val(self.names[item],val=NoneValue)) + sep)
        fp.write(str(self.get_val(self.names[-1],val=NoneValue)) + "\n")

def read_general_file(filename,sep=";",key=None,out_type="dict",ldebug=False):
    '''
    Reads a table file.
    First line is the header which are used to name the attributes.
    Returns a dictionary, with the key found in the first column when key=None. Otherwise the given key will be used.
    Each dictionary contains a common class, with attributes given in the header and the values given in the table.  
    '''

    # Read input file.
    if (ldebug):
        print("Reading input file: ",filename)
    dat = my_sys.my_readfile(filename)

    # Find header
    headline = -1
    for linenum in range(len(dat)):
        # Look whether the first column contains a ! or a #. When it does, then it is comment. Otherwise it is the header.
        fields = dat[linenum].split(sep)
        if (len(fields) == 0):
            continue
        if ((len(fields) == 1) and (len(fields[0].lstrip()) == 0)):
            continue
        if not ("!" in fields[0].lstrip() or "#" in fields[0].lstrip()):
            # Header found.
            headline = linenum
            number_of_fields = len(fields)
            # Put header in a list
            header = []
            for item in range(len(fields)):
                header.append(fields[item])
            break
            
    if (headline < 0):
        raise MyError("Input file " + filename + " is not correct.",\
                       "Header is not found.")

    # Store information from file in the General class
    if (out_type == 'dict'):
        out = {}
    elif(out_type == 'list'):
        out = []
    else:
        raise MyError("Reading input file " + filename + " goes wrong.",\
                      "Wrong value for argument out_type.")

    # Make a standard empty General class for this file.
    G = General()

    # Make a list of all attributes
    key_item = 1000000
    for item in range(0,len(header)):
        if (len(str(header[item]))!=0):
            if (type(out) == dict):
                if (key == None and item == 0):
                    # This is the key of the dictionary
                    key_item = item
                    continue
                else:
                    if (str(header[item]) == key):
                        key_item = item
                        continue
            G.add_item(header[item],None)

    if (type(out) == dict and key_item == 1000000):
        raise MyError("Error in reading file: " + filename,\
                      " Given key for dictionay is not found. ",\
                      " Available keys are: ",\
                      " ".join(G.names))

    for linenum in range(headline+1,len(dat)):
        lineinfo = G.copy()
       
        # Split line into fields
        fields = dat[linenum].split(sep)
        
        if (len(fields) == 0):
            # Empty line
            continue
        
        # Check whether we have the same amount of data fields as header fields.
        if (len(fields) != number_of_fields):
            # Try to fix the problem, that field separator is used when it i splaced between quotes.
            fields = merge_fields(fields)
            if (len(fields) != number_of_fields):
                print(fields)
                raise MyError("Input file " + filename + " is not correct. ",\
                           "Number of headers fields found: " + str(number_of_fields),\
                           ". Number of data fields found: " + str(len(fields)) + " on linenumber: " +str(linenum+1)+".")

        for item in range(0,number_of_fields):
            try:
                if (len(str(header[item]))==0):
                    # Header is empty, so no data. Skip this column.
                    continue
                elif (type(out) == dict and item == key_item):
                    dictkey =  fields[item]
                else:
                    lineinfo.set_val(header[item],fields[item])   
            except AttributeError:
                raise MyError("Input file " + filename + " is not correct.",\
                       " Error on linenumber: " + str(linenum+1) + " and column number: " + str(item+1)+".")

        if (type(out) == dict):
            # Add line to the output dictionary
            try:
                qq = out[dictkey]
                # Key does exist, so information is overwritten.
                raise MyError("READ_GENERAL_FILE: Double entry found in file: " + filename + ".",\
                              " Key which is double: " + str(dictkey))
            except KeyError:
                # Key does not exist.
                if (len(dictkey) > 0):
                    out[dictkey] = lineinfo
        elif (type(out) == list):
            # Add line to the output list
            out.append(lineinfo)

    # Return the list with General objects.
    return out


def merge_fields(fields):
    '''
    Try to merge fields which are accidental split. The case here is that a field start with a quote 
    but does not ends with a quote.
    '''
    delnum = []
    for ifield in range(len(fields)):                
        # Merge fields which are separated by accident.
        if (len(fields[ifield])==0):
            # Empty field
            continue
        elif (fields[ifield].count("\"")%2 == 0):
            # This field has begin and end quote, so skip
            continue

        if (fields[ifield][0] == "\"" and fields[ifield][-1] != "\""):
            # Look for the end of the string (the quotes).
            try:
                fields[ifield] += fields[ifield+1]
            except IndexError:
                continue
            delnum.append(ifield+1)
            ilocal = ifield + 1
            try:
                condit = (fields[ilocal][-1] != "\"")
            except IndexError:
                condit = True

            while (condit):
                ilocal += 1
                #print "SIT1 ",fields[0],ilocal,len(fields)
                try: 
                    fields[ifield] += fields[ilocal]
                except IndexError:
                    # In this case we have the last element of the list
                    break
                delnum.append(ilocal)          
                try:
                    condit = (fields[ilocal][-1] != "\"")
                except IndexError:
                    condit = True

        else:
             # We have a problem because the number of quotes in the fields in not even and begin and end has a quote.
             # So make the next field start of next search only when this field is not deleted yet.
             sp=":"
             if not (ifield in delnum):
                 try:
                     #print "SIT2 ",sp+str(fields[0])+sp+str(ifield+1)+sp+str(len(fields))+sp+str(fields[ifield])+sp+str(fields[ifield][0])+sp+str(fields[ifield][-1])+sp
                     fields[ifield+1] = fields[ifield] + fields[ifield+1]
                     delnum.append(ifield)
                 except IndexError:
                     pass
    
    # Delete all double fields.
    for i in range(len(delnum)-1,-1,-1):
        del fields[delnum[i]]
     
    return fields
