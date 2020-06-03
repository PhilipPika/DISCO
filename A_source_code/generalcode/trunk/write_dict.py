# ******************************************************
## Revision "$LastChangedDate: 2018-06-01 15:05:44 +0200 (Fri, 01 Jun 2018) $"
## Date "$LastChangedRevision: 1 $"
## Author "$LastChangedBy: arthurbeusen $"
## URL "$HeadURL: https://pbl.sliksvn.com/generalcode/trunk/write_dict.py $"
## Copyright 2017, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
 
import sys
from error import *
from my_sys import my_print

def write_dict(fp_in=None,filename=None,dic=None,headerkey="",sep=" ",lkey=1,ldebug=1,fatal=0):
    '''
    This function writes a dictionairy to an output file
    Arguments:
    filename: name of output file
    dic: name of the dictionairy. 
         A element of the dictionary may contain list, own type of python
         (like number or string) or a own class. In the last case there must 
         be a write function!
    headerkey: the name of the key of the header. When headerkey is empty, no header will be printed.
    sep: Seperator between the fields in the output file
    lkey: 0 (no key is written) or 1 (key is included in the output file).
    ldebug: 0 (no debug information to screen), 1 (debug information is written to screen).
    fatal: 1 (a raise is given in case of an error), 0 (a warning is given in case of an error, programme continues).
           In case of serious error, the programme will stop, independent of fatal!
    '''
    if (dic == None): raise MyError("Dictionary must be given.")
    try:
        if (ldebug == 1): print("Start of function write_dict.")
        try:
            if (fp_in == None):
                fp = open(filename,"w")
                if (ldebug == 1): print("Output file " + filename + " is openend.")
            else:
                fp = fp_in
        except:
            raise MyError("Opening error of file: " + filename)

        # Teller gives the number of lines that will be written to the output file.
        teller = 0 
        # First write the header to the output file.
        if (not headerkey == ""):
            # Check whether the giver headerkey exists.
            if headerkey not in dic:
                my_print("WRITE_DICT: Given headerkey: " + headerkey + " is not found.")
            else:
                if (type(dic[headerkey]) == type([])):
                    if (ldebug == 1): print("Dictionary element is of type list.")
                else:
                    if (ldebug == 1): print("Dictionary element is of type which has a write function.")
                
                # Write key element
                if (lkey == 1): fp.write(headerkey + sep)
                # Write rest of the header
                if (type(dic[headerkey]) == type([])):
                    # In case of list
                    for veld in range(0,len(dic[headerkey])-1):
                        fp.write(str(dic[headerkey][veld]) + sep)
                    fp.write(str(dic[headerkey][len(dic[headerkey])-1]) + "\n")
                else:
                    # In case of non list type.
                    fp.write(str(dic[headerkey])+"\n")
                if (ldebug == 1): 
                    print("The header of the file is written to file.")
                
                # Add number of lines written to output file.    
                teller +=1
                
        # Sort the dictionary so the output is always the same.
        sortkey = []
        for key in list(dic.keys()):
            # Do not take the header in the sorted list
            if not (str(key).upper() == str(headerkey).upper()):
                sortkey.append(key)
        sortkey.sort()
        
        if (ldebug): print("Dictionary is sorted") 

        # Write all elements to output file:
        for item in range(0,len(sortkey)):
            key = sortkey[item]
            # Write key element
            if (lkey == 1): fp.write(str(key) + sep)
            # Write rest of the element
            if (type(dic[key]) == type([])):
                # Dictionary met lists.
                for veld in range(0,len(dic[key])-1):
                    fp.write(str(dic[key][veld]) + sep)
                fp.write(str(dic[key][len(dic[key])-1]) + "\n")
            else:
                # In case of non list type.
                dic[key].write(fp,sep)
            teller +=1
        # Close output file
        if (fp_in == None): 
            fp.close()
            del fp
        
        if (ldebug == 1):
            if (filename == None): 
                print("There are "+str(teller)+" lines written to file "+ fp.name)
            else:    
                print("There are "+str(teller)+" lines written to file "+ filename)
    except MyError as val:
        val.write()
        raise MyError()       
    except:
        print("WRITE_DICT: Something goes wrong with error type: ", sys.exc_info()[0])
        if (filename == None): 
            if not fp.closed:
                print("There are "+str(teller)+" lines written to file "+ fp.name)
                filename = fp.name
            else:    
                print("There are "+str(teller)+" lines written to file.")
                filename = "not_known"
                
        else:    
            print("There are "+str(teller)+" lines written to file "+ filename)
        if not fp.closed: fp.close()
        del fp
        raise MyError("ERROR in WRITE_DICT for file: "+filename)

