# ******************************************************
## Revision "$LastChangedDate: 2018-11-20 11:45:12 +0100 (di, 20 nov 2018) $"
## Date "$LastChangedRevision: 8 $"
## Author "$LastChangedBy: laurianevilmin $"
## URL "$HeadURL: https://pbl.sliksvn.com/dgnm/core/general_func.py $"
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

# Import general python modules
import os
import sys
import traceback
import pickle
from scipy import *
import time
import reactions #LV 19-11-2018

# Import own general modules
try:
    from error import *
except ImportError:
    # Set the general path for the own python modules
    import general_path
    from error import *

def calc_o2_sat(temperature):
    '''
    Calculates O2 saturation concentration in mgO2.L-1, based on equation developed by
    American Public Health Authority (Cox, 2003 Science of the Total Environment)
    '''
    tempK = temperature + 273.15
    c0 = -139.34411
    c1 = 1.575701 * 1.0e5
    c2 = -6.642308 * 1.0e7
    c3 = 1.243800 * 1.0e10
    c4 = -8.621949 * 1.0e11
    ln_csat = c0 + c1 / tempK + c2 / math.pow(tempK,2.0) + c3 / math.pow(tempK,3.0) + c4 / math.pow(tempK,4.0)
    #print("csat  " + str(math.exp(ln_csat)))
    return math.exp(ln_csat)


def start_process(p,jobs,max_number_of_processes):
    '''
    P: proces handler.
    Start the current process and continue only when the number of processes running is smaller than
    the specified max_number_of_processes    
    '''
    # Start the process.
    p.start()

    # The following while loop checks the number of processes which are active.
    # When more than the defined number of cpu are active, it stops with starting new jobs until
    # there are cpu cores available.
    while (len(jobs) >= max_number_of_processes):
        # This sleep is needed, otherwise this while loop consumes all cpu time. 
        time.sleep(0.05)
        # Delete all jobs in the job list which are finished.
        for j in range(len(jobs)-1,-1,-1):
            if not (jobs[j].is_alive()):
                del jobs[j]


def get_number_of_cores(params,number_of_cells):
    '''
    Assign the number of cpu's needed for this amount of cells in the calculation packet.
    '''

    if (number_of_cells > 3 * params.minimal_number_of_cells):
       return 3
    elif (number_of_cells > 2 * params.minimal_number_of_cells):
       return 2
    else:
       return 1

def stryear(year):
    '''
    In order to avoid too long names in the case of monthly simulations...
    '''
    
    if (isinstance(year, float)):
        strtime = '%.3f' % round(year,3)
    else:
        strtime = str(year)
    return strtime

def path_conversion(path):
    '''
    When a path (direcory) is used in a system call (DOS or Linux) then the path
    should have the right slashes (Linux) or backslashes (Windows).
    A conversion is performed depending on the platform.
    '''

    # Get current directory from the system.
    currdir = os.getcwd()
    if ("/" in currdir):
        # Change all "\" into "/"
        outpath = path.replace("\\","/")
    elif ("\\" in currdir):
        # Change all "/" into "\\"
        outpath = path.replace("/","\\")
    else:
        outpath = path
    return outpath

    
def avg_within_range_2dim(domain,values,begin,end):
    '''
    Find all the values which are needed to cover the range begin to end.
    The calculate the average as the sum of the values divided by the number of values.
    domain must be sorted!
    '''
    if (len(domain) != len(values)):
        raise MyError("avg_within_range_2dim: domain and values do not have same length.",\
                      "Domain: "+str(len(domain))+" Values: "+str(len(values)))

    # Get start and end number of part of the list1 that is needed.
    ibegin,iend = find_within_range(domain,begin,end)
    # Extract them in an new list
    selected_numbers = values[ibegin:iend+1]
    # Calculate the average
    return sum(selected_numbers)/float(len(selected_numbers))

def sum_within_range_2dim(domain,values,begin,end):
    '''
    Find all the values which are needed to cover the range begin to end.
    The calculate the sum of the selected values.
    domain must be sorted!
    '''
    if (len(domain) != len(values)):
        raise MyError("sum_within_range_2dim: domain and values do not have same length.",\
                      "Domain: "+str(len(domain))+" Values: "+str(len(values)))

    # Get start and end number of part of the list1 that is needed.
    ibegin,iend = find_within_range(domain,begin,end)
    # Extract them in an new list
    #selected_numbers = values[ibegin:iend]
    # Calculate the sum
    return sum(values[ibegin:iend+1])

def avg_within_range(list1,begin,end):
    '''
    Find all the values which are needed to cover the range begin to end.
    The calculate the average as the sum of the values divided by the number of values.
    list1 must be sorted!
    '''
    # Get start and end number of part of the list1 that is needed.
    ibegin,iend = find_within_range(list1,begin,end)
    # Extract them in an new list
    selected_numbers = list1[ibegin:iend+1]
    # Calculate the average
    return sum(selected_numbers)/float(len(selected_numbers))

def sum_within_range(list1,begin,end):
    '''
    Find all the values which are needed to cover the range begin to end.
    The calculate the sum of the selected values.
    list1 must be sorted!
    '''
    # Get start and end number of part of the list1 that is needed.
    ibegin,iend = find_within_range(list1,begin,end)

    # Calculate the sum
    return sum(list1[ibegin:iend+1])

def find_within_range(list1,begin,end):
    '''
    Find all the values which are needed to cover the range begin to end.
    This means that the last element before begin is taken into acount and the the first after end is included.
    list1 must be sorted!
    '''
    ibegin = 0
    iend = len(list1) 
    for item in range(1,len(list1)):
        if (float(list1[item]) < begin):
            ibegin = item
    for item in range(len(list1)-1,-1,-1):
        if (float(list1[item]) > end):
            iend = item
    return ibegin,iend
        
def pickleLoader(pklFile):
    try:
        while True:
            yield pickle.load(pklFile)
    except EOFError:
        pass

def get_content(pklFile):
    '''
    Get the content of the file (one line)..
    '''
    fp = open(pklFile,"rb")
    out = pickle.load(fp)
    fp.close()
    return out

def get_content_rm(pklFile):
    '''
    Get the content of the file (one line) and remove the file.
    '''
    fp = open(pklFile,"rb")
    out = pickle.load(fp)
    fp.close()
    os.remove(pklFile)
    return out

def check_header(filename,species,year):
    # Check the startup file
    fp = open(filename, 'rb')
    print("Open " + filename)
    # Read header of file, whether this file belongs to this situation
    header=pickle.load(fp)
    # Check year of startup file
    if (abs(header[0] - year) > 0.01):
        #raise MyError("File " + filename + " does not have the year.",\
        #              "Year found in file: "+str(header[0]),\
        #              "Year needed in file: "+str(year))
        print("WARNING: File " + filename + " does not have the year.",\
              "Year found in file: "+str(header[0]),\
              "Year needed in file: "+str(year)) 
    # Check whether the speciesnames and order are equal
    for item in range(len(species)):
        try:
            if (header[item+1] != species[item].get_name()):
                raise MyError("File " + filename + " does not have the same species or same order.",\
                              "Species found in file: "+str(header[item]),\
                              "Species needed in file: "+str(species[item-1].get_name()))
        except:
            raise MyError("File " + filename + " does not have the same species or same order.",\
                          "Last species found in file: "+str(header[item]),\
                          "Species needed in file: "+str(species[item].get_name()))

    fp.close()

def read_header(fp=None,filename=None):
    # Read header of output file. Read year and name of all the species.
    if (fp == None and filename != None):
        fp_loc = open(filename, 'rb')
    elif (fp != None and filename == None):
        # File is already open.
        fp_loc = fp
        pass
    else:
        raise MyError("To read a header, you must specify a filepointer or filename")

    # Read header
    list = pickle.load(fp_loc)

    # Close file if it was opened here
    if fp == None:
        fp_loc.close()

    # Get time
    year = list[0]
    # Remove year from specienames list
    list.pop(0)
  
    return year,list

def make_header(filename,species,year):
    # Make output file. Start with year and name of all the species.
    fp = open(filename, 'wb')
    line = [year]
    for item in range(len(species)):
        line.append(species[item].get_name())
    pickle.dump(line,fp, -1)
    fp.close()

def make_header_budget(filename,species,proc,year):
    '''
    Make output file with process rates. 
    year, load, outflow, change in volume, particle exchanges, processes
    '''
    fp = open(filename, 'wb')
    line = [year]
    for ispec in range(len(species)):
        line.append(species[ispec].get_name()+"_loadIN")
        line.append(species[ispec].get_name()+"_load_src")
        line.append(species[ispec].get_name()+"_load_up")
        line.append(species[ispec].get_name()+"_load_hw")
        line.append(species[ispec].get_name()+"_loadOUT")
    for iproc in range(len(proc)):
        line.append(proc[iproc].get_val("name"))
    pickle.dump(line,fp, -1)
    fp.close()

def initialize_budget(params, species, proc):
    '''
    Create a list of empty budget items (species' loads and process rates) for every order
    '''
    xbud = []

    for iorder in range(params.norder):
        xbud.append([0.]*(5*len(species)+len(proc)))

    return xbud

def add_budget_load(xbud, iorder, ispec, term, value):
    '''
    Add load terms in budget for a specified order
    '''
    if (term == "loadIN"):
        xbud[iorder][5*ispec] += value
    elif (term == "load_src"):
        xbud[iorder][5*ispec+1] += value
    elif (term == "load_up"):
        xbud[iorder][5*ispec+2] += value
    elif (term == "load_hw"):
        xbud[iorder][5*ispec+3] += value
    elif (term == "loadOUT"):
        xbud[iorder][5*ispec+4] += value
        
def add_budget_procs(species, xbud, iorder, nrivers, dt, proc_rates):
    '''
    Add process rates in budget for a specified order
    '''
    for iproc in range(len(proc_rates)):
        xbud[iorder][5*len(species)+iproc] += nrivers * proc_rates[iproc] * dt
    

def get_val_dict(dict,key):
    try:
        return dict[key]
    except KeyError:
        return None

def find_spec(list,name,lfatal=0):
    ''' 
    Returns the index of the list where the name is found
    '''
    for i in range(len(list)):
        lab = list[i].get_name()
        if (lab == name): return i
    if (lfatal):
        raise MyError("Species " + name + " is not found in the list.")
    else:
        print("Species " + name + " is not found in the list.")
        
def is_spec(list,name): 
    ''' 
    Returns 1 if the species is defined, 0 otherwise
    '''
    found = False
    i = 0
    while ((found == False) & (i < len(list))):
        lab = list[i].get_name()
        if (lab == name): found = True
        i = i + 1
    return found

def find_func(list,name):
    ''' 
    Returns the index of the list where the name is found.
    '''
    for i in range(len(list)):
        lab = list[i].get_val("name")
        if (lab == name): return i
    
    raise MyError("Process " + name + " is not found in the list.")

def find_func_dict(list):
    ''' 
    Makes a dictionary for name of the functions and the index in the list. 
    Dictionary returs the index of the list where the name is found.
    '''
    dict = {}
    for i in range(len(list)):
        lab = list[i].get_val("name")
        dict[lab] = i
    return dict
    
#def get_conc(list):
def get_amount(list):
    ''' 
    Returns a list with the concentration.
    List contains Spec items.
    '''
    out = []
    for i in range(len(list)):
        #out.append(list[i].get_conc())
        out.append(list[i].get_amount())
    
    return out

def MM(conc,half_sat):
    '''
    Michaelis-Menten equation
    '''
    if ((half_sat + conc) > 0.):
        return (conc / (half_sat + conc))
    else:
        return 0.

def fT(topt, sigma, temperature):
    '''
    Temperature dependancy for biogeochemical processes (cf. RIVE)
    '''
    return math.exp(-(topt - temperature)**2 / sigma**2)

def light_extinct(ss_conc, eta_ss, depth):
    '''
    Reduction factor of photosynthesis due to 
    suspended sediments in the water column (light extinction)
    '''
    lext = 1.0
    if ((eta_ss * ss_conc * depth) > 0.0):
        #print "eta_ss  " + str(eta_ss) + " -- ss_conc  " + str(ss_conc) + " -- depth  " + str(depth)
        lext = 1 - math.exp((-1) * eta_ss * ss_conc * depth)
        lext /= eta_ss * ss_conc * depth
    return lext

def calculate_vel(Q, width, depth):
    '''
    Returns flow velocity in km/yr
    Q is in km3/yr, width and depth are in m
    '''
    vel = 0.
    if (width * depth > 0.):
        vel = Q * 1.e6 / (width * depth)
    return vel

def get_years(dirname,pre_text=None):
    '''
    Read all objects on the directory dirname and returns the list of all years present.
    The pre_text is used for directories which have a text before the year. So if pre_text
    is given than all directory names are returned which have the form "pre_text"float(...).
    Otherwise all directory names are returned which are a float.
    '''
    # Read all the directory on the input directory and convert this to different forms.
    years = os.listdir(dirname)
    # Throw away all objects which are not a year or a directory
    for item in range(len(years)-1,-1,-1):
        if (not os.path.isdir(os.path.join(dirname,years[item]))):
            # Throw away all objects which not a directory
            del years[item]
        else:
            # Not a year
            try:
                if (pre_text == None): 
                    qq = float(years[item])
                else:
                    if (years[item].startswith(pre_text)):
                        qq = float(years[item][len(pre_text):])
                        # Change the directory name without pre_text.
                        years[item] = years[item][len(pre_text):]
                    else:
                        del years[item]
            except:
                del years[item]

    # Sort the years.
    years.sort()

    return years

if __name__ == "__main__":

    # Set the general path for the own python modules
    import general_path

    # Import general python modules
    import os
    import sys
    import traceback
    import pickle
 
    # Import own general modules
    from error import *

    domain = [1,2,3,4,5,6,7,8,9]
    values = [1,2,3,4,5,6,7,8,9]
    begin = 2.00000000001
    end = 9.5
    print(domain)
    print("Start and end index: ",find_within_range(domain,begin,end)," Value expected: (1,9)")
    ibegin,iend = find_within_range(domain,begin,end)
    print("Length of part of domain: ",min(len(domain)-1,iend) - ibegin + 1," Value expected: 8")
    print("Sum of part of domain: ",sum_within_range(domain,begin,end)," Value expected: 44")
    print("Avg of part of domain: ",avg_within_range(domain,begin,end)," Value expected: 5.5")
    print("Sum of part of values: ",sum_within_range_2dim(domain,values,begin,end)," Value expected: 44")
    print("Avg of part of values: ",avg_within_range_2dim(domain,values,begin,end)," Value expected: 5.5")
    print("SUM: ",sum_within_range(domain,2.1,2.1)," Value expected: 5")
    begin = 2.00000000001
    end = 3.5
    print("Start and end index: ",find_within_range(domain,begin,end)," Value expected: (1,3)")
    ibegin,iend = find_within_range(domain,begin,end)
    print("Length of part of domain: ",min(len(domain)-1,iend) - ibegin + 1," Value expected: 3")
    print("Sum of part of domain: ",sum_within_range(domain,begin,end)," Value expected: 9")
    print("Avg of part of domain: ",avg_within_range(domain,begin,end)," Value expected: 3.0")
    print("Sum of part of values: ",sum_within_range_2dim(domain,values,begin,end)," Value expected: 9")
    print("Avg of part of values: ",avg_within_range_2dim(domain,values,begin,end)," Value expected: 3.0")
    begin = -2.00000000001
    end = 3.5
    print("Start and end index: ",find_within_range(domain,begin,end)," Value expected: (0,3)")
    ibegin,iend = find_within_range(domain,begin,end)
    print("Length of part of domain: ",min(len(domain)-1,iend) - ibegin + 1," Value expected: 4")
    print("Sum of part of domain: ",sum_within_range(domain,begin,end)," Value expected: 10")
    print("Avg of part of domain: ",avg_within_range(domain,begin,end)," Value expected: 2.5")
    print("Sum of part of values of same begin and end: ",sum_within_range_2dim(domain,values,2.1,2.1)," Value expected: 5")
    print("Avg of part of values of same begin and end: ",avg_within_range_2dim(domain,values,2.1,2.1)," Value expected: 2.5")
    values = [5,2,3,4,5,6,7,8,9]
    print("Sum of part of values: ",sum_within_range_2dim(domain,values,begin,end)," Value expected: 14")
    print("Avg of part of values: ",avg_within_range_2dim(domain,values,begin,end)," Value expected: 3.5")
    # Test error handling
    try:
         print(sum_within_range_2dim(domain,values[1:],begin,end))
    except MyError:
         print("Error handling went okay.")
    except:
         print("Error handling of sum_within_range_2dim went wrong! ERROR")
    

