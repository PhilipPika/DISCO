def read_parameters_ncinifile(inifile):
    '''
    Reads general and variable parameters in an inifile
    Stores it and returns in a dictionary form
    '''
    general_params = general_class.General() #dictionary
    src_params = []                          #list of dictionaries
    var_params = []                          #list of dictionaries
    
    check_file(inifile)
  
    # Open file with all parameters
    fp = open(inifile,'r')

    # Determine all the information per line
    type_option = None
    for line in fp.readlines():
        option = [elem.strip() for elem in line.split('#')[0].split('=')]
        if (len(option[0]) < 1):
            # Empty line
            continue
        
        # Check whether the first element is start of a new class 
        # (first character is a '[' )
        if option[0].startswith('['):
            # Start of a new class
            type_option = option[0].strip().upper()
            print('\n' + type_option)
            if (type_option == '[SOURCE]'):
                src_params.append(general_class.General())
            elif (type_option == '[SPECIE]'):
                var_params.append(general_class.General())
            elif (type_option != '[GENERAL]'):
                raise MyError("File '%s' contains errors. Class %s is not known." % (inifile,option[0]))

        else:
            # Here information from a class is given
            print(option[0] + ' = ' + option[1])
            try:
                if (type_option == '[GENERAL]'):
                    general_params.add_item(option[0],option[1])
                elif (type_option == '[SOURCE]'):
                    src_params[-1].add_item(option[0],option[1])
                elif (type_option == '[SPECIE]'):
                    var_params[-1].add_item(option[0],option[1])
            except IndexError:
                raise MyError("File '%s' contains errors." % inifile,\
                              "Can not handle line: %s" % line)
    # Close input file
    fp.close()

    # Check if all variables have a name
    # 2DI: other checks... src_params...
    print('\n')
    for ivar in range(len(var_params)):
        if not('name' in var_params[ivar].get_attrib()):
            raise MyError("File '%s' contains errors. No name was given to species nb %d." % (inifile,ivar+1))

    #return general_params, var_params
    return general_params, src_params, var_params
