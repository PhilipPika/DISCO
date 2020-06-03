# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import general_class

def do():
    arguments = general_class.General()
    arguments.add_item("vol",None)
    arguments.add_item("dvoldt",None)
    arguments.add_item("globrad",None) 
    arguments.add_item("temperature",None)
    arguments.add_item("windspeed",None)
    arguments.add_item("flow_velocity", None) #WJ201708
    arguments.add_item("width",None) 
    arguments.add_item("depth",None) 
    arguments.add_item("slope",None) 
    arguments.add_item("discharge", None)
    arguments.add_item("length", None)
    arguments.add_item("residence_time", None)
    arguments.add_item("area", None)
    arguments.add_item("icell", None)
    arguments.add_item("next_cell", None)
    return arguments
