# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************

import os

def do(cmd_str, tasktxt, nodetxt, errtxt, stdtxt, jobtxt, outputdir):
  basename = "sbatch.sh"
  full_filename = os.path.join(outputdir,basename)
  with open(full_filename, 'w') as fp:
    fp.write('#!/bin/bash \n')
    fp.write('#SBATCH '+tasktxt+'\n')
    fp.write('#SBATCH '+nodetxt+'\n')
    fp.write('#SBATCH '+errtxt+'\n')
    fp.write('#SBATCH '+stdtxt+'\n')
    fp.write('#SBATCH '+jobtxt+'\n')
    fp.write(cmd_str)