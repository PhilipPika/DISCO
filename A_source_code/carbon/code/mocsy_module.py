# ******************************************************
## Copyright 2019, PBL Netherlands Environmental Assessment Agency and Utrecht University.
## Reuse permitted under Gnu Public License, GPL v3.
# ******************************************************
'''
imports mocsy module
'''
try: 
  import mocsy
except ModuleNotFoundError:
  import os
  import sys
  root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
  if "/" in root:
      # working in linux OS
      path = os.path.join(root,"libs","mocsy","linux")
  if (os.path.exists(path)):
      sys.path.insert(3, path)
