def do(nc_data):
  dumlist = list(nc_data.variables.keys())
  dumlist.remove('time')
  try:
      dumlist.remove('lat')
  except(ValueError):
      dumlist.remove('latitude')
  try:
      dumlist.remove('lon')
  except(ValueError):
      dumlist.remove('longitude')
  varname = dumlist[0]
  return varname