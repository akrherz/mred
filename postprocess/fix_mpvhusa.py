# Fix standard_name

import netCDF3, glob

files = glob.glob("final/*/mpvhusa*")
for file in files:
  nc = netCDF3.Dataset( file, 'a')
  print nc.variables['mpvhusa'].standard_name 
  nc.close()
