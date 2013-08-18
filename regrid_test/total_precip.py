
import netCDF3, glob, numpy

for file in glob.glob("*.nc"):
  nc = netCDF3.Dataset(file, 'r')
  print file, numpy.sum( nc.variables['pr'][:] )
  nc.close()
