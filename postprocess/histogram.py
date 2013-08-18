# Actually compute the histograms, yummy
# Dump out to netCDF

import os
import numpy
import netCDF3

years = range(1982, 2009)
days = ['1111', '1112', '1113', '1114', '1115', 
        '1121', '1122', '1123', '1124', '1125',
        '1129', '1130', '1201', '1202', '1203']

# Figure out our precipitation bins
# Floor: 0.25 mm/hour
# Max:   75   mm/hour  ??
# Interval: 0.25 mm/hour

bins = numpy.arange( 0.25 / 3600.0, 75.0 / 3600.0, 0.25 / 3600.)

output = netCDF3.Dataset("mred_precip_histogram.nc", 'w')
output.createDimension("bins", len(bins) -1 )
output.createDimension("runid", len(years) * len(days) )

data = output.createVariable("count", numpy.float32, ("runid", "bins") )
data.long_name = "Grid cell count"

ncbins = output.createVariable("bins", numpy.float32, ("bins") )
ncbins.long_name = "Precipitation Bins"
ncbins.units = "kg m-2 s-1"
ncbins[:] = bins[:-1]

cnt = 0
for year in years:
  for day in days:
    ar = numpy.zeros( None, 'f')
    for box in range(6):
      fp = "final.prechist/box%s_pr_IMM5_%s%s03_CFS01.nc" %(box, year, day)
      if not os.path.isfile(fp):
        continue
      nc = netCDF3.Dataset(fp, 'r')
      ar = numpy.append( ar, nc.variables['pr'][:] )
      nc.close()
    hist, edges = numpy.histogram(ar, bins)
    print numpy.shape(hist)
    data[cnt,:] = hist
    cnt += 1

output.close()
