# Actually compute the histograms, yummy
# Dump out to netCDF

import os
import numpy
import netCDF3
import mx.DateTime

years = range(1982, 2009)
days = ['1111', '1112', '1113', '1114', '1115', 
        '1121', '1122', '1123', '1124', '1125',
        '1129', '1130', '1201', '1202', '1203']

# Figure out our precipitation bins
# Floor: 0.25 mm/hour
# Max:   75   mm/hour  ??
# Interval: 0.25 mm/hour

bins = numpy.arange( 0.25 / 3600.0, 75.0 / 3600.0, 0.25 / 3600.)

output = netCDF3.Dataset("mred_precip_histogram_monthly.nc", 'w')
output.createDimension("bins", len(bins) -1 )
output.createDimension("runid", len(years) * len(days) )
output.createDimension("month", 4 )
output.createDimension("latlonbox", 6 )

month = output.createVariable("month", numpy.float32, ("month"))
month.units = "months since 2000-01-01"
month[:] = range(4)

latlonbox = output.createVariable("latlonbox", numpy.float32, ("latlonbox"))
latlonbox[:] = range(6)

data = output.createVariable("count", numpy.float32, ("runid", "latlonbox", "month", "bins") )
data.long_name = "Grid cell count"

ncbins = output.createVariable("bins", numpy.float32, ("bins") )
ncbins.long_name = "Precipitation Bins"
ncbins.units = "kg m-2 s-1"
ncbins[:] = bins[:-1]

cnt = 0
for year in years:
  for day in days:
    ts = mx.DateTime.strptime("%s%s" % (year,day), "%Y%m%d")
    for box in range(6):
      fp = "final.prechist/box%s_pr_IMM5_%s%s03_CFS01.nc" %(box, year, day)
      if not os.path.isfile(fp):
        continue
      nc = netCDF3.Dataset(fp, 'r')
      for month in range(1,5):
        ts0 = ts + mx.DateTime.RelativeDateTime(years=1,month=month,day=1)
        ts1 = ts0 + mx.DateTime.RelativeDateTime(months=1)
        offset0 = int((ts0 - ts).hours / 3.0 + 1.)
        offset1 = int((ts1 - ts).hours / 3.0 + 2.)
        #print ts, ts0, ts1, offset0, offset1
        ar = nc.variables['pr'][offset0:offset1] 
        hist, edges = numpy.histogram(ar, bins)
        data[cnt,box,month-1,:] = hist
      nc.close()
    cnt += 1

output.close()
