# Print out an ASCII table of which files we have and don't have :)

import os, mx.DateTime

years = range(1982, 2009)
days = ['1111', '1112', '1113', '1114', '1115', 
        '1121', '1122', '1123', '1124', '1125',
        '1129', '1130', '1201', '1202', '1203']

print "YEAR",
for day in days:
  print "%s" % (day,) ,
print

thres = mx.DateTime.DateTime(2010,3,1)

for year in years:
  print "%s" % (year,) ,
  for day in days:
    runid = "%s%s00" % (year, day)
    fp = "final/%s/cuhusa_IMM5_%s%s03_CFS01.nc" % (runid, year, day)
    if os.path.isfile( fp ):
      ts = mx.DateTime.DateTime(1970,1,1) + mx.DateTime.RelativeDateTime(seconds = os.stat(fp)[-2]) - mx.DateTime.RelativeDateTime(hours=5)
      if ts > thres:
        print " XX ",
      else:
        print " OO ",
    else:
      print "    ",
  print
