# Generate a histogram of 3 hourly precipitation

# 
boxes = [
 '-123.0,-121.0,37.0,38.0',
 '-122.0,-120.0,36.0,37.0',
 '-121.7,-119.0,35.0,36.0',
 '-120.7,-117.0,34.0,35.0',
 '-118.7,-116.0,33.0,34.0',
 '-117.5,-116.0,32.5,33.0',
]

import os, glob, sys

for dir in glob.glob("final/*"):
  for fp in glob.glob(dir+"/pr_*.nc"):
    bogus, runid, nc = fp.split("/")
    newdir = "final.prechist/"
    #if not os.path.isdir(newdir):
    #  os.mkdir( newdir )
    for i in range(len(boxes)):
      #cmd = "cdo sellonlatbox,%s -histcount,0.00000500,0.00000826,0.00001309,0.00002007,0.00002993,0.00004363,0.00006364,0.00009318,0.00013925,0.00022762,0.00037651,0.0005,inf %s %s/hist_box%s_%s" % (
      cmd = "cdo sellonlatbox,%s %s %s/box%s_%s" % (
          boxes[i], fp, newdir, i, nc)
      print cmd
      os.system( cmd )
      print
