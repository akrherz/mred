# Remap the CFS netcdf files to MRED grid

import os, glob

for dir in glob.glob("final/2008111500"):
  for fp in glob.glob(dir+"/*.nc"):
    bogus, runid, nc = fp.split("/")
    newdir = "final.monavg/%s" % (runid,)
    if not os.path.isdir(newdir):
      os.mkdir( newdir )

    #cmd = "cdo remapbil,cdo_mredgrid.txt %s final.mredgrid/%s/%s" % (
    #      fp, runid, nc)
    #os.system( cmd )
    cmd = "cdo monavg %s final.monavg/%s/mon_%s" % (
          fp, runid, nc)
    print cmd
    os.system( cmd )
