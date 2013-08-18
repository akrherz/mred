# Something simple to print out ranges of variables and such

import sys, os, glob
import netCDF3, numpy

def find_grid(nc):
    lats = nc.variables['lat']
    lons = nc.variables['lon']
    findlat = 41.99
    findlon = -93.62
    for i in range(len(lats[:])):
        if lats[i] > findlat:
            gridy = i
            break
    for i in range(len(lons[:])):
        if lons[i] > findlon:
            gridx = i
            break
    return gridx, gridy


def main():
    runid = sys.argv[1]
    os.chdir("final/%s" % (runid,))
    files = glob.glob("*.nc")
    gridx = 0
    gridy = 0
    print "%-8s %-12s %14s %14s %14s" % ("VARIABLE", "UNITS", "MINIMUM",
         "AVERAGE", "MAXIMUM") 
    for file in files:
        nc = netCDF3.Dataset(file, 'r')
        if gridx == 0:
            gridx, gridy = find_grid( nc )

        vname = file.split("_")[0]
        data = nc.variables[vname][:,gridy,gridx]
        print "%-8s %-12s %14.6f %14.6f %14.6f" % (vname, 
             nc.variables[vname].units, numpy.min( data ),
             numpy.average( data ), numpy.max( data ) )

main()
