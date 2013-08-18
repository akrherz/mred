# Something simple to print out ranges of variables and such

import sys, os, glob
import netCDF3, numpy

def find_grid(nc):
    """
    This is method, which is called with an argument of nc, which is
    a netcdf object
    """
    # Query two variables out of the netcdf object
    lats = nc.variables['lat']
    lons = nc.variables['lon']
    # I am looking for this
    findlat = 43.99
    findlon = -99.62
    # Loop over lats 
    for i in range(len(lats[:])):
        if lats[i] > findlat:
            gridy = i
            break
    # Loop over lons
    for i in range(len(lons[:])):
        if lons[i] > findlon:
            gridx = i
            break
    # return results
    return gridx, gridy


def main():
    """
    My main method
    """
    # Command line argument number one
    vname = sys.argv[1]
    # Change directory
    os.chdir("final")
    gridx = 0
    gridy = 0
    # print out a header
    print "%-10s %-8s %-12s %14s %14s %14s" % ("RUNID", "VARIABLE", "UNITS", "MINIMUM",
         "AVERAGE", "MAXIMUM") 
    # Look for all directories in the final/ folder
    for runid in glob.glob("*"):
        # Change to that directory
        os.chdir("/mnt/tera11/mred/postprocess/final/"+ runid)
        # open the netcdf file I care about
        file = "%s_IMM5_%s3_CFS01.nc" % (vname, runid[:-1])
        nc = netCDF3.Dataset(file, 'r')
        if gridx == 0:
            gridx, gridy = find_grid( nc )

        # pull out the data
        data = nc.variables[vname][:,gridy,gridx]
        # Check if max value in the data is larger than 1
        if numpy.max( data ) > 1:
            # Print out some information
            print "%s %-8s %-12s %14.6f %14.6f %14.6f" % (runid, vname, 
             nc.variables[vname].units, numpy.min( data ),
             numpy.average( data ), numpy.max( data ) )

# Call main()
main()
