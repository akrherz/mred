# All in one processor for ISUMM5 output to MRED output specification
# Simply call with a argument for the CFS run
# daryl herzmann, 15 Jul 2009
# 2 Apr 2010 - Fix bug with surface runoff calculations

# Standard Library Python imports
import sys, glob, os
# Special Python libraries that need to be installed
import mx.DateTime, mm5_class, Numeric
import netCDF3


run = sys.argv[1] # 1982111100
ts0 = mx.DateTime.strptime(run[:8], "%Y%m%d")
tsend = ts0 + mx.DateTime.RelativeDateTime(years=1,month=5,day=1)

BASEDIR = "/mnt/tera11/mred"
DATADIR = "%s/postprocess/data/%s/" % (BASEDIR, run)

# Simple reference for filenames
NCOUT, MMOUT = range(2)
datasources = ["NCOUT_DOMAIN1_??_interp.nc", "MMOUTP_DOMAIN1_??_interp.nc"]

# Dictionary of Lookup values for how variables are to be retrieved from
# the ISUMM5 produced output files. Each entry has the following attrs
#  source: Either NCOUT or MMOUT to specify which file to use
#  vname : Variable name in the respective file to use
#  lname : Output NetCDF long_name
#  units : Units attribute to use when writting the output netcdf
#  sname : Output NetCDF standard_name attribute
#  coord : Optional setting for NetCDF Variable Attribute coordinates
#  cellm : Optional setting if cell_methods were used
#  quo   : Optionally set if you need the value divided by something
#  positive: Optional if necessary
outputvars = {
 'tas': {'source': MMOUT,
         'vname' : 't2',
         'units' : 'K', 
         'coord' : 'lon lat height',
         'cellm' : 'time: Instantaneous',
         'lname' : 'Surface Air Temperature',
         'sname' : 'air_temperature'},
 'tasmax': {'units'  : 'K',
            'source' : NCOUT,
            'vname'  : 'tmax',
            'coord' : 'lon lat height',
            'cellm' : 'time: maximum (interval: 3 hours)',
            'lname' : 'Maximum Surface Air Temperature',
            'sname'  : 'air_temperature'},
 'tasmin': {'units'  : 'K', 
            'source' : NCOUT,
            'vname'  : 'tmin',
            'coord' : 'lon lat height', 
            'lname' : 'Minimum Surface Air Temperature',
            'cellm' : 'time: minimum (interval: 3 hours)',
            'sname'  : 'air_temperature'},
 'rsds': {'source' : NCOUT,
          'vname'  : 'gsw1',
          'units'  : 'W m-2',
             'coord' : "lon lat",
          'cellm'  : 'time: mean (interval: 3 hours)',
          'lname'  : 'Surface Downwelling Shortwave Radiation',
          'sname'  : 'surface_downwelling_shortwave_flux_in_air'},
 'uas' : {'source' : MMOUT, 
          'vname'  : 'u10',
          'units'  : 'm s-1',
          'coord'  : 'lon lat height',
          'cellm'  : 'time: Instantaneous',
          'lname'  : 'Zonal Surface Wind Speed',
          'sname'  : 'eastward_wind'},
 'vas' : {'source' : MMOUT, 
          'vname'  : 'v10',
          'units'  : 'm s-1', 
          'coord'  : 'lon lat height', 
          'cellm'  : 'time: Instantaneous',
          'lname'  : 'Meridional Surface Wind Speed',
          'sname'  : 'northward_wind'},
 'ps'  : {'source' : MMOUT, 
          'vname'  : 'psfc',
          'units'  : 'Pa',
             'coord' : "lon lat",
          'cellm'  : 'time: Instantaneous',
          'lname'  : 'Surface Pressure',
          'sname'  : 'surface_air_pressure'},
 'psl' : {'source' : NCOUT, 
          'vname'  : 'pmsl', 
          'units'  : 'Pa',
             'coord' : "lon lat",
          'cellm'  : 'time: Instantaneous',
          'lname'  : 'Sea Level Pressure',
          'sname'  : 'air_pressure_at_sea_level'},
 'zg850' : {'source' : MMOUT, 
            'vname'  : 'h', 
            'vindex' : 2, #!important, as it supports 3d var
            'units'  : 'm',
            'cellm'  : 'time: Instantaneous',
            'lname'  : '850 hPa Geopotential Height',
            'coord'  : 'lon lat plev',
            'sname'  : 'geopotential_height'},
 'zg500' : {'source' : MMOUT,
            'vname'  : 'h',
            'vindex' : 3, #!important, as it supports 3d var
            'units'  : 'm',
            'cellm'  : 'time: Instantaneous',
            'lname'  : '500 hPa Geopotential Height',
            'coord'  : 'lon lat plev' ,
            'sname'  : 'geopotential_height'},
 'zg200' : {'source' : MMOUT, 
            'vname'  : 'h', 
            'vindex' : 4, #!important, as it supports 3d var
            'units'  : 'm', 
            'cellm'  : 'time: Instantaneous',
            'lname'  : '200 hPa Geopotential Height',
            'coord'  : 'lon lat plev' ,
            'sname'  : 'geopotential_height'},
 'snd'   : {'source' : MMOUT, 
            'vname'  : 'snowh',
            'units'  : 'm', 
             'coord' : "lon lat",
            'cellm'  : 'time: Instantaneous',
            'lname'  : 'Snow Depth',
            'sname'  : 'surface_snow_thickness'},
 'hfls'  : {'units'  : 'W m-2', 
            'source' : NCOUT, 
            'vname'  : 'qfx1',
            'coord'  : "lon lat",
             'positive': 'up',
            'cellm'  : 'time: average (interval: 3 hours)',
            'lname'  : 'Surface Latent Heat Flux',
            'sname'  : 'surface_upward_latent_heat_flux'},
 'hfss'  : {'units'  : 'W m-2', 
            'source' : NCOUT, 
            'vname'  : 'hfx1', 
            'coord'  : "lon lat",
             'positive': 'up',
            'cellm'  : 'time: average (interval: 3 hours)',
            'cellm'  : 'time: average (interval: 3 hours)',
            'lname'  : 'Surface Sensible Heat Flux',
            'sname'  : 'surface_upward_sensible_heat_flux'},
 'rlds'  : {'units'  : 'W m-2', 
            'source' : NCOUT, 
            'vname'  : 'glw1',
            'coord'  : "lon lat",
            'cellm'  : 'time: average (interval: 3 hours)',
            'lname'  : 'Surface Downwelling Longwave Radiation',
            'sname'  : 'surface_downwelling_longwave_flux_in_air'},
 'prw'   : {'units'  : 'kg m-2', 
            'source' : NCOUT, 
            'vname'  : 'tpw',
             'coord' : "lon lat",
            'cellm'  : 'time: Instantaneous',
            'lname'  : 'Precipitable Water',
            'sname'  : 'atmosphere_water_vapor_content'},
 'cuhusa' : {'units' : 'kg m-1 s-1', 
             'source': NCOUT, 
             'quo'   : 10800.0,
             'coord' : "lon lat",
             'vname' : 'uqvflux',
             'positive': 'east',
             'cellm' : 'time: average (interval: 3 hours)',
             'lname' : 'Vertically integrated zonal moisture flux',
             'sname' :
   'eastward_atmosphere_water_vapor_transport_across_unit_distance'},
 'cvhusa' : {'units' : 'kg m-1 s-1', 
             'source': NCOUT, 
             'quo'   : 10800.0,
             'coord' : "lon lat",
             'vname' : 'vqvflux',
             'positive': 'north',
             'cellm' : 'time: average (interval: 3 hours)',
             'lname' : 'Vertically integrated meridional moisture flux',
             'sname' : 
   'northward_atmosphere_water_vapor_transport_across_unit_distance'},
 'mpuhusa': {'units' : 'm s-1', 
             'source': NCOUT, 
	     'positive': 'east',
             'vname' : 'uqv850',
             'coord' : "lon lat plev",
             'cellm' : 'time: Instantaneous',
             'lname' : '850 hPa zonal flux of specific humidity',
             'sname' : 'product_of_eastward_wind_and_specific_humidity'},
 'mpvhusa': {'units' : 'm s-1', 
             'source': NCOUT, 
             'vname' : 'vqv850',
	     'positive': 'north',
             'coord' : "lon lat plev",
             'cellm' : 'time: Instantaneous',
             'lname' : '850 hPa meridional flux of specific humidity',
             'sname' : 'product_of_northward_wind_and_specific_humidity'},
}

def generate_filename(varname):
    """
    Return the filename that should be used given a variable name
    """
    return "%s_IMM5_%s3_CFS01.nc" % (varname, run[:-1])

def special_soilw():
    """
    Need a special function to compute soil moisture content kg m-2
    """
    os.chdir( DATADIR )
    fname = generate_filename("mrso")
    nc = netCDF3.Dataset( fname , 'w')
    nc_setup( nc , 'mrso')

    myvar = nc.createVariable('mrso', 'f', ('time','lat','lon'))
    myvar.units = 'kg m-2'
    myvar.standard_name = 'soil_moisture_content'
    myvar.long_name = "Total Soil Moisture Content"
    myvar.cell_methods = "time: instantaneious"
    myvar.coordinates = "lon lat"
    myvar._FillValue = 1.e20
    myvar.original_name = '((soil_w_1*0.1)+(soil_w_2*0.3)+(soil_w_3*0.6)+soil_w_4)/1000.0/2.0'

    files = glob.glob( datasources[ MMOUT ])
    files.sort() # Make sure time lines up this way!
    tcounter = 0
    for file in files:
        tnc = netCDF3.Dataset( file )
        # layer 1, m-3 m-3  * depth of layer [m]
        l1 = tnc.variables['soil_w_1'][:] * 0.1
        # layer 2
        l2 = tnc.variables['soil_w_2'][:] * 0.3
        # layer 3
        l3 = tnc.variables['soil_w_3'][:] * 0.6
        # layer 4
        l4 = tnc.variables['soil_w_4'][:] * 1.0
        # Total depth of water in 2m, m  1mm = 1kg
        tot = l1 + l2 + l3 + l4
        tsteps = l1.shape[0]
        myvar[tcounter:(tcounter+tsteps)] = (tot / 1000.0 / 2.0).astype('f')
        tcounter += tsteps
        del(tnc)
    nc.close()

def special_runoff2():
    """
    Special function to compute surface runoff
    """
    os.chdir( DATADIR )
    fname = generate_filename("mrros")
    nc = netCDF3.Dataset( fname , 'w')
    nc_setup( nc , 'mrros')

    myvar = nc.createVariable('mrros', 'f', ('time','lat','lon'))
    myvar.units = 'kg m-2 s-1'
    myvar.standard_name = 'surface_runoff_flux'
    myvar.long_name = "Surface Runoff"
    myvar.coordinates = "lon lat"
    myvar.cell_methods = "time: average (interval: 3 hours)"
    myvar._FillValue = 1.e20
    myvar.original_name = 'sfcrnoff'

    files = glob.glob( datasources[ MMOUT ])
    files.sort() # Make sure time lines up this way!
    tcounter = 0
    for file in files:
        tnc = netCDF3.Dataset( file )
        # Surface runoff
        surface = tnc.variables['sfcrnoff'][:]
        tsteps = surface.shape[0]
        for i in range(tsteps):
            if i == 0 and tcounter == 0:
                s0 = surface[i]
            else:
                s0 = surface[i] - sold
            sold = surface[i]
            myvar[tcounter+i] = s0 / 10800.0
        tcounter += tsteps
        del(tnc)

    nc.close()


def special_runoff():
    """
    Special function to compute total runoff, actually easy
    """
    os.chdir( DATADIR )
    fname = generate_filename("mrro")
    nc = netCDF3.Dataset( fname , 'w')
    nc_setup( nc , 'mrro')

    myvar = nc.createVariable('mrro', 'f', ('time','lat','lon'))
    myvar.units = 'kg m-2 s-1'
    myvar.standard_name = 'runoff_flux'
    myvar.coordinates = "lon lat"
    myvar.long_name = "Surface and Subsurface Runoff"
    myvar.cell_methods = "time: average (interval: 3 hours)"
    myvar._FillValue = 1.e20
    myvar.original_name = 'sfcrnoff+ugdrnoff'

    files = glob.glob( datasources[ MMOUT ])
    files.sort() # Make sure time lines up this way!
    tcounter = 0
    for file in files:
        tnc = netCDF3.Dataset( file )
        # Surface runoff
        surface = tnc.variables['sfcrnoff'][:]
        # Subsurface runoff
        subsurface = tnc.variables['ugdrnoff'][:]
        tsteps = surface.shape[0]
        for i in range(tsteps):
            if i == 0 and tcounter == 0:
                s0 = surface[i]
                ss0 = subsurface[i]
            else:
                s0 = surface[i] - sold
                ss0 = subsurface[i] - ssold
            sold = surface[i]
            ssold = subsurface[i]
            myvar[tcounter+i] = (s0 + ss0) / 10800.0
        tcounter += tsteps
        del(tnc)

    nc.close()

def special_spechumidity():
    """
    Special function to compute Specific Humidity
    """
    os.chdir( DATADIR )
    fname = generate_filename("huss")
    nc = netCDF3.Dataset( fname , 'w')
    nc_setup( nc , 'huss')

    myvar = nc.createVariable('huss', 'f', ('time','lat','lon'))
    myvar.units = 'kg kg-1'
    myvar.standard_name = 'specific_humidity'
    myvar.long_name = "Surface Specific Humidity"
    myvar.cell_methods = "time: instantaneious"
    myvar._FillValue = 1.e20
    myvar.coordinates = "lon lat height"
    myvar.original_name = 'q2/(1+q2)'

    files = glob.glob( datasources[ MMOUT ])
    files.sort() # Make sure time lines up this way!
    tcounter = 0
    for file in files:
        tnc = netCDF3.Dataset( file )
        # 2m mixing ratio kg kg-1
        q2 = tnc.variables['q2'][:]
        tsteps = q2.shape[0]
        myvar[tcounter:(tcounter+tsteps)] = (q2 / (1.0 + q2)).astype('f')
        tcounter += tsteps
        del(tnc)

      # Compute Dew Point
      #d2 = t2[i] / (1+ 0.000425 * t2[i] * -(Numeric.log10(rh[i]/100.0)) )
      # Compute Saturation vapor pressure
      #pws = Numeric.exp( 77.3450 + (0.0057 * d2) - (7235 / d2)) / Numeric.power(d2,8.2)  
      #sh = 0.62198 * pws / (100000 + pws)

    nc.close()

def special_precip():
    """
    Need a special function to compute precip
    """
    os.chdir( DATADIR )
    fname = generate_filename("pr")
    nc = netCDF3.Dataset( fname , 'w')
    nc_setup( nc , 'pr')

    myvar = nc.createVariable('pr', 'f', ('time','lat','lon'))
    myvar.units = 'kg m-2 s-1'
    myvar.standard_name = 'precipitation_flux'
    myvar.long_name = "Precipitation"
    myvar._FillValue = 1.e20
    myvar.coordinates = "lon lat"
    myvar.cell_methods = "time: average (interval: 3 hours)"
    myvar.original_name = 'raincon+rainnon'
    myvar.history = 'v2 code 20090806'

    files = glob.glob( datasources[ MMOUT ])
    files.sort() # Make sure time lines up this way!
    tcounter = 0
    for file in files:
        tnc = netCDF3.Dataset( file )
        # Values are in cm
        non = tnc.variables['rain_non']
        tsteps = non.shape[0]
        for i in range(tsteps):
            con = tnc.variables['rain_con'][i]
            non = tnc.variables['rain_non'][i]
            tot = non + con
            #print tcounter, max( max( tot ) )
            # Write out! Convert to kg m-2 and divide by 10800 secs
            # Input is cm need to x10 to get to mm == kg m^-2
            myvar[tcounter] = (tot * 10.0 / 10800.0).astype('f')
            tcounter += 1
        del(tnc)

    nc.close()

def finalize_step():
    """
    This is going to be the most complicated, have to generate the necessary
    output files with everything set A-OK
    """
    os.chdir( DATADIR )
    # Loop over each variable
    for vname in outputvars.keys():
        meta = outputvars[ vname ]
        # Create the NetCDF File
        fname = generate_filename(vname)
        nc = netCDF3.Dataset( fname , 'w')
        nc_setup( nc , vname)
        # Now we create the variable of interest
        myvar = nc.createVariable(vname, 'f', ('time','lat','lon'))
        # Assign Units
        myvar.units = meta['units']
        # Assign Standard Name
        myvar.standard_name = meta['sname']
        # Assign Long Name
        myvar.long_name = meta["lname"]
        # Fill Value
        myvar._FillValue = 1.e20
        # Coordinate Singleton
        if meta.has_key('coord'):
            myvar.coordinates = meta['coord']
        # Cell Methods
        if meta.has_key('cellm'):
            myvar.cell_methods = meta['cellm']
        # Positive
        if meta.has_key('positive'):
            myvar.positive = meta['positive']
        # Original Name of the variable in the MM5 file
        myvar.original_name = meta['vname']

        # Lets get data already!
        files = glob.glob( datasources[ meta['source'] ])
        files.sort() # Make sure time lines up this way!
        tcounter = 0
        for file in files:
            tnc = netCDF3.Dataset( file )
            # Get the variable from the NetCDF File
            tvar = tnc.variables[ meta['vname'] ][:,:,:]
            if meta.has_key('vindex'):
                tvar = tnc.variables[ meta['vname'] ][:,meta['vindex'],:,:]
            # Figure out the time dimension length, always first
            tlen = tvar.shape[0]
            print '%s from tdx %s to %s from %s sample: %.2f' % (vname, tcounter,
                  tcounter+tlen, file, tvar[0,10,10] )
            if meta.has_key('quo'):
                tvar = tvar / meta['quo']
            myvar[tcounter:(tcounter+tlen)] = tvar
            tcounter += tlen
            del(tnc)

        # Done with var
        nc.close()

def nc_setup(nc, varname):
    """
    Does the standard stuff to a new NetCDF file
    """
    # Requirements
    nc.institution   = "Iowa State University, Ames, IA, USA"
    nc.source        = "MM5 (2009): atmosphere: MM5v3.6.3 non-hydrostatic, split-explicit; sea ice: Noah; land: Noah"
    nc.project_id    = "MRED"
    nc.table_id      = "Table 2"
    nc.realization   = 1
    nc.forcing_data  = "CFS01"
    nc.experiment_id = "CFS Seasonal Run %s" % (run,)
    # Optional
    nc.Conventions   = 'CF-1.0'
    nc.contact       = "Daryl Herzmann, akrherz@iastate.edu, 515-294-5978"
    nc.history       = "%s Generated" % (mx.DateTime.now().strftime("%d %B %Y"),)
    nc.comment       = "Runs processed on derecho@ISU, output processed on mred@ISU"
    nc.title         = "ISU MM5 model output prepared for MRED using CFS input"

    # Setup Dimensions
    tsteps = int( (tsend - ts0).days * 8 ) # 3 hourly
    nc.createDimension('time', tsteps)
    nc.createDimension('lat', 66)
    nc.createDimension('lon', 155)
    nc.createDimension('bnds', 2)

    # Generate the coordinate variables
    tm = nc.createVariable('time', 'd', ('time',) )
    tm.units = "days since %s 00:00:00.0" % (ts0.strftime("%Y-%m-%d"),)
    tm.calendar = "gregorian"
    tm.long_name = "time"
    tm.standard_name = "time"
    tm.axis = "T"
    tm[:] = Numeric.arange(0.125, (tsend - ts0).days + 0.125, 0.125)

    # Compute the time_bnds variable, somewhat yucky
    if varname in ['cuhusa', 'cvhusa', 'hfls', 'hfss', 'mrro',
                   'pr', 'psl', 'rsds', 'rlds', 'mrros']:
        tm.bounds = "time_bnds"

        tb = nc.createVariable('time_bnds', 'd', ('time','bnds') )
        val = Numeric.zeros( (tsteps,2), 'd' )
        val[:,0] = Numeric.arange(0, (tsend - ts0).days , 0.125)
        val[:,1] = Numeric.arange(0.125, (tsend - ts0).days + 0.125, 0.125)
        tb[:] = val

    lat = nc.createVariable('lat', 'd', ('lat',) )
    lat.units = "degrees_north"
    lat.long_name = "latitude"
    lat.standard_name = "latitude"
    lat.axis = "Y"
    lat[:] = Numeric.arange(24.75,49.5,0.375)

    lon = nc.createVariable('lon', 'd', ('lon',) )
    lon.units = "degrees_east"
    lon.long_name = "longitude"
    lon.standard_name = "longitude"
    lon.axis = "X"
    lon[:] = Numeric.arange(-124.75,-66.625,0.375)

    # Generate singletons
    if varname in ['huss','tas','tasmax', 'tasmin','uas','vas']:
        height = nc.createVariable('height', 'd')
        height.long_name = "height"
        height.standard_name = "height"
        height.units = "m"
        height.positive = "up"
        height.axis = "Z"
        if varname in ['huss','tas','tasmax','tasmin']:
            height[:] = 2.
        elif varname in ['uas','vas']:
            height[:] = 10.

    if varname in ['zg850','zg500','zg200','mpuhusa','mpvhusa']:
        plev = nc.createVariable('plev', 'd')
        plev.long_name = "pressure"
        plev.standard_name = "air_pressure"
        plev.units = "Pa"
        plev.positive = "down"
        plev.axis = "Z"
        if varname in ['zg850','mpuhusa','mpvhusa']:
            plev[:] = 85000.
        elif varname in ['zg500',]:
            plev[:] = 50000.
        elif varname in ['zg200',]:
            plev[:] = 20000.


def extract_times(mm5file):
    """
    Function to return an array of mx.DateTime instances for the timestamps
    in this MM5 file
    """
    mm5 = mm5_class.mm5(mm5file)
    # Sample for first timestamp with variable u
    # Requires a modification to mm5_class.py to read header
    tstr = mm5.get_field_list()['u']['time'][:13]
    ts0 = mx.DateTime.strptime(tstr, "%Y-%m-%d_%H")
    interval = mm5.timeincrement
    ts1 = ts0 + mx.DateTime.RelativeDateTime(
                seconds=(interval * (mm5.tsteps-1) ))
    taxis = []
    now = ts0
    while (now <= ts1):
        taxis.append( now )
        now += mx.DateTime.RelativeDateTime(seconds=interval)
    del(mm5)
    return taxis

def interpb_step():
    os.chdir(DATADIR)
    # Figure out a list of MMOUT files
    files = glob.glob("MMOUT_DOMAIN1_??")
    files.sort()
    # Move us to interpb
    os.chdir("%s/INTERPB" % (BASEDIR,))
    for file in files:
        # We don't wish to convert this file, it is just along for the ride
        if file == "MMOUT_DOMAIN1_00":
            continue
        # Figure out time axis
        taxis = extract_times(DATADIR+file)

        # Setup variable substitution values
        vars = {}
        vars['mm5file'] = DATADIR+file
        vars['syear'] = taxis[0].year
        vars['smonth'] = taxis[0].month
        vars['sday'] = taxis[0].day
        vars['shour'] = taxis[0].hour
        vars['eyear'] = taxis[-1].year
        vars['emonth'] = taxis[-1].month
        vars['eday'] = taxis[-1].day
        vars['ehour'] = taxis[-1].hour

        # Edit the namelist.input for interb
        data = open('namelist.tpl', 'r').read()
        out = open('namelist.input', 'w')
        out.write( data % vars )
        out.close()

        # Run interb for each file
        print "Running INTERPB for %s [%s - %s]" % (file,
              taxis[0].strftime("%Y-%m-%d %H"), taxis[-1].strftime("%Y-%m-%d %H"))
        os.system("./interpb >& interpb.log")

        # Move output file to right location
        os.rename("MMOUTP_DOMAIN1", DATADIR + file.replace("UT", "UTP"))
        # Cleanup
        os.system("rm -rf FILE_*")

def convert_to_netcdf_step():
    """
    Conversion of MM5 Format files to NetCDF
    """
    # Change directory
    os.chdir( DATADIR )
    # Look for any MMOUTP and NCOUT files
    files = glob.glob("NCOUT_DOMAIN1_??")
    files = files + glob.glob("MMOUTP_DOMAIN1_??")
    files.sort()
    # Loop over the files
    for file in files:
        if file == "NCOUT_DOMAIN1_00" or file == "MMOUTP_DOMAIN1_00":
            continue
        # Figure out how many timesteps there are.
        mm5 = mm5_class.mm5(file)
        cmd = "archiver %s 0 %s" % (file, mm5.tsteps)
        print "Converting %s to NetCDF %s tsteps" % (file, mm5.tsteps)
        si,so = os.popen4( cmd )
        a = so.read() # Necessary to keep things blocking?
        if not os.path.isfile( file+".nc" ):
            print "FAIL!", file
            print a
            sys.exit()

def convert_to_ll_step():
    """
    Conversion of NetCDF files in lambert to latlon regular
    """
    # Change directory
    os.chdir( DATADIR )
    # Look for any MMOUTP and NCOUT NC files
    files = glob.glob("NCOUT_DOMAIN1_??.nc")
    files = files + glob.glob("MMOUTP_DOMAIN1_??.nc")
    files.sort()
    # Loop over the files
    for file in files:
        cmd = "proj2ll %s 49.125 24.75 0.375 -67 -124.75 0.375" % (file,)
        print "Convert %s to LatLon Coords" % (file,)
        si,so = os.popen4( cmd )
        a = so.read() # Necessary to keep things blocking?
        if not os.path.isfile( file.replace(".nc", "_interp.nc") ):
            print "FAIL!", file
            print a
            sys.exit()
        # Remove intermediate netcdf file
        os.remove( file )

def cleanup():
    """
    And finally, we need to cleanup and compress things down
    """
    os.chdir( DATADIR )
    # Make sure final resting place dir exists
    finaldir = "../../final/%s" % (run,)
    if not os.path.isdir( finaldir ):
        os.mkdir( finaldir )
    # Move output files to necessary location
    for file in glob.glob("*_IMM5_*.nc"):
        os.rename(file, finaldir+"/"+file)
    # Remove any netcdf files
    os.system("rm *.nc")
    # Remove Interpolated data
    os.system("rm MMOUTP*")

# Begin Processing! 

# Interpolate the MMOUT files from Sigma to Pressure Levels
# via MM5's INTERPB
interpb_step()

# Convert MM5 Format NCOUT and MMOUT files to NetCDF
convert_to_netcdf_step()

# Convert these NetCDF Files from Lambert to Lat/lon WGS84
convert_to_ll_step()

# Process the variables from the NetCDF files into 1 per variable
finalize_step()
# Special logic is needed for these variables
special_precip()
special_spechumidity()
special_runoff()
special_runoff2()
special_soilw()

# Finally Cleanup
cleanup()
