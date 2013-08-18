# Generate a land/sea mask for the MRED domain
import netCDF3
import mm5_class
import mx.DateTime

mm5 = mm5_class.mm5('TERRAIN_DOMAIN1')
land = mm5.get_field('landmask', 0)

# 1,143,208
data = land['values']

lats =  mm5.get_field('latitdot',0)['values']
lons = mm5.get_field('longidot',0)['values']

nc = netCDF3.Dataset('LANDSEA_IMM5.nc', 'w')
nc.institution   = "Iowa State University, Ames, IA, USA"
nc.source        = "MM5 (2009): atmosphere: MM5v3.6.3 non-hydrostatic, split-explicit; sea ice: Noah; land: Noah"
nc.project_id    = "MRED"
nc.table_id      = "Table 2"
nc.realization   = 1
nc.forcing_data  = "CFS01"

# Optional
nc.Conventions   = 'CF-1.0'
nc.contact       = "Daryl Herzmann, akrherz@iastate.edu, 515-294-5978"
nc.history       = "%s Generated" % (mx.DateTime.now().strftime("%d %B %Y"),)
nc.comment       = "Runs processed on derecho@ISU, output processed on mred@ISU"
nc.title         = "ISU MM5 model output prepared for MRED using CFS input"

nc.createDimension('y', 143)
nc.createDimension('x', 208)

lat = nc.createVariable('lat', 'd', ('y','x') )
lat.units = "degrees_north"
lat.long_name = "latitude"
lat.standard_name = "latitude"
lat.axis = "Y"
lat[:] = lats

lon = nc.createVariable('lon', 'd', ('y', 'x',) )
lon.units = "degrees_east"
lon.long_name = "longitude"
lon.standard_name = "longitude"
lon.axis = "X"
lon[:] = lons


lsea = nc.createVariable('landmask', 'd', ('y', 'x'))
lsea.long_name = "land mask"
lsea.standard_name = "land mask"
lsea[:] = data[0,:,:]

nc.close()
mm5.close()   