from __future__ import division
import numpy as np
import netCDF4
from os import path as op


def create_test_file(log=False):

    # _________________________________________________________________________
    # put everything in netCDF file.

    # create netCDF file
    nc = netCDF4.Dataset('test.nc', mode='w')
    

    # create global attributes
    nc.history = 'Created by '+op.basename(__file__)
    nc.title    = 'Test file for script `depth_average_nc.py`'

    # create dimensions
    nc.createDimension('time', 3)
    nc.createDimension('getmGrid3D_getm_1', 10)
    nc.createDimension('getmGrid3D_getm_2', 15)
    nc.createDimension('getmGrid3D_getm_3', 5)
    nc.createDimension('getmGrid2D_getm_1', 10)
    nc.createDimension('getmGrid2D_getm_2', 15)


    VARS = list()

    var = dict()
    var['Name']                   = 'time'
    var['Nctype']                 = 'float'
    var['attrs']                  = dict()
    var['Dims']                   = ['time']
    var['attrs']['long_name']     = 'time'

    VARS.append(var)


    var = dict()
    var['Name']                   = 'water_depth_at_soil_surface'
    var['Nctype']                 = 'float'
    var['_FillValue']             = -1.e+30
    var['attrs']                  = dict()
    var['Dims']                   = ['time', 'getmGrid2D_getm_2', 'getmGrid2D_getm_1']
    var['attrs']['long_name']     = 'water_depth_at_soil_surface'
    var['attrs']['units']         = 'm'
    var['attrs']['missing_value'] = -1.e+30


    VARS.append(var)


    var = dict()
    var['Name']                   = 'concentration_of_SPM_in_water_001'
    var['Nctype']                 = 'float'
    var['_FillValue']             = -1.e+30
    var['attrs']                  = dict()
    var['Dims']                   = ['time', 'getmGrid3D_getm_3', 'getmGrid3D_getm_2', 'getmGrid3D_getm_1']
    var['attrs']['long_name']     = 'concentration_of_SPM_in_water_001'
    var['attrs']['units']         = 'mg/l'
    var['attrs']['missing_value'] = -1.e+30


    VARS.append(var)


    var = dict()
    var['Name']                   = 'concentration_of_SPM_in_water_002'
    var['Nctype']                 = 'float'
    var['_FillValue']             = -1.e+30
    var['attrs']                  = dict()
    var['Dims']                   = ['time', 'getmGrid3D_getm_3', 'getmGrid3D_getm_2', 'getmGrid3D_getm_1']
    var['attrs']['long_name']     = 'concentration_of_SPM_in_water_002'
    var['attrs']['units']         = 'mg/l'
    var['attrs']['missing_value'] = -1.e+30


    VARS.append(var)

    var = dict()
    var['Name']                   = 'getmGrid3D_getm_layer'
    var['Nctype']                 = 'float'
    var['_FillValue']             = -1.e+30
    var['attrs']                  = dict()
    var['Dims']                   = ['getmGrid3D_getm_3', 'getmGrid3D_getm_2', 'getmGrid3D_getm_1']
    var['attrs']['long_name']     = 'getmGrid3D_getm_layer'
    var['attrs']['info']          = 'Initial distance from water surface to the middle of the layer. Always negative'
    var['attrs']['units']         = 'mg'
    var['attrs']['missing_value'] = -1.e+30


    VARS.append(var)

    # create variables in netcdf file
    for var in VARS:
        fv = False  #default
        if 'attrs' in var.keys():
            if 'missing_value' in var['attrs'].keys():
                fv = var['attrs']['missing_value']
            elif '_FillValue' in var['attrs'].keys():
                fv = var['attrs']['_FillValue']

        dims = tuple()
        if 'Dims' in var.keys():
            dims = tuple(var['Dims'])

        if log:
            print '\tadd: var <{0}>, type <{1}>, dims <{2}>, fv <{3}>'.format(var['Name'], var['Nctype'], dims, fv)
        nc_var = nc.createVariable(var['Name'], var['Nctype'], dimensions=dims, fill_value=fv)
        
        if 'attrs' in var.keys():
            for attr_name, attr_val in var['attrs'].iteritems():
                nc_var.setncattr(attr_name, attr_val)

    

    # save values...
    nc.variables['time'][:] = np.array([3600*i for i in xrange(3)])

    nc.variables['water_depth_at_soil_surface'][:] = 5  # 5m depth everywhere

    nc.variables['concentration_of_SPM_in_water_001'][:, 0, :, :] = 10  # concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_001'][:, 1, :, :] = 9   # concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_001'][:, 2, :, :] = 8   # concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_001'][:, 3, :, :] = 7   # concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_001'][:, 4, :, :] = 6   # concentration at all t,y,x


    nc.variables['concentration_of_SPM_in_water_002'][:, 0, :, :] = 2*10  # concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_002'][:, 1, :, :] = 2*9   # concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_002'][:, 2, :, :] = 2*8   # concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_002'][:, 3, :, :] = 2*7   # concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_002'][:, 4, :, :] = 2*6   # concentration at all t,y,x

    nc.variables['getmGrid3D_getm_layer'][0, :, :] = -4.5   # depth of the layer middle in [m] downside negative
    nc.variables['getmGrid3D_getm_layer'][1, :, :] = -3.5   # depth of the layer middle in [m] downside negative
    nc.variables['getmGrid3D_getm_layer'][2, :, :] = -2.5   # depth of the layer middle in [m] downside negative
    nc.variables['getmGrid3D_getm_layer'][3, :, :] = -1.5   # depth of the layer middle in [m] downside negative
    nc.variables['getmGrid3D_getm_layer'][4, :, :] = -0.5   # depth of the layer middle in [m] downside negative


    # close the file
    nc.close()
    if log:
        print 'File created successfully'

if __name__ == '__main__':
    create_test_file(log=True)
