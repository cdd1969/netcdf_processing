from __future__ import division
import numpy as np
import netCDF4
from os import path as op


def create_test_file(fname='test.nc', layer_height='eq', log=False, testmode=False):
    '''
        Args:
        -----

        layer_height(str):
            param to control the thickness of vertical layers
            in this example we have 5 vertical layers

            'eq' - all 5 layers are 1m
            'noneq' - 0.1, 0.3, 0.6, 1.5, 2.5 m
    '''
    if layer_height == 'eq':
        predefined_layer_depth = [-4.5, -3.5, -2.5, -1.5, -0.5]
    elif layer_height == 'noneq':
        predefined_layer_depth = [-4.95, -4.75, -4.3, -3.25, -1.25]



    # _________________________________________________________________________
    # put everything in netCDF file.

    # create netCDF file
    nc = netCDF4.Dataset(fname, mode='w')
    

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
    var['attrs']['coordinates']   = 'lat lon'


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
    var['attrs']['coordinates']   = 'layer lat lon'


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

    nc.variables['concentration_of_SPM_in_water_001'][:, 0, :, :] = 10  # [mg/l]=[g/m3] concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_001'][:, 1, :, :] = 9   # [mg/l]=[g/m3] concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_001'][:, 2, :, :] = 8   # [mg/l]=[g/m3] concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_001'][:, 3, :, :] = 7   # [mg/l]=[g/m3] concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_001'][:, 4, :, :] = 6   # [mg/l]=[g/m3] concentration at all t,y,x


    nc.variables['concentration_of_SPM_in_water_002'][:, 0, :, :] = 2*10  # [mg/l]=[g/m3] concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_002'][:, 1, :, :] = 2*9   # [mg/l]=[g/m3] concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_002'][:, 2, :, :] = 2*8   # [mg/l]=[g/m3] concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_002'][:, 3, :, :] = 2*7   # [mg/l]=[g/m3] concentration at all t,y,x
    nc.variables['concentration_of_SPM_in_water_002'][:, 4, :, :] = 2*6   # [mg/l]=[g/m3] concentration at all t,y,x

    nc.variables['getmGrid3D_getm_layer'][0, :, :] = predefined_layer_depth[0]   # depth of the layer middle in [m] downside negative
    nc.variables['getmGrid3D_getm_layer'][1, :, :] = predefined_layer_depth[1]   # depth of the layer middle in [m] downside negative
    nc.variables['getmGrid3D_getm_layer'][2, :, :] = predefined_layer_depth[2]   # depth of the layer middle in [m] downside negative
    nc.variables['getmGrid3D_getm_layer'][3, :, :] = predefined_layer_depth[3]   # depth of the layer middle in [m] downside negative
    nc.variables['getmGrid3D_getm_layer'][4, :, :] = predefined_layer_depth[4]   # depth of the layer middle in [m] downside negative



    # the integration results should be....
    ''' for one cell (y,x) vertical column (5,1,1) (z-size,y-size,x-size)

    relative thickness:

        [0.2, 0.2, 0.2, 0.2, 0.2]

    mass_spm001 [g/m2]:
        [10*0.2*5, 9*0.2*5, 8*0.2*5, 7*0.2*5, 6*0.2*5]=
       =[10, 9, 8, 7, 6]
    
    mass_spm001_depth_ave [g/m2]:
        sum(10, 9, 8, 7, 6)/5 = 8

    mass_spm001 [g/m2]:
       [20, 18, 16, 14, 12]
    
    mass_spm002_depth_ave [g/m2]:
        sum(20, 18, 16, 14, 12)/5 = 16

    mass_summ_depth_ave [g/m2]:
        8+16 = 24

    '''
    if not testmode:
        # close the file
        nc.close()
        if log:
            print 'File created successfully'
    
    else:
        return nc



if __name__ == '__main__':
    create_test_file(log=True)
