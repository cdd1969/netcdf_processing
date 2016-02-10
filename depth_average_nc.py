#!/usr/bin/env python
import netCDF4 as netcdf
import numpy as np
import click
import os.path

__author__ = 'nikolai chernikov, <nikolai.chernikov.gmail.com>'
__date__ = 'January 2016'
__info__ = """
For UI info run this script from command-line:
    $ python depth_average_nc.py --help

This code in script is grouped into two large parts:
    1) the actual logic that is behind processing (back-end): 2 functions:
            caclulate_relative_layer_thickness()
            create_depth_averaged_nc()
    2) commadline user-interface (front-end): 1 function:
            run()

See doc strings of corresponding routines for more details. Be aware of
assumptions introduces within these routines.
"""


def caclulate_relative_layer_thickness(layer_depth, water_depth, include_time=False, soil_surface='first'):
    '''Calculate relative layer thickness for mutlidimensional array
    (t, z, y, x) at timestep t=0

    Assumptions:
        1) z=0 is the closest to soil-surface layer, where `z` is the index
           of the vertical dimension of the matrix

        2) relative layer thickness is not time dependent !!!

    Args:
    -----
        layer_depth (3D numpy array):
            3D-Array with (z, y, x) dimensions; represents layer depth below
            water surface (0, level) at center of each z,y,x cell.
            Units are shared with `water_depth` array. Always negative.
        
        water_depth (3D numpy array):
            3D-Array with (time, y, x) dimensions; represents water depth at soil
            surface. Units are shared with `layer_depth` array. Always positive.

        include_time (Optional[bool]):
            Flag to control the shape of the output array. It does not affect the data
            that is stored in output array: it is completely identical.
            ---
            If `True` - the shape will be (time, z, y, x)
            If `False` - the shape will be (z, y, x)  (DEFAULT)

        soil_surface (Optional['first'|'last']):
            String flag that tells which `layer` should be considered closest to soil
            surface. The `layer` is the index of the z-dimension in `layer_depth` array.
            ---
            If 'first' - z=0 is considered to be closest to soil surface layer, i.e.
                         layer_depth(t, 0, y, x) - is the near-bottom layer at t,y,x
            If 'last'  - z=-1 is considered to be closest to soil surface layer, i.e.
                         layer_depth(t, -1, y, x) - is the near-bottom layer at t,y,x

    Return:
    -------
        layer_relthickness (3D|4D numpy array):
            3D|4D-Array with (z, y, x)|(t, z, y, x) dimensions; represents relative layer thickness.
            Is dimensionless and is constant within simulation period
    '''
    # >>> Get dimensions
    z, y, x = layer_depth.shape
    t, y, x = water_depth.shape
    # >>> Allocate memory for `layear_thickness` array, initialize it. This array represents layer thickness at timestep t=0. Values are always positive
    layer_thickness = np.empty((z, y, x), dtype=float)
    # >>> Allocate memory for `relative_thcikness` array, initialize it. This array represents relaitve layer thickness at timestep t=0 with respect to total water-depth. Values are always positive, dimensionless
    layer_relthickness = np.empty((z, y, x), dtype=float)
    
    if soil_surface == 'first':
        # >>> Layer thickness of the near-bottom layer
        layer_thickness[0, :, :] = (water_depth[0, :, :] - (-layer_depth[0, :, :]) ) * 2.0
        # >>> Layer thickness of the rest layers
        for k in xrange(1, z):
            layer_thickness[k, :, :] = (water_depth[0, :, :] - (-layer_depth[k, :, :]) - layer_thickness[k-1, :, :]) * 2.0
        # >>> Now calculate relative layer thickness
        for k in xrange(z):
            layer_relthickness[k, :, :] = abs(layer_thickness[k, :, :] / water_depth[0, :, :])
    
    # >>> Finally return the result
    if include_time is False:
        return layer_relthickness
    else:
        # add new axis at 0-index position
        # and repeat the array (z,y,x) `t` times along 0-index position axis
        return layer_relthickness.reshape(1, z, y, x).repeat(t, 0)



def create_depth_averaged_nc(nc_in,
        nc_out='depth_averaged_spm.nc',
        var_list=['concentration_of_SPM_in_water_001', 'concentration_of_SPM_in_water_002'],
        z_dimname='getmGrid3D_getm_3',
        layers=(),
        waterdepth_varname='water_depth_at_soil_surface',
        layerdepth_varname='getmGrid3D_getm_layer',
        log=True):
    ''' Function reads variable (or list of variables) from netcdf datafile,
    and performs depth-averaging. It is done in two steps:
        1)  multiplying data-value at given cell at given timestep
            by dimensionless factor - relative thickness of corresponding
            layer.
        2)  summing the mutliplication result along given z-axis

    Creates new netcdf file and stores the results there.

    Args:
    -----
        nc_in (str):
            absolute name of the netcdf datafile to be read

        nc_out (Optional[str]):
            absolute name of the result necdf file to be created. By default,
            creates a file "depth_averaged_spm.nc" within current working dir

        var_list (Optional[list(str)]):
            list of the names of the variables within file `nc_in` which will
            be processed. The averaging will be performed independently to
            each of the variables in this list; The result will be saved in
            the output file `nc_out` under name "varname_averaged", where
            "varname" is the name from this parameter `var_list`.
            If this list contains more than one variable, the sum of these
            averaged values will be saved as new variable with name
            "SUM_averaged"

        z_dimname (Optional[str]):
            The name of the dimension, along which the depth-averaging will be
            performed. By default will set this param to "getmGrid3D_getm_3".

        layers (Optional[list(int)]):
            List of two integers, indicating the indexes (0-indexed) of the
            layers to be averaged. This is useful to do averaging only within
            certain part of layers. By default `layers` is empty tuple (), meaning
            that all layers will be taken into account. Consider example below:
            -----------------------------------------------------------------
                We have an 3d array with velocity magnitudes (x, y, z),
                with actual shape of (10, 10, 5).
                We want to get averaged values along z-axis.
                Task 1: get averaging considering all z-layers
                    `layers` = () (default)
                Task 2: get averaging only of two z-layers indexed 0, 1
                    `layers` = [0, 1]
                Task 3: get averaging of three middle z-layers indexed 2, 3, 4
                    `layers` = [2, 4]
        
        waterdepth_varname (Optional[str]):
            name of the variable that represents water depth at soil surface.
            Units must be shared with `layerdepth_varname`. Always positive.
            3D-Array with (time, y, x) dimensions;
            By default will use "water_depth_at_soil_surface"
        
        layerdepth_varname (Optional[str]):
            name of the variable in `nc_in` netcdf file, that represents
            layer depth below water surface at the element center. Units
            must be shared with `waterdepth_varname`. Always negative.
            4D-Array with (time, z, y, x) dimensions;
            By default will use "getmGrid3D_getm_layer"

        log (Optional[bool]):
            flag to print additional info in console, while processing

    Return:
    -------
        the new file with processed results is created under `nc_out` path

    '''
    savedData = list()
    vnames = list()

    nc    = netcdf.Dataset(nc_in , mode='r')
    ncout = netcdf.Dataset(nc_out, mode='w', format='NETCDF4_CLASSIC')
    ncout.history = 'generated by script <{0}>'.format(__name__)


    #>>> Get the dimensions of the variable of interest...
    original_dims_names = [str(d) for d in nc.variables[var_list[0]].dimensions]
    original_dims_sizes = list(nc.variables[var_list[0]].shape)
    z_dim_index = original_dims_names.index(z_dimname)

    #>>> Get the dimensions of the averaged-variable of interest...
    averaged_dims_names = [d for d in original_dims_names if d != z_dimname]
    averaged_dims_sizes = tuple([nc.dimensions[d].__len__() for d in averaged_dims_names])

    #>>> Create them in netcdf out-file
    for dim_n, dim_s in zip(original_dims_names, original_dims_sizes):
        if dim_n not in ncout.dimensions.keys():
	    if log: print 'Creating dimension:', dim_n, 'of size', dim_s
            ncout.createDimension(dim_n, size=dim_s)

    #>>> Decide which layers I'd like to consider for averaging
    if not layers:
        # if empty list will integrate all layers
        l1 = 0
        l2 = original_dims_sizes[z_dim_index]-1
    else:
        # or will integrate within given layer-range
        l1, l2 = layers[0], layers[1]


    #>>> Now get the relative layer thickness
    layer_relthickness = caclulate_relative_layer_thickness(nc.variables[layerdepth_varname][:], nc.variables[waterdepth_varname][:], include_time=True)

    selected_layer_relthickness = np.take(layer_relthickness, np.arange( l1, l2+1, 1), axis=z_dim_index)
    
    # >>> Save into netcdf file
    if log: print 'creating variable:', 'layer_relative_thickness'
    newvar = ncout.createVariable('layer_relative_thickness', float, dimensions=nc.variables[layerdepth_varname].dimensions)
    newvar.setncattr('units', 'dimensionless')
    newvar.setncattr('info', 'Variable is generated automatically during script execution. See function `caclulate_relative_layer_thickness()` in script `'+__name__+'`')
    newvar[:] = layer_relthickness[0, ...]

    # >>> Continue with variables of interest
    if log: print u'Reading file: {2}. Calculating depth averaged data for layer range {0}:{1}'.format(l1, l2, nc_in)
    for v in var_list:
        if log: print u'Working with variable `{0}`'.format(v)
        var = nc.variables[v]

        name  = v+'_averaged'
        dType = var.datatype
        fv    = var._FillValue
        units = var.units if 'units' in var.ncattrs() else 'unknown'

        #>>> Copy original variable
        original_var = ncout.createVariable(v, dType, dimensions=original_dims_names, fill_value=fv)
        for attr_n in var.ncattrs():
            original_var.setncattr(attr_n, var.getncattr(attr_n))
        original_var[:] = var[:]

        #>>> Now make sure that coordinate-variables are saved
        #    1) It could be the name of the dimension
        #    2) It could be stored within `coordinates` attribute
        possible_coord_var_list = list(var.dimensions)
        if 'coordinates' in var.ncattrs():
            coords = var.coordinates
            if isinstance(coords, (list, tuple)):
                possible_coord_var_list += list(coords)
            elif isinstance(coords, (str, unicode)):
                possible_coord_var_list += coords.split()
        
        #>>> Now check if these coord-vars already exist. If not - save them!
        for dim_name in possible_coord_var_list:
            if dim_name in nc.variables.keys() and dim_name not in ncout.variables.keys():
                #>>> if conditions are met > copy our coordinate-variable
                if log: print '\tcopying coordinate variable:', dim_name
                try:
                    fv = nc.variables[dim_name]._FillValue
                    coord_var = ncout.createVariable(dim_name, nc.variables[dim_name].datatype, dimensions=nc.variables[dim_name].dimensions, fill_value=fv)
                except:
                    coord_var = ncout.createVariable(dim_name, nc.variables[dim_name].datatype, dimensions=nc.variables[dim_name].dimensions)

                for attr_n in nc.variables[dim_name].ncattrs():
                    coord_var.setncattr(attr_n, nc.variables[dim_name].getncattr(attr_n))
                coord_var[:] = nc.variables[dim_name][:]


        # create depth averaging
        data = var[:]
        if log: print '\toriginal data shape:', data.shape
        selected_data = np.take(data, np.arange( l1, l2+1, 1), axis=z_dim_index)
        
        if log: print '\tselected data shape:', selected_data.shape
        # careful here!
        #   selected_layer_relthickness.shape = (time, z-selected, y, x)
        #   selected_data.shape = (time, z-selected, y, x)
        averaged_data = np.sum(selected_layer_relthickness * selected_data, axis=z_dim_index)

        if log: print '\taveraged data shape:', averaged_data.shape
        if log: print '\tdepth averaging >>> ok'
        


        if log: print '\tcreating variable:', name
        newvar = ncout.createVariable(name, dType, dimensions=averaged_dims_names, fill_value=fv)
        newvar.setncattr('units', units)
        newvar.setncattr('original_var_name', v)
        newvar.setncattr('layers_averaged', [l1, l2])
        newvar.setncattr('averaged_along', z_dimname)
        newvar.setncattr('info', 'this variable has been generated by averaging values from <original_var_name> within layers <layers_depth_averaged> including both')

        newvar[:] = averaged_data
        savedData.append(averaged_data)
        vnames.append(name)

    if len(var_list) > 1:
        if log: print 'SUM_averaged:', name
        newvar = ncout.createVariable(u'SUM_averaged', dType, dimensions=averaged_dims_names, fill_value=var._FillValue)
        newvar.setncattr('units', units)
        newvar.setncattr('vars_to_sum', v)
        newvar.setncattr('averaged_along', z_dimname)
        newvar.setncattr('layers_averaged', [l1, l2])
        newvar.setncattr('info', 'this variable has been generated by summing data from variables <vars_to_sum>')
        newvar[:] = np.sum(savedData)

    nc.close()
    ncout.close()
    #print 'Finished: <{0}> created successfully.'.format(nc_out)














# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------- ALL LOGIC IS CODED ABOVE ---------------------------------------------------
# ---------------------------- BELOW IS THE UI CODE -----------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------







@click.command(short_help='putin')
@click.argument('nc_in', type=click.Path(exists=True, dir_okay=False), metavar='nc_in'
    )
@click.option('--nc_out', '-o', type=click.Path(exists=False, dir_okay=False), default='depth_averaged.nc',
                help='Name of the output netcdf file with results. Default: `depth_averaged.nc`'
    )
@click.option('--layers', '-l', type=click.IntRange(min=0), default=None, nargs=2,
    help='Two integers, indicating the indexes (0-indexed) of layers to be averaged for the given `z_dimname` axis. This is useful to do averaging within certain layers (i.e 3 near-bed layers or 5 top-layers). Default: all layers of given `z_dimname` will be considered'
    )
@click.option('--varname', '-v', type=click.STRING, multiple=True,
            default=('concentration_of_SPM_in_water_001', 'concentration_of_SPM_in_water_002'),
    help='Name of the variable within `nc_in` to be processed. Multiple variables can be passed with additional `-v` prefix for each. In example `-v name1 -v name2`. Default: `-v concentration_of_SPM_in_water_001 -v concentration_of_SPM_in_water_002`'
    )
@click.option('--z_dimname', '-z', type=click.STRING, default='getmGrid3D_getm_3',
    help='Name of the dimension, along which the depth-averaging will be performed. Default: `getmGrid3D_getm_3`'
    )
@click.option('--waterdepth_varname', '--wd', 'waterdepth_varname', type=click.STRING, default='water_depth_at_soil_surface',
    help='Name of the variable that represents water depth at soil surface. Units must be shared with `layerdepth_varname`. Always positive. 3D-Array with (time, y, x) dimensions. Default: `water_depth_at_soil_surface`'
    )
@click.option('--layerdepth_varname', '--ld', 'layerdepth_varname', type=click.STRING, default='getmGrid3D_getm_layer',
    help='Name of the variable that represents layer depth below water surface at the element center. Units must be shared with `waterdepth_varname`. Always negative. 4D-Array with (time, z, y, x) dimensions. Default `getmGrid3D_getm_layer`'
    )
@click.option('--verbose', is_flag=True, default=False,
    help='Flag to print additional output during processing.'
    )
def run(nc_in, nc_out, layers, varname, z_dimname, waterdepth_varname, layerdepth_varname, verbose):
    ''' INFO: Reads variable (or list of variables) from netcdf file, and performs depth-averaging.
    Results are stored within newly created netcdf file.
    
    DESCRIPTION: Depth averaging is done in two steps:
     1)  calculate RLT - dimensionless relative layer thickness.

     2)  multiply data-value at given cell at given timestep by RLT.
    
     3)  summ the mutliplication results along given z-axis
    
    EXAMPLES: Lets assume we have netcdf file <gfsen.nc> and we want to average spm-data
    over the depth. Within this file we have two spm variables <spm_c1> and <spm_c2> both of
    them are 4D with dimensions (time, z_layer, y, x) of shape (20, 10, 50, 100). The file also
    contains waterdepth information stored in 3D variable <wd> with dimensions (time, y, x) and
    the layerdepth variable <ld> with dimensions (time, z_layer, y, x).
    
    ---Problem 1---

    Generate <z_layer>-averaged <spm_c1>, averaged over all 10 layers, and store output in file <out1.nc>

    ---Solution 1---
    
    $ python depth_average_nc.py gfsen.nc -o out1.nc -v spm_c1 -z z_layer --waterdepth_varname wd --layerdepth_varname ld

    ---Problem 2---

    Generate <z_layer>-averaged <spm_c1>, <spm_c2> averaged over layers [2, 3, 4, 5]. Store output in file <out2.nc>
    
    ---Solution 2---
    
    $ python depth_average_nc.py gfsen.nc -o out2.nc -v spm_c1 -v spm_c2 -l 2 5 -z z_layer --waterdepth_varname wd --layerdepth_varname ld
    '''
    try:
        nc    = netcdf.Dataset(nc_in , mode='r')
    except Exception, err:
        raise click.BadParameter('( {1} ) Can not read NetCDF file {0}'.format(nc_in, err), param_hint=['nc_in'])
    if waterdepth_varname not in nc.variables.keys():
        raise click.BadParameter('Variable `{1}` does not exist in file {0}'.format(nc_in, waterdepth_varname), param_hint=['--waterdepth_varname'])
    if layerdepth_varname not in nc.variables.keys():
        raise click.BadParameter('Variable `{1}` does not exist in file {0}'.format(nc_in, layerdepth_varname), param_hint=['--layerdepth_varname'])
    for v in varname:
        if v not in nc.variables.keys():
            raise click.BadParameter('Variable `{1}` does not exist in file {0}'.format(nc_in, v), param_hint=['--varname'])
        if z_dimname not in nc.variables[v].dimensions:
            raise click.BadParameter('Dimension `{1}` doesnot exist in variable {v}'.format(nc_in, z_dimname), param_hint=['--z_dimname'])
    if layers:
        if layers[1]+1 > nc.dimensions[z_dimname]:
            msg = 'Layer range-index `{0}` exceeds maximum size {1} of dimension `{2}`'.format(layers[1], nc.dimensions[z_dimname], z_dimname)
            raise click.BadParameter(msg, param_hint=['--layers'])
        
        if layers[0] > layers[1]:
            msg = 'Layer lower range-index `{0}` is greater than upper range-index `{1}`'.format(layers[0], layers[1])
            raise click.BadParameter(msg, param_hint=['--layers'])
    nc.close()

    if os.path.exists(nc_out):
        click.confirm('File `{0}` already exists. Do you want to overwrite it?'.format(nc_out), abort=True)

    create_depth_averaged_nc(nc_in,
        nc_out=nc_out,
        var_list=varname,
        z_dimname=z_dimname,
        layers=layers,
        waterdepth_varname=waterdepth_varname,
        layerdepth_varname=layerdepth_varname,
        log=verbose)

    click.echo(click.style('Finished: <{0}> created successfully.'.format(nc_out), fg='green', bold=True))

if __name__ == '__main__':
    run()
