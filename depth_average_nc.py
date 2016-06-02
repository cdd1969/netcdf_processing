#!/usr/bin/env python
from __future__ import division
import netCDF4 as netcdf
import numpy as np
import click
import os.path
import glob


__author__ = 'nikolai chernikov, <nikolai.chernikov.gmail.com>'
__date__ = 'January 2016'
__info__ = """
For UI info run this script from command-line:
    $ python depth_average_nc.py --help

This code in script is grouped into two large parts:
    1) the actual logic that is behind processing (back-end): 2 functions:
            caclulate_relative_layer_thickness()
            create_depth_averaged_nc()
            copy_nc_var() - a very general procedure to copy nc-variable
    2) commadline user-interface (front-end): 1 function:
            run()

See doc strings of corresponding routines for more details. Be aware of
assumptions introduces within these routines.
"""


def caclulate_relative_layer_thickness(layer_depth, water_depth, include_time=False):
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

    Return:
    -------
        layer_relthickness (3D|4D numpy array):
            3D|4D-Array with (z, y, x)|(t, z, y, x) dimensions; represents relative layer thickness.
            Is dimensionless and is constant within simulation period
    '''
    # >>> Get dimensions
    z, y, x = layer_depth.shape
    t, y, x = water_depth.shape

    # >>> Check water_depth to be valid
    nonzero_wd_y, nonzero_wd_x = water_depth[0, :, :].nonzero()  # See ISSUE 3
    if len(nonzero_wd_y) < 1 or len(nonzero_wd_x) < 1:
        msg = u'{0}WARNING! All entries of the `water_depth` at the initial timestep are MASKED or ZERO. Will not be able to properly calculate relative layer thickness. Proceeding anyway.{0}'.format('\n'+'-'*50+'\n')
        click.echo(click.style(msg, fg='red', bold=True))
        #print msg


    # >>> Get optional mask
    if isinstance(layer_depth, np.ma.MaskedArray):
        mask = layer_depth.mask
        # >>> Allocate memory for `layear_thickness` array, initialize it. This array represents layer thickness at timestep t=0. Values are always positive
        layer_thickness = np.ma.array(np.empty((z, y, x), dtype=float), mask=mask)
        # >>> Allocate memory for `relative_thcikness` array, initialize it. This array represents relaitve layer thickness at timestep t=0 with respect to total water-depth. Values are always positive, dimensionless
        layer_relthickness = np.ma.array(np.empty((z, y, x), dtype=float), mask=mask)
    elif isinstance(water_depth, np.ma.MaskedArray):
        # add new axis at 0-index position
        # and repeat the array (z,y,x) `t` times along 0-index position axis
        mask = water_depth[0, :, :].mask
        mask = mask.reshape(1, y, x).repeat(z, 0)
        layer_thickness = np.ma.array(np.empty((z, y, x), dtype=float), mask=mask)
        layer_relthickness = np.ma.array(np.empty((z, y, x), dtype=float), mask=mask)
    else:
        layer_thickness = np.empty((z, y, x), dtype=float)
        layer_relthickness = np.empty((z, y, x), dtype=float)
    


    # >>> Layer thickness of the near-bottom layer
    layer_thickness[0, :, :] = (water_depth[0, :, :] - (-layer_depth[0, :, :]) ) * 2.0
    # >>> Layer thickness of the rest layers
    for k in xrange(1, z):
        #layer_thickness[k, :, :] = (water_depth[0, :, :] - (-layer_depth[k, :, :]) - layer_thickness[k-1, :, :]) * 2.0
        layer_thickness[k, :, :] = (water_depth[0, :, :] - (-layer_depth[k, :, :]) - np.sum(layer_thickness[0:k, :, :], axis=0)) * 2.0
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
        nc_out='out.nc',
        var_list=['concentration_of_SPM_in_water_001', 'concentration_of_SPM_in_water_002'],
        z_dimname='getmGrid3D_getm_3',
        layers=(),
        waterdepth_varname='water_depth_at_soil_surface',
        layerdepth_varname='getmGrid3D_getm_layer',
        coord_attr='coordinates',
        copy_vars=(),
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

        nc_out (Optional[str or None]):
            absolute name of the result necdf file to be created. By default,
            creates a file "out.nc" within current working dir. If `None`
            no newfile is created, but the result is appended to the file,
            which is given as `nc_in`.

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
            3D-Array with (z, y, x) dimensions;
            By default will use "getmGrid3D_getm_layer"

        coord_attr (Otional[str]):
            name of the attribute of a varibles from `var_list` that
            describes the coordinates. For example (cdl):
                spm001.coordinates = 'level lat lon';

        copy_vars (Optional[Tuple[str]]):
            list with names of the variables that will be copied from the input
            netcdf file to the output without processing. This is done at the very
            end of the processing script.
            Default: empty tuple => no additional variables will be copied.

        log (Optional[bool]):
            flag to print additional info in console, while processing

    Return:
    -------
        the new file with processed results is created under `nc_out` path

    '''
    savedData = list()
    vnames = list()

    # >>> If `nc_out` has been passed, create a new file (if not - append)
    if nc_out:

        nc    = netcdf.Dataset(nc_in , mode='r')
        ncout = netcdf.Dataset(nc_out, mode='w', format='NETCDF4_CLASSIC')
        ncout.history = 'generated by script <{0}>'.format(__file__)
    # >>> If `nc_out` has not been passed, append to existing file
    else:
        ncout = netcdf.Dataset(nc_in , mode='a')
        nc = ncout

    #>>> Get the dimensions of the variable of interest...
    original_dims_names = [str(d) for d in nc.variables[var_list[0]].dimensions]
    original_dims_sizes = list(nc.variables[var_list[0]].shape)
    z_dim_index = original_dims_names.index(z_dimname)

    #>>> Get the dimensions of the averaged-variable of interest...
    averaged_dims_names = [d for d in original_dims_names if d != z_dimname]
    averaged_dims_sizes = tuple([nc.dimensions[d].__len__() for d in averaged_dims_names])

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
    if log: print u'Calculating relative layer thickness array of shape {0}, where {1} is the initial number of z-layers'.format(layer_relthickness.shape, layer_relthickness.shape[1])
    
    selected_layer_relthickness = np.take(layer_relthickness, np.arange( l1, l2+1, 1), axis=z_dim_index)
    if log: print u'Selected (with z-layers {1}) relative layer thickness array has shape {0}'.format(selected_layer_relthickness.shape, np.arange(l1, l2+1, 1))
    
    nonzero_cells = selected_layer_relthickness[0, 0, :, :].nonzero()  # will return a 2D array of (j, i) indices of an unmasked non-zero cells
    if log: print u'Found non-zero cells: {0}'.format(nonzero_cells)
    if len(nonzero_cells[0]) > 0 and len(nonzero_cells[1]) > 0:  #length is not zero. See ISSUE #3
        valid_cell_ji = (nonzero_cells[0][0], nonzero_cells[1][0])  # see ISSUE #2
    else:
        valid_cell_ji = (0, 0)
    if log:
        print u'Relative layer thickness of all z-layers at valid cell {0} (y, x) is:'.format(valid_cell_ji)
        for l_ in xrange(layer_relthickness.shape[1]):
            print u'\t layer {0} >>> {1}'.format(l_, layer_relthickness[0, l_, valid_cell_ji[0], valid_cell_ji[1]])
    
    rel_thick_factor = 1. / selected_layer_relthickness[0, :, valid_cell_ji[0], valid_cell_ji[1]].sum()  # see ISSUE #1 , #2
    selected_layer_relthickness = selected_layer_relthickness * rel_thick_factor
    if log:
        print u'Relative layer thickness of the selected {0} z-layers (considering multiplication factor {2}) is: {1}'.format(np.arange(l1, l2+1, 1), selected_layer_relthickness[0, :, valid_cell_ji[0], valid_cell_ji[1]], rel_thick_factor)

    # >>> Continue with variables of interest
    if log: print u'Reading file: {2}. Calculating depth averaged data for layer range {0}:{1}'.format(l1, l2, nc_in)
    
    for v in var_list:
        if log: print u'Working with variable `{0}`'.format(v)
        var = nc.variables[v]

        # >>> If `nc_out` has been passed, copy origina var and coord-vars
        if nc_out:
            #>>> Copy original variable
            if log: print u'\tCopying original variable `{0}` and dependencies...'.format(v)
            copy_nc_var(nc, ncout, v, coord_attr=coord_attr, log=log, indent='\t\t')
            
        # create depth averaging
        if log: print u'\tProcessing depth_averaging...'
        data = var[:]
        if log: print '\t\toriginal data shape:', data.shape
        selected_data = np.take(data, np.arange( l1, l2+1, 1), axis=z_dim_index)
        if log: print '\t\tselected data shape:', selected_data.shape
        # careful here!
        #   selected_layer_relthickness.shape = (time, z-selected, y, x)
        #   selected_data.shape = (time, z-selected, y, x)
        averaged_data = np.sum(selected_layer_relthickness * selected_data, axis=z_dim_index)

        if log: print '\t\taveraged data shape:', averaged_data.shape
        if log: print '\t\tdepth averaging >>> ok'

        name  = v+'_averaged'
        units = var.units if 'units' in var.ncattrs() else 'unknown'
        fv = var._FillValue if '_FillValue' in nc.variables[v].ncattrs() else None
        if log: print '\t\tCreating variable:', name
        newvar = ncout.createVariable(name, var.datatype, dimensions=averaged_dims_names, fill_value=fv)
        newvar.setncattr('units', units)
        newvar.setncattr('original_var_name', v)
        newvar.setncattr('layers_averaged', range(l1, l2+1))
        newvar.setncattr('averaged_along', z_dimname)
        newvar.setncattr('averaged_along', z_dimname)
        newvar.setncattr('info', 'generated by averaging values from <original_var_name> within layers <layers_depth_averaged>')
        if coord_attr in var.ncattrs():
            new_coords = var.getncattr(coord_attr).split()[1::]  # get coordinate attribute and get rid of the first dimension (z)
            newvar.setncattr(coord_attr, ' '.join(new_coords))

        newvar[:] = averaged_data
        savedData.append(averaged_data)
        vnames.append(name)

    if len(var_list) > 1:
        if log: print 'Creating variable: `SUM_averaged`; from :', var_list
        newvar = ncout.createVariable(u'SUM_averaged', var.datatype, dimensions=averaged_dims_names, fill_value=fv)
        newvar.setncattr('units', units)
        newvar.setncattr('vars_to_sum', ' '.join(var_list))
        newvar.setncattr('averaged_along', z_dimname)
        newvar.setncattr('layers_averaged', range(l1, l2+1))
        newvar.setncattr('info', 'this variable has been generated by summing data from variables <vars_to_sum>')
        if coord_attr in var.ncattrs():
            newvar.setncattr(coord_attr, ' '.join(new_coords))
        newvar[:] = sum(savedData)
        if isinstance(newvar[:], np.ma.MaskedArray):
            newvar[:].set_fill_value(fv)


    # >>> Save Relative thickness into netcdf file
    if log: print 'Creating variable:', '`layer_relative_thickness`'
    newvar = ncout.createVariable('layer_relative_thickness', float, dimensions=nc.variables[layerdepth_varname].dimensions)
    newvar.setncattr('units', '')
    newvar.setncattr('long_name', 'Dimensionless initial relative thickness of the layer')
    newvar.setncattr('info', 'Variable is generated automatically with script `'+__file__+'`')
    if coord_attr in var.ncattrs():
        newvar.setncattr(coord_attr, var.getncattr(coord_attr))
    newvar[:] = layer_relthickness[0, ...]

    # >>> Finally copy the additional variables
    if copy_vars and nc_out:
        if log: print u'Copying additional variables from the list {0} ...'.format(copy_vars)
        for copy_varname in copy_vars:
            copy_nc_var(nc, ncout, copy_varname, coord_attr=coord_attr, log=log, indent='')
            
    if nc_out:
        nc.close()
    ncout.close()



def copy_nc_var(nc_in, nc_out, varname, coord_attr='coordinates', log=False, indent='\t'):
    ''' Procedure. Copy variable from one netcdf file to another, including
    necessary dimensions and so-called "coordinate-variables". Coordinate
    variables are variables that meet following conditions:
        - have the same name as dimensions of `varname` variable
        - are declared in `coord_attr` attribute of the `varname` variable

    Args:
    -----
        nc_in (netCDF4.Dataset):
            opened instance of netcdf (copy from there)

        nc_out (netCDF4.Dataset):
            opened instance of netcdf (copy there)
        
        varname (str):
            variable name to copy

        coord_attr (str):
            name of the attribute of the variable, where
            the information about its coordinates is strored.

            i.e
            my_var.getncattr(coord_attr) = "level lat lon"

        log (bool):
            flag to print out more info

        indent (str):
            indent for nice output when log=True
    
    Return:
    -------
        This is a procedure, it does not return anything. It copies
        a variable `varname` from `nc_in` to `nc_out`
    '''
    if varname not in nc_in.variables.keys():
        if log: print u'Can not copy variable `{0}` from {1} to {2}. Variable not found in the input netcdf file'.format(varname, nc_in.filepath(), nc_out.filepath())
        return

    if varname in nc_out.variables.keys():
        if log: print u'Can not copy variable `{0}` from {1} to {2}. Variable already exists in the output netcdf file'.format(varname, nc_in.filepath(), nc_out.filepath())
        return

    var = nc_in.variables[varname]

    # add requred dimensions
    for d_name in var.dimensions:
        if d_name not in nc_out.dimensions.keys():
            nc_out.createDimension(d_name, size=len(nc_in.dimensions[d_name]))
    
    # add variable, copy attributes, copy data
    if log: print indent+u'Copying variable `{0}`'.format(varname)
    var_copy = nc_out.createVariable(varname,
        nc_in.variables[varname].datatype,
        dimensions=nc_in.variables[varname].dimensions,
        fill_value=nc_in.variables[varname]._FillValue if '_FillValue' in nc_in.variables[varname].ncattrs() else None)
    for attr_n in var.ncattrs():
        var_copy.setncattr(attr_n, var.getncattr(attr_n))
    var_copy[:] = var[:]
   
    # search for additional variables
    if log: print indent+'Searching for possible coordinate-vars of variable `{0}`'.format(varname)
    additional_var_list = list(var.dimensions)
    if coord_attr in var.ncattrs():
        coords = var.getncattr(coord_attr)
        if isinstance(coords, (str, unicode)):
            additional_var_list += coords.split()

    # now add additional variables
    for v_name in additional_var_list:
        if log: print indent+u'\tFound `{1}` coordinate-var of variable `{0}`'.format(varname, v_name),
        if v_name not in nc_in.variables.keys():
            if log: print u'... skipping (not found in `nc_in`)'
            continue
        if v_name in nc_out.variables.keys():
            if log: print u'... skipping (already exist)'
            continue
        if log: print u'... adding'
        # entering recursion
        copy_nc_var(nc_in, nc_out, v_name, coord_attr=coord_attr, log=log, indent=indent+'\t')












# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------- ALL LOGIC IS CODED ABOVE ---------------------------------------------------
# ---------------------------- BELOW IS THE UI CODE -----------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------







@click.command()
#@click.argument('nc_in', type=click.Path(exists=True, dir_okay=False), metavar='nc_in'
@click.argument('nc_in', nargs=-1
    )
@click.option('--nc_out', '-o', type=click.Path(exists=False, dir_okay=False), default='out.nc',
                help='Name (or basename if more than one input file is given) of the output netcdf file(-s) with results. Default: `out.nc`'
    )
@click.option('-a', '--append', is_flag=True, default=False,
        help='Flag to append result to the existing file `nc_in` If this option is used `--nc_out` is ignored.'
    )
@click.option('--layers', '-l', type=click.IntRange(min=0), default=None, nargs=2,
    help='Two integers, indicating the indexes (0-indexed) of layers to be averaged for the given `z_dimname` axis. This is useful to do averaging within certain layers (i.e 3 near-bed layers or 5 top-layers). Default: all layers of given `z_dimname` will be considered'
    )
@click.option('--varname', '-v', type=click.STRING, multiple=True,
            default=('concentration_of_SPM_in_water_001', 'concentration_of_SPM_in_water_002'),
    help='Name of the variable within `nc_in` to be processed. Multiple variables can be passed with additional `-v` prefix for each. For example `-v name1 -v name2`. Default: `-v concentration_of_SPM_in_water_001 -v concentration_of_SPM_in_water_002`'
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
@click.option('--coord_attr', '--ca', 'coord_attr', type=click.STRING, default='coordinates',
    help='Name of the coordinate attribute of netcdf variables given with `-v`. Default `coordinates`'
    )
@click.option('--copy_vars', '-c', type=click.STRING, multiple=True,
            default=(),
    help='Copy variable(s) from the input to the output file without processing. Pass the name(s) of the variable(s) to be copied. Multiple variables can be passed with additional `-c` prefix for each. For example `-c name1 -c name2`.'
    )
@click.option('--verbose', is_flag=True, default=False,
    help='Flag to print additional output during processing.'
    )
def run(nc_in, nc_out, append, layers, varname, z_dimname, waterdepth_varname, layerdepth_varname, coord_attr, copy_vars, verbose):
    ''' INFO: Reads variable(-s) from netcdf file(-s), and performs depth-averaging.
    Results are stored within newly created netcdf file(-s) or appended to the input file(-s).
    Additionally can copy variables without processing.
    
    DESCRIPTION: Depth averaging is done in 3 steps:

      1)  calculate RLT - dimensionless relative layer thickness.

      2)  multiply data-value at given cell at given timestep by RLT.
    
      3)  summ the mutliplication results along given z-axis
    
    EXAMPLES: Lets assume we have netcdf file <gfsen.nc> and we want to average spm-data
    over the depth. Within this file we have two spm variables <spm_c1> and <spm_c2> both of
    them are 4D with dimensions (time, z_layer, y, x) of shape (20, 10, 50, 100). The file also
    contains waterdepth information stored in 3D variable <wd> with dimensions (time, y, x) and
    the layerdepth variable <ld> with dimensions (z_layer, y, x).
    
    <<< Problem 1:
    Generate <z_layer>-averaged <spm_c1>, averaged over all 10 layers, and store output in file <out1.nc>

    >>> Solution 1:
    $ python depth_average_nc.py gfsen.nc -o out1.nc -v spm_c1 -z z_layer --wd wd --lv ld

    <<< Problem 2:
    Generate <z_layer>-averaged <spm_c1>, <spm_c2> averaged over layers [2, 3, 4, 5]. Store output in file <out2.nc>
    
    >>> Solution 2:
    $ python depth_average_nc.py gfsen.nc -o out2.nc -v spm_c1 -v spm_c2 -l 2 5 -z z_layer --wv wd --lv ld
    
    <<< Problem 3:
    Append existing file same data as in "Problem 2"
    
    >>> Solution 3:
    $ python depth_average_nc.py gfsen.nc -a -v spm_c1 -v spm_c2 -l 2 5 -z z_layer --wv wd --lv ld

    <<< Problem 4:
    Solve "Problem 2" for all nectdf file in current directory. Store output in directory "output" under name "out.001.nc", "out.002.nc", etc.
    
    >>> Solution 4:
    $ python depth_average_nc.py *.nc -o output/out.nc -v spm_c1 -v spm_c2 -l 2 5 -z z_layer --wv wd --lv ld
    '''


    for index, fname_in in enumerate(nc_in):
        # generate nc_out names
        if len(nc_in) > 1:
            name, ext = os.path.splitext(nc_out)
            fname_out = name+'.{0:02d}'.format(index)+ext
        else:
            fname_out = nc_out
        try:
            nc    = netcdf.Dataset(fname_in , mode='r')
        except Exception, err:
            raise click.BadParameter('( {1} ) Can not read NetCDF file {0}'.format(fname_in, err), param_hint=['nc_in'])
        if waterdepth_varname not in nc.variables.keys():
            raise click.BadParameter('Variable `{1}` does not exist in file {0}'.format(fname_in, waterdepth_varname), param_hint=['--waterdepth_varname'])
        if layerdepth_varname not in nc.variables.keys():
            raise click.BadParameter('Variable `{1}` does not exist in file {0}'.format(fname_in, layerdepth_varname), param_hint=['--layerdepth_varname'])
        for v in varname:
            if v not in nc.variables.keys():
                raise click.BadParameter('Variable `{1}` does not exist in file {0}'.format(fname_in, v), param_hint=['--varname'])
            if z_dimname not in nc.variables[v].dimensions:
                raise click.BadParameter('Dimension `{1}` doesnot exist in variable {v}'.format(fname_in, z_dimname), param_hint=['--z_dimname'])
        if layers:
            if layers[1]+1 > nc.dimensions[z_dimname]:
                msg = 'Layer range-index `{0}` exceeds maximum size {1} of dimension `{2}`'.format(layers[1], nc.dimensions[z_dimname], z_dimname)
                raise click.BadParameter(msg, param_hint=['--layers'])
            
            if layers[0] > layers[1]:
                msg = 'Layer lower range-index `{0}` is greater than upper range-index `{1}`'.format(layers[0], layers[1])
                raise click.BadParameter(msg, param_hint=['--layers'])
        nc.close()

        click.echo(click.style('Processing: <{0}>'.format(fname_in), fg='yellow', bold=True))
        if append is False:
            if os.path.exists(fname_out):
                click.confirm(click.style('File `{0}` already exists. Do you want to overwrite it?'.format(fname_out), fg='red'), abort=True)
        else:
            fname_out = None
        create_depth_averaged_nc(fname_in,
            nc_out=fname_out,
            var_list=varname,
            z_dimname=z_dimname,
            layers=layers,
            waterdepth_varname=waterdepth_varname,
            layerdepth_varname=layerdepth_varname,
            coord_attr=coord_attr,
            copy_vars=copy_vars,
            log=verbose)
        
        if append:
            click.echo(click.style('Finished: <{0}> appended successfully.'.format(fname_in), fg='green', bold=True))
        else:
            click.echo(click.style('Finished: <{0}> created successfully.'.format(fname_out), fg='green', bold=True))

if __name__ == '__main__':
    run()
