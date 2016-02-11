''' This is the xecution file for built-in tests.
Run this file in package mode:

$ cd ../..
$ python -m netcdf_processing.test.run -v
'''

import unittest
from create_test_netcdf_file import create_test_file
from ..depth_average_nc import create_depth_averaged_nc
import os
import netCDF4


class DepthAverageScriptBaseTestCase(unittest.TestCase):
    '''Base class for all Script tests. '''

    def assertItemsAlmostEqual(self, list1, list2, tol):
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(a, b, tol)
    
    def basetest_relative_layer_thickness(self, var, correct_relth, tol=7):
        self.assertItemsAlmostEqual(var[:, 0, 0], correct_relth, tol)
        self.assertItemsAlmostEqual(var[:, 4, 2], correct_relth, tol)
        self.assertItemsAlmostEqual(var[:, 6, 8], correct_relth, tol)

    def basetest_depth_averaged_spm001(self, var, correct_val, tol=7):
        self.assertAlmostEqual(var[0, 0, 0], correct_val, tol)
        self.assertAlmostEqual(var[1, 4, 2], correct_val, tol)
        self.assertAlmostEqual(var[2, 6, 8], correct_val, tol)
        self.assertAlmostEqual(var[:].sum(), correct_val*var.shape[0]*var.shape[1]*var.shape[2], tol)
    
    def basetest_depth_averaged_spm002(self, var, correct_val, tol=7):
        self.assertAlmostEqual(var[0, 0, 0], correct_val, tol)
        self.assertAlmostEqual(var[1, 4, 2], correct_val, tol)
        self.assertAlmostEqual(var[2, 6, 8], correct_val, tol)
        self.assertAlmostEqual(var[:].sum(), correct_val*var.shape[0]*var.shape[1]*var.shape[2], tol)

    def basetest_depth_averaged_summ_spm001_spm002(self, var, correct_val, tol=7):
        self.assertAlmostEqual(var[0, 0, 0], correct_val, tol)
        self.assertAlmostEqual(var[1, 4, 2], correct_val, tol)
        self.assertAlmostEqual(var[2, 6, 8], correct_val, tol)
        self.assertAlmostEqual(var[:].sum(), correct_val*var.shape[0]*var.shape[1]*var.shape[2], tol)




class DepthAverageScriptTest_EqualLayers(DepthAverageScriptBaseTestCase):
    def setUp(self):
        self.ncin_fname = 'in_test_equal.nc'
        self.ncout_fname = 'out_test_equal.nc'
        create_test_file(fname=self.ncin_fname, layer_height='eq')
        create_depth_averaged_nc(self.ncin_fname, self.ncout_fname, log=False)
        self.nc = netCDF4.Dataset(self.ncout_fname, 'r')

    def tearDown(self):
        self.nc.close()
        os.remove(self.ncin_fname)
        os.remove(self.ncout_fname)


    def test_relative_layer_thickness(self):
        var = self.nc.variables['layer_relative_thickness']
        self.basetest_relative_layer_thickness(var, [0.2, 0.2, 0.2, 0.2, 0.2])

    def test_depth_averaged_spm001(self):
        var = self.nc.variables['concentration_of_SPM_in_water_001_averaged']
        self.basetest_depth_averaged_spm001(var, 8)
    
    def test_depth_averaged_spm002(self):
        var = self.nc.variables['concentration_of_SPM_in_water_002_averaged']
        self.basetest_depth_averaged_spm002(var, 16)

    def test_depth_averaged_summ_spm001_spm002(self):
        var = self.nc.variables['SUM_averaged']
        self.basetest_depth_averaged_summ_spm001_spm002(var, 24)



class DepthAverageScriptTest_NonEqualLayers(DepthAverageScriptBaseTestCase):
    def setUp(self):
        self.ncin_fname = 'in_test_nonequal.nc'
        self.ncout_fname = 'out_test_nonequal.nc'
        create_test_file(fname=self.ncin_fname, layer_height='noneq')
        create_depth_averaged_nc(self.ncin_fname, self.ncout_fname, log=False)
        self.nc = netCDF4.Dataset(self.ncout_fname, 'r')

    def tearDown(self):
        self.nc.close()
        os.remove(self.ncin_fname)
        os.remove(self.ncout_fname)

    def test_relative_layer_thickness(self):
        var = self.nc.variables['layer_relative_thickness']
        self.basetest_relative_layer_thickness(var, [0.02, 0.06, 0.12, 0.3, 0.5])

    def test_depth_averaged_spm001(self):
        var = self.nc.variables['concentration_of_SPM_in_water_001_averaged']
        self.basetest_depth_averaged_spm001(var, 6.8)
    
    def test_depth_averaged_spm002(self):
        var = self.nc.variables['concentration_of_SPM_in_water_002_averaged']
        self.basetest_depth_averaged_spm002(var, 13.6)

    def test_depth_averaged_summ_spm001_spm002(self):
        var = self.nc.variables['SUM_averaged']
        self.basetest_depth_averaged_summ_spm001_spm002(var, 20.4)


if __name__ == '__main__':


    unittest.main()