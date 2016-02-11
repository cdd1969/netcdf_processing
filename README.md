##About
This folder contains two scripts for processing netcdf files:

1. dim_mean_nc.py
2. depth_average_nc.py


**dim_mean_nc.py** - calculate of the mean value along selected dimension (with respect only to the number of layers) in netcdf-variable.

**depth_average_nc.py** - calculate depth-averaged values along selected depth-dimension (real depth-values are used for integration)

## Installation
Download the source code, no specific installation is needed:

```sh
$ git clone https://github.com/cdd1969/netcdf_processing.git scripts
```
## Testing
**Note:** *testing requires the write-right in the `scripts/test` directory, because some netcdf files are generated during testing (they are deleted after tests are done)*

It is recommended to run test cases before using the script on a real data. To run tests type:

```sh
$ python -m scripts.test.run -v
```
or alternatively run tests with displaying the coverage rate:

```sh
$ coverage run --source scripts -m scripts.test.run
$ coverage report
```

## Running
Scripts have command line user interface.

```sh
$ cd scripts
$ python depth_average_nc.py --help
```

##Documentation
1. run the script with `--help` option
2. read docstrings in source