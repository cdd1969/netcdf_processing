##About
This folder contains two scripts for processing netcdf files:

1. dim_mean_nc.py
2. depth_average_nc.py


**dim_mean_nc.py** - use for calculating mean value along selected dimension (with respect only to the number of layers) in netcdf-variable.

**depth_average_nc.py** - use for calculate depth-averaged values along selected depth-dimension (real depth-values are used for integration)

## Installation
Download the source code, no specific installation is needed:

```sh
$ git clone https://github.com/cdd1969/netcdf_processing.git scripts
```
## Testing
It is recommended afterwards to run test cases:

```sh
$ python -m scripts.test.run -v
```
or alternatively to test with displaying the coverage rate:

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