[![licence](https://img.shields.io/badge/licence-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# WhereIsMyPlume

### Description

_Objective:_ This code uses measured ocean currents data to predict the trajectory of virtual drifters or a plume of released dye. We assume the ocean current measurements are obtained by shipboard [acoustic Doppler current profiler](https://en.wikipedia.org/wiki/Acoustic_Doppler_current_profiler). Predictng drifter trajectories only based on the ADCP data is a computationally challenging task. This is becase the ADCP data are sparse both in space and time, and they are time series of one-dimensional profile of the ocean current that probed only a column underneath the moving ship. In addition, some other factors can drastically affect the trajectory prediction, such as the wind (for surface drifters), the uncertainty of the measurement devices, and the rapid change of the ocean dynamics, such as the tidal effects.

_Method:_ This code reconstructs a full spatio-temporal velocity field of the ocean currents using the Gaussian process prediction. The predicted velocity field is then used to trace virtual drifter trajectories, either by deterministic or stochastic models. Using Kalman filtering method, the trajectory predictions of virtual drifters are further enhanced by the assimilation of actual drifters that are floating around near the field of study.

### Requirements

* _Matlab:_ version 2017a or later with the following Matlab toolboxes are required:

  * Statistics toolbox
  * Econometrics toolbox
  * Financial toolbox
  * Mapping toolbox

* _Python_: (optional, you may not need this step) To reproduce the output plots in python, the Jupyter noteook and the following python packages are required:

  * numpy
  * scipy
  * matplotlib
  * basemap
  * matlab engine
  
  You can install these packages, for instance, with anaconda by
  
      sudo conda install -c conda-forge jupyter numpy scipy matplotlib basemap
  
  Also, install Matlab engine for Python by
  
      cd matlabroot/extern/engines/python
      sudo python setup.py install
  
  In the above, `matlabroot` is the directory where Matlab is installed.

### Download the code

Clone the repository:

    git clone https://github.com/ameli/WhereIsMyPlume.git

### Download the data

The cloned repository contains the directory [`data`](https://github.com/ameli/WhereIsMyPlume/tree/master/data) which is an empty folder. Download the [actual data files](http://transport.me.berkeley.edu/thredds/fileServer/trajectories/FieldExperiment-2018/data.tar.gz) (650 MB) to replace with the empty `data` directory in the repository. You may do so either manually or by the following commands:

    cd WhereIsMyPlume
    wget http://transport.me.berkeley.edu/thredds/fileServer/trajectories/FieldExperiment-2018/data.tar.gz
    gunzip data.tar.gz
    tar -xvf data.tar

These data were assimilated during two-weeks field experiment on the offshore of the south coast of Martha's Vineyard, MA, in August 2018. The data directory contain:

* `/data/adcp`: The raw ADCP files from Teledyne RDI workhorse instrument. These include both `LTA` and `STA` file types.
* `/data/wind`: NAM wind forecast data at 10m above the ground as netCDF4 file types.
* `/data/drifters`: Raw data of the GPS coordinates of 45 surface and subsurface code drifters deployed during the field experiments. These include both `xlsx` and `csv` formats.
* `/data/Drifters.mat`: Processed data of the drifters and ready to be used in Matlab.
* `/data/ShipTrack.mat`: Processed GPS coordinates of the ship trajectory during the field experiment.

### Run the code

1. If you are running Matlab on a server, start Matlab with:

        matlab -nodekstop -nosplash

2. The main script of the code is [`WhereIsMyPlume.m`](https://github.com/ameli/WhereIsMyPlume/blob/master/WhereIsMyPlume.m) and accepts an input configuration file to configure how the code should run. Some sample configuration files can be found in [`config`](https://github.com/ameli/WhereIsMyPlume/tree/master/config) directory. As an example, run the code with the configuration given in [`/config/config_subsurface.m`](https://github.com/ameli/WhereIsMyPlume/blob/master/config/config_subsurface.m) by

        cd WhereIsMyPlume
        WhereIsMyPlume 'config_subsurface'

3. (Optional) For publication quality plots, run the jupyter notebook script in [`python`](https://github.com/ameli/WhereIsMyPlume/tree/master/python) directory. This should be done after the Matlab code is executed which generates the output data in [`/python/data`](https://github.com/ameli/WhereIsMyPlume/tree/master/python/data). These data will be used by the python script to reproduce the plots.

### Output of the code

The results of the code are written into two directories:

* ['/log'](https://github.com/ameli/WhereIsMyPlume/tree/master/log): Two `mat` files will be written to the directory this directory, namely:
  * `/log/InputLog.mat`: The code reads the raw data (ADCP, wind, drifters) and stores them as Matlab's struct arrays in this file. On the second run of the code, if this file still exists, the code will no longer reads the raw data again, rather, directly uses this file directory as an input. If this file is deleted, the code reads the raw data again, which could be time consuming.
  * `/log/OutputLog.mat`: All arrays generated by the code will be written as a Matlab struct array in this file.
* [`/output`](https://github.com/ameli/WhereIsMyPlume/tree/master/output): The plots will be written to this directory. The type and number of output plots depend on the settings in the configuration file. For instance, they might be the trajectory predictions of drifters, probability density function of the stochastic model of trajectories, velocity column profile of the ADCP data, etc. 

### Configuration file

The user can interact with the code's seetings by editing the configuration files as desired. Samples of the configuration files in [`/config`](https://github.com/ameli/WhereIsMyPlume/tree/master/config) are as follows:

* [`/config/surface`](https://github.com/ameli/WhereIsMyPlume/blob/master/config/config_surface.m): example of settings for surface drifters without assimilating the drifter data.
* [`/config/subsurface`](https://github.com/ameli/WhereIsMyPlume/blob/master/config/config_subsurface.m): example of settings for subsurface drifters without assimilating the drifter data.
* [`/config/kalman_surface`](https://github.com/ameli/WhereIsMyPlume/blob/master/config/config_kalman_surface.m): example of settings for surface drifters with assimilating the drifter data using Kalman filter.
* [`/config/kalman_subsubsurface`](https://github.com/ameli/WhereIsMyPlume/blob/master/config/config_kalman_subsubsurface.m): example of settings for subsurface drifters with assimilating the drifter data using Kalman filter.

### Directory Structure

| Directory/File | Purpose |
| -------------- | ------- |
| [`/config`](https://github.com/ameli/WhereIsMyPlume/tree/master/config) | Examples of configuration files are stored here. A configuration file is used as an argument to the main script. |
| [`/data`](https://github.com/ameli/WhereIsMyPlume/tree/master/data) | Includes adcp data, ship tracks, drifter data, and wind data. The user data can be any other directory, provided that the data directory is set in config file. |
| [`/doc`](https://github.com/ameli/WhereIsMyPlume/tree/master/doc) | Documentation files, such as output images. |
| [`/log`](https://github.com/ameli/WhereIsMyPlume/tree/master/log) | Both input and output `mat` files are stored here. If an input log file exists, the code will use it in the next, otherwise, the code generates a new log file. |
| [`/output`](https://github.com/ameli/WhereIsMyPlume/tree/master/output) | Stores all plots as an output of the code. |
| [`/python`](https://github.com/ameli/WhereIsMyPlume/tree/master/python) | Python scripts to reproduce the plots with a better quality (than Matlab) for publication purposes only. |
| [`/src`](https://github.com/ameli/WhereIsMyPlume/tree/master/src) | Source code `m` files. |
| [`/test`](https://github.com/ameli/WhereIsMyPlume/tree/master/test) | Scripts to be run by continuous integration, Travis. |
| [`/utilities`](https://github.com/ameli/WhereIsMyPlume/tree/master/utilities) | Directory containing scripts that are not a part of the package, but useful, in particular, for data files conversion and plots. |
| [`WhereIsMyPlume.m`](https://github.com/ameli/WhereIsMyPlume/blob/master/WhereIsMyPlume.m) | Main script of the package. This script calls the codes in `/src` directory. |

### Acknowledgements

* This work was supported by National Science Foundation award number [1520825](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1520825).
* The script `/src/rdradcp.m`, which reads the raw ADCP data files, is written by [Rich Pawlowicz](http://www.eoas.ubc.ca/~rich/).

### How to cite

The publications will be updated.
