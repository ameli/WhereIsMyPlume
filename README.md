# WhereIsMyPlume

* [Homepage](homepage)
* [![licence](https://img.shields.io/badge/licence-MIT-blue.svg)](https://opensource.org/licenses/MIT)

### Description

You are on a boat in the middle of nowhere ocean. You lost the treasure map which was in a bottle, and now it is floating somewhere. Fortunately, your boat comes with with an [acoustic Doppler current profiler](https://en.wikipedia.org/wiki/Acoustic_Doppler_current_profiler). Using the ocean current measurements, this code guides you to the location of the drifting bottle.

### Requirements

Matlab 2017a and later with the following toolbox are required:

* Statistics toolbox
* Econometrics toolbox
* Financial toolbox
* Map toolbox

### Download the code

1. Clone the repository:

        git clone https://github.com/ameli/WhereIsMyPlume.git

2. The cloned repository contains the directory `data` which is an empty folder. Download the [actual data files](http://transport.me.berkeley.edu/thredds/fileServer/trajectories/FieldExperiment-2018/data.tar.gz) (650 MB) to replace with empty `data` directory in the repository. You may do so manually or by the following commands:

        cd WhereIsMyPlume
        wget http://transport.me.berkeley.edu/thredds/fileServer/trajectories/FieldExperiment-2018/data.tar.gz
        gunzip data.tar.gz
        tar -xvf data.tar

### Run the code

1. If running Matlab on a server, start Matlab with:

        matlab -nodekstop -nosplash

2. The main script of the code is `WhereIsMyPlume.m` and accepts a configuration file to configure how the code should run. The configuration. Some sample configuration files can be found in `config` directory. As an example, run the code with the configuration given in `/config/config_subsurface.m` file by

        cd WhereIsMyPlume
        WhereIsMyPlume 'config_subsurface'

   After the code is finished, the resulted `mat` files will be written to the directory `log`. The plots will be written to the directory `output`.

4. For publication quality plots, run the jupyter notebook script in `python` directory. This should be done after the matlab code is executed which generates the output data in `/python/data`. These data will be used by the python script to reproduce the plots.

### Directory Structure

| Directory/File | Purpose |
| -------------- | ------- |
| `/config` | Examples of configuration files are stored here. A configuration file is used as an argument to the main script. |
| `/data` | Includes adcp data, ship tracks, drifter data, and wind data. The user data can be any other directory, provided that the data directory is set in config file. |
| `/doc` | Documentation files, such as output images. |
| `/log` | Both input and output `mat` files are stored here. If an input log file exists, the code will use it in the next, otherwise, the code generates a new log file. |
| `/output` | Stores all plots as an output of the code. |
| `/python` | Python scripts to reproduce the plots with a better quality (than Matlab) for publication purposes only. |
| `/src` | Source code `m` files. |
| `/test` | Scripts to be run by continuous integration, Travis. |
| `/utilities` | Directory containing scripts that are not a part of the package, but useful, in particular, for data files conversion and plots. |
| `WhereIsMyPlume.m` | Main script (runner) of the package. This script runs the codes in `/src` directory. |


### Acknowledgements

* This work was supported by National Science Foundation award number [1520825](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1520825).
* The script `/src/rdradcp.m`, which reads the raw ADCP data files, is written by [Rich Pawlowicz](http://www.eoas.ubc.ca/~rich/).

### How to cite

Coming.
