This directory is initially empty after cloning the repository. After running the Matlab code, this directory stores user input and output log files, which are:

1. `InputLog.mat`: This file stores input data as Matlab's `mat` file. The raw data files of ADCP, drifters, ship tracks, wind data (which are in the directory) come in various file formats. When the code is run for the first time, these raw data are read and written into `InputLog.mat`. By this, the next time you run the code, the raw data will not be needed to be read again. However, to re-read the raw data files again (such as using new configuration settings), simply delete the `InputLog.mat` file in this directory in the future runs.
2. `OutputLog.mat`: All results of the computations will be stored in this file.
