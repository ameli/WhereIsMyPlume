This directory is initially empty when the repository is cloned. Download the [actual data files](http://transport.me.berkeley.edu/thredds/fileServer/trajectories/FieldExperiment-2018/data.tar.gz) (650 MB) to replace with the empty `data` in this directory. You may do so either manually or by the following commands:

    cd WhereIsMyPlume
    wget http://transport.me.berkeley.edu/thredds/fileServer/trajectories/FieldExperiment-2018/data.tar.gz
    gunzip data.tar.gz
    tar -xvf data.tar

After downloading the data, the contents of this directory should include:

1. `/adcp` folder containing ADCP data with subfolders `LTA`, `STA`, and `MAT`.
2. `/drifters` folder containing `DrifterIds.txt`, and subfolders `csv` and `xlsx`.
3. `/wind` folder containing `nc` files for NAM wind forecasts.
4. `Drifters.mat` file.
5. `ShipTrack.mat` file.
