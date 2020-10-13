% This script is used for continuous integration with Travis.

% Add the root directory to the path
[CurrentDirectory,~,~] = fileparts(mfilename('fullpath'));
[RootDirectory,~,~] = fileparts(CurrentDirectory);
addpath(RootDirectory);

% Main script
MainScript = 'WhereIsMyPlume';

% A config filename for test
ConfigFilename = 'config_subsurface';

% Run a test
feval(MainScript,ConfigFilename);
