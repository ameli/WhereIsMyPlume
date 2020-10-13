% -------------------------------------------------------------------------
%         filename:  WhereIsMyPlume.m 
%           author:  Siavash Ameli
%             date:  2018/08/11
%          project:  NSF ALPHA, full study field experiment, August 2018
% -------------------------------------------------------------------------

% =================
% Where Is My Plume
% =================

function WhereIsMyPlume(ConfigFilename)

    % This function sets the user's hardcoded configuration.
    % and calls the Manager (main function).

    clear

    % Project directory
    [ProjectDirectory,~,~] = fileparts(mfilename('fullpath'));

    % Add the script path
    addpath('config');
    addpath('src');

    % Config files
    if ~exist('ConfigFilename','var')

        % Set a default config filename
        % ConfigFilename = 'config_kalman_surface';
        % ConfigFilename = 'config_kalman_subsurface';
        % ConfigFilename = 'config_surface';
        ConfigFilename = 'config_subsurface';
    end

    % Run the main code
    Config = feval(ConfigFilename,ProjectDirectory);
    Manager(Config);

end
