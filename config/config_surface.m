function Config = config_surface(ProjectDirectory)

    % Directory where the ADCP files are stored. This should be full path
    Config.DataDirectory = fullfile(ProjectDirectory,'data','adcp','LTA');

    % ADCP filename (*.LTA or *.STA format)
    Config.Filename = '*.LTA';

    % Time Adjustment shift
    Config.TimeShift = 7 + 53/60;  % In hours, This is added to Data.Time which will be in UTC time    % CHANGED: 53/60 to 60/60
                                   % The adcp.nav_stime seem to be in UTC time. So no shift will be added here.
                                   % The adcp.mtime seem to be 7:53' behind UTC time.

    % Cut time array to begining and end of numerical integration
    Config.ExtraTimeCutMarginIndex = 0; %200;         % Margin before and after the begining and end of computation time CHANGED: 0 to 100
                                                % The computation time is from start to end of numerical integration
                                                % The whole data is cut within a margin from start and end of computation time.

    % Plume Release Start Time and End Time, initial depth release
    % Time format [yyyy,mm,dd,hh,mm,ss];
    % Config.PlumeRelease.StartReleaseTime = [2018, 08, 16, 11, 27, 19];      % Surface release
    % Config.PlumeRelease.EndReleaseTime   = [2018, 08, 16, 11, 49, 49];

    Config.PlumeRelease.StartReleaseTime = [2018, 08, 16, 12, 00, 00];        % Surface release
    Config.PlumeRelease.EndReleaseTime   = [2018, 08, 16, 12, 00, 01];

    % Config.PlumeRelease.StartReleaseTime = [2018, 08, 16, 11, 40, 00];      % Surface release
    % Config.PlumeRelease.EndReleaseTime   = [2018, 08, 16, 11, 40, 01];

    % Config.PlumeRelease.StartReleaseTime = [2018, 08, 16, 08, 50, 00];      % Surface release
    % Config.PlumeRelease.EndReleaseTime   = [2018, 08, 16, 08, 50, 01];

    % Config.PlumeRelease.StartReleaseTime = [2018, 08, 14, 13, 36, 00];    % Subsurface release
    % Config.PlumeRelease.EndReleaseTime   = [2018, 08, 14, 14, 52, 00];
    Config.PlumeRelease.ReleaseDepth = 0;   % Meters. Positive downward. (3m surface, 15m subsurface)  % CHANGED: 3 to 2.88
    Config.PlumeRelease.ComputeStreakLine = true;

    % Initial time of plume track
    % Update position from last seen
    Config.PlumeLastSeen.UseUpdatedPosition = true;                                           % If true, set [B1] and either of [B2], [B4] or [B5-B6]. If false, it uses position at end of plume release
    Config.PlumeLastSeen.LastObservedTime = [2018, 08, 16, 12, 00, 01];                       % [B1]
    Config.PlumeLastSeen.UseDrifterPositionAtObservedTimeForPlume = false;                    % [B2] If true, it uses drifter poisition, set [B3]
    Config.PlumeLastSeen.DrifterIdsToUseForPositionAtObservedTimeForPlume = [25];             % [B3] Set if [B2] is true
    Config.PlumeLastSeen.UseShipPositionAtObservedTimeForPlume = false;                       % [B4] If true, it uses the ship positon at observed time
    % Config.PlumeLastSeen.LastObservedLongitude = [-70.555,-70.570];                            % [B5] Manually set position of plume in [B5] and [B6]
    % Config.PlumeLastSeen.LastObservedLatitude = [41.0945,41.100];                              % [B6] Manually set position of plume in [B5] and [B6]
    % Config.PlumeLastSeen.LastObservedLongitude = linspace(-70.555,-70.570,5);                            % [B5] Manually set position of plume in [B5] and [B6]    % used for subsurfaces
    % Config.PlumeLastSeen.LastObservedLatitude = linspace(41.0945,41.100,5);                              % [B6] Manually set position of plume in [B5] and [B6]
    % Config.PlumeLastSeen.LastObservedLongitude = linspace(-70.555,-70.585,7);                            % [B5] Manually set position of plume in [B5] and [B6]    % Used for subsurfaces
    % Config.PlumeLastSeen.LastObservedLatitude = linspace(41.0945,41.1055,7);                             % [B6] Manually set position of plume in [B5] and [B6]
    % Config.PlumeLastSeen.LastObservedLongitude = linspace(-70.555,-70.585,7);                            % [B5] Manually set position of plume in [B5] and [B6]    % Used for surfaces
    % Config.PlumeLastSeen.LastObservedLatitude = linspace(41.0945,41.088,7);                              % [B6] Manually set position of plume in [B5] and [B6]
    Config.PlumeLastSeen.LastObservedLongitude = linspace(-70.555,-70.555,7);                              % [B5] Manually set position of plume in [B5] and [B6]    % Used for surfaces
    Config.PlumeLastSeen.LastObservedLatitude = linspace(41.0915,41.11,7);                                 % [B6] Manually set position of plume in [B5] and [B6]

    % Final time of plume track
    % Inquiry Time (Current time)
    Config.PlumeInquiry.UseLastAvailableTime = false;                     % If false, set [A1]
    Config.PlumeInquiry.InquiryTime = [2018, 08, 17, 12, 00, 00];         % [A1] Surface release
    % Config.PlumeInquiry.InquiryTime = [2018, 08, 16, 18, 00, 00];         % [A1] Surface release  % Test with shorter integration time for model selection
    % Config.PlumeInquiry.InquiryTime = [2018, 08, 14, 18, 00, 00];       % [A1] Subsurface release

    % Update Depth from last seen
    Config.PlumeLastSeen.UseReleaseDepth = true;                          % If false, set [C1]. If true, it use the same depth at the plume release
    Config.PlumeLastSeen.CurrentDepth = 15;                               % [C1]

    % Velocities
    Config.ReadFromMatFiles = false;                                           % Uses Matlab mat files that exported from Teledyne's velocity software. Velocity is in earth coordinate
    Config.MatFileDirectory = fullfile(ProjectDirectory,'data','adcp','MAT');  % If Config.ReadFromMatFiles is true, set the directory of MAT files.
    Config.FillMissingAlongTimeMovingAverageWindow = 10;                       % If 0, it does not fill missing. If nonzero, if uses moving average to fill NaN
    Config.SmoothingWindowAlongVerticalRange = 8;                              % Effective if window is larger than 1. CHANGED: 3 to 1
    Config.SmoothingDurationAlongTime = [0 0 0 0 15 0];                        % Moving average smoothing in the format of: yyyy mm dd HH MM SS
    Config.OceanVelocityUnitToMeterPerSecond = 1;                              % ocean velocity in rdradcp.m is converted to m/s
    Config.ShipVelocityUnitToMeterPerSecond = 1;                               % ship velocity in rdradcp.m is converted to m/s
    Config.OceanVelocityIsRelativeToShip = true; 
    Config.ShipVelocityScale = 1;                                              % Scale ship velocity for velocity-dependent error due to ship movement
    Config.DiscardOceanVelocitiesBelowIntensitySpike = true;
    Config.NumberOfStreakLineParticles = 30;

    % Preprocess Velocity
    Config.PreprocessVelocity.FirstSmoothOverTimeThenInterpolateInDepth = true;          % if false, does the opposie, ie, interpolate in depth then smooth over time
    Config.PreprocessVelocity.PolynomialOrderForRegression = 2;                          % Used only in extrapolation when release depth is less than 2.88m, 
    Config.PreprocessVelocity.NumberOfPointsForRegression = 4;                           % Used only in extrapolation when release depth is less than 2.88m.

    % Windage
    Config.Wind.IncludeWindData = true;
    Config.Wind.WindFilename = fullfile(ProjectDirectory,'data','wind','nam.conusnest.hiresf.20180816.t06z.nc');
    Config.Wind.WindDataOriginTime = '1970-01-01 00:00:00';
    Config.Wind.WindageCoefficient = 0.0019;                      % CHANGED: 0.0019 to 0.005

    % Stochastic modeling
    Config.StochasticModel.UseStochasticModel = false;           % false uses deterministic model
    Config.StochasticModel.NumberOfSamples = 100;                % Num Monte-Carlo draws from Gaussian distributions for Brownian motion
    Config.StochasticModel.NumberOfRefinementSteps = 2;          % Number of evaulations between each dt, where dt is distance of two data times
    Config.StochasticModel.UseParallelSimulation = false;        % Note: do not use true, since it does not produce centered Gaussian sampling
    Config.StochasticModel.UseFirstOrderMarkovModel = false;     % If false, uses zero-th order markon model
    Config.StochasticModel.VelocitySTD = 14e-2;                  % In m/s. If negative is given, it uses analytic formula from beam covariances

    % Kriging
    Config.Kriging.UseKriging = true;                                         % If false, interpolates velocity with interp1 and assumes velocity v(x,t) = v(t) everywhere
    Config.Kriging.UseDualKriging = true;                                     % If false, uses primal krignign. Dual krigning is much faster than primal kriging
    Config.Kriging.SpatialCorrelationKernelType = 'matern-3/2';               % 'gaussian', 'matern-1/2', 'matern-3/2', and 'matern-5/3'
    Config.Kriging.SpatialBasisFunctionPolynomialDegree = 0;                  % Polynomial degree for only spatial variables x and y.
    Config.Kriging.ScaleLengthInSpatialBasisFunctions = 0.001;                % Distances in meters make large numbers. We scale them for better conditioned matrices. 0.001 converts m to km.
    Config.Kriging.TemporalBasisFunctionPeriodicity = 60*60*12 + 60*25;       % A tide cycle, 12 hours and 25 minutes in unit of seconds
    Config.Kriging.SplineStifnessInDualKriging = 10;                          % Spline stiffness weight
    Config.Kriging.UseFutureTimesInKriging = true;                            % If true, uses past and future ship data, if 'false', only uses past data till current time
    Config.Kriging.TemporalDecorrelationScale  = [0 0 0 00 40 00];            % In the unit of yyyy,mm,dd,HH,MM,SS
    Config.Kriging.TemporalDecorrelationCutoff = [0 0 0 03 00 00];            % In the unit of yyyy,mm,dd,HH,MM,SS
    Config.Kriging.SpatialDecorrelationScale = 2000;                          % meters
    Config.Kriging.SpatialDecorrelationCutoff = 10000;                         % meters
    Config.Kriging.ShipTimeIndexSkipping = 2; %1;                                 % Index. 1 reads all ship positions. 2 reads every 2 index, and so on.

    % Kalman Filter
    Config.KalmanFilter.UseKalmanFilter = false;                                                % If false, Kalman filter is not applied. If true, make sure [K1] below is also true.
    Config.KalmanFilter.DoTheUpdateStep = true;                                                 % [K1] If false, only the prediction step of KF is performed, thus Kalman filter is not really used.
    Config.KalmanFilter.ObservationStandardDeviation = 0;                                       % Error of position of drifters. Used in matrix R. in unit of meters
    Config.KalmanFilter.DriftersFileData = fullfile(ProjectDirectory,'data','Drifters.mat');    % Drifter data to be used for observation
    Config.KalmanFilter.ListOfControlDrifters = [24,25,28,31,32,33,34,35,36,37,38,39,40,41,42,43];
    % Config.KalmanFilter.ListOfControlDrifters = [24,28,32,34,36,38,40,42];
    % Config.KalmanFilter.ListOfControlDrifters = [17,18,19,20];
    Config.KalmanFilter.ObservationTimeIntervalsOfDrifters = [0 0 0 0 5 0];                    % every ten minutes, in the unit of yyyy,mm,dd,HH,MM,SS

    % Input Log file (matlab mat file that stores cancatenation of all input ADCP data)
    % If this file is deleted, all ADCP files are read again, otherwise this file is used.
    % The content of this file stores the 'Data' struct.
    Config.InputLogFilename = fullfile(ProjectDirectory,'log','InputLog.mat');

    % If UseInputLogFileFromPreviousRun is set to true, the Data struct will be read from the InputLog.mat file
    % providing if this file exits. This saves time on next runs, since reading all adcp raw
    % files might be time consuming. 
    Config.WriteInputDataToInputLogFile = true;
    Config.UseInputLogFileFromPreviousRun = true;

    % Output Log file (matlab mat file that stores previously ran computations)
    % If this file is deleted, all computation in the next run will start over.
    % The content of this file stores the 'Output' struct.
    Config.OutputLogFilename = fullfile(ProjectDirectory,'log','OutputLog.mat');

    % If UseOutputLogFileFromPreviousRun set to true, the 'StreakLine' trajectories are read from the OutputLog.mat file
    % providing this file exists. This saves time to not repeat the computation of strakline.
    % It is advised to set this option to false.
    Config.WriteOutputDataToOutputLogFile = true;
    Config.UseOutputLogFileFromPreviousRun = false;

    % Plots
    Config.Plots.ExtraPlots = false;                   % Extra plots of velocities
    Config.Plots.FontName = 'Latin Modern Math';       % Check available fonts with 'listfonts'. Example: 'Latin Modern Math', 'CMU Serif'
    Config.Plots.FiguresDirectory = fullfile(ProjectDirectory,'output');
    Config.Plots.FiguresFormat = '.pdf';
    % Config.Plots.DriftersToPlot = [24,25,28,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46];    % If empty [], it plots all drifters
    Config.Plots.DriftersToPlot = [24,25,28,31,32,33,34,35,36,37,38,39,40,41,42,43];               % If empty [], it plots all drifters
    % Config.Plots.DriftersToPlot = [];                                                            % Plots all drifters
    % Config.Plots.DriftersToPlot = [17,18,19,20,21];
    % Config.Plots.DriftersToPlot = [17,18,19,20];

    % Plots for python (external code)
    Config.SaveForPythonPlots = true;
    Config.SaveForPythonDirectory = fullfile(ProjectDirectory,'python','data');

end
