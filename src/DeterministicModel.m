classdef DeterministicModel

    methods(Static)

    % =============================
    % Trajectory Initial Conditions
    % =============================

    function [TrajectoryStartTime,TrajectoryStartLongitudes,TrajectoryStartLatitudes] = TrajectoryInitialConditions( ...
        Config,Data,Output,InquiryTime)

        % Convert datetimes to times (second from 0000 Jan 00)
        DayToSecond = 24.0 * 60.0 * 60.0;

        % Initial Position
        if Config.PlumeLastSeen.UseUpdatedPosition == true

            % Last Observed time
            TrajectoryStartTime = datenum(Config.PlumeLastSeen.LastObservedTime) * DayToSecond;

            % Use a drifter position for initial conditon
            if Config.PlumeLastSeen.UseDrifterPositionAtObservedTimeForPlume == true

                % Load drifter data
                DriftersFileData = load(Config.KalmanFilter.DriftersFileData);
                Drifters = DriftersFileData.Drifters;

                % Initilaize an array of trajecotries
                NumberOfInitialTrajectories = length(Config.PlumeLastSeen.DrifterIdsToUseForPositionAtObservedTimeForPlume);
                TrajectoryStartLongitudes = zeros(1,NumberOfInitialTrajectories);
                TrajectoryStartLatitudes = zeros(1,NumberOfInitialTrajectories);
                for TracerIterator = 1:NumberOfInitialTrajectories

                    % Find drifter Id
                    DrifterId = -1;
                    for i = 1:size(Drifters,2)
                        if Drifters(i).Id == Config.PlumeLastSeen.DrifterIdsToUseForPositionAtObservedTimeForPlume(TracerIterator)
                            DrifterId = i;
                            break
                        end
                    end
                    if DrifterId == -1
                        error('Cannot find drifter Id: %d',Config.PlumeLastSeen.DrifterIdsToUseForPositionAtObservedTimeForPlume(TracerIterator))
                    end

                    % Drifter times
                    DrifterTimes = Drifters(DrifterId).Time * DayToSecond;
                    if (DrifterTimes(1) > TrajectoryStartTime) || (DrifterTimes(end) < TrajectoryStartTime)
                        error('Drifter %d time is from %s to %s but trajectory start time is %s.', ...
                        DrifterId,datestr(DrifterTimes(1)/DayToSecond),datestr(DrifterTimes(end)/DayToSecond),datestr(TrajectoryStartTime/DayToSecond));
                    end
                    if (DrifterTimes(1) > InquiryTime) || (DrifterTimes(end) < InquiryTime)
                        error('Drifter %d time is from %s to %s but inquiry start time is %s.', ...
                        DrifterId,datestr(DrifterTimes(1)/DayToSecond),datestr(DrifterTimes(end)/DayToSecond),datestr(InquiryTime/DayToSecond));
                    end

                    % Interrpolate drifter positon
                    TrajectoryStartLongitudes(TracerIterator) = interp1(DrifterTimes,Drifters(DrifterId).Longitude,TrajectoryStartTime);
                    TrajectoryStartLatitudes(TracerIterator) = interp1(DrifterTimes,Drifters(DrifterId).Latitude,TrajectoryStartTime);
                end

            % Use the ship position for initial conditon
            elseif Config.PlumeLastSeen.UseShipPositionAtObservedTimeForPlume == true

                % Interpolate ship position
                TrajectoryStartLongitudes = interp1(Data.Time,Data.Ship.Longitude,TrajectoryStartTime,'pchip');
                TrajectoryStartLatitudes = interp1(Data.Time,Data.Ship.Latitude,TrajectoryStartTime,'pchip');

            % Use a custom position for initial conditon
            else

                % Number of trajectories
                NumberOfInitialTrajectories = length(Config.PlumeLastSeen.LastObservedLongitude);
                TrajectoryStartLongitudes = zeros(1,NumberOfInitialTrajectories);
                TrajectoryStartLatitudes = zeros(1,NumberOfInitialTrajectories);

                for TracerIterator = 1:NumberOfInitialTrajectories

                    % Use user custom input
                    TrajectoryStartLongitudes(TracerIterator) = Config.PlumeLastSeen.LastObservedLongitude(TracerIterator);
                    TrajectoryStartLatitudes(TracerIterator) = Config.PlumeLastSeen.LastObservedLatitude(TracerIterator);

                end
            end

        else

            % Use center of initial streakline
            TrajectoryStartTime = Output.PlumeInitialStreakLine.EndReleaseTime;
            TrajectoryStartLongitudes = Output.PlumeInitialStreakLine.CenterLongitudeAtEndReleaseTime;
            TrajectoryStartLatitudes = Output.PlumeInitialStreakLine.CenterLatitudeAtEndReleaseTime;

        end

    end

    % ============================
    % Control Drifter Trajectories
    % ============================

    function [DrifterControlTimes,ControlTrajectoryLongitudes,ControlTrajectoryLatitudes] = ControlDrifterTrajectories( ...
        Config, ...
        TrajectoryStartTime,InquiryTime,ListOfControlDrifterIds)

        % Number of control drifters: D
        % Number of discrete times from start to end: T
        %
        % Output:
        % ControlTrajectoryLongitudes: T*D matrix
        % ControlTrajectoryLatitudes:  T*D matrix

        % Time step to discretize drifters
        DayToSecond = 60 * 60 * 24;
        TimeStepAsDate = Config.KalmanFilter.ObservationTimeIntervalsOfDrifters;    % Duration of time intervals to observe drifter positions. Here every ten minutes
        TimeStep = DayToSecond * datenum(TimeStepAsDate);   % convert to the unit of seconds

        % Array of discrete times from start to end
        NumTimes = floor(0.5 + (InquiryTime - TrajectoryStartTime) / TimeStep) + 1;
        DrifterControlTimes = linspace(TrajectoryStartTime,InquiryTime,NumTimes);

        % Initialize arrays
        NumDrifters = length(ListOfControlDrifterIds);
        ControlTrajectoryLongitudes = zeros(NumTimes,NumDrifters);
        ControlTrajectoryLatitudes = zeros(NumTimes,NumDrifters);

        % Load drifter data
        DriftersFileData = load(Config.KalmanFilter.DriftersFileData);
        Drifters = DriftersFileData.Drifters;

        % Iterate over each drifter
        DrifterCounter = 1;
        for DrifterIdToSearch = ListOfControlDrifterIds

            % Find drifter Id
            DrifterIdFound = -1;
            for i = 1:size(Drifters,2)
                if Drifters(i).Id == DrifterIdToSearch
                    DrifterIdFound = i;
                    break
                end
            end
            if DrifterIdFound == -1
                error('Cannot find drifter Id: %d',Config.PlumeLastSeen.DrifterIdToUseForPositionAtObservedTimeForPlume)
            end

            % Drifter times
            DrifterTimes = Drifters(DrifterIdFound).Time * DayToSecond;
            if (DrifterTimes(1) > TrajectoryStartTime) || (DrifterTimes(end) < TrajectoryStartTime)
                error('Drifter %d time is from %s to %s but trajectory start time is %s.', ...
                DrifterIdFound,datestr(DrifterTimes(1)/DayToSecond),datestr(DrifterTimes(end)/DayToSecond),datestr(TrajectoryStartTime/DayToSecond));
            end
            if (DrifterTimes(1) > InquiryTime) || (DrifterTimes(end) < InquiryTime)
                error('Drifter %d time is from %s to %s but inquiry start time is %s.', ...
                DrifterIdFound,datestr(DrifterTimes(1)/DayToSecond),datestr(DrifterTimes(end)/DayToSecond),datestr(InquiryTime/DayToSecond));
            end

            % Interrpolate drifter positon
            ControlTrajectoryLongitudes(:,DrifterCounter) = interp1(DrifterTimes,Drifters(DrifterIdFound).Longitude,DrifterControlTimes)';
            ControlTrajectoryLatitudes(:,DrifterCounter) = interp1(DrifterTimes,Drifters(DrifterIdFound).Latitude,DrifterControlTimes)';

            DrifterCounter = DrifterCounter + 1;
        end

    end

    % ================
    % Get Inquiry Time
    % ================

    function InquiryTime = GetInquiryTime(Config,Data)

        % Convert datetimes to times (second from 0000 Jan 00)
        DayToSecond = 24.0 * 60.0 * 60.0;
    
        % Inquiry time
        if Config.PlumeInquiry.UseLastAvailableTime == true
            InquiryTime = Data.Time(end);
        else
            InquiryTime = datenum(Config.PlumeInquiry.InquiryTime) * DayToSecond;
        end

        % Check inquity time is larger than data time
        if InquiryTime > Data.Time(end) - 10*eps
            error('Inquiry time: %s is larger then the last data time: %s.', ...
            datestr(InquiryTime/DayToSecond),datestr(Data.Time(end)))
        end

    end

    % ==========================
    % Find Plume At Inquiry Time
    % ==========================

    function Output = FindPlumeAtInquiryTime(Data,Config,Output)

        disp('Computing deterministic trajectories ...')

        % Convert datetimes to times (second from 0000 Jan 00)
        DayToSecond = 24.0 * 60.0 * 60.0;

        % Plume Depth
        if Config.PlumeLastSeen.UseReleaseDepth == true
            PlumeDepth = Config.PlumeRelease.ReleaseDepth;
        else
            PlumeDepth = Config.PlumeLastSeen.CurrentDepth;
        end

        % Inquiry time
        InquiryTime = DeterministicModel.GetInquiryTime(Config,Data);

        %  Trajectory Initial conditons
        [TrajectoryStartTime,TrajectoriesStartLongitudes,TrajectoriesStartLatitudes] = DeterministicModel.TrajectoryInitialConditions( ...
            Config,Data,Output,InquiryTime);

        % Control Drifter Trajectories
        ListOfControlDrifterIds = Config.KalmanFilter.ListOfControlDrifters;
        [ControlDriftersTimes,ControlDriftersLongitude,ControlDriftersLatitude] = DeterministicModel.ControlDrifterTrajectories( ...
            Config, ...
            TrajectoryStartTime,InquiryTime,ListOfControlDrifterIds);

        % Specify an origin for coordinate system
        OriginLongitude = Data.Ship.Longitude(1);
        OriginLatitude = Data.Ship.Latitude(1);

        % Convert lon/lat to XY for trajectory, ship and control drifters
        [TrajectoriesStartX,TrajectoriesStartY] = Utilities.ConvertLonLatToXY(OriginLongitude,OriginLatitude,TrajectoriesStartLongitudes,TrajectoriesStartLatitudes);
        [ControlDriftersX,ControlDriftersY] = Utilities.ConvertLonLatToXY(OriginLongitude,OriginLatitude,ControlDriftersLongitude,ControlDriftersLatitude);
        [ShipX,ShipY] = Utilities.ConvertLonLatToXY(OriginLongitude,OriginLatitude,Data.Ship.Longitude',Data.Ship.Latitude');
        ShipCoordinates = [ShipX,ShipY]; % in meters

        % Preprocess ocean velocity data (interpolate at plume depth, add wind, smooth over time and extract mean and noise)
        if Config.PreprocessVelocity.FirstSmoothOverTimeThenInterpolateInDepth == true
            % First smoooth over time, then interpoate over depth
            OceanVelocityData = PreprocessVelocityData.Process_FirstSmoothOverTimeThenInterpolateInDepth(Config,Data,PlumeDepth);  % Better to use this one
        else
            % First interpoate over depth, then smooth over time
            OceanVelocityData = PreprocessVelocityData.Process_FirstInterpolateInDepthThenSmoothOverTime(Config,Data,PlumeDepth);
        end

        Ocean_VelAtDepth_u       = OceanVelocityData.Ocean_VelAtDepth_u;
        Ocean_VelAtDepth_v       = OceanVelocityData.Ocean_VelAtDepth_v;
        Ocean_VelAtDepth_u_Mean  = OceanVelocityData.Ocean_VelAtDepth_u_Mean;
        Ocean_VelAtDepth_v_Mean  = OceanVelocityData.Ocean_VelAtDepth_v_Mean;
        Ocean_VelAtDepth_u_Noise = OceanVelocityData.Ocean_VelAtDepth_u_Noise;
        Ocean_VelAtDepth_v_Noise = OceanVelocityData.Ocean_VelAtDepth_v_Noise;

        % if Config.KalmanFilter.UseKalmanFilter == true  % TODO
            % Standard deviation of velocity data
            [StandardDeviation_U,StandardDeviation_V,VelocityStandardDeviationAvergae] = PreprocessVelocityData.EstimateStandardDeviationOfVelocity(Ocean_VelAtDepth_u_Noise,Ocean_VelAtDepth_v_Noise);
        % end

        % Choose how to obtain ocean velocity field (interpolation, kriging, dual kriging)
        if Config.Kriging.UseKriging == true

            % Decorrelation settings
            TemporalDecorrelationScale  = DayToSecond * datenum(Config.Kriging.TemporalDecorrelationScale);    % in the unit of seconds
            TemporalDecorrelationCutoff = DayToSecond * datenum(Config.Kriging.TemporalDecorrelationCutoff);   % in the unit of seconds
            SpatialDecorrelationScale = Config.Kriging.SpatialDecorrelationScale;                              % meters
            SpatialDecorrelationCutoff = Config.Kriging.SpatialDecorrelationCutoff;                            % Index (integer)
            Skip = Config.Kriging.ShipTimeIndexSkipping;

            % Basis functions at cluster of points
            [F,BasisPolynomialOrder,ScaleLength] = BasisFunctions.GenerateBasisFunctions(Config,ShipCoordinates(1:Skip:end,:),Data.Time(1:Skip:end));

            % Covariance of cluster of ship points with itself. Method 1: using user-defined covariance parameters
            DataCovariance = CovarianceFunction.CovarianceOfDataPointsWithItself( ...
                Config, ...
                ShipCoordinates(1:Skip:end,:),Data.Time(1:Skip:end), ...
                TemporalDecorrelationScale,TemporalDecorrelationCutoff, ...
                SpatialDecorrelationScale,SpatialDecorrelationCutoff);

            % Covariance of cluster of data points with itself. Method 2:  using Bayesian Inference
            % [DataCovariance2,CovarianceParameters] = ModelSelection.CovarianceOfDataPointsWithItself( ...
            %     Config, ...
            %     F, ...
            %     ShipCoordinates(1:Skip:end,:),Data.Time(1:Skip:end), ...
            %     Ocean_VelAtDepth_u(1:Skip:end),Ocean_VelAtDepth_v(1:Skip:end));

            % Coefficients of basis functions for a given velocity data
            beta_U = BasisFunctions.SolveCoefficientsOfBasisFunctions(DataCovariance,F,Ocean_VelAtDepth_u_Mean(1:Skip:end)');
            beta_V = BasisFunctions.SolveCoefficientsOfBasisFunctions(DataCovariance,F,Ocean_VelAtDepth_v_Mean(1:Skip:end)'); 

            if Config.Kriging.UseDualKriging == false

                % Use Primal Kriging
                OceanVelocity = @(t,X) GaussianProcessRegression.KrigVelocityAtInquiryPoint( ...
                    Config, ...
                    X,t, ...
                    ShipCoordinates(1:Skip:end,:),Data.Time(1:Skip:end),...
                    DataCovariance,F, ...
                    TemporalDecorrelationScale,TemporalDecorrelationCutoff, ...
                    SpatialDecorrelationScale,SpatialDecorrelationCutoff, ...
                    Ocean_VelAtDepth_u_Mean(1:Skip:end),Ocean_VelAtDepth_v_Mean(1:Skip:end));

            else

                % Column vectors b and d to be used for dual kriging
                [b_u,b_v,d_u,d_v] = GaussianProcessRegression.DualKrigCoefficients( ...
                    Config, ...
                    ShipCoordinates(1:Skip:end,:),Data.Time(1:Skip:end),...
                    DataCovariance,F, ...
                    Ocean_VelAtDepth_u_Mean(1:Skip:end),Ocean_VelAtDepth_v_Mean(1:Skip:end));

                % Use Dual kriging
                OceanVelocity = @(t,X) GaussianProcessRegression.DualKrigVelocityAtInquiryPoint( ...
                    Config, ...
                    X,t, ...
                    ShipCoordinates(1:Skip:end,:),Data.Time(1:Skip:end), ...
                    b_u,b_v,d_u,d_v, ...
                    TemporalDecorrelationScale,TemporalDecorrelationCutoff, ...
                    SpatialDecorrelationScale,SpatialDecorrelationCutoff);

            end
        else
            % Do not use Krigning. Function handle to interpolate velocity in time
            OceanVelocity =@(t,x) [interp1(Data.Time,Ocean_VelAtDepth_u_Mean,t);interp1(Data.Time,Ocean_VelAtDepth_v_Mean,t)];
        end

        % Trace trajectory (choose either kalman filter or straightforward integration)
        if Config.KalmanFilter.UseKalmanFilter == false

            % Do not use Kalmand filter. So, the drifter data will not correct the solution
            opts = odeset('RelTol',1e-5);
            [TrajectoriesT,TrajectoryCoordinates] = ...
                ode45(OceanVelocity,[TrajectoryStartTime,InquiryTime],[TrajectoriesStartX,TrajectoriesStartY],opts);

            % Separate X and Y in array
            NumberOfTrajectories = length(TrajectoriesStartX);
            TrajectoriesX = TrajectoryCoordinates(:,1:NumberOfTrajectories);
            TrajectoriesY = TrajectoryCoordinates(:,NumberOfTrajectories+1:end);

        else

            % Kalman filter should be used when Kriging is enabled
            if Config.Kriging.UseKriging == false
                error('When Kalman filter is used, the Krigning (either primal or dual) should be enabled.')
            end
        
            % Hybrid Kalman filter
            [TrajectoriesT, ...
                TrajectoriesXStandardDeviation,TrajectoriesYStandardDeviation, ...
                TrajectoriesX,TrajectoriesY, ...
                ControlTrajectoriesX,ControlTrajectoriesY] = ...
            KalmanFilter.HybridKalmanFilter( ...
                Config, ...
                TrajectoriesStartX,TrajectoriesStartY, ...
                ControlDriftersTimes,ControlDriftersX,ControlDriftersY, ...
                OceanVelocity, ...
                BasisPolynomialOrder,ScaleLength,beta_U,beta_V, ...
                VelocityStandardDeviationAverage,SpatialDecorrelationScale,SpatialDecorrelationCutoff);

        end

        % Convert back XY to LonLat
        [TrajectoriesLongitude,TrajectoriesLatitude] = Utilities.ConvertXYToLonLat( ...
            OriginLongitude,OriginLatitude, ...
            TrajectoriesX,TrajectoriesY);

        % If Kalmand filter was used, convert control trajectories too
        if Config.KalmanFilter.UseKalmanFilter == true
            [ControlTrajectoriesLongitude,ControlTrajectoriesLatitude] = Utilities.ConvertXYToLonLat( ...
                OriginLongitude,OriginLatitude, ...
                ControlTrajectoriesX,ControlTrajectoriesY);
        end

        % Output: trajectory
        Output.PlumeCenterTrajectory.Longitudes = TrajectoriesLongitude;
        Output.PlumeCenterTrajectory.Latitudes = TrajectoriesLatitude;
        Output.PlumeCenterTrajectory.Times = TrajectoriesT;

        % Ship trajectory during plume track
        Output.ShipTrajectory.Longitude = interp1(Data.Time,Data.Ship.Longitude,TrajectoriesT,'pchip');
        Output.ShipTrajectory.Latitude = interp1(Data.Time,Data.Ship.Latitude,TrajectoriesT,'pchip');

        % Output: first and last point on a trajectory
        TrajectoryIndex = 1;
        Output.PlumeCenterTrajectory.FirstLongitude = TrajectoriesLongitude(1,TrajectoryIndex);
        Output.PlumeCenterTrajectory.FirstLatitude  = TrajectoriesLatitude(1,TrajectoryIndex);
        Output.PlumeCenterTrajectory.LastLongitude  = TrajectoriesLongitude(end,TrajectoryIndex);
        Output.PlumeCenterTrajectory.LastLatitude   = TrajectoriesLatitude(end,TrajectoryIndex);
        Output.PlumeCenterTrajectory.LastTime       = TrajectoriesT(end);

        if Config.KalmanFilter.UseKalmanFilter == true
            Output.ControlDrifters.Time       = ControlDriftersTimes;
            Output.ControlDrifters.Longitudes = ControlTrajectoriesLongitude;
            Output.ControlDrifters.Latitudes  = ControlTrajectoriesLatitude;
        end

    end

    % =====================
    % Plot Plume Trajectory
    % =====================

    function PlotPlumeTrajectory(Data,Config,Output)

        fig = figure();
        oceanColor = [.6 .8 .9];
        ax = axesm('MapProjection','mercator');
        setm(ax,'FFaceColor',oceanColor);

        % Ship whole trajectory
        % h1 = plot(Data.Ship.Longitude,Data.Ship.Latitude,'-o','color','black','DisplayName','Ship trajectory');
        h1 = plotm(Data.Ship.Latitude,Data.Ship.Longitude,'-o','color','black','DisplayName','Ship trajectory');
        hold on
         
        % Plume center trajectory
        NumberOfTrajectories = size(Output.PlumeCenterTrajectory.Longitudes,2);
        for i = 1:NumberOfTrajectories
            h2 = plotm(Output.PlumeCenterTrajectory.Latitudes(:,i),Output.PlumeCenterTrajectory.Longitudes(:,i),'color','blue','linewidth',2,'DisplayName','Center plume trajectory');
            plotm(Output.PlumeCenterTrajectory.Latitudes(1,i),Output.PlumeCenterTrajectory.Longitudes(1,i),'o','MarkerEdgeColor','blue','MarkerFaceColor','blue','MarkerSize',8);
            plotm(Output.PlumeCenterTrajectory.Latitudes(end,i),Output.PlumeCenterTrajectory.Longitudes(end,i),'v','MarkerEdgeColor','blue','MarkerFaceColor','blue','MarkerSize',8);
            textm(Output.PlumeCenterTrajectory.Latitudes(1,i),Output.PlumeCenterTrajectory.Longitudes(1,i),datestr(Output.PlumeCenterTrajectory.Times(1)/(24*3600),' mm/dd HH:MM'));
            textm(Output.PlumeCenterTrajectory.Latitudes(end,i),Output.PlumeCenterTrajectory.Longitudes(end,i),datestr(Output.PlumeCenterTrajectory.Times(end)/(24*3600),' mm/dd HH:MM'));
        end

        % Control drifter trajectories
        if isfield(Output,'ControlDrifters')
            NumberOfTrajectories = size(Output.ControlDrifters.Longitudes,2);
            for i = 1:NumberOfTrajectories
                h22 = plotm(Output.ControlDrifters.Latitudes(:,i),Output.ControlDrifters.Longitudes(:,i),'color','magenta','linewidth',2,'DisplayName','Virtual drifter trajectory');
                plotm(Output.ControlDrifters.Latitudes(1,i),Output.ControlDrifters.Longitudes(1,i),'o','MarkerEdgeColor','magenta','MarkerFaceColor','magenta','MarkerSize',8);
                plotm(Output.ControlDrifters.Latitudes(end,i),Output.ControlDrifters.Longitudes(end,i),'v','MarkerEdgeColor','magenta','MarkerFaceColor','magenta','MarkerSize',8);
                textm(Output.ControlDrifters.Latitudes(1,i),Output.ControlDrifters.Longitudes(1,i),datestr(Output.PlumeCenterTrajectory.Times(1)/(24*3600),' mm/dd HH:MM'));
                textm(Output.ControlDrifters.Latitudes(end,i),Output.ControlDrifters.Longitudes(end,i),datestr(Output.PlumeCenterTrajectory.Times(end)/(24*3600),' mm/dd HH:MM'));
            end
        end

        % Ship trajectory during current plume track
        h3 = plotm(Output.ShipTrajectory.Latitude,Output.ShipTrajectory.Longitude,'color','green','linewidth',2,'DisplayName','Ship trajectory during plume track');
        plotm(Output.ShipTrajectory.Latitude(1),Output.ShipTrajectory.Longitude(1),'o','MarkerEdgeColor','green','MarkerFaceColor','green','MarkerSize',8);
        plotm(Output.ShipTrajectory.Latitude(end),Output.ShipTrajectory.Longitude(end),'v','MarkerEdgeColor','green','MarkerFaceColor','green','MarkerSize',8);
        textm(Output.ShipTrajectory.Latitude(1),Output.ShipTrajectory.Longitude(1),datestr(Output.PlumeCenterTrajectory.Times(1)/(24*3600),' mm/dd HH:MM'));
        textm(Output.ShipTrajectory.Latitude(end),Output.ShipTrajectory.Longitude(end),datestr(Output.PlumeCenterTrajectory.Times(end)/(24*3600),' mm/dd HH:MM'));

        % Translated StreakLine Position to the observed time
        PlumeObservedStreakLineLongitude = Output.PlumeCenterTrajectory.FirstLongitude - Output.PlumeInitialStreakLine.CenterLongitudeAtEndReleaseTime + Output.PlumeInitialStreakLine.LongitudeAtEndReleaseTime;
        PlumeObservedStreakLineLatitude = Output.PlumeCenterTrajectory.FirstLatitude - Output.PlumeInitialStreakLine.CenterLatitudeAtEndReleaseTime + Output.PlumeInitialStreakLine.LatitudeAtEndReleaseTime;

        % plot translated plume streakline
        h4 = plotm(PlumeObservedStreakLineLatitude,PlumeObservedStreakLineLongitude,'color',[1,0.7,0.7],'linewidth',2,'DisplayName','Plume streakline at observed time');

        % Translated StreakLine Position to the inquiry time
        PlumeCurrentStreakLineLongitude = Output.PlumeCenterTrajectory.LastLongitude - Output.PlumeInitialStreakLine.CenterLongitudeAtEndReleaseTime + Output.PlumeInitialStreakLine.LongitudeAtEndReleaseTime;
        PlumeCurrentStreakLineLatitude = Output.PlumeCenterTrajectory.LastLatitude - Output.PlumeInitialStreakLine.CenterLatitudeAtEndReleaseTime + Output.PlumeInitialStreakLine.LatitudeAtEndReleaseTime;

        % plot translated plume streakline
        h5 = plotm(PlumeCurrentStreakLineLatitude,PlumeCurrentStreakLineLongitude,'color','red','linewidth',2,'DisplayName','Plume streakline at inquiry time');

        % Plot Drifters
        DriftersFileData = load(Config.KalmanFilter.DriftersFileData);
        Drifters = DriftersFileData.Drifters;
        DayToSecond = 24.0 * 3600.0;
        h_drift = [];
        Counter = 1;
        for i = 1:size(Drifters,2)

            % Filter Drifters (empty DriftersToPlot will plot all)
            DrifterId = Drifters(i).Id;
            if isnan(DrifterId)
                continue
            end
            if ~isempty(Config.Plots.DriftersToPlot) && all(DrifterId ~= Config.Plots.DriftersToPlot)
                continue
            end

            DrifterTime = Drifters(i).Time * DayToSecond;
            [UniqueDrifterTime,Index] = unique(DrifterTime);   % Remove duplicates
            DrifterColor = [0.45,0.45,0.45];

            if DrifterTime(1) <= Data.Time(1) && DrifterTime(end) >= Data.Time(end)

                fprintf('Plotting drifter: %d\n',DrifterId)

                % Get Drifter info
                DrifterLongitude = Drifters(i).Longitude;
                DrifterLatitude = Drifters(i).Latitude;

                % Interpolate lat and lon to Data.Time
                TrackingTimes = Output.PlumeCenterTrajectory.Times;
                % Longitude = interp1(UniqueDrifterTime,DrifterLongitude(Index),Data.Time,'pchip');
                % Latitude = interp1(UniqueDrifterTime,DrifterLatitude(Index),Data.Time,'pchip');
                Longitude = interp1(UniqueDrifterTime,DrifterLongitude(Index),TrackingTimes,'pchip');
                Latitude = interp1(UniqueDrifterTime,DrifterLatitude(Index),TrackingTimes,'pchip');

                h_drift(Counter) = plotm(Latitude,Longitude,'color',DrifterColor,'DisplayName',strcat('Drifter ',num2str(DrifterId)));
                hold on
                plotm(Latitude(1),Longitude(1),'o','MarkerEdgeColor','black','MarkerFaceColor',DrifterColor,'MarkerSize',4,'DisplayName','')
                plotm(Latitude(end),Longitude(end),'v','MarkerEdgeColor','black','MarkerFaceColor',DrifterColor,'MarkerSize',4,'DisplayName','')
                textm(Latitude(floor(end/2)),Longitude(floor(end/2)),num2str(DrifterId),'HorizontalAlignment','center','FontSize',8)
                Counter = Counter + 1;

            end

        end

        % Determine Border of map
        MinLon = min([PlumeObservedStreakLineLongitude(:);PlumeCurrentStreakLineLongitude(:);Output.ShipTrajectory.Longitude(:);Data.Ship.Longitude(:)]);
        MaxLon = max([PlumeObservedStreakLineLongitude(:);PlumeCurrentStreakLineLongitude(:);Output.ShipTrajectory.Longitude(:);Data.Ship.Longitude(:)]);
        MinLat = min([PlumeObservedStreakLineLatitude(:);PlumeCurrentStreakLineLatitude(:);Output.ShipTrajectory.Latitude(:);Data.Ship.Latitude(:)]);
        MaxLat = max([PlumeObservedStreakLineLatitude(:);PlumeCurrentStreakLineLatitude(:);Output.ShipTrajectory.Latitude(:);Data.Ship.Latitude(:)]);

        LatExt = MaxLat - MinLat;
        LonExt = MaxLon - MinLon;
        Ext = max(LatExt,LonExt);
        Ratio = 0.1;

        LatLim = [MinLat-Ratio*Ext,MaxLat+Ratio*Ext];
        LonLim = [MinLon-Ratio*Ext,MaxLon+Ratio*Ext];

        % xlim([MinLon,MaxLon])
        % ylim([MinLat,MaxLat])

        setm(gca,'MapLatLimit',LatLim,'MapLonLimit',LonLim)
        % setm(gca,'FLatLimit',LatLim,'FLonLimit',LonLim)
        % axis('equal')

        xlabel('Longitude(deg)')
        ylabel('Latitude (deg)')
        title('Plume trajectory')

        setm(gca,'MLineLocation',0.005,'MLabelLocation',0.005,'MLabelRound',-3,'MeridianLabel','on','MLabelParallel','south');
        setm(gca,'PLineLocation',0.005,'PLabelLocation',0.005,'PLabelRound',-3,'ParallelLabel','on','PLabelMeridian','west');

        tightmap;    % Tighten up the map area shown
        framem on;   % Turn on the black frame
        gridm on;    % Turn on grid lines
        mlabel on;
        plabel on;
        showaxes;
        grid off;
        axis off;

        scaleruler('unit','km')
        setm(handlem('scaleruler1'),'MajorTick',0:0.5:1,'MinorTick',0:0.1:0.5,'TickDir','down')

        if isfield(Output,'ControlDrifters')
            legend([h1,h2,h22,h3,h4,h5])
        else
            legend([h1,h2,h3,h4,h5])
        end
        
        % Font
        set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

        % Save
        saveas(gcf,fullfile(Config.Plots.FiguresDirectory,strcat('DeterministicModelTrajectory',Config.Plots.FiguresFormat)))

    end

    % ---------------------
    % End of static methods

    end
end
