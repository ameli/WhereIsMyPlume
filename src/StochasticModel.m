classdef StochasticModel

    methods(Static)

    % ===========================
    % Compute Velocity Covariance
    % ===========================

    function VelocityCovariance = ComputeVelocityCovariance(Data,Config,VelocityStandardDeviation)

    % References:
    % 
    % 1. Hummon and Firing, 2002
    %    A Direct Comparison of Two RDI Shipboard ADCPs: A 75-kHz Ocean Surveyor and a 150-kHz Narrow Band
    % 
    % 2. Teledyn RD Instrument (RDI)
    %    WorkHouse manual, commands and output data format
    % 
    % 3. Gilcoto, et al. 2009
    %    Robust Estimations of Current Velocities with Four-Beam Broadband ADCPs
    %
    % Notes:
    %
    % 1. Mode:
    %    RDI has several profiling modes (see RDI manual, section WM)
    %    The adcp.config.prof_mode = 1, which means the profile mode 1 is used.
    %    Mode 1 is used for dynamic ocean velocities where an abrupt change of velocity might be 
    %    or high ship speed might be expected. It has a higher range of velocity meaurement but
    %    error is also higher. Mode 1 is characterized by code WB, and it can be WB0 or WB1.
    %    According to RDI manual (section WM), for 300 hertz frezuency, the mode WB0 is used.
    %
    % 2. Bandwidth:
    %    It can be 100%, 25%, 6.25% etc.
    %    For WB0, the bandwidth is wide (= 25%).
    %
    % 3. Cycles:
    %    This is 100 / Bandwidth and is equal to tau * f. Here f is frequency and tau is the time
    %    of each bit (see Hummon and Firing, appendix).
    %
    % 4. Ambiguity velocity:
    %    For Mode 1, the V_am is high, hence cause more error. According to RDI manual (section WV)
    %    the max V_am is 700 cm/s.
    %
    % Setting Velocity Standard Deviation manually (without above calculations)
    % See: http://www.teledynemarine.com/Lists/Downloads/mariner_datasheet_hr.pdf
    %      Also, ADCP Coordinate Transformation Formulas and Calculations


    SpeedOfSound = 1498.0;   % In water, m/s

    % For profile mode 1, with 300 hertz, the configuration is WB0, with bandwidth 25%
    Bandwidth = 12.0;                    % 25% for WB0, 12% for WB1 (see workhouse manual, section WB)
    Cycles = 100.0 / Bandwidth;          % 4 cycles (see Hummon and Firing 2002, appendix). Cycle = Frequency*tau (see Hummon, paragraph between equation A3 to A4)
    Frequency = Data.Beam.Frequency;     % 300 hertz
    LagLength = Data.Beam.LagLength;     % 7
    CellSize = Data.Beam.CellSize;       % 1 meter (see STA log files for WS100, meaning cell size 100cm)
    BeamAngle = Data.Beam.Angle;         % 20 degrees (in radian)

    % Ambiguity velocity
    % AmbiguityVelocity = SpeedOfSound / (4.0 * LagLength * Cycles)   % 13.375 m/s, See Hummon et al, Equation A3. 

    % % RDI's max ambiguity velocity for WB0 (see warehouse manual, section WV)
    % if AmbiguityVelocity > 7.0
    %     AmbiguityVelocity = 7.0;
    % end
    AmbiguityVelocity = 3.3;   % According to setting in STA log files, this is set to WV330 (330 cm/s). Note, the actual V_a might be different.

    % Coefficient
    SaftyFactor = 1.5;
    PI = 3.14159265;
    Coefficient = (((SaftyFactor * AmbiguityVelocity) / (2.0 * PI))^2) * ((Cycles * SpeedOfSound * cos(BeamAngle)) / (Frequency * CellSize));

    % Beam covariance (sigma^2) for each beam, at each depth and each time. Dimension is [D,4,T]. See Gilcoto, et al. 2009.
    BeamCovariance = Coefficient * (1.0./(Data.Beam.Correlation.^2) - 1.0);

    % Beam covariances for each beam (= sigma_i^2), with dimension of [D,T], D is depth and T is time.
    S1 = squeeze(BeamCovariance(:,1,:));    % Beam 1 and 3 are in opposite directions
    S2 = squeeze(BeamCovariance(:,2,:));    % Beam 2 and 4 are in opposite directions
    S3 = squeeze(BeamCovariance(:,3,:));
    S4 = squeeze(BeamCovariance(:,4,:));

    % Velocity covariance, dimension is [D,T,2,2]
    VelocityCovariance = ones(size(BeamCovariance,1),size(BeamCovariance,3),2,2);

    % Iterate over time
    for t = 1:size(BeamCovariance,3)

        % Matrix to convert ship coordinate to earth coordinate
        theta = Data.Ship.Heading(t);   % In radian, clockwise from noth

        % Coordinate x-y is as follow. y is in the direction of heading of ship, and x is lateral to the right side of ship
        % We want to counter-clockwise rotate x-y to be on the east-north. The angle between these to are heading (angle between north and y, or east and x).
        Q = [[cos(theta),sin(theta)];[-sin(theta),cos(theta)]];     % Counter-clockwise rotation to take heading back to north, to counteract clockwise heading.

        % Iterate over depth
        for d = 1:size(BeamCovariance,1)

            % Velocity covariances (on ship coordinate, not earth coordinate), Dimension is [2,2]
            if Config.StochasticModel.VelocitySTD > 0
                VelocityCovarianceInShipCoordinate = eye(2) * (Config.StochasticModel.VelocitySTD)^2;   % Standard deviation in m/s
                % VelocityCovarianceInShipCoordinate = eye(2) * (VelocityStandardDeviation)^2;   % Standard deviation in m/s
            else
                VelocityCovarianceInShipCoordinate = [[S1(d,t)+S2(d,t),0];[0,S3(d,t)+S4(d,t)]] ./ (4.0 * (sin(BeamAngle))^2);
            end

            % Convert velocity covariance from ship coordinate to earth coordinate. If X_earth = Q*X_ship, then cov(X_earth) = Q*cov(X_ship)*Q'.
            VelocityCovariance(d,t,:,:) = Q * VelocityCovarianceInShipCoordinate * Q';
        end
    end

    % Modify velocity covariance based on number of pings per ensemble
    VelocityCovariance = VelocityCovariance ./ Data.Beam.PingsPerEnsemble;

    % Plot results
    if Config.Plots.ExtraPlots == true
        StochasticModel.PlotVelocityCovariance(Data,Config,VelocityCovariance);
    end

    end

    % ========================
    % Plot Velocity Covariance
    % ========================

    function PlotVelocityCovariance(Data,Config,VelocityCovariance)

    figure()

    % Find max bathymetry
    NonZero = find(all(isnan(Data.Ocean.AbsoluteVel_u),2));
    if isempty(NonZero)
        MaxBathymetry = size(Data.Ocean.AbsoluteVel_u,1);
    else
        MaxBathymetry = NonZero(1);
    end

    SecondToDay = 1.0 / (24.0 * 3600.0);
    Datetime = Data.Time * SecondToDay;

    subplot(2,2,1)
    pcolor(Datetime,Data.Ocean.Range(1:MaxBathymetry,1),VelocityCovariance(1:MaxBathymetry,:,1,1))
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    shading interp;
    h1 = colorbar;
    set(gca,'Ydir','reverse');
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Ranges (m)');
    ylabel(h1,'(m^2/s^2)');
    title('East-East Variance');
    set(gca,'XGrid','on')
    set(gca,'color',[0.8,0.52,0.24])

    subplot(2,2,2)
    pcolor(Datetime,Data.Ocean.Range(1:MaxBathymetry,1),VelocityCovariance(1:MaxBathymetry,:,1,2))
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    shading interp;
    h2 = colorbar;
    set(gca,'Ydir','reverse');
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Ranges (m)');
    ylabel(h2,'(m^2/s^2)');
    title('East-North Covariance');
    set(gca,'XGrid','on')
    set(gca,'color',[0.8,0.52,0.24])

    subplot(2,2,3)
    pcolor(Datetime,Data.Ocean.Range(1:MaxBathymetry,1),VelocityCovariance(1:MaxBathymetry,:,2,1))
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    shading interp;
    h3 = colorbar;
    set(gca,'Ydir','reverse');
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Ranges (m)');
    ylabel(h3,'(m^2/s^2)');
    title('North-East Covariance');
    set(gca,'XGrid','on')
    set(gca,'color',[0.8,0.52,0.24])

    subplot(2,2,4)
    pcolor(Datetime,Data.Ocean.Range(1:MaxBathymetry,1),VelocityCovariance(1:MaxBathymetry,:,2,2))
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    shading interp;
    h4 = colorbar;
    set(gca,'Ydir','reverse');
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Ranges (m)');
    ylabel(h4,'(m^2/s^2)');
    title('North-North Variance');
    set(gca,'XGrid','on')
    set(gca,'color',[0.8,0.52,0.24])

    % Font
    set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

    % Save
    saveas(gcf,fullfile(Config.Plots.FiguresDirectory,strcat('VelocityCovariance',Config.Plots.FiguresFormat)))

    end

    % ================================
    % Compute AutoCorrelation Function
    % ================================

    function [AutoCovariance,AutoCorrelation,Lags] = ComputeAutoCorrelationFunction(Vel_u,Vel_v,Time,NumTimeLags,Config)

    % Ouput: AutoCorrelation of size (2,2,NumTimeLags)
    % AutoCorrelation(i,j,:) is the autocorrelation between i-th and j-th components of velocity.
    % AutoCorrelation(i,j,1) at initial time is the identity tensor.
    %
    % Negatine NumTimeLag will use all length of time series

    % Default time lag
    if NumTimeLags < 0
        NumTimeLags = size(Vel_u,2) - 1;
    end

    % Cut length
    if NumTimeLags >= size(Vel_u,2)
        NumTimeLags = size(Vel_u,2) - 1;
    end

    % Velocity
    Vel = [Vel_u(:),Vel_v(:)];

    % Time difference
    dt = Time(2) - Time(1);

    Dimension = 2;
    AutoCovariance = zeros(NumTimeLags,Dimension,Dimension);

    for i = 1:Dimension
        for j = 1:Dimension
            for LagIndex = 0:NumTimeLags-1

                % Assuming non-periodic signals
                Signal1 = squeeze(Vel(1:end-LagIndex,i));
                Signal2 = squeeze(Vel(1+LagIndex:end,j));
                TimeWindow = Time(1:end-LagIndex);
                Normalization = TimeWindow(end) - TimeWindow(1);     % Unbiased
                % Normalization = Time(end) - Time(1);                   % Biased

                % Assuming Periodic signal (results will be the same as using FFT)
                % Signal1 = squeeze(Vel(:,i))';
                % Signal2 = [squeeze(Vel(1+LagIndex:end,j))',squeeze(Vel(1:LagIndex,j))'];
                % TimeWindow = Time(:);
                % Normalization = TimeWindow(end) - TimeWindow(1);

                % Use whole signal length
                Mean1 = mean(Vel(:,i));
                Mean2 = mean(Vel(:,j));
                Std1 = std(Vel(:,i));
                Std2 = std(Vel(:,j));

                % Use only the portion of signal in convolution
                % Mean1 = mean(Signal1);
                % Mean2 = mean(Signal2);
                % Std1 = std(Signal1);
                % Std2 = std(Signal2);

                % AutoCovariance(LagIndex+1,i,j) = sum((Signal1 - Mean1).*(Signal2 - Mean2)) * dt / (Normalization);
                AutoCovariance(LagIndex+1,i,j) = trapz(TimeWindow,((Signal1 - Mean1).*(Signal2 - Mean2))) / (Normalization);

            end
        end
    end

    % Symmetrize AutoCovariance, so that alpha(tau) = alpha(-tau)
    % for LagIndex = 1:NumTimeLags-1
    %     for i = 1:Dimension
    %         for j = i:Dimension
    %             if i ~= j
    %                 AutoCovariance(LagIndex+1,i,j) = 0.5 * (AutoCovariance(LagIndex+1,i,j) + AutoCovariance(LagIndex+1,j,i));
    %                 AutoCovariance(LagIndex+1,j,i) = AutoCovariance(LagIndex+1,i,j);
    %             end
    %         end
    %     end
    % end

    % Normalize autocovariance to obtain autocorrelation. Method I (divide by diagonals)
    % Covariance = squeeze(AutoCovariance(1,:,:));
    % Cinv = diag(1./sqrt(diag(Covariance)));
    % AutoCorrelation = zeros(NumTimeLags,Dimension,Dimension);
    % for LagIndex = 1:NumTimeLags
    %     AutoCorrelation(LagIndex,:,:) = Cinv * squeeze(AutoCovariance(LagIndex,:,:)) * Cinv;
    % end

    % Normalize autocovariance to obtain autocorrelation. Method II
    Covariance = squeeze(AutoCovariance(1,:,:));
    C_inv = inv(Covariance);
    C_inv_sqrt = sqrtm(C_inv);
    AutoCorrelation = zeros(NumTimeLags,Dimension,Dimension);
    for LagIndex = 1:NumTimeLags
        AutoCorrelation(LagIndex,:,:) = C_inv * squeeze(AutoCovariance(LagIndex,:,:));                          % ACF will not be symmetric
        % AutoCorrelation(LagIndex,:,:) = C_inv_sqrt * squeeze(AutoCovariance(LagIndex,:,:)) * C_inv_sqrt;      % ACF will be symmetric
    end

    % Lags in seconds
    Lags = Time(1:NumTimeLags) - Time(1);  % in seconds

    if Config.Plots.ExtraPlots == true

        % Plot velocity
        figure();
        hold on
        TimeInDays = Time / (24*3600);
        plot(TimeInDays,Vel_u,'DisplayName','u')
        plot(TimeInDays,Vel_v,'DisplayName','v')
        xticks(linspace(TimeInDays(1),TimeInDays(end),4))
        datetick('x','mm/dd HH:MM','keepticks','keeplimits')
        xtickangle(30)
        xlabel('Time')
        ylabel('Velocity (m s^{-1})')
        grid on
        legend()
        xlim([TimeInDays(1),TimeInDays(end)])
        title('Velocity')
        set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

        % Save
        saveas(gcf,fullfile(Config.Plots.FiguresDirectory,strcat('VelocityTimeSeries',Config.Plots.FiguresFormat)))

        % Plot ACF
        figure();
        hold on
        LagsScaled = Lags / (60*60); % In hours
        plot(LagsScaled,squeeze(AutoCorrelation(:,1,1)),'DisplayName','\alpha_{11}')
        plot(LagsScaled,squeeze(AutoCorrelation(:,1,2)),'DisplayName','\alpha_{12}')
        plot(LagsScaled,squeeze(AutoCorrelation(:,2,1)),'DisplayName','\alpha_{21}')
        plot(LagsScaled,squeeze(AutoCorrelation(:,2,2)),'DisplayName','\alpha_{22}')
        legend()
        xlim([0,LagsScaled(end)])
        xlabel('Lags (hours)')
        ylabel('\alpha')
        title('Autocorrelation of components of velocity')
        ax = gca;
        ax.YGrid = 'on';
        set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

        % Save
        saveas(gcf,fullfile(Config.Plots.FiguresDirectory,strcat('AutoCorrelation',Config.Plots.FiguresFormat)))
    end

    end

    % =============================
    % Compute Lagrangian Time Scale
    % =============================

    function T = ComputeLagrangianTimeScale(AutoCorrelation,Lags)

    % T is (2,2) tensor, where each component (i,j) is computed from the autocorrelation AutoCorrelation(:,i,j).
    % Here we use normalized time lag integral.

    Dimension = size(AutoCorrelation,2);
    T = zeros(Dimension,Dimension);

    for i = 1:Dimension
        for j = 1:Dimension
            T(i,j) = trapz(Lags,squeeze(AutoCorrelation(:,i,j)));
            % T(i,j) = sum(squeeze(AutoCorrelation(:,i,j)));
        end
    end

    % Simplified version by making it just diagonal
    % T = eye(Dimension) * mean(T(:));

    end

    % ===========================================
    % Stochastic Plume Trajectory At Inquiry Time
    % ===========================================

    function [Output,ErrorCode] = StochasticPlumeTrajectoryAtInquiryTime(Data,Config,Output)

    disp('Computing stochastic trajectories ...')
    ErrorCode = 0;

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

    % Plume Depth
    if Config.PlumeLastSeen.UseReleaseDepth == true
        PlumeDepth = Config.PlumeRelease.ReleaseDepth;
    else
        PlumeDepth = Config.PlumeLastSeen.CurrentDepth;
    end

    % Initial Position
    if Config.PlumeLastSeen.UseUpdatedPosition == true

        % Last Observed time
        TrajectoryStartTime = datenum(Config.PlumeLastSeen.LastObservedTime) * DayToSecond;

        % Use a drifter position for initial conditon
        if Config.PlumeLastSeen.UseDrifterPositionAtObservedTimeForPlume == true

            % Load drifter data
            DriftersFileData = load(Config.KalmanFilter.DriftersFileData);
            Drifters = DriftersFileData.Drifters;

            % Find drifter Id
            DrifterId = -1;
            for i = 1:size(Drifters,2)
                if Drifters(i).Id == Config.PlumeLastSeen.DrifterIdsToUseForPositionAtObservedTimeForPlume
                    DrifterId = i;
                    break
                end
            end
            if DrifterId == -1
                error('Cannot find drifter Id: %d',Config.PlumeLastSeen.DrifterIdToUseForPositionAtObservedTimeForPlume)
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
            TrajectoryStartLongitude = interp1(DrifterTimes,Drifters(DrifterId).Longitude,TrajectoryStartTime);
            TrajectoryStartLatitude = interp1(DrifterTimes,Drifters(DrifterId).Latitude,TrajectoryStartTime);

        % Use the ship position for initial conditon
        elseif Config.PlumeLastSeen.UseShipPositionAtObservedTimeForPlume == true

            % Interpolate ship position
            TrajectoryStartLongitude = interp1(Data.Time,Data.Ship.Longitude,TrajectoryStartTime,'pchip');
            TrajectoryStartLatitude = interp1(Data.Time,Data.Ship.Latitude,TrajectoryStartTime,'pchip');

        % Use a custom position for initial conditon
        else

            % Use user custom input
            TrajectoryStartLongitude = Config.PlumeLastSeen.LastObservedLongitude;
            TrajectoryStartLatitude = Config.PlumeLastSeen.LastObservedLatitude;

        end

    else

        % Use center of initial streakline
        TrajectoryStartTime = Output.PlumeInitialStreakLine.EndReleaseTime;
        TrajectoryStartLongitude = Output.PlumeInitialStreakLine.CenterLongitudeAtEndReleaseTime;
        TrajectoryStartLatitude = Output.PlumeInitialStreakLine.CenterLatitudeAtEndReleaseTime;

    end

    % Specify an origin for coordinate system
    OriginLongitude = Data.Ship.Longitude(1);
    OriginLatitude = Data.Ship.Latitude(1);

    % Convert lon/lat to XY
    [TrajectoryStartX,TrajectoryStartY] = Utilities.ConvertLonLatToXY(OriginLongitude,OriginLatitude,TrajectoryStartLongitude,TrajectoryStartLatitude);

    % Preprocess ocean velocity data (interpolate at plume depth, add wind, smooth over time and extract mean and noise)
    % OceanVelocityData = PreprocessVelocityData.Process_FirstInterpolateInDepthThenSmoothOverTime(Config,Data,PlumeDepth);
    OceanVelocityData = PreprocessVelocityData.Process_FirstSmoothOverTimeThenInterpolateInDepth(Config,Data,PlumeDepth);

    Ocean_VelAtDepth_u       = OceanVelocityData.Ocean_VelAtDepth_u;
    Ocean_VelAtDepth_v       = OceanVelocityData.Ocean_VelAtDepth_v;
    Ocean_VelAtDepth_u_Mean  = OceanVelocityData.Ocean_VelAtDepth_u_Mean;
    Ocean_VelAtDepth_v_Mean  = OceanVelocityData.Ocean_VelAtDepth_v_Mean;
    Ocean_VelAtDepth_u_Noise = OceanVelocityData.Ocean_VelAtDepth_u_Noise;
    Ocean_VelAtDepth_v_Noise = OceanVelocityData.Ocean_VelAtDepth_v_Noise;

    % Standard deviation of velocity data
    VelocityStandardDeviation = PreprocessVelocityData.EstimateStandardDeviationOfVelocity(Ocean_VelAtDepth_u_Noise,Ocean_VelAtDepth_v_Noise);

    % Autocorrelation Tensor
    NumTimeLags = floor(length(Data.Time)/3);   % Use a quarter of total dataset length so that diagonals of ACF be positive
    [AutoCovariance,AutoCorrelation,Lags] = StochasticModel.ComputeAutoCorrelationFunction(Ocean_VelAtDepth_u_Noise,Ocean_VelAtDepth_v_Noise,Data.Time,NumTimeLags,Config);
    % Covariance = squeeze(AutoCovariance(1,:,:));
    Covariance = (Config.StochasticModel.VelocitySTD)^2 * eye(2);
    fprintf('Covariance:\n')
    disp(Covariance)

    % Integral Time Scale (choose either of methods below)

    % Method 1: Use Auto-Correlation
    % T = StochasticModel.ComputeLagrangianTimeScale(AutoCorrelation,Lags);

    % Method 2: Use Auto-Covariance
    % P = StochasticModel.ComputeLagrangianTimeScale(AutoCovariance,Lags);
    % T = inv(Covariance)*P;

    % Method 3: Guess a constant decorrelation time
    % T = eye(2) * 3600*1;  % One hour in seconds
    T = eye(2) * 60*60*0.05;  % three minues duration (in the unit of seconds)

    % Check if P is positive definite
    % Eigenvalues_P = eig(P+P');
    % if any(Eigenvalues_P < 0)
    %     fprintf('P is not positive definite.\n')
    %     disp(P)
    %     disp(Eigenvalues_P)
    %     ErrorCode = 1;
    %     return
    % end

    Dimension = 2;
    A = inv(T)';

    % Check if T is positive definite, otherwise B^2 = C*A + A'*C will not be positive definite, and B can not be real matrix.
    % Eigenvalues_T = eig(T+T');
    % if any(Eigenvalues_T < 0)
    %     fprintf('Lagrangian decorrelation length tensor T is not positive definite.\n')
    %     fprintf('To solve this issue, descrease time lag length used in computation of autocorrelation function.\n')
    %     fprintf('Eigenvalues of T: %f, %f\n',Eigenvalues_T(1),Eigenvalues_T(2));
    %     disp('Lagrangian decorrelation length T: ')
    %     disp(T)
    %     ErrorCode = 1;
    %     return
    % end

    % Check if real part of eigenvalues of A are positive
    Eigenvalues_A = eig(A);
    if any(real(Eigenvalues_A) < 0)
        fprintf('A has eigenvalues with negative real part.\n')
        disp(Eigenvalues_A)
        ErrorCode = 1;
        return
    end

    % Interpolate velocity covariance at depth, Array size: (Time,Direction,Direction)
    VelocityCovariance = StochasticModel.ComputeVelocityCovariance(Data,Config,VelocityStandardDeviation);
    VelocityCovarianceAtDepth = zeros(size(VelocityCovariance,2),2,2);
    for i = 1:Dimension
        for j = 1:Dimension
            % Interpolates along rows for each column of 2D array (depths). Definitely needs extrapolation when plume depth is near ocean surface, which is out of Range
            VelocityCovarianceAtDepth(:,i,j) = interp1(Data.Ocean.Range,VelocityCovariance(:,:,i,j),PlumeDepth,'pchip','extrap');
        end
    end

    % Function handle to interpolate velocity in time
    V = @(t,x) [interp1(Data.Time,Ocean_VelAtDepth_u_Mean,t);interp1(Data.Time,Ocean_VelAtDepth_v_Mean,t)];   % V is deterministic velocity vector, V_i = E[u_i]
    C = @(t) squeeze(interp1(Data.Time,VelocityCovarianceAtDepth,t));     % C is covariance tensor of velocity fluctuations, C_ij = E[u_i,u_j]

    % Check if B can be well defined. B*B' should be positive definite.
    Initial_C = C(Data.Time(1));
    Initial_B2 = A*Initial_C + Initial_C*A';
    [~,Initial_Eigenvalues] = eig(Initial_B2);
    if any(diag(Initial_Eigenvalues) < 0)
        disp('B2 is not positive definite');
        disp(diag(Initial_Eigenvalues))
        ErrorCode = 1;
        return
    end

    % Diffusion tensor
    function B_ = B(t)
        C_ = C(t);
        B2 = A*C_ + C_*A';
        B_ = sqrtm(B2);
    end

    NumBrownians = Dimension;

    % Choose to use zero-th or first order Markov model
    if Config.StochasticModel.UseFirstOrderMarkovModel == true

        % Drift Rate
        DriftRate = @(t,X) [V(t,X(1:2));zeros(NumBrownians,1)] + [[zeros(Dimension),eye(Dimension)];[zeros(Dimension),-A]] * X;    % First order Markov process
        % DriftRate = @(t,X) [V(t,X(1:2));zeros(NumBrownians,1)];                                                                  % No Brownian motion

        % Diffusion Rate
        DiffusionRate = @(t,X) [zeros(Dimension);B(t)];                                   % First order Markov process
        % DiffusionRate = @(t,X) [zeros(Dimension);zeros(NumBrownians)];                  % No Brownian motion

        % Start state
        StartState = [TrajectoryStartX(:);TrajectoryStartY(:);zeros(NumBrownians,1)];     % For First order Markov process

    else

        % Drift rate
        DriftRate = @(t,X) V(t,X);                                                % Zero-th order Markov process

        % Diffusion rate
        delta_t = Data.Time(2) - Data.Time(1);                                   % In second
        DiffusionRate = @(t,X) sqrtm(delta_t*Covariance);                        % Zero-th order Markov process

        % Start state
        StartState = [TrajectoryStartX;TrajectoryStartY];                         % Zero-th order Markov process

    end

    % SDE
    SDE = sde(DriftRate,DiffusionRate,'StartTime',TrajectoryStartTime,'StartState',StartState,'Correlation',eye(Dimension));

    % Simulate settings
    IntegrationDuration = InquiryTime - TrajectoryStartTime;
    dt = Data.Time(2) - Data.Time(1);                                   % In second
    NumPeriods = int64(floor((IntegrationDuration / dt)));              % N: Number of times. Total integration time T = N * dt.
    NumTrials = Config.StochasticModel.NumberOfSamples;                 % Number of Monte-Carlo sample draws from Gaussian distribution
    NumSteps = Config.StochasticModel.NumberOfRefinementSteps;          % Number of integration refinements inside each dt
    StateVectorSize = length(StartState);

    % Simulate
    rng('default');
    disp('Stochastic simulation started ...')

    % Run parallel or serial
    if Config.StochasticModel.UseParallelSimulation == true
        %Run in parallel
        RandomSamples = randn(NumPeriods*NumSteps,NumBrownians,NumTrials);
        StochasticTrajectoryCoordinates = zeros(NumPeriods+1,StateVectorSize,NumTrials);
        StochasticTrajectoryTimes = zeros(NumPeriods+1,NumTrials);
        parfor i = 1:NumTrials
            [StochasticTrajectoriesCoordinates(:,:,i),StochasticTrajectoriesTimes(:,i)] = simByEuler(SDE,NumPeriods,'DeltaTime',dt,'nSteps',NumSteps,'Antithetic',true,'Z',RandomSamples(:,:,i));
        end
        StochasticTrajectoriesTimes = StochasticTrajectoriesTimes(:,1);

    else
        % Run in serial
        [StochasticTrajectoriesCoordinates,StochasticTrajectoriesTimes] = simByEuler(SDE,NumPeriods,'DeltaTime',dt,'nTrials',NumTrials,'nSteps',NumSteps,'Antithetic',true);

    end
    disp('Stochastic simulation finished.')

    % Convert back XY to LonLat
    StochasticTrajectoriesCoordinates_X = squeeze(StochasticTrajectoriesCoordinates(:,1,:));
    StochasticTrajectoriesCoordinates_Y = squeeze(StochasticTrajectoriesCoordinates(:,2,:));
    [StochasticTrajectoriesLongitudes,StochasticTrajectoriesLatitudes] = Utilities.ConvertXYToLonLat(OriginLongitude,OriginLatitude,StochasticTrajectoriesCoordinates_X,StochasticTrajectoriesCoordinates_Y);

    Output.StochasticPlumeCenterTrajectories.Longitudes = squeeze(StochasticTrajectoriesLongitudes);
    Output.StochasticPlumeCenterTrajectories.Latitudes = squeeze(StochasticTrajectoriesLatitudes);
    Output.StochasticPlumeCenterTrajectories.Times = StochasticTrajectoriesTimes;

    % Deterministic trajectory
    opts = odeset('RelTol',1e-7);
    [TrajectoryTimes,TrajectoryCoordinates] = ...
        ode45(V,[TrajectoryStartTime,InquiryTime],[TrajectoryStartX,TrajectoryStartY],opts);

     % Convert back XY to LonLat
     [TrajectoryLongitudes,TrajectoryLatitudes] = Utilities.ConvertXYToLonLat(OriginLongitude,OriginLatitude,TrajectoryCoordinates(:,1),TrajectoryCoordinates(:,2));

    % Output: trajectory
    Output.PlumeCenterTrajectory.Longitudes = TrajectoryLongitudes;
    Output.PlumeCenterTrajectory.Latitudes = TrajectoryLatitudes;
    Output.PlumeCenterTrajectory.Times = TrajectoryTimes;

    % Ship trajectory during plume track
    Output.ShipTrajectory.Longitude = interp1(Data.Time,Data.Ship.Longitude,TrajectoryTimes,'pchip');
    Output.ShipTrajectory.Latitude = interp1(Data.Time,Data.Ship.Latitude,TrajectoryTimes,'pchip');

    % Output: last point on trajectory
    Output.PlumeCenterTrajectory.FirstLongitude = TrajectoryLongitudes(1);
    Output.PlumeCenterTrajectory.FirstLatitude = TrajectoryLatitudes(1);
    Output.PlumeCenterTrajectory.LastLongitude = TrajectoryLongitudes(end);
    Output.PlumeCenterTrajectory.LastLatitude = TrajectoryLatitudes(end);
    Output.PlumeCenterTrajectory.LastTime = TrajectoryTimes(end);

    % Kernel density estimation of end point of stochastic trajectories
    [PDF_ValuesOnGrid,PDF_LongitudeGrid,PDF_LatitudeGrid,PDF_XGrid,PDF_YGrid] = ProbabilityDensityFunctions.Estimate2DKernelDensityFunction(StochasticTrajectoriesLongitudes,StochasticTrajectoriesLatitudes);
    Output.StochasticPDF.PDF_ValuesOnGrid = PDF_ValuesOnGrid;
    Output.StochasticPDF.LongitudeGrid = PDF_LongitudeGrid;
    Output.StochasticPDF.LatitudeGrid = PDF_LatitudeGrid;

    % Compute Complementrary Cumulative Distribution Function
    CCDF_ValuesOnGrid = ProbabilityDensityFunctions.ComplementaryCumulativeDistributionFunction(PDF_XGrid,PDF_YGrid,PDF_ValuesOnGrid);
    Output.StochasticPDF.CCDF_ValuesOnGrid = CCDF_ValuesOnGrid;

    end

    % =====================
    % Plot Plume Trajectory
    % =====================

    function PlotPlumeTrajectory(Data,Config,Output)

    fig = figure;
    OceanColor = [.6 .8 .9];
    ax = axesm('MapProjection','mercator');
    setm(ax,'FFaceColor',OceanColor);

    % Ship whole trajectory
    % h1 = plot(Data.Ship.Longitude,Data.Ship.Latitude,'-o','color','black','DisplayName','Ship trajectory');
    h1 = plotm(Data.Ship.Latitude,Data.Ship.Longitude,'-o','color','black','linewidth',0.5,'MarkerSize',1,'DisplayName','Ship trajectory');
    hold on

    % Plot stochastic trajectories
    NumStochasticTrajectories = size(Output.StochasticPlumeCenterTrajectories.Longitudes,2);
    StochasticTrajectoryColor = [0.7 0.7 1];
    StochasticPointsColor = [0.65 0.65 1];
    for Trial = 1:NumStochasticTrajectories

        % Stochastic trajectory
        h6 = plotm(squeeze(Output.StochasticPlumeCenterTrajectories.Latitudes(:,Trial)), ...
        squeeze(Output.StochasticPlumeCenterTrajectories.Longitudes(:,Trial)), ...
        '-','color',StochasticTrajectoryColor,'MarkerEdgeColor',StochasticTrajectoryColor,'MarkerFaceColor',StochasticTrajectoryColor,'MarkerSize',2);

        % End of stochastic trajectory
        plotm(squeeze(Output.StochasticPlumeCenterTrajectories.Latitudes(end,Trial)), ...
        squeeze(Output.StochasticPlumeCenterTrajectories.Longitudes(end,Trial)), ...
        'o','MarkerEdgeColor',StochasticPointsColor,'MarkerFaceColor',StochasticPointsColor,'MarkerSize',2);

        if Trial == NumStochasticTrajectories
            h6.set('DisplayName','Stochastic trajectories');
        end
    end

    % Plot Stochastic PDF (in unit of km^{-2})
    pcolorm(Output.StochasticPDF.LatitudeGrid,Output.StochasticPDF.LongitudeGrid,Output.StochasticPDF.PDF_ValuesOnGrid*1e6,'AlphaDataMapping','scaled','AlphaData',Output.StochasticPDF.PDF_ValuesOnGrid,'FaceAlpha','flat')
    cmap = parula(1024);
    colormap(cmap);
    c = colorbar;
    c.Label.String = 'Probability density of plume (km^{-2})';

    % Plot Stochastic CCDF
    [C,hc] = contourm(Output.StochasticPDF.LatitudeGrid,Output.StochasticPDF.LongitudeGrid,Output.StochasticPDF.CCDF_ValuesOnGrid,linspace(0.1,0.9,9),'--');
    hct = clabelm(C,hc);
    set(hct,'BackgroundColor','none','FontSize',9)

    % Set back the gcf, since the previous command (contourm) loses it and makes next error commands
    figure(fig);
    
    % Translated StreakLine Position to the observed time
    PlumeObservedStreakLineLongitude = Output.PlumeCenterTrajectory.FirstLongitude - Output.PlumeInitialStreakLine.CenterLongitudeAtEndReleaseTime + Output.PlumeInitialStreakLine.LongitudeAtEndReleaseTime;
    PlumeObservedStreakLineLatitude = Output.PlumeCenterTrajectory.FirstLatitude - Output.PlumeInitialStreakLine.CenterLatitudeAtEndReleaseTime + Output.PlumeInitialStreakLine.LatitudeAtEndReleaseTime;

    % plot translated plume streakline
    h4 = plotm(PlumeObservedStreakLineLatitude,PlumeObservedStreakLineLongitude,'color',[1,0.6,0.6],'linewidth',2,'DisplayName','Plume streakline at observed time');

    % Translated StreakLine Position to the inquiry time
    PlumeCurrentStreakLineLongitude = Output.PlumeCenterTrajectory.LastLongitude - Output.PlumeInitialStreakLine.CenterLongitudeAtEndReleaseTime + Output.PlumeInitialStreakLine.LongitudeAtEndReleaseTime;
    PlumeCurrentStreakLineLatitude = Output.PlumeCenterTrajectory.LastLatitude - Output.PlumeInitialStreakLine.CenterLatitudeAtEndReleaseTime + Output.PlumeInitialStreakLine.LatitudeAtEndReleaseTime;

    % plot translated plume streakline
    h5 = plotm(PlumeCurrentStreakLineLatitude,PlumeCurrentStreakLineLongitude,'color','red','linewidth',2,'DisplayName','Plume streakline at inquiry time');

    % Ship trajectory during plume release
    ShipTrajectoryColor1 = [0.06,0.65,0.95];
    h70 = plotm(Output.PlumeInitialStreakLine.LatitudeAtStartReleaseTime,Output.PlumeInitialStreakLine.LongitudeAtStartReleaseTime,'-','color',ShipTrajectoryColor1,'linewidth',1.5,'DisplayName','Ship trajectory during plume release');
    h71 = plotm(Output.PlumeInitialStreakLine.LatitudeAtStartReleaseTime(1),Output.PlumeInitialStreakLine.LongitudeAtStartReleaseTime(1),'o','MarkerEdgeColor','black','MarkerFaceColor',ShipTrajectoryColor1,'MarkerSize',7,'DisplayName','Ship at start of plume release');
    h72 = plotm(Output.PlumeInitialStreakLine.LatitudeAtStartReleaseTime(end),Output.PlumeInitialStreakLine.LongitudeAtStartReleaseTime(end),'v','MarkerEdgeColor','black','MarkerFaceColor',ShipTrajectoryColor1,'MarkerSize',7,'DisplayName','Ship at end of plume release');

    % Ship trajectory during plume track
    ShipTrajectoryColor2 = [0.06,0.65,0.06];
    h30 = plotm(Output.ShipTrajectory.Latitude,Output.ShipTrajectory.Longitude,'-','color',ShipTrajectoryColor2,'linewidth',1.5,'DisplayName','Ship trajectory during plume track');
    h31 = plotm(Output.ShipTrajectory.Latitude(1),Output.ShipTrajectory.Longitude(1),'o','MarkerEdgeColor','black','MarkerFaceColor',ShipTrajectoryColor2,'MarkerSize',7,'DisplayName','Ship at start of plume track');
    h32 = plotm(Output.ShipTrajectory.Latitude(end),Output.ShipTrajectory.Longitude(end),'v','MarkerEdgeColor','black','MarkerFaceColor',ShipTrajectoryColor2,'MarkerSize',7,'DisplayName','Ship at end of plume track');
    textm(Output.ShipTrajectory.Latitude(1),Output.ShipTrajectory.Longitude(1),strcat(' Ship at: ',datestr(Output.PlumeCenterTrajectory.Times(1)/(24*3600),' yyyy/mm/dd HH:MM')));
    textm(Output.ShipTrajectory.Latitude(end),Output.ShipTrajectory.Longitude(end),strcat(' Ship at: ',datestr(Output.PlumeCenterTrajectory.Times(end)/(24*3600),' yyyy/mm/dd HH:MM')));

    % Plume center trajectory
    PlumeCenterColor = [0.3,0.3,1];
    h20 = plotm(Output.PlumeCenterTrajectory.Latitudes,Output.PlumeCenterTrajectory.Longitudes,'color',PlumeCenterColor,'linewidth',2,'DisplayName','Center plume trajectory');
    h21 = plotm(Output.PlumeCenterTrajectory.Latitudes(1),Output.PlumeCenterTrajectory.Longitudes(1),'o','MarkerEdgeColor','black','MarkerFaceColor',PlumeCenterColor,'MarkerSize',7,'DisplayName','Start of plume track');
    h22 = plotm(Output.PlumeCenterTrajectory.Latitudes(end),Output.PlumeCenterTrajectory.Longitudes(end),'v','MarkerEdgeColor','black','MarkerFaceColor',PlumeCenterColor,'MarkerSize',7,'DisplayName','End of plume track');
    textm(Output.PlumeCenterTrajectory.Latitudes(1),Output.PlumeCenterTrajectory.Longitudes(1),strcat(' Plume at: ',datestr(Output.PlumeCenterTrajectory.Times(1)/(24*3600),' yyyy/mm/dd HH:MM')));
    textm(Output.PlumeCenterTrajectory.Latitudes(end),Output.PlumeCenterTrajectory.Longitudes(end),strcat(' Plume at: ',datestr(Output.PlumeCenterTrajectory.Times(end)/(24*3600),' yyyy/mm/dd HH:MM')));

    % Plot Drifters
    DriftersFileData = load('/home/sia/work/alpha/WhereIsMyPlume/Drifters.mat');
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

        if DrifterTime(1) <= Data.Time(1) & DrifterTime(end) >= Data.Time(end)

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
    MinLon = min([PlumeObservedStreakLineLongitude(:);PlumeCurrentStreakLineLongitude(:);Output.ShipTrajectory.Longitude(:);Data.Ship.Longitude(:);Output.StochasticPDF.LongitudeGrid(1,1)]);
    MaxLon = max([PlumeObservedStreakLineLongitude(:);PlumeCurrentStreakLineLongitude(:);Output.ShipTrajectory.Longitude(:);Data.Ship.Longitude(:);Output.StochasticPDF.LongitudeGrid(1,end)]);
    MinLat = min([PlumeObservedStreakLineLatitude(:);PlumeCurrentStreakLineLatitude(:);Output.ShipTrajectory.Latitude(:);Data.Ship.Latitude(:);Output.StochasticPDF.LatitudeGrid(1,1)]);
    MaxLat = max([PlumeObservedStreakLineLatitude(:);PlumeCurrentStreakLineLatitude(:);Output.ShipTrajectory.Latitude(:);Data.Ship.Latitude(:);Output.StochasticPDF.LatitudeGrid(end,1)]);

    % MinLon = min([PlumeObservedStreakLineLongitude(:);PlumeCurrentStreakLineLongitude(:);Output.ShipTrajectory.Longitude(:);Data.Ship.Longitude(:)]);
    % MaxLon = max([PlumeObservedStreakLineLongitude(:);PlumeCurrentStreakLineLongitude(:);Output.ShipTrajectory.Longitude(:);Data.Ship.Longitude(:)]);
    % MinLat = min([PlumeObservedStreakLineLatitude(:);PlumeCurrentStreakLineLatitude(:);Output.ShipTrajectory.Latitude(:);Data.Ship.Latitude(:)]);
    % MaxLat = max([PlumeObservedStreakLineLatitude(:);PlumeCurrentStreakLineLatitude(:);Output.ShipTrajectory.Latitude(:);Data.Ship.Latitude(:)]);

    LatExt = MaxLat - MinLat;
    LonExt = MaxLon - MinLon;
    Ext = max(LatExt,LonExt);
    Ratio = 0.0;

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

    setm(gca,'MLineLocation',0.01,'MLabelLocation',0.01,'MLabelRound',-2,'MeridianLabel','on','MLabelParallel','south');
    setm(gca,'PLineLocation',0.01,'PLabelLocation',0.01,'PLabelRound',-2,'ParallelLabel','on','PLabelMeridian','west');

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

    legend([h1,h20,h21,h22,h30,h31,h32,h4,h5,h6,h70,h71,h72])

    % Font
    set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

    % Save
    saveas(gcf,fullfile(Config.Plots.FiguresDirectory,strcat('StochasticTrajectories',Config.Plots.FiguresFormat)))

    end

    % ---------------------
    % End of static methods

    end
end
