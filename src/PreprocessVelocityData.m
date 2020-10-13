classdef PreprocessVelocityData

    methods(Static)

    % ========================================================================
    % Process Velocity Data - First Interpolate In Depth Then Smooth Over Time
    % ========================================================================

    function OceanVelocityData = Process_FirstInterpolateInDepthThenSmoothOverTime(Config,Data,PlumeDepth)

        % This function, first interpolates at depth, then smooth over time.

        % Check Depth to interpolate or extrapolate column profile of velocity at plume depth
        if PlumeDepth > Data.Ocean.Range(end)
            % Requested depth is not in the data range. 
            error('Plume depth: %f is not in the data range from %f(m) to %f(m).',PlumeDepth,Data.Ocean.Range(1),Data.Ocean.Range(end))
        elseif PlumeDepth < Data.Ocean.Range(1)
            % Extrapolate column of velocity profile at plume depth
            % [Ocean_VelAtDepth_u,Ocean_VelAtDepth_v] = PreprocessVelocityData.ExtrapolateVelocityNearSurface(Config,Data,Data.Ocean.AbsoluteVel_u,Data.Ocean.AbsoluteVel_v,PlumeDepth);
            [Ocean_VelAtDepth_u,Ocean_VelAtDepth_v] = PreprocessVelocityData.ExtrapolateVelocityNearSurface_ZeroShearAtOrigin(Config,Data,Data.Ocean.AbsoluteVel_u,Data.Ocean.AbsoluteVel_v,PlumeDepth);
        else
            % Interpolate column of velocity profile at plume depth
            Ocean_VelAtDepth_u = interp1(Data.Ocean.Range,Data.Ocean.AbsoluteVel_u,PlumeDepth,'linear','extrap');
            Ocean_VelAtDepth_v = interp1(Data.Ocean.Range,Data.Ocean.AbsoluteVel_v,PlumeDepth,'linear','extrap');
        end

        % Include wind to the ocean velocity
        if Config.Wind.IncludeWindData == true
            [Ocean_VelAtDepth_u,Ocean_VelAtDepth_v] = PreprocessVelocityData.IncludeWindToOceanVelocityData(Config,Data,Ocean_VelAtDepth_u,Ocean_VelAtDepth_v,PlumeDepth);
        end

        % Smooth velocity along time
        [Ocean_VelAtDepth_u_Mean,Ocean_VelAtDepth_v_Mean,Ocean_VelAtDepth_u_Noise,Ocean_VelAtDepth_v_Noise] = ...
        PreprocessVelocityData.SmoothVelocityAlongTime(Config,Data,Ocean_VelAtDepth_u,Ocean_VelAtDepth_v);

        % Store all outputs in a struct
        OceanVelocityData.Ocean_VelAtDepth_u = Ocean_VelAtDepth_u;
        OceanVelocityData.Ocean_VelAtDepth_v = Ocean_VelAtDepth_v;
        OceanVelocityData.Ocean_VelAtDepth_u_Mean = Ocean_VelAtDepth_u_Mean;
        OceanVelocityData.Ocean_VelAtDepth_v_Mean = Ocean_VelAtDepth_v_Mean;
        OceanVelocityData.Ocean_VelAtDepth_u_Noise = Ocean_VelAtDepth_u_Noise;
        OceanVelocityData.Ocean_VelAtDepth_v_Noise = Ocean_VelAtDepth_v_Noise;

    end

    % ========================================================================
    % Process Velocity Data - First Smooth Over Time Then Interpolate In Depth
    % ========================================================================

    function OceanVelocityData = Process_FirstSmoothOverTimeThenInterpolateInDepth(Config,Data,PlumeDepth)

        % This function, first, smooths over time, them interpolates at depth

        % Smooth velocity along time
        [Ocean_Vel_u_Mean,Ocean_Vel_v_Mean,Ocean_Vel_u_Noise,Ocean_Vel_v_Noise] = ...
        PreprocessVelocityData.SmoothVelocityAlongTime(Config,Data,Data.Ocean.AbsoluteVel_u,Data.Ocean.AbsoluteVel_v);

        % Check Depth to interpolate or extrapolate column profile of velocity at plume depth
        if PlumeDepth > Data.Ocean.Range(end)
            % Requested depth is not in the data range. 
            error('Plume depth: %f is not in the data range from %f(m) to %f(m).',PlumeDepth,Data.Ocean.Range(1),Data.Ocean.Range(end))
        elseif PlumeDepth < Data.Ocean.Range(1)
            % Extrapolate column of velocity profile at plume depth
            [Ocean_VelAtDepth_u,Ocean_VelAtDepth_v] = PreprocessVelocityData.ExtrapolateVelocityNearSurface(Config,Data,Data.Ocean.AbsoluteVel_u,Data.Ocean.AbsoluteVel_v,PlumeDepth);

            % [Ocean_VelAtDepth_u_Mean,Ocean_VelAtDepth_v_Mean] = PreprocessVelocityData.ExtrapolateVelocityNearSurface(Config,Data,Ocean_Vel_u_Mean,Ocean_Vel_v_Mean,PlumeDepth);
            [Ocean_VelAtDepth_u_Mean,Ocean_VelAtDepth_v_Mean] = PreprocessVelocityData.ExtrapolateVelocityNearSurface_ZeroShearAtOrigin(Config,Data,Ocean_Vel_u_Mean,Ocean_Vel_v_Mean,PlumeDepth);

            Ocean_VelAtDepth_u_Noise = Ocean_VelAtDepth_u - Ocean_VelAtDepth_u_Mean;
            Ocean_VelAtDepth_v_Noise = Ocean_VelAtDepth_v - Ocean_VelAtDepth_v_Mean;

        else
            % Interpolate column of velocity profile at plume depth
            Ocean_VelAtDepth_u = interp1(Data.Ocean.Range,Data.Ocean.AbsoluteVel_u,PlumeDepth,'linear','extrap');
            Ocean_VelAtDepth_v = interp1(Data.Ocean.Range,Data.Ocean.AbsoluteVel_v,PlumeDepth,'linear','extrap');

            Ocean_VelAtDepth_u_Mean = interp1(Data.Ocean.Range,Ocean_Vel_u_Mean,PlumeDepth,'linear','extrap');
            Ocean_VelAtDepth_v_Mean = interp1(Data.Ocean.Range,Ocean_Vel_v_Mean,PlumeDepth,'linear','extrap');

            Ocean_VelAtDepth_u_Noise = interp1(Data.Ocean.Range,Ocean_Vel_u_Noise,PlumeDepth,'linear','extrap');
            Ocean_VelAtDepth_v_Noise = interp1(Data.Ocean.Range,Ocean_Vel_v_Noise,PlumeDepth,'linear','extrap');
        end

        % Include wind to the ocean velocity
        if Config.Wind.IncludeWindData == true
            [Ocean_VelAtDepth_u_Mean,Ocean_VelAtDepth_v_Mean] = PreprocessVelocityData.IncludeWindToOceanVelocityData(Config,Data,Ocean_VelAtDepth_u_Mean,Ocean_VelAtDepth_v_Mean,PlumeDepth);
        end


        % Store all outputs in a struct
        OceanVelocityData.Ocean_VelAtDepth_u = Ocean_VelAtDepth_u;
        OceanVelocityData.Ocean_VelAtDepth_v = Ocean_VelAtDepth_v;
        OceanVelocityData.Ocean_VelAtDepth_u_Mean = Ocean_VelAtDepth_u_Mean;
        OceanVelocityData.Ocean_VelAtDepth_v_Mean = Ocean_VelAtDepth_v_Mean;
        OceanVelocityData.Ocean_VelAtDepth_u_Noise = Ocean_VelAtDepth_u_Noise;
        OceanVelocityData.Ocean_VelAtDepth_v_Noise = Ocean_VelAtDepth_v_Noise;

    end

    % ===================================
    % Include Wind To Ocean Velocity Data
    % ===================================

    function [Ocean_VelAtDepth_u,Ocean_VelAtDepth_v] = IncludeWindToOceanVelocityData(Config,Data,Ocean_VelAtDepth_u,Ocean_VelAtDepth_v,PlumeDepth)

        % Check if depth is non-zero
        if PlumeDepth < 3.0
            Ocean_VelAtDepth_u = Ocean_VelAtDepth_u + Config.Wind.WindageCoefficient * Data.Wind.OnShipTrack.EastWind;
            Ocean_VelAtDepth_v = Ocean_VelAtDepth_v + Config.Wind.WindageCoefficient * Data.Wind.OnShipTrack.NorthWind;
        else
            error('The wind can only be added to ocean surface data where the depth is zero. Current depth: %f',PlumeDepth);
            error('Set Config.Wind.IncludeWindData to false.');
        end

    end

    % =================================
    % Extrapolate Velocity Near Surface
    % =================================

    function [Ocean_VelAtDepth_u,Ocean_VelAtDepth_v] = ExtrapolateVelocityNearSurface(Config,Data,Ocean_Vel_u,Ocean_Vel_v,PlumeDepth)

        % Polynomial of order N-1, Uses the first N points in Data.Ocean.Range

        % Use first N data points in range
        N = 2;  % Polynomial of order N-1

        % Matrix of coefficients for solving a linear system
        Vandermond = vander([Data.Ocean.Range(1:N)']);

        % Solving the system
        InvVandermond = inv(Vandermond);

        % Known values
        Ocean_VelAtDepth_u = zeros(1,size(Data.Ocean.AbsoluteVel_u,2));
        Ocean_VelAtDepth_v = zeros(1,size(Data.Ocean.AbsoluteVel_u,2));

        Powers = N-1:-1:0;

        for i = 1:N
            Ocean_VelAtDepth_u = Ocean_VelAtDepth_u + PlumeDepth^Powers(i) * InvVandermond(i,:) * Ocean_Vel_u(1:N,:);
            Ocean_VelAtDepth_v = Ocean_VelAtDepth_v + PlumeDepth^Powers(i) * InvVandermond(i,:) * Ocean_Vel_v(1:N,:);
        end

    end

    % ========================================================
    % Extrapolate Velocity Near Surface - Zero Shear At Origin
    % ========================================================

    function [Ocean_VelAtDepth_u,Ocean_VelAtDepth_v] = ExtrapolateVelocityNearSurface_ZeroShearAtOrigin(Config,Data,Ocean_Vel_u,Ocean_Vel_v,PlumeDepth)

        % Polynomial of order N, zero-slope at origin. Uses N+1 points the first N points in Data.Ocean.Range, and a slope condition at surface

        % Use first N data points in range
        N = 2;  % Polynomial of order N

        % Matrix of coefficients for solving a linear system
        Vandermond = vander([0,Data.Ocean.Range(1:N)']);

        % We want the polynomial to have zero slope at origin (zero-shear flow at ocean surface). So, we remove end-1 column, corresponding to the term x^1
        Vandermond = Vandermond(2:end,[1:end-2,end]);

        % Solving the system
        InvVandermond = inv(Vandermond);

        % Known values
        Ocean_VelAtDepth_u = zeros(1,size(Data.Ocean.AbsoluteVel_u,2));
        Ocean_VelAtDepth_v = zeros(1,size(Data.Ocean.AbsoluteVel_u,2));

        Powers = N:-1:1;
        Powers(end) = 0;

        % for i = 1:N
        %     Ocean_VelAtDepth_u = Ocean_VelAtDepth_u + PlumeDepth^Powers(i) * InvVandermond(i,:) * Ocean_Vel_u(1:N,:);
        %     Ocean_VelAtDepth_v = Ocean_VelAtDepth_v + PlumeDepth^Powers(i) * InvVandermond(i,:) * Ocean_Vel_v(1:N,:);
        % end

        x1 = Data.Ocean.Range(1);
        x2 = Data.Ocean.Range(2);

        % Quadratic from first two points and zero slope at origin
        % Ocean_VelAtDepth_u = (-x2^2*Ocean_Vel_u(1,:) + x1^2*Ocean_Vel_u(2,:))/(x1^2-x2^2);
        % Ocean_VelAtDepth_v = (-x2^2*Ocean_Vel_v(1,:) + x1^2*Ocean_Vel_v(2,:))/(x1^2-x2^2);

        % Same value as first point
        % Ocean_VelAtDepth_u = Ocean_Vel_u(1,:);
        % Ocean_VelAtDepth_v = Ocean_Vel_v(1,:);

        % Linear from first two points
        % Ocean_VelAtDepth_u = (-x2*Ocean_Vel_u(1,:) + x1*Ocean_Vel_u(2,:))/(x1-x2);
        % Ocean_VelAtDepth_v = (-x2*Ocean_Vel_v(1,:) + x1*Ocean_Vel_v(2,:))/(x1-x2);

        % Least square fitting a line (first order) from the first N points
        NumPoints = Config.PreprocessVelocity.NumberOfPointsForRegression;
        PolynomialDegree = Config.PreprocessVelocity.PolynomialOrderForRegression;
        A = vander(Data.Ocean.Range(1:NumPoints));
        A = A(:,end-PolynomialDegree:end);

        % Make polynomial to have zero slope at origin
        % A = A(:,[1:end-2,end]);

        C = inv(A'*A)*A';

        Ocean_VelAtDepth_u = C(end,:)*Ocean_Vel_u(1:NumPoints,:);
        Ocean_VelAtDepth_v = C(end,:)*Ocean_Vel_v(1:NumPoints,:);

    end

    % ==========================
    % Smooth Velocity Along Time
    % ==========================

    function [Ocean_VelAtDepth_u_Mean,Ocean_VelAtDepth_v_Mean,Ocean_VelAtDepth_u_Noise,Ocean_VelAtDepth_v_Noise] = ...
        SmoothVelocityAlongTime( ...
        Config,Data, ...
        Ocean_VelAtDepth_u, ...
        Ocean_VelAtDepth_v)

        % Mean of the data is found by Gaussian convolution filter with prespecified window length.
        % The residual of original data with the mean is called noise.

        DayToSecond = 24 * 60 * 60;

        % Moving average filter for velocity over time
        delta_t = Data.Time(2) - Data.Time(1);                                             % (in unit of seconds)
        MeanWindow = DayToSecond * datenum(Config.SmoothingDurationAlongTime);             % Filtering window duration (in unit of seconds)
        MeanWindowLength = floor(MeanWindow / delta_t);

        if MeanWindowLength < 1
            error('Moving average filter window is less than one.')
        end

        % Ocean_VelAtDepth_u_Mean = movmean(Ocean_VelAtDepth_u,MeanWindowLength,2);          % Filters fluctuation over time window
        % Ocean_VelAtDepth_v_Mean = movmean(Ocean_VelAtDepth_v,MeanWindowLength,2);

        Ocean_VelAtDepth_u_Mean = smoothdata(Ocean_VelAtDepth_u,2,'gaussian',MeanWindowLength);          % Filters fluctuation over time window
        Ocean_VelAtDepth_v_Mean = smoothdata(Ocean_VelAtDepth_v,2,'gaussian',MeanWindowLength);

        Ocean_VelAtDepth_u_Noise = Ocean_VelAtDepth_u - Ocean_VelAtDepth_u_Mean;           % Fluctuations of velocities
        Ocean_VelAtDepth_v_Noise = Ocean_VelAtDepth_v - Ocean_VelAtDepth_v_Mean;

    end


    % =======================================
    % Estimate Standard Deviation of Velocity
    % =======================================
    
    function [StandardDeviation_U,StandardDeviation_V,VelocityStandardDeviationAverage] = EstimateStandardDeviationOfVelocity( ...
        Ocean_VelAtDepth_u_Noise, ...
        Ocean_VelAtDepth_v_Noise)

        % Uses Ljung-Box Q-test to test if the velocity data is close to a Gaussian distribution.
        % lbqtest needs Econometrix toolbox to be nstalled.

        % Ljung-Box Q-test
        h_u = lbqtest(Ocean_VelAtDepth_u_Noise);
        h_v = lbqtest(Ocean_VelAtDepth_v_Noise);

        if (h_u == 0)
            % error('Ljung-Box Q-text on east velocity is negative, meaning that the noise of east velocity is not stationary Gaussian.');
            warning('Ljung-Box Q-text on east velocity is negative, meaning that the noise of east velocity is not stationary Gaussian.');
        elseif (h_v == 0)
            % error('Ljung-Box Q-text on north velocity is negative, meaning that the noise of north velocity is not stationary Gaussian.');
            warning('Ljung-Box Q-text on north velocity is negative, meaning that the noise of north velocity is not stationary Gaussian.');
        end

        % Standard deviation of noise assuming noise of u and v are stationary gaussian
        StandardDeviation_U = std(Ocean_VelAtDepth_u_Noise);
        StandardDeviation_V = std(Ocean_VelAtDepth_v_Noise);

        % Average standard deviations for both east and north velocities
        VelocityStandardDeviationAverage = 0.5 * (StandardDeviation_U + StandardDeviation_V);

    end

    % ---------------------
    % End of Static methods

    end
end
