classdef GaussianProcessRegression

    methods(Static)

    % ==============================
    % Krig Velocity At Inquiry Point
    % ==============================

    function Vel_InquiryPoints = KrigVelocityAtInquiryPoint( ...
        Config, ...
        X_InquiryPoints,T_InquiryPoint, ...
        X_DataPoints,T_DataPoints, ...
        DataCovariance,F, ...
        TemporalDecorrelationScale,TemporalDecorrelationCutoff, ...
        SpatialDecorrelationScale,SpatialDecorrelationCutoff, ...
        U_DataPoints,V_DataPoints)

        % References:
        % Wackernagel, Multivariate geostatistis
        % 1. Chapter 3:  Simple kriging, assumes mean is constant (stationary) and known
        % 2. Chapter 11: Ordinary Kriging (OK), assumes mean is constant (stationary), but not known. We solve for mean as well.
        % 3. Chapter 37: External drift, assumes mean is not constant (non-stationary), and is a constant plus a linear drift term
        % 4. Chapter 38: Universal kriging (UK) (also known as Kriging with Trend, KT), assumes mean is a function of some known basis functions with unknown coefficidents
        % 5. Chapter 39: Translation Invariant Drift, a special KT with basis functions that are translation invariant.
        %
        % What I am using here is Kriging with drift (KT) with translation invariant basis functions.
        % 
        % U_DataPoints and V_DataPoints are Nx1 shape matrices (a column vector)

        % Number of data points
        N = size(X_DataPoints,1);

        % Cutoff future time (choose either of two options below)
        if Config.Kriging.UseFutureTimesInKriging == true
            % uses past and future data
            CurrentTimeIndex = N;
        else
            % Do not use future times. Cut off the array till the current time
            CurrentTimeIndex = find(T_DataPoints > T_InquiryPoint,1,'first');
        end

        % Number of Inquiry points (this is now a 1D array of size 2M)
        M = length(X_InquiryPoints) / 2;

        % Reshape inquiry points to be M*2
        X_InquiryPoints = reshape(X_InquiryPoints,M,2);

        % Covariance between inquiry point and the cluster of points
        k = CovarianceFunction.CovarianceOfInquiryPointsWithDataPoints( ...
            Config, ...
            X_InquiryPoints,T_InquiryPoint, ...
            X_DataPoints(1:CurrentTimeIndex,:),T_DataPoints(1:CurrentTimeIndex), ...
            TemporalDecorrelationScale,TemporalDecorrelationCutoff, ...
            SpatialDecorrelationScale,SpatialDecorrelationCutoff);

        % Basis function at inquiry point
        f = BasisFunctions.GenerateBasisFunctions(Config,X_InquiryPoints,T_InquiryPoint);

        % Number of basis functions
        l = size(F,2);

        % Kriging matrices
        C = [DataCovariance(1:CurrentTimeIndex,1:CurrentTimeIndex),F(1:CurrentTimeIndex,:);F(1:CurrentTimeIndex,:)',zeros(l,l)];
        D = [k;f'];

        % Solve for weights
        Solution = C\D;

        % Extract weights from solution. Solution = [Weights;LagrangianMultiplier].
        Weights = Solution(1:CurrentTimeIndex,:);

        % Interpolation
        U_InquiryPoints = Weights'*U_DataPoints(1:CurrentTimeIndex)';
        V_InquiryPoints = Weights'*V_DataPoints(1:CurrentTimeIndex)';

        Vel_InquiryPoints = [U_InquiryPoints;V_InquiryPoints];

    end

    % ======================
    % Dual Krig Coefficients
    % ======================

    function [b_u,b_v,d_u,d_v] = DualKrigCoefficients( ...
        Config, ...
        X_DataPoints,T_DataPoints, ...
        DataCovariance,F, ...
        U_DataPoints,V_DataPoints)

        % References:
        % Wackernagel, Multivariate geostatistis
        % See chapter 39: Translation Invariant Drift, in particual page 313-314 for Dual Kriging and page 314-315 for spline.

        % Number of data points
        N = size(X_DataPoints,1);

        % Cutoff future time (choose either of two options below)
        if Config.Kriging.UseFutureTimesInKriging == true
            % uses past and future data
            CurrentTimeIndex = N;
        else
            % Do not use future times. Cut off the array till the current time
            % CurrentTimeIndex = find(T_DataPoints > T_InquiryPoint,1,'first');

            error('In the current implementation, when dual kriging is enabled, all past and future data of the ship measurements should be used.')
        end

        % Number of basis functions
        m = size(F,2);

        % Kriging matrices
        StiffnessFactor = Config.Kriging.SplineStifnessInDualKriging;      % Spline stiffness weight
        Lambda = StiffnessFactor*eye(CurrentTimeIndex,CurrentTimeIndex);   % Spline smoothing, added as a diagonal matrix to the covariance K
        C = [Lambda+DataCovariance(1:CurrentTimeIndex,1:CurrentTimeIndex),F(1:CurrentTimeIndex,:);F(1:CurrentTimeIndex,:)',zeros(m,m)];
        Z_u = [U_DataPoints';zeros(m,1)];
        Z_v = [V_DataPoints';zeros(m,1)];

        % Solve for column vector bd = [b,d] of both b and d
        bd_u = C\Z_u;
        bd_v = C\Z_v;

        % Extract b and d from column vector bd
        b_u = bd_u(1:CurrentTimeIndex);
        d_u = bd_u(CurrentTimeIndex+1:end);

        b_v = bd_v(1:CurrentTimeIndex);
        d_v = bd_v(CurrentTimeIndex+1:end);

    end

    % ===================================
    % Dual Krig Velocity At Inquiry Point
    % ===================================

    function Vel_InquiryPoints = DualKrigVelocityAtInquiryPoint( ...
        Config, ...
        X_InquiryPoints,T_InquiryPoint, ...
        X_DataPoints,T_DataPoints, ...
        b_u,b_v,d_u,d_v, ...
        TemporalDecorrelationScale,TemporalDecorrelationCutoff, ...
        SpatialDecorrelationScale,SpatialDecorrelationCutoff)

        % References:
        % Wackernagel, Multivariate geostatistis
        % Chapter 39, Translation invariant drift

        % Number of data points
        N = size(X_DataPoints,1);

        % Cutoff future time (choose either of two options below)
        if Config.Kriging.UseFutureTimesInKriging == true
            % uses past and future data
            CurrentTimeIndex = N;
        else
            % Do not use future times. Cut off the array till the current time
            CurrentTimeIndex = find(T_DataPoints > T_InquiryPoint,1,'first');
        end

        % Number of Inquiry points (this is now a 1D array of size 2M)
        M = length(X_InquiryPoints) / 2;

        % Reshape inquiry points to be M*2
        X_InquiryPoints = reshape(X_InquiryPoints,M,2);

        % Covariance betwene inquiry point and the cluster of points
        k = CovarianceFunction.CovarianceOfInquiryPointsWithDataPoints( ...
            Config, ...
            X_InquiryPoints,T_InquiryPoint, ...
            X_DataPoints(1:CurrentTimeIndex,:),T_DataPoints(1:CurrentTimeIndex), ...
            TemporalDecorrelationScale,TemporalDecorrelationCutoff, ...
            SpatialDecorrelationScale,SpatialDecorrelationCutoff);
        k = k(1:CurrentTimeIndex,:);

        % Basis function at inquiry point
        f = BasisFunctions.GenerateBasisFunctions(Config,X_InquiryPoints,T_InquiryPoint);

        % Interpolation
        U_InquiryPoints = k'*b_u + f*d_u;
        V_InquiryPoints = k'*b_v + f*d_v;

        Vel_InquiryPoints = [U_InquiryPoints;V_InquiryPoints];

    end

    % ---------------------
    % End of Static methods

    end
end
