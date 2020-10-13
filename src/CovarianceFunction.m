classdef CovarianceFunction

    methods(Static)

    % ====================
    % Temporal Correlation
    % ====================

    function TemporalCorrelation = TemporalCorrelation( ...
        T_InquiryPoint, ...
        T_DataPoints, ...
        TemporalDecorrelationScale, ...
        TemporalDecorrelationCutoff)

        % Format of t
        % t = [t_1,...t_n]'
        % t is in the unit of seconds
        % TemporalDecorrelationScale is scalar in unit of seconds

        % Distance
        TemporalDistance = T_DataPoints(:) - repmat(T_InquiryPoint,length(T_DataPoints),1);

        % Normalize to decorrelation scale
        ScaledTemporalDistance = TemporalDistance / TemporalDecorrelationScale;

        % Correlation kernel (using exponential decay)
        TemporalCorrelation = exp(-abs(ScaledTemporalDistance));

        % Zero-out distance larger than cutoff
        TemporalCorrelation(TemporalDistance > TemporalDecorrelationCutoff) = 0.0;

    end

    % ======================================================
    % Spatial Correlation of Inquiry Points With Data Points
    % ======================================================

    function SpatialCorrelation = SpatialCorrelationOfInquiryPointsWithDataPoints( ...
        Config, ...
        X_InquiryPoints, ...
        X_DataPoints, ...
        SpatialDecorrelationScale, ...
        SpatialDecorrelationCutoff)

        % Format of points X: (n x 2) matrix, where n is the number of points
        % First column is x coordinates in meters and second column is y coordinates in meters.
        %
        % M inquiry points. X_InquiryPoints is M*2 matrix
        % N data points. X_DataPoints is N*2 matrix
        %
        % X_InquiryPoints = [Point_1_x, Point_1_y
        %                    Point_2_x, Point_2_y
        %                    ...
        %                    Point_M_x, Point_M_y];
        %
        % X_DataPoints =    [Point_1_x, Point_1_y
        %                    Point_2_x, Point_2_y
        %                    ...
        %                    Point_N_x, Point_N_y];
        %
        % SpatialCorrelation is N*M matrix.

        % Number of data points
        N = size(X_DataPoints,1);

        % Number of inquiry points (here, it is the same as the number of data points, N)
        M = size(X_InquiryPoints,1);

        % Difference along x axis (N*M matrix)
        Difference_x = repmat(X_DataPoints(:,1),1,M) - repmat(X_InquiryPoints(:,1)',N,1);

        % Difference along y axis (N*M matrix)
        Difference_y = repmat(X_DataPoints(:,2),1,M) - repmat(X_InquiryPoints(:,2)',N,1);

        % Spatial distance (N*M matrix)
        SpatialDistance = sqrt(Difference_x.^2+Difference_y.^2);

        % Normalize to decorrelation scale
        ScaledSpatialDistance = SpatialDistance / SpatialDecorrelationScale;

        % Choose correlation kernel type
        if strcmp(Config.Kriging.SpatialCorrelationKernelType,'gaussian')

            % Gaussian
            SpatialCorrelation = exp(-0.5*ScaledSpatialDistance.^2);

        elseif strcmp(Config.Kriging.SpatialCorrelationKernelType,'matern-1/2')

            % Matern class, nu = 1/2
            SpatialCorrelation = exp(-ScaledSpatialDistance);

        elseif strcmp(Config.Kriging.SpatialCorrelationKernelType,'matern-3/2')

            % Matern class, nu = 3/2
            SpatialCorrelation = (1 + sqrt(3)*ScaledSpatialDistance) .* exp(-sqrt(3)*ScaledSpatialDistance);

        elseif strcmp(Config.Kriging.SpatialCorrelationKernelType,'matern-5/3')

            % Matern class, nu = 5/3
            SpatialCorrelation = (1 + sqrt(5)*ScaledSpatialDistance + (5/3)*ScaledSpatialDistance) .* exp(-sqrt(5)*ScaledSpatialDistance);

        else
            error('Valid kernels are "gaussian", "matern-1/2", "matern-3/2", and "matern-5/3".');
        end
 
        % Zero-out distance larger than cutoff
        SpatialCorrelation(SpatialDistance > SpatialDecorrelationCutoff) = 0.0;
        % SpatialCorrelation(SpatialCorrelation < 0.3) = 0.3;

    end

    % =====================================
    % Covariance Of Data Points With Itself
    % =====================================

    function DataCovariance = CovarianceOfDataPointsWithItself( ...
        Config, ...
        X_DataPoints,T_DataPoints, ...
        TemporalDecorrelationScale,TemporalDecorrelationCutoff, ...
        SpatialDecorrelationScale,SpatialDecorrelationCutoff)

        % Spatio-temporal covariance between all two mutual points in a list of points.
        % The spatio-temporal covariance is the product of temporal correlation and spatial correlation.
        % Output is NxN matrix.

        % Number of points
        N = size(X_DataPoints,1);
        DataCovariance = zeros(N,N);

        for i = 1:N

            % Temporal correlation
            PointTemporalCorrelation = CovarianceFunction.TemporalCorrelation( ...
                T_DataPoints(i),T_DataPoints(i:end),TemporalDecorrelationScale,TemporalDecorrelationCutoff);

            % Spatial correlation
            PointSpatialCorrelation = CovarianceFunction.SpatialCorrelationOfInquiryPointsWithDataPoints( ...
                Config, ...
                X_DataPoints(i,:),X_DataPoints(i:end,:), ...
                SpatialDecorrelationScale,SpatialDecorrelationCutoff);

            % Spatio-temporal covariance
            PointSpatioTemporalCovariance = PointTemporalCorrelation .* PointSpatialCorrelation;
            % PointSpatioTemporalCovariance = PointTemporalCorrelation;
            % PointSpatioTemporalCovariance = PointSpatialCorrelation;

            % Store to covariance matrix
            DataCovariance(i,i:end) = PointSpatioTemporalCovariance';
            DataCovariance(i:end,i) = PointSpatioTemporalCovariance;

        end

    end

    % =============================================
    % Covariance Of Inquiry Points With Data Points
    % =============================================

    function InquiryPointsCovariance = CovarianceOfInquiryPointsWithDataPoints( ...
        Config, ...
        X_InquiryPoints,T_InquiryPoint, ...
        X_DataPoints,T_DataPoints, ...
        TemporalDecorrelationScale,TemporalDecorrelationCutoff, ...
        SpatialDecorrelationScale,SpatialDecorrelationCutoff)

        % Spatio-temporal covariance between a single point (an inquiry point) and a list of points.
        % The spatio-temporal covariance is the product of temporal correlation and spatial correlation
        % Output is Nx1 shape matrix (a column vector)

        % Temporal correlation
        TemporalCorrelation = CovarianceFunction.TemporalCorrelation(T_InquiryPoint,T_DataPoints,TemporalDecorrelationScale,TemporalDecorrelationCutoff);

        % Spatial correlation
        SpatialCorrelation = CovarianceFunction.SpatialCorrelationOfInquiryPointsWithDataPoints( ...
            Config, ...
            X_InquiryPoints,X_DataPoints, ...
            SpatialDecorrelationScale,SpatialDecorrelationCutoff);

        % Spatio-temporal covariance
        InquiryPointsCovariance = TemporalCorrelation .* SpatialCorrelation;
        % InquiryPointsCovariance = TemporalCorrelation;

    end

    % ========================================
    % Covariance of Inquity Points With Itself
    % ========================================

    function Q = CovarianceOfInquiryPointsWithItself( ...
        Config, ...
        X_InquiryPoints, ...
        VelocityStandardDeviation, ...
        SpatialDecorrelationScale,SpatialDecorrelationCutoff)

        % The inquiry points are all having the same time. So there is no temporal correlation between the
        % inquiry points. 

        % Number of inquiry points
        M = size(X_InquiryPoints,1);
        Q = zeros(M,M);

        for i = 1:M

            % Spatial correlation
            PointSpatialCorrelation = CovarianceFunction.SpatialCorrelationOfInquiryPointsWithDataPoints( ...
                Config, ...
                X_InquiryPoints(i,:),X_InquiryPoints(i:end,:),SpatialDecorrelationScale,SpatialDecorrelationCutoff);

            % Covariance
            PointSpatialCovariance = VelocityStandardDeviation^2 * PointSpatialCorrelation;

            % Store to covariance matrix
            Q(i,i:end) = PointSpatialCovariance';
            Q(i:end,i) = PointSpatialCovariance;

        end

    end

    % ---------------------
    % End of Static methods

    end
end
