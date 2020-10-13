classdef BasisFunctions

    methods(Static)

    % ===============
    % Basis Functions
    % ===============

    function [F,BasisPolynomialOrder,ScaleLength] = GenerateBasisFunctions(Config,X_DataPoints,T_DataPoints)

        % Used for Kriging with Trend (KT). Basis functions here are translation invariant, hence KT is intrinsic.
        % Inrtrinsic basis fucntions (IRF-k) are from exponential family, i.e. either polynomials, sin, cos funcitons, or exp.
        % For the time variable, I use sin,cos functions to capture the periodic nature of velicity changes within the tidal cycles.
        % For the space variables, I use a linear model (model order k = 1).
        %
        % basis functions = b_0*sin(Omega*t) +  b_1*cos(Omega*t)a_0 + a_0*1 + a_1*x + a_2*y + a_3*x*y + a_4*x^2 + a_5*y^2 + ...
        %
        % Degree of spatial basis functions:
        %
        % deg = 0;    % basis functions are: 1
        % deg = 1;    % basis functions are: 1,x,y
        % deg = 2;    % basis functions are: 1,x,y,x*y,x^2,y^2
        % deg = 3;    % basis functions are: 1,x,y,x*y,x^2,y^2,x^2*y,x*y^2,x^3,y^3

        % Number of data points
        N = size(X_DataPoints,1);

        % Number of temporal basis functions
        l_t = 2;    % basis functions are sin(omega*t) and cos(omega*t)

        % Order of spatial polynomial model.
        deg = Config.Kriging.SpatialBasisFunctionPolynomialDegree;

        % Rename deg for output
        BasisPolynomialOrder = deg;

        % Number of spatial basis functions
        l_x = (deg+1)*(deg+2) / 2;

        % Total number of basis functions
        l = l_t + l_x;

        % Initializing output
        F = zeros(N,l);

        % Temporal period, one tidal cycle (12 hours and 25 minutes)
        TemporalPeriod = Config.Kriging.TemporalBasisFunctionPeriodicity;   % A tide cycle, 12 hours and 25 minutes in unit of seconds
        Omega = 2*pi/TemporalPeriod;

        % Spatial basis functions
        ColumnCounter = 0;
        ScaleLength = Config.Kriging.ScaleLengthInSpatialBasisFunctions;
        for i_x = 0:deg
            for i_y = 0:deg-i_x
                ColumnCounter = ColumnCounter + 1;

                % x^i_x * y^i_y
                F(:,ColumnCounter) = (ScaleLength*X_DataPoints(:,1)).^i_x .* (ScaleLength*X_DataPoints(:,2)).^i_y;
            end
        end

        % Temporal basis functions
        F(:,ColumnCounter+1) = sin(Omega * T_DataPoints(:));
        F(:,ColumnCounter+2) = cos(Omega * T_DataPoints(:));

    end

    % =====================================
    % Solve Coefficients Of Basis Functions
    % =====================================

    function beta = SolveCoefficientsOfBasisFunctions(...
        DataCovariance,F,Z)

        % F*beta = data
        % F is N*m matrix, N is number of data points, m is number of basis functions
        % beta is coefficients of basis functions, m*1 column vector
        % Z is velocity Data (either east or north), a N*1 column vector.
        % F*beta = Z
        %
        % Note:
        % beta is independent of variance. So, instead of DataCovariance, the correlation can also be used.

        % Cinv = inv(DataCovariance);          % This  might be computationally expensive
        % F_inv = inv(F'*Cinv*F)*F'*Cinv;

        % We find Y = C_inv*F (or Y' = F'*C_inv) where C is DataCovariance
        Y = DataCovariance \ F;
        F_inv = inv(F'*Y)*Y';         % This is inverse of F with respect to Mahalanobis norm
        beta = F_inv*Z;

    end

    % ===========================
    % Jacobian of Basis Functions
    % ===========================

    function J = JacobianOfBasisFunctions( ...
        BasisPolynomialOrder, ...
        X_InquiryPoints, ...
        ScaleLength, ...
        beta_U,beta_V)

        % J = [du/dx, du/dy
        %      dv/dx, dv/dy]

        % Dimension of state for each tracer
        Dimension = size(X_InquiryPoints,2);
        if Dimension ~= 2
            error('Dimension should be 2.')
        end

        % Number of tracers
        M = size(X_InquiryPoints,1);

        % Initialize Jacobian
        J_dudx = zeros(M,M);
        J_dudy = zeros(M,M);
        J_dvdx = zeros(M,M);
        J_dvdy = zeros(M,M);

        % Rename basis polynomial order
        deg = BasisPolynomialOrder;

        % Number of spatial basis functions
        l_x = (deg+1)*(deg+2) / 2;

        % Derivarive with respect to x
        i = 0;
        for i_x = 0:deg
            for i_y = 0:deg-i_x

                % Update diagonal entry counter
                i = i + 1;

                % du/dx: beta_U * (i_x)*x^(i_x-1) * y^(i_y)
                J_dudx = J_dudx +  diag(beta_U(i) * ((i_x)*(ScaleLength*X_InquiryPoints(:,1)).^(i_x-1)) .* ((ScaleLength*X_InquiryPoints(:,2)).^(i_y))) * ScaleLength;

                % dv/dx: beta_V * (i_x)*x^(i_x-1) * y^(i_y)
                J_dvdx = J_dvdx +  diag(beta_V(i) * ((i_x)*(ScaleLength*X_InquiryPoints(:,1)).^(i_x-1)) .* ((ScaleLength*X_InquiryPoints(:,2)).^(i_y))) * ScaleLength;

                % du/dy: beta_U * x^(i_x) * (i_y)*y^(i_y-1)
                J_dudy = J_dudy +  diag(beta_U(i) * ((ScaleLength*X_InquiryPoints(:,1)).^(i_x)) .* ((i_y) * (ScaleLength*X_InquiryPoints(:,2)).^(i_y-1))) * ScaleLength;

                % dv/dy: beta_V * x^(i_x) * (i_y)*y^(i_y-1)
                J_dvdy = J_dvdy +  diag(beta_V(i) * ((ScaleLength*X_InquiryPoints(:,1)).^(i_x)) .* ((i_y) * (ScaleLength*X_InquiryPoints(:,2)).^(i_y-1))) * ScaleLength;
            end
        end

        % Concatenate
        J = [J_dudx,J_dudy
             J_dvdx,J_dvdy];

    end

    % ---------------------
    % End of Static methods

    end
end
