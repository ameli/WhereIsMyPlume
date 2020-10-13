classdef ModelSelection

    methods(Static)

    % ================================
    % Compute Spatio-Temporal Distance
    % ================================

    function [SpatioTemporalDistanceMatrix,TimeDifferenceMatrix] = ComputeSpatioTemporalDistanceMatrix( ...
        X_DataPoints,T_DataPoints, ...
        U_DataPoints,V_DataPoints, ...
        X_InquiryPoints,T_InquiryPoints)

        % To find correlation of data points with itself, call this function with T_InquiryPoints with
        % the same array of T_DataPoints, also for X_InquiryPoints use X_DataPoints.

        % Smooth velocities U and V to get tidal trends in velocity data and filter out fluctuations due to change of location of ship
        dT_DataPoints = mean(diff(T_DataPoints));                         % difference of two consecutive time points in data
        DayToSecond = 24*3600;
        MeanWindow = datenum([0 0 0 4 0 0]) * DayToSecond;                % In format of yyyy mm dd HH MM SS, SETTING
        MeanWindowLength = floor(MeanWindow / dT_DataPoints + 0.5);
        U_DataPoints_Mean = movmean(U_DataPoints,MeanWindowLength);
        V_DataPoints_Mean = movmean(V_DataPoints,MeanWindowLength);

        % U and V at inquiry points (using the same U and V as data points, but on inquiry times)
        if all(T_DataPoints == T_InquiryPoints)
            % Use the same data
            U_InquiryPoints_Mean = U_DataPoints_Mean;
            V_InquiryPoints_Mean = V_DataPoints_Mean;
        else
            % Interpolate
            U_InquiryPoints_Mean = interp1(T_DataPoints,U_DataPoints_Mean,T_InquiryPoints);
            V_InquityPoints_Mean = interp1(T_DataPoints,V_DataPoints_Mean,T_InquiryPoints);
        end

        % Number of data points (observed points)
        N = size(X_DataPoints,1);

        % Number of inquiry points (unobserved points)
        M = size(X_InquiryPoints,1);

        % Concatenated form of data points (observed points)
        X_obs = repmat(X_DataPoints(:,1),1,M);         % Observed X (N*M matrix)
        Y_obs = repmat(X_DataPoints(:,2),1,M);         % Observed Y (N*M matrix)
        T_obs = repmat(T_DataPoints',1,M);             % Observed T (N*M matrix)
        U_obs = repmat(U_DataPoints_Mean',1,M);        % Observed U (N*M matrix)
        V_obs = repmat(V_DataPoints_Mean',1,M);        % Observed V (N*M matrix)

        % Concatenated form of inquiry points (unobserved points)
        X_inq = repmat(X_InquiryPoints(:,1)',N,1);     % Unobserved X (N*M matrix)
        Y_inq = repmat(X_InquiryPoints(:,2)',N,1);     % Unobserved Y (N*M matrix)
        T_inq = repmat(T_InquiryPoints,N,1);           % Unobserved T (N*M matrix)
        U_inq = repmat(U_InquiryPoints_Mean,N,1);      % Unobserved U (N*M matrix)
        V_inq = repmat(V_InquiryPoints_Mean,N,1);      % Unobserved V (N*M matrix)

        % Difference of x and y, and t
        dX = X_inq - X_obs;
        dY = Y_inq - Y_obs;
        dT = T_inq - T_obs;

        % Using average of velocities between data points and inquiry points
        U = (U_inq + U_obs) * 0.5;
        V = (V_inq + V_obs) * 0.5;

        % Spatio-temporal differences
        FrozenFieldFlag = 0;                           % boolean, to turn on and off the "frozen covariance" field, SETTING
        dXT = dX - U.*dT*FrozenFieldFlag;
        dYT = dY - V.*dT*FrozenFieldFlag;

        % Spatio-temporal distance between data points and inquiry points 
        SpatioTemporalDistanceMatrix = sqrt(dXT.^2 + dYT.^2);           % Matrix of N*M size
        TimeDifferenceMatrix = dT;

        % Plot
        PlotFlag = 0;
        if PlotFlag
            f = figure('visible','off');
            imagesc(SpatioTemporalDistanceMatrix)
            shading interp
            colorbar()
            colormap hsv
            % caxis([0,6e4])
            title('Spatio-temporal distance')
            xlabel('Unobserved (Inquiry) points')
            ylabel('Observed (Data) points')
            % axis equal
            axis square
            saveas(f,'./plots/DistanceMatrix')

        end

    end

    % ==========================
    % Compute Correlation Matrix
    % ==========================

    function K = ComputeCorrelationMatrix( ...
        DistanceMatrix,TimeDifferenceMatrix, ...
        Hyperparameters)

        % DistanceMatrix:
        %   N*M matrix of the spatial, or spatio-temporal distance between two group of N points and M points.
        %   The distance can be spatial, meaning Euclidean distance between two points (x1,y1) and (x2,y2), or
        %   spatio-temporal distance between (x1,y1,t1) and (x2,y2,t2). See SpatioTemporalDistance function on how
        %   this distance is defined.
        %                 
        % CorrelationScale:
        %   Also known as the rho paramter in Matern class. Normalization of distance. 
        %
        % CorrelationSmoothness:
        %   Also known as the parameter nu in Mattern covariance function. A covariance (here, correlation) of degree nu
        %   is a ceil(nu-1) differentible function. For example if nu = 1.5, the function is first order differentiable.
        %
        % K:
        %   Output of this function, N*M matrix. Each entry corresponds to the the correlation of emtries of DitanceMatrix.
        %   This function computes correlation, not covariance. By multuplying correlation with variance (sigma^2), the
        %   covariance can be achieved.
        %
        % Note:
        % The Matern correlation function with large CorrelationSmoothness parameter is numerically unstable, particularly
        % for small distances (where distance is zero), becase it encounters multiplication of zero by infinity.
        % On the other hand, for CorrelationSmoothness parameter larger than 25~30, the Matern class approaches the
        % square exponential correlaiton (Gaussian function) which happends exactly at CorrelationSmoothness = infinity.
        % So, to avoid the numerical instability, at CorrelationSmoothness larger than 30, we use Gaussian function instead.

        % Hyperparameters
        tau = Hyperparameters(1);    % Temporal correlation scale, in hours
        rho = Hyperparameters(2);    % Spatial correlation scale, in Km
        alpha = Hyperparameters(3);  % Temporal smoothness
        nu = Hyperparameters(4);     % Spatial smoothness
        beta = Hyperparameters(5);   % Separability
        delta = Hyperparameters(6);

        % Unit conversion
        HourToSecond = 3600;
        KmToMeter = 1000;
        tau = tau * HourToSecond;  % in hours
        rho = rho * KmToMeter;  % In meters

        t = abs(TimeDifferenceMatrix) ./ tau;
        x = DistanceMatrix ./ rho;

        % Hyperparameters
        d = 2;           % Dimension of the spatial space
        % alpha = 0.5;     % Temporal process smoothness
        % beta = 1;        % In [0,1]. Zero means spatio-temporal kernel is separable. 1 means they spatio-temporal kernel is highly non-separable.
        % delta = 0;       % non-negative
        % nu = 50;       % Spatial process smoothness

        % Gneiting correlation function
        psi = (1+t.^(2*alpha)).^beta;
        % psi = (1+t.^(2)./(2*beta)).^beta;
        s = x ./ sqrt(psi);

        MaxSpatialSmoothness = 30;
        if nu == 0.5
            K = (1./psi.^(delta+d/2)) .* exp(-s);
        elseif nu == 1.5
            K = (1./psi.^(delta+d/2)) .* (1.0 + sqrt(3.0)*s).*exp(-sqrt(3.0)*s);
        elseif nu == 2.5
            K = (1./psi.^(delta+d/2)) .* (1.0 + sqrt(5.0)*s + (5.0/3.0)*(s.^2)).*exp(-sqrt(5.0)*s);
        elseif nu < MaxSpatialSmoothness
            K = (1./psi.^(delta+d/2)) .* ((2^(1-nu))/gamma(nu)) .* (sqrt(2*nu)*s).^nu .* besselk(nu,sqrt(2*nu)*s);
            % K = (1./psi.^(delta+d/2)) .* exp(-s.^nu);
            % K = (1 + t) ./ ((1 + t).^2 + x.^2).^(1.5);  % Cressie-Huang

            % Avoid NaN when distance is zero (since BesselK goes to infinity)
            K(x == 0) = 1;

            % Check K for nan and inf
            if any(isnan(K))
                error('Correlation has nan values.')
            elseif any(isinf(K))
                error('Correlation has inf values.')
            end
        else
            % Use Gaussian correlation function when SmoothnessCorrelation is larger than 30 (we assume Smoothness is infinity)
            K = (1./psi.^(delta+d/2)) .* exp(-0.5*s.^2);
        end

        % Plot
        PlotFlag = 0;
        if PlotFlag
            figure()
            imagesc(K)
            shading interp
            colorbar()
            colormap hsv
            % caxis([0,6e4])
            title('K')
            xlabel('Unobserved (Inquiry) points')
            ylabel('Observed (Data) points')
            % axis equal
            axis square
        end

    end

    % ===================
    % Find Optimal Nugget
    % ===================

    function [Optimal_eta,Optimal_sigma,Optimal_sigma0] = FindOptimalNugget(K,X,z)

        % Trace Estimation Utilities
        TraceEstimationUtilities.UseEigenvaluesMethod = true;
        K_eigenvalues = abs(eig(K));
        TraceEstimationUtilities.EigenvaluesMethodUtilities.K_eigenvalues = K_eigenvalues;

        % Use eigenvalues method to estimate trace
        Interval_eta = [1e-3,1e3];
        Results = LikelihoodEstimation.FindZeroOfLogLikelihoodFirstDerivative(z,X,K,TraceEstimationUtilities,Interval_eta);

        % Output
        Optimal_eta = Results.eta;
        Optimal_sigma = Results.sigma;
        Optimal_sigma0 = Results.sigma0;

        % Plot first derivative of likelihood function
        PlotFlag = 0;
        if PlotFlag
            LikelihoodEstimation.PlotLogLikelihoodFirstDerivative(z,X,K,TraceEstimationUtilities,Optimal_eta)
        end

    end

    % ========================
    % Hyperparameters LogPrior
    % ========================

    function LogPrior = HyperparametersLogPrior(...
            LowerBoundHyperparameters, ...
            UpperBoundHyperparameters, ...
            Hyperparameters);

        Large = 1e10;

        NumberOfHyperparameters = length(Hyperparameters);

        LogPrior = 0;
        for i = 1:NumberOfHyperparameters
            if (Hyperparameters(i) < LowerBoundHyperparameters(i)) || (Hyperparameters(i) > UpperBoundHyperparameters(i))
                LogPrior = LogPrior + Large;
            end
        end

    end

    % ======================
    % Log Posterior Function
    % ======================

    function LogPosterior = LogPosteriorFunction( ...
        DistanceMatrix,TimeDifferenceMatrix,X,z, ...
        LowerBoundHyperparameters, ...
        UpperBoundHyperparameters, ...
        Hyperparameters)
        
        % beta can be computed either by covariance (K*Variance) or correlation (K). It is independent of variance.
        % Note: beta and F are scaled (see Config.ScaleLengthInSpatialBasisFunction = S)
        % the columns 3:end of X are scaled by scale and columns 3:end of beta are scaled by the inverse of scale.
        % However, the product X*beta remains inact. Thus, as long as beta is used together with X, everything is consistent.

        % LogPrior of hyperparameters
        LogPrior = ModelSelection.HyperparametersLogPrior(...
            LowerBoundHyperparameters, ...
            UpperBoundHyperparameters, ...
            Hyperparameters);

        % If prior is unbounded, terminate the evaluation of posterior
        % Large = 1e10;
        % if LogPrior == Large
        %     LogPosterior = Large;
        %     return
        % end

        % Correlation
        K = ModelSelection.ComputeCorrelationMatrix( ...
            DistanceMatrix,TimeDifferenceMatrix, ...
            Hyperparameters);

        % Find optimal nugget
        [eta,sigma,sigma0] = ModelSelection.FindOptimalNugget(K,X,z);
        sigma2 = sigma^2;
        sigma02 = sigma0^2;
        log_eta = log10(eta);

        % n: number of observed data points, m: number of basis functions
        [n,m] = size(X);

        if ~isinf(eta)

            % Correlation with noise term
            Kn = K + eta.*eye(size(K));

            % Compute variance
            % beta = Kriging.SolveCoefficientsOfBasisFunctions(Kn,X,z);
            % deviation = z - X*beta;                                               % n*1 column vector
            % Prior for Variance of z using Inverse Gamma distribution
            % b0 = 1;                                 % Just a guess
            % a0 = (sigma2 / b0) + 1;       % Mean of inverse gamma distribution is b0/(a0-1). Solving for "a0" yields.
            % sigma = (2*b0 + (deviation'*(K \ deviation))) / (n-m + 2*(a0+1));    % a scalar
            % sigma2 = (deviation'*(Kn \ deviation)) / (n-m);    % a scalar

            % Inverse of correlation (K)
            B = X' * (Kn \ X);

            % Posterior probability
            A0 = -0.5*(n-m)*(1+log(sigma2));
            Kn_eigenvalues = abs(eig(Kn));
            A1 = -0.5*sum(log(Kn_eigenvalues));
            A2 = -0.5*log(det(B));
            LogPosterior = A0 + A1 + A2;

        else

            % eta is infinity. Instead of Kn, use covariance matrix C  = sigma02 * I
            B = X'*X;
            A0 = -0.5*(n-m)*(1+log(sigma02));
            A1 = 0.0;
            A2 = -0.5*log(det(B));
            LogPosterior = A0 + A1 + A2;

        end

        % Add prior to posterior
        % LogPosterior = LogPosterior + LogPrior;

        % Hyperparameters
        tau = Hyperparameters(1);
        rho = Hyperparameters(2);
        alpha = Hyperparameters(3);
        nu = Hyperparameters(4);
        beta = Hyperparameters(5);
        delta = Hyperparameters(6);

        fprintf('tau: %0.4f,\t rho: %0.4f,\t alpha: %0.4f,\t nu:%0.4f,\t beta: %0.4f,\t delta: %0.4f,\t LogEta: %0.4f,\t sigma0: %0.4f,\t sigma: %0.4f,\t Posterior: %0.16f\n', ...
            tau,rho,alpha,nu,beta,delta,log_eta,sigma0,sigma,LogPosterior);

    end

    % ==============================
    % Maximize Posterior Probability
    % ==============================

    function CovarianceParameters = MaximizePosteriorProbability( ...
        DistanceMatrix,TimeDifferenceMatrix,F,Z)

        % Guess search parameters
        Guess_tau = 1;           % Temporal correlation scale, In hours
        Guess_rho = 1;           % Spatial correlation scale, In Km
        Guess_alpha = 0.5;       % Temporal smoothness, Between 0 and 1
        Guess_nu = 2;            % Spatial smoothness, Greater than zero
        Guess_beta = 0.5;        % Spatiotemporal separability, Between 0 and 1
        Guess_delta = 1;
        GuessParameters = [Guess_tau,Guess_rho,Guess_alpha,Guess_nu,Guess_beta,Guess_delta];
        % LowerBoundParameters = [0.01,0.1,0.0,0.5,0];
        LowerBoundParameters = [0.01,0.1,0,0.5,0,0];
        % UpperBoundParameters = [5,5,1,31,1];
        UpperBoundParameters = [100,100,1,40,1,100];

        % Functional to maximize
        LogPosteriorPartialFunction = @(Parameters) -ModelSelection.LogPosteriorFunction( ...
            DistanceMatrix,TimeDifferenceMatrix,F,Z, ...
            LowerBoundParameters,UpperBoundParameters, ...
            Parameters);

        % Optimize hyperparameters. choose between ooption 1 (local search) and option 2  (global search)

        % Option 1: local search methods
        Options = optimoptions(@fminunc,'UseParallel',true,'MaxIterations',5000); 
        % OptimalParameters = fminsearch(LogPosteriorPartialFunction,GuessParameters);
        % OptimalParameters = fminunc(LogPosteriorPartialFunction,GuessParameters);
        OptimalParameters = fmincon(LogPosteriorPartialFunction,GuessParameters,[],[],[],[],LowerBoundParameters,UpperBoundParameters,[],Options);  % GOOD

        % Option 2:
        % OptimalParameters = simulannealbnd(LogPosteriorPartialFunction,GuessParameters,LowerBoundParameters,UpperBoundParameters,Options);
        % OptimalParameters = patternsearch(LogPosteriorPartialFunction,GuessParameters,[],[],[],[],LowerBoundParameters,UpperBoundParameters,optimoptions('patternsearch','UseParallel',true));
        % OptimalParameters = surrogateopt(LogPosteriorPartialFunction,LowerBoundParameters,UpperBoundParameters,optimoptions('surrogateopt','Display','off','UseParallel',true,'MaxFunctionEvaluations',1000));
        % OptimalParameters = particleswarm(LogPosteriorPartialFunction,length(LowerBoundParameters),LowerBoundParameters,UpperBoundParameters,optimoptions('particleswarm','UseParallel',true,'SwarmSize',100,'HybridFcn',@fmincon));

        % GlobalSearch (this is not parallel. For parallel, use multistart)
        % gs = GlobalSearch('FunctionTolerance',1e-4,'NumTrialPoints',1000);
        % problem = createOptimProblem('fmincon','x0',GuessParameters,'objective',LogPosteriorPartialFunction,'lb',LowerBoundParameters,'ub',UpperBoundParameters);
        % OptimalParameters = run(gs,problem)

        % Multistart
        % Options = optimoptions(@fmincon,'Algorithm','sqp');
        % problem = createOptimProblem('fmincon','objective',LogPosteriorPartialFunction,'x0',GuessParameters,'lb',LowerBoundParameters,'ub',UpperBoundParameters,'options',Options);
        % ms = MultiStart('UseParallel',true');
        % StartPoints = 20;
        % OptimalParameters = run(ms,problem,StartPoints);

        disp(OptimalParameters)

        % Get the optimal value of the function
        OptimialValue = ModelSelection.LogPosteriorFunction( ...
            DistanceMatrix,TimeDifferenceMatrix,F,Z, ...
            LowerBoundParameters,UpperBoundParameters, ...
            OptimalParameters);

    end

    % ==========================
    % Plot Posterior Probability
    % ==========================

    function CovarianceParameters = PlotPosteriorProbability( ...
        DistanceMatrix,TimeDifferenceMatrix,F, ...
        Z,Z_Variance_Nugget)

        % DistanceMatrix: N*M (where M is here N) of distancxe between all mutual data points (x1,y1,t1) and (x2,y2,t2)
        % F: design matrix. N*l matrix of basis functions
        % Z: N*1 matrix. velocity data at observed location. This is either U_DataPoints_Mean or V_DataPoints_Mean
        % Z_Variance_Nugget: variance of Z (either U or V velocity). This variance is the nugget effect.

        % Scale
        Scales = logspace(2,4,20);
        % Scales = [2000];

        % Smoothness
        % CorrelationSmoothnesses = linspace(0.5,4,4);
        % CorrelationSmoothnesses = [0.5,1,5];
        CorrelationSmoothnesses = [1.5];

        Z_Variance_Nugget = 0.02^2;

        % Prior variance
        % Z_Variance_Prior = Z_Variance_Nugget * logspace(-2,2,10);
        % Z_Variance_Prior = logspace(-5,-1,10);
        NuggetLogRatios = linspace(0.5,1.5,20);

        % LogPosterior = zeros(length(Scales),length(CorrelationSmoothnesses),length(Z_Variance_Prior));
        LogPosterior = zeros(length(Scales),length(CorrelationSmoothnesses),length(NuggetLogRatios));
        Nuggets = zeros(size(LogPosterior));

        % Loop through parameters
        for i = 1:length(Scales)
            for j = 1:length(CorrelationSmoothnesses)
                % for k = 1:length(Z_Variance_Prior)
                for k = 1:length(NuggetLogRatios)

                    % Covariance Parameters (these paramters have to be optimized)
                    % CovarianceParameters.CorrelationScale = Scales(i);
                    % CovarianceParameters.CorrelationSmoothness = CorrelationSmoothnesses(j);
                    % CovarianceParameters.Variance = Z_Variance_Prior(k);

                    % Parameters to optimize as a vector argument
                    CorrelationScale = Scales(i);
                    CorrelationSmoothness = CorrelationSmoothnesses(j);
                    NuggetLogRatio = NuggetLogRatios(k);
                    CovarianceParameters = [CorrelationScale,CorrelationSmoothness,NuggetLogRatio];

                    % Posterior
                    % LogPosterior(i,j,k) = ModelSelection.LogPosteriorFunction(DistanceMatrix,TimeDifferenceMatrix,F,Z,Z_Variance_Nugget,CovarianceParameters);
                    [LogPosterior(i,j,k),Nuggets(i,j,k)] = ModelSelection.LogPosteriorFunction(DistanceMatrix,TimeDifferenceMatrix,F,Z,CovarianceParameters);
                    % fprintf('i: %d, j: %d, \t Scale: %0.2f,\t Smoothness: %0.2f,\t Posterior: %0.2f\n',i,j,Scales(i),CorrelationSmoothnesses(j),LogPosterior(i,j))
                end
            end
        end

        Colors = jet(size(LogPosterior,1));
        for i = 1:size(LogPosterior,1)
            for j = 1:size(LogPosterior,2)
                % semilogx(Scales,LogPosterior(:,j,1),'-o','Displayname',sprintf('nu = %0.2f',CorrelationSmoothnesses(j)))
                % semilogx(Z_Variance_Prior,squeeze(LogPosterior(i,j,:)),'-o','Displayname',sprintf('rho = %0.2f',Scales(i)),'color',Colors(i))
                plot(NuggetLogRatios,squeeze(LogPosterior(i,j,:)),'-o','color',Colors(i,:),'Displayname',sprintf('nu = %0.2f',Scales(i)))
                hold on
            % plot(CorrelationSmoothnesses,LogPosterior(1,:),'-o')
            end
        end
        legend()
        title('Log Posterior')


        figure()
        for i = 1:size(LogPosterior,1)
            for j = 1:size(LogPosterior,2)
                % semilogx(Scales,LogPosterior(:,j,1),'-o','Displayname',sprintf('nu = %0.2f',CorrelationSmoothnesses(j)))
                % semilogx(Z_Variance_Prior,squeeze(LogPosterior(i,j,:)),'-o','Displayname',sprintf('rho = %0.2f',Scales(i)),'color',Colors(i))
                plot(NuggetLogRatios,sqrt(squeeze(Nuggets(i,j,:))),'-o','color',Colors(i,:),'Displayname',sprintf('nu=%0.2f',Scales(i)))
                hold on
            % plot(CorrelationSmoothnesses,LogPosterior(1,:),'-o')
            end
        end
        legend()
        title('Nuggets')

        % Pack output
        % CovarianceParameters.CorrelationSmoothness = CorrelationSmoothness;
        % CovarianceParameters.CorrelationScale = CorrelationScale;
        % CovarianceParameters.Variance = Variance;

    end

    % =====================================
    % Covariance Of Data Points With Itself
    % =====================================
    
    function [U_DataCovariance,V_DataCovariance,U_CovarianceParameters,V_CovarianceParameters] = CovarianceOfDataPointsWithItself( ...
        Config, ...
        F, ...
        X_DataPoints,T_DataPoints, ...
        U_DataPoints,V_DataPoints, ...
        U_CovarianceParameters,V_CovarianceParameters)

        % Uses maximum a posterior method to find best correlation parameters, then constructs the correlation
        % and covariance matrix.

        % Spatio-temporal distance matrix
        [SpatioTemporalDistanceMatrix,TimeDifferenceMatrix] = ModelSelection.ComputeSpatioTemporalDistanceMatrix( ...
            X_DataPoints,T_DataPoints, ...
            U_DataPoints,V_DataPoints, ...
            X_DataPoints,T_DataPoints);

        % Find optimal covariance parameters
        if ~exist('U_CovarianceParameters','var')
            U_CovarianceParameters = ModelSelection.MaximizePosteriorProbability( ...
                SpatioTemporalDistanceMatrix,TimeDifferenceMatrix,F, ...
                U_DataPoints');
        end

        if ~exist('V_CovarianceParameters','var')
            V_CovarianceParameters = ModelSelection.MaximizePosteriorProbability( ...
                SpatioTemporalDistanceMatrix,TimeDifferenceMatrix,F, ...
                V_DataPoints');
        end

        % Build correlation with the new parameters for both data of U and V fields
        U_CorrelationMatrix = ModelSelection.ComputeCorrelationMatrix(SpatioTemporalDistanceMatrix,U_CovarianceParameters);
        U_DataCovariance = U_CorrelationMatrix * U_CovarianceParameters.Variance; 

        V_CorrelationMatrix = ModelSelection.ComputeCorrelationMatrix(SpatioTemporalDistanceMatrix,V_CovarianceParameters);
        V_DataCovariance = V_CorrelationMatrix * V_CovarianceParameters.Variance; 

    end
        
    % ---------------------
    % End of Static methods

    end

end
