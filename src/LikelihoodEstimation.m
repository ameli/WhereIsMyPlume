classdef LikelihoodEstimation

    methods(Static)

    % ================================
    % Log Likelihood With Sigma Sigma0
    % ================================

    function lp = LogLikelihoodWithSigmaSigma0(z,X,K,TraceEstimationUtilities,SignSwitch,sigma,sigma0)
        
        % Here we use direct parameter, sigma and sigma0
        % 
        % SignSwitch chnages the sign of the output from lp to -lp. When True, this is used to
        % minimizing (instad of maximizing) the negative of log-likelihood function.

        % Covariance
        n = size(K,1); 
        I = eye(n);
        S = (sigma^2)*K + (sigma0^2)*I;

        % Compute log det (S)
        if TraceEstimationUtilities.UseEigenvaluesMethod == true
            % Use eigenvalues of K to estimate trace
            K_eigenvalues = TraceEstimationUtilities.EigenvaluesMethodUtilities.K_eigenvalues;
            LogDet_S = sum(log((sigma^2)*K_eigenvalues + (sigma0^2)));
        else
            % Use Cholesky factorization
            L = chol(S,'lower');
            Diag_L = diag(L);
            LogDet_L = real(sum(log(Diag_L)));
            LogDet_S = 2.0*LogDet_L;
        end

        % Compute log det (X.T*Sinv*X)
        Y = S \ X;
        w = S \ z;

        XtSinvX = X'*Y;
        LogDet_XtSinvX = log(det(XtSinvX));

        Binv = inv(XtSinvX);
        YBinvYt = Y*Binv*Y';

        % Log likelihood
        lp = -0.5*LogDet_S -0.5*LogDet_XtSinvX -0.5*dot(z,w-YBinvYt*z);

        % If lp is used in scipy.optimize.minimize, change the sign to optain the minimum of -lp
        if SignSwitch == true
            lp = -lp;
        end

    end

    % ===============================
    % Log Likelihood First Derivative
    % ===============================

    function dlp_deta =  LogLikelihoodFirstDerivative(z,X,K,TraceEstimationUtilities,LogEta)

        % lp is the log likelihood probability
        % dlp_deta is d(lp)/d(eta), is the derivative of lp with respect to eta when the optimal
        % value of sigma is subtituted in the likelihood function per given eta.

        % Change LogEta to eta
        if isinf(LogEta) && (LogEta < 0)
            eta = 0.0;
        else
            eta = 10.0^LogEta;
        end

        % Correlation
        I = eye(size(K,1));
        Kn = K + eta*I;
   
        % Compute Kn_inv*X and Kn_inv*z
        Y = Kn \ X;
        w = Kn \ z;

        [n,m] = size(X);

        % Splitting M into M1 and M2. Here, we compute M2
        B = X'*Y;
        Binv = inv(B);
        Ytz = Y'*z;
        Binv_Ytz = Binv*Ytz;
        Y_Binv_Ytz = Y*Binv_Ytz;
        Mz = w - Y_Binv_Ytz;

        % Traces
        if TraceEstimationUtilities.UseEigenvaluesMethod == true
            K_eigenvalues = TraceEstimationUtilities.EigenvaluesMethodUtilities.K_eigenvalues;
            TraceKninv = sum(1.0./(K_eigenvalues + eta));
        else
            TraceKninv = TraceEstimation.EstimateTrace(TraceEstimationUtilities,eta);
            % TraceKninv = TraceEstimation.ComputeTraceOfInverse(Kn)    % Use direct method without interpolation, Test
        end

        YtY = Y'*Y;
        TraceBinvYtY = trace(Binv*YtY);
        TraceM = TraceKninv - TraceBinvYtY;

        % Derivative of log likelihood
        zMz = dot(z,Mz);
        zM2z = dot(Mz,Mz);
        sigma02 = zMz/(n-m);
        dlp_deta = -0.5*(TraceM - zM2z/sigma02);

    end

    % ================================
    % Log Likelihood Second Derivative
    % ================================

    function d2lp_deta2 = LogLikelihoodSecondDerivative(z,X,K,TraceEstimationUtilities,eta)

        % The second derivative of lp is computed as a function of only eta. Here, we
        % substituted optimal value of sigma, which istself is a function of eta.

        % Correlation
        I = eye(size(K,1));
        Kn = K + eta*I;

        Y = Kn \ X;
        V = Kn \ Y;
        w = Kn \ z;

        [n,m] = size(X);

        % Splitting M
        B = X'*Y;
        Binv = inv(B);
        Ytz = Y'*z;
        Binv_Ytz = Binv*Ytz;
        Y_Binv_Ytz = Y*Binv_Ytz;
        Mz = w - Y_Binv_Ytz;

        % Trace of M
        if TraceEstimationUtilities.UseEigenvaluesMethod == true
            % Use eigenvalues method
            K_eigenvalues = TraceEstimationUtilities.EigenvaluesMethodUtilities.K_eigenvalues;
            TraceKninv = sum(1.0/(K_eigenvalues + eta));
        else
            % Use interpolation method
            TraceKninv = TraceEstimation.EstimateTrace(TraceEstimationUtilities,eta);
        end

        YtY = Y'*Y;
        A = Binv*YtY;
        TraceA = trace(A);
        TraceM = TraceKninv - TraceA;

        % Trace of M^2
        if TraceEstimationUtilities.UseEigenvaluesMethod == true
            K_eigenvalues = TraceEstimationUtilities.EigenvaluesMethodUtilities.K_eigenvalues;
            TraceKn2inv = sum(1.0./((K_eigenvalues + eta).^2));
        else
            Kn2 = Kn*Kn;
            TraceKn2inv = TraceEstimation.ComputeTraceOfInverse(Kn2);
        end

        YtV = Y'*V;
        C = Binv*YtV;
        TraceC = trace(C);
        AA = A*A;
        TraceAA = trace(AA);
        TraceM2 = TraceKn2inv - 2.0*TraceC + TraceAA;

        % Find z.T * M^3 * z
        YtMz = Y'*Mz;
        Binv_YtMz = Binv*YtMz;
        Y_Binv_YtMz = Y*Binv_YtMz;

        v = Kn \ Mz;
        MMz = v - Y_Binv_YtMz;

        % Second derivative (only at the location of zero first derivative)
        zMz = dot(z,Mz);
        zM2z = dot(Mz,Mz);
        zM3z = dot(Mz,MMz);
        sigma02 = zMz / (n-m);
        d2lp_deta2 = (0.5/sigma02)*((TraceM2/(n-m) + (TraceM/(n-m))^2) * zMz - 2.0*zM3z);

    end

    % ======================================
    % Find Zero of Log Likelihood Derivative
    % ======================================

    function Results = FindZeroOfLogLikelihoodFirstDerivative(z,X,K,TraceEstimationUtilities,Interval_eta)

        % root finding of the derivative of lp.
        % The log likelihood function is implicitly a function of eta. We have substituted the
        % value of optimal sigma, which itself is a function of eta.

        % ------------------
        % Find Optimal Sigma
        % ------------------

        function sigma = FindOptimalSigma(z,X,K,eta)
        
            % Based on a given eta, finds optimal sigma

            I = eye(size(K,2));
            Kn = K +  eta*I;

            Y = Kn \ X;
            w = Kn \ z;

            [n,m] = size(X);
            B = X'*Y;
            Binv = inv(B);
            Ytz = Y'*z;
            v = Y*Binv*Ytz;
            sigma2 = dot(z,w-v) / (n-m);
            sigma = sqrt(sigma2);

        end

        % -------------------
        % Find Optimal Sigma0
        % -------------------

        function sigma0 = FindOptimalSigma0(z,X)

            % When eta is very large, we assume sigma is zero. Thus, sigma0 is computed by this function.

            [n,m] = size(X);
            B = X'*X;
            Binv = inv(B);
            Xtz = X'*z;
            v = X*Binv*Xtz;
            sigma02 = dot(z,z-v) / (n-m);
            sigma0 = sqrt(sigma02);

        end

        % -----------------

        % Find an interval that the function changes sign before finding its root (known as bracketing the function)
        LogEta_Start = log10(Interval_eta(1));
        LogEta_End = log10(Interval_eta(2));

        % Partial function with minus to make maximization to a minimization
        LogLikelihoodFirstDerivative_PartialFunction = @(LogEta) ...
                LikelihoodEstimation.LogLikelihoodFirstDerivative( ...
                z,X,K,TraceEstimationUtilities,LogEta);

        Tolerance = 1e-6;     % SETTING
        Options = optimset('TolX',Tolerance);
        NumTrials = 4;    % SETTING
        FirstLogEtaGuess = 0;
        LogEtaGuess = FirstLogEtaGuess;
        GuessIncrement = 1;

        for i = [1:NumTrials,-1:-1:-NumTrials]

            % Attempt to find a root
            [Root,FuncValue,ExitFlag,OutputStruct] = fzero(LogLikelihoodFirstDerivative_PartialFunction,LogEtaGuess,Options);

            % Check success of finding the root
            if ExitFlag == 1

                % Check second derivative
                d2lp_deta2 = LikelihoodEstimation.LogLikelihoodSecondDerivative(z,X,K,TraceEstimationUtilities,Root);
                if d2lp_deta2 <= 0
                    break
                end

            else
                % No root found. Change initial guess
                LogEtaGuess = FirstLogEtaGuess + i*GuessIncrement;

            end
        end

        % Initial points
        Bracket = [LogEta_Start,LogEta_End];
        % BracketFound,Bracket,BracketValues = RootFinding.FindIntervalWithSignChange(LogLikelihoodFirstDerivative_PartialFunction,Bracket,NumTrials,args=(),)

        % if BracketFound
        if (ExitFlag == 1)
            eta = 10^Root;
            sigma = FindOptimalSigma(z,X,K,eta);
            sigma0 = sqrt(eta) * sigma;

        elseif ExitFlag == -5

            % Function is singular. Assume solution is either in eta = 0 or eta = inf
            % Case 1: assume eta = 0
            Case1_eta = 0;
            Case1_sigma0 = 0;
            Case1_sigma = FindOptimalSigma(z,X,K,Case1_eta);
            SignSwitch = false;
            Case1_lp = LikelihoodEstimation.LogLikelihoodWithSigmaSigma0(z,X,K,TraceEstimationUtilities,SignSwitch,Case1_sigma,Case1_sigma0);
            
            % Case 2: assume eta = inf
            Case2_sigma = 0;
            Case2_sigma0 = FindOptimalSigma0(z,X);
            SignSwitch = false;
            Case2_lp = LikelihoodEstimation.LogLikelihoodWithSigmaSigma0(z,X,K,TraceEstimationUtilities,SignSwitch,Case2_sigma,Case2_sigma0);

            if Case1_lp > Case2_lp
                sigma = Case1_sigma;
                sigma0 = Case1_sigma0;
                eta = 0;
            else
                sigma = Case2_sigma;
                sigma0 = Case2_sigma0;
                eta = inf;
            end

        else
            % Bracket with sign change was not found.

            % Evaluate the function in intervals
            log_eta_left = Bracket(1);
            log_eta_right = Bracket(2);
            dlp_deta_left = LogLikelihoodFirstDerivative_PartialFunction(log_eta_left);
            dlp_deta_right = LogLikelihoodFirstDerivative_PartialFunction(log_eta_right);
            dlp_deta_root = LogLikelihoodFirstDerivative_PartialFunction(Root);

            % Second derivative of log likelihood at eta = zero, using either of the two methods below
            eta_zero = 0.0;
            d2lp_deta2_zero_eta = LikelihoodEstimation.LogLikelihoodSecondDerivative(z,X,K,TraceEstimationUtilities,eta_zero);   % Method 1 directly from analytical equation
            % dlp_deta_zero_eta = LikelihoodEstimation.LogLikelihoodFirstDerivative(z,X,K,TraceEstimationUtilities,log10(eta_zero))
            % d2lp_deta2_zero_eta = (dlp_deta_lowest_eta - dlp_deta_zero_eta) / eta_lowest        % Method 2 usng forward differencing from first derivative

            % fprintf('dL/deta   at eta = 0.0\t %0.2f\n',dlp_deta_zero_eta)
            fprintf('dL/deta   at log_eta = %0.2f\t %0.2f\n',log_eta_left,dlp_deta_left)
            fprintf('dL/deta   at log_eta = %0.2f\t %0.16f\n',log_eta_right,dlp_deta_right)
            fprintf('d2L/deta2 at eta = 0.0\t %0.2f\n',d2lp_deta2_zero_eta)

            % No sign change. Can not find a root
            if ((dlp_deta_left <= 0) && (dlp_deta_right >= 0)) || ((dlp_deta_left >= 0) && (dlp_deta_right <= 0))
                error('Sign change detected, but fzero did not find the root.')

            elseif (dlp_deta_left > 0) && (dlp_deta_right > 0)
                if d2lp_deta2_zero_eta > 0
                    eta = 0.0;
                else
                    eta = inf;
                end

            elseif (dlp_deta_left < 0) && (dlp_deta_right < 0)
                if d2lp_deta2_zero_eta < 0
                    eta = 0.0;
                else
                    eta = inf;
                end
            end

            % Find sigma and sigma0
            if eta == 0
                sigma0 = 0;
                sigma = FindOptimalSigma(z,X,K,eta);
            elseif isinf(eta)
                sigma = 0;
                sigma0 = FindOptimalSigma0(z,X);
            else
                error('eta must be zero or inf at this point.')
            end
        end

        % Output struct
        Results.sigma = sigma;
        Results.sigma0 = sigma0;
        Results.eta = eta;
        
    end

    % ==================================
    % Compute Bounds of First Derivative
    % ==================================

    function [dlp_deta_UpperBound,dlp_deta_LowerBound] = ComputeBoundsOfFirstDerivative(X,K,eta)
        
        % Upper and lower bound.

        [n,m] = size(X);
        Eigenvalue_smallest = abs(eigs(K,1,'smallestabs'));
        Eigenvalue_largest = abs(eigs(K,1,'largestabs'));
        dlp_deta_UpperBound = 0.5*(n-m)*(1./(eta+Eigenvalue_smallest) - 1./(eta+Eigenvalue_largest));
        dlp_deta_LowerBound = -dlp_deta_UpperBound;

    end

    % =====================================
    % Compute Asymptote of First Derivative
    % =====================================

    function [Asymptote_1_order,Asymptote_2_order,Roots_1,Roots_2] = ComputeAsymptoteOfFirstDerivative(z,X,K,eta)

        % Computes first and second order asymptote to the first derivative of log marginal likelihood function.

        % Initialize output
        Asymptote_1_order = zeros(length(eta),1);
        Asymptote_2_order = zeros(length(eta),1);

        [n,m] = size(X);
        I = eye(n);
        Im = eye(m);
        Q = X*inv(X'*X)*X';
        R = I - Q;
        N  = K*R;
        N2 = N*N;
        N3 = N2*N;
        N4 = N3*N;

        mtrN = trace(N)/(n-m);
        mtrN2 = trace(N2)/(n-m);
        
        A0 = -R*(mtrN*I - N);
        A1 =  R*(mtrN*N + mtrN2*I - 2*N2);
        A2 = -R*(mtrN*N2 + mtrN2*N - 2*N3);
        A3 =  R*(mtrN2*N2 - N4);
        
        zRz = dot(z,R*z);
        z_Rnorm = sqrt(zRz);
        zc = z / z_Rnorm;
        
        a0 = dot(zc,A0*zc);
        a1 = dot(zc,A1*zc);
        a2 = dot(zc,A2*zc);
        a3 = dot(zc,A3*zc);
        
        for i = 1:length(eta)

            Asymptote_1_order(i) = (-0.5*(n-m))*(a0 + a1/eta(i))/eta(i)^2;
            Asymptote_2_order(i) = (-0.5*(n-m))*(a0 + a1/eta(i) + a2/eta(i)^2 + a3/eta(i)^3)/eta(i)^2;

        end

        % Roots
        Polynomial_1 = [a0,a1];
        Polynomial_2 = [a0,a1,a2,a3];

        Roots_1 = roots(Polynomial_1);
        Roots_2 = roots(Polynomial_2);

        % Remove complex roots
        Roots_2 = sort(real(Roots_2(abs(imag(Roots_2)) < 1e-10)))

        disp('Asymptote roots:')
        disp(Roots_1)
        disp(Roots_2)

    end

    % ====================================
    % Plot Log Likelihood First Derivative
    % ====================================

    function PlotLogLikelihoodFirstDerivative(z,X,K,TraceEstimationUtilities,Optimal_eta)

        % Plots the derivative of log likelihood as a function of eta.
        % Also it shows where the optimal eta is, which is the location
        % where the derivative is zero.

        disp('Plot first derivative ...')

        if (~isnan(Optimal_eta)) && (Optimal_eta ~= 0) && (~isinf(Optimal_eta))
            PlotOptimal_eta = true;
        else
            PlotOptimal_eta = false;
        end

        % Specify which portion of eta array be high resolution for plotting in the inset axes
        LogEtaStart = -3;
        LogEtaEnd = 3;

        if PlotOptimal_eta
            LogEtaStartHighRes = floor(log10(Optimal_eta));
            LogEtaEndHighRes = LogEtaStartHighRes + 3;

            % Arrays of low and high resolutions of eta
            eta_HighRes = logspace(LogEtaStartHighRes,LogEtaEndHighRes,100);
            eta_LowRes_left = logspace(LogEtaStart,LogEtaStartHighRes,50);
            eta_LowRes_right = logspace(LogEtaEndHighRes,LogEtaEnd,20);

            % array of eta as a mix of low and high res
            if LogEtaEndHighRes >= LogEtaEnd
                eta = [eta_LowRes_left,eta_HighRes];
            else
                eta = [eta_LowRes_left,eta_HighRes,eta_LowRes_right];
            end

        else
            eta = logspace(LogEtaStart,LogEtaEnd,100);
        end

        % Compute derivative of L
        dlp_deta = zeros(length(eta),1);
        for i  = 1:length(eta)
            dlp_deta(i) = LikelihoodEstimation.LogLikelihoodFirstDerivative(z,X,K,TraceEstimationUtilities,log10(eta(i)));
        end

        % Compute upper and lower bound of derivative
        [dlp_deta_UpperBound,dlp_deta_LowerBound] = LikelihoodEstimation.ComputeBoundsOfFirstDerivative(X,K,eta);

        % Compute asymptote of first derivative, using both first and second order approximation
        try
            % eta_HighRes migh not be defined, depending on PlotOptimal_eta
            x = eta_HighRes;
        catch
            x = logspace(1,LogEtaEnd,100);
        end
        [dlp_deta_Asymptote_1,dlp_deta_Asymptote_2,Roots_1,Roots_2] = LikelihoodEstimation.ComputeAsymptoteOfFirstDerivative(z,X,K,x);

        % Main plot
        figure()
        semilogx(eta,dlp_deta_UpperBound,'--','color','black','DisplayName','Upper bound')
        hold on
        semilogx(eta,dlp_deta_LowerBound,'-.','color','black','DisplayName','Lower bound')
        semilogx(eta,dlp_deta,'color','black','DisplayName','Exact')
        if PlotOptimal_eta
            semilogx(Optimal_eta,0,'o','markersize',4,'color','black','DisplayName','Root')
        end

        % Min of plot limit
        MaxPlot = max(dlp_deta);
        MaxPlotLim = ceil(abs(MaxPlot)/10.0)*10.0*sign(MaxPlot);

        MinPlotLim1 = -100;
        yticks([MinPlotLim1,0,MaxPlotLim])
        ylim([MinPlotLim1,MaxPlotLim])
        xlim([eta(1),eta(end)])
        xlabel('$\eta$','Interpreter','latex')
        ylabel('$\mathrm{d} \ell_{\hat{\sigma}^2(\eta)}(\eta)/\mathrm{d} \eta$','Interpreter','latex')
        title('Derivative of Log Marginal Likelihood Function')
        grid()
        legend('Location','southwest')

        % Inset plot
        if PlotOptimal_eta

            axes('Position',[0.43,0.39,0.45,0.5])
            semilogx(eta,abs(dlp_deta_UpperBound),'--','color','black','DisplayName','Upper bound')
            hold on
            semilogx(eta,abs(dlp_deta_LowerBound),'-.','color','black','DisplayName','Lower bound')
            semilogx(x,dlp_deta_Asymptote_1,'DisplayName','1st order asymptote','color','red')
            semilogx(x,dlp_deta_Asymptote_2,'DisplayName','2nd order asymptote','color','green')
            semilogx(eta_HighRes,dlp_deta(length(eta_LowRes_left)+1:length(eta_LowRes_left)+length(eta_HighRes)),'color','black','DisplayName','Exact')
            semilogx(Optimal_eta,0,'o','markersize',6,'color','black','DisplayName',sprintf('Exact root at eta = 10^{%0.2f}',log10(Optimal_eta)))
            semilogx(Roots_1(end),0,'o','markersize',6,'color','red','DisplayName',sprintf('Approximated root at eta_1 = 10^{%0.2f}',log10(Roots_1(end))))
            semilogx(Roots_2(end),0,'o','markersize',6,'color','green','DisplayName',sprintf('Approximated root at eta_2 = 10^{%0.2f}',log10(Roots_2(end))))
            % set_xlim([eta_HighRes[0],eta_HighRes[-1]])

            % Find suitable range for plot limits
            MinPlot = abs(min(dlp_deta))
            MinPlotBase = 10^floor(log10(abs(MinPlot)))
            % MinPlotLim = numpy.ceil(MinPlot/MinPlotBase)*MinPlotBase
            MinPlotLim = ceil(MinPlot/MinPlotBase + 1.0)*MinPlotBase
            ylim([-MinPlotLim,MinPlotLim])
            set(gca,'YGrid','on')
            yticks([-abs(MinPlotLim),0,abs(MinPlotLim)])

            text(Optimal_eta*10^0.05,MinPlotLim*0.05,'\eta','horizontalalignment','left','verticalalignment','bottom','fontsize',10)
            text(Roots_1(end)*10^0.05,MinPlotLim*0.05,'\eta_1','horizontalalignment','left','verticalalignment','bottom','fontsize',10)
            text(Roots_2(end)*10^0.05,MinPlotLim*0.05,'\eta_2','horizontalalignment','left','verticalalignment','bottom','fontsize',10)

        end

        legend()

    end

    % ---------------------
    % End of Static methods

    end

end
