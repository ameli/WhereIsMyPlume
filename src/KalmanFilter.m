classdef KalmanFilter

    methods(Static)

    % ====================
    % Hybrid Kalman Filter
    % ====================

    function [InquiryTrajectoriesAtControlTimesT, ...
        InquiryTrajectoriesXStandardDeviationAtControlTimes,InquiryTrajectoriesYStandardDeviationAtControlTimes, ...
        InquiryTrajectoriesAtControlTimesX,InquiryTrajectoriesAtControlTimesY, ...
        ControlTrajectoriesAtControlTimesX,ControlTrajectoriesAtControlTimesY] = ...
        HybridKalmanFilter( ...
        Config, ...
        TrajectoriesStartX,TrajectoriesStartY, ...
        ControlDriftersTimes,ControlDriftersX,ControlDriftersY, ...
        OceanVelocity, ...
        BasisPolynomialOrder,ScaleLength,a_u,a_v, ...
        VelocityStandardDeviationAverage,SpatialDecorrelationScale,SpatialDecorrelationCutoff)
    
        % Number of control drifters
        NumberOfControlDrifters = size(ControlDriftersX,2);

        % Number of inquiry points (non drifter trajectories)
        NumberOfInquiryPoints = size(TrajectoriesStartX,2);

        % Number of all trajectories (both inquiry and control drifters combined)
        NumberOfAllTrajectories = NumberOfInquiryPoints + NumberOfControlDrifters;

        % Shorthand variable names
        M = NumberOfInquiryPoints;
        N = NumberOfAllTrajectories;

        % Dimension
        Dimension = 2;

        % Size of State vector
        StateVectorSize = Dimension * NumberOfAllTrajectories;

        % Number of drifter control times
        NumberOfDrifterControlTimes = length(ControlDriftersTimes);

        % Array of solutions (inquiry trajectories)
        InquiryTrajectoriesAtControlTimesT = zeros(NumberOfDrifterControlTimes,1);
        InquiryTrajectoriesAtControlTimesX = zeros(NumberOfDrifterControlTimes,NumberOfInquiryPoints);
        InquiryTrajectoriesAtControlTimesY = zeros(NumberOfDrifterControlTimes,NumberOfInquiryPoints);
        InquiryTrajectoriesXStandardDeviationAtControlTimes = zeros(NumberOfDrifterControlTimes,NumberOfInquiryPoints);
        InquiryTrajectoriesYStandardDeviationAtControlTimes = zeros(NumberOfDrifterControlTimes,NumberOfInquiryPoints);

        % Array of solutions (control trajectory, when covariance R is zero, these are the same as ControlDritfter trajectories)
        ControlTrajectoriesAtControlTimesX = zeros(NumberOfDrifterControlTimes,NumberOfControlDrifters);
        ControlTrajectoriesAtControlTimesY = zeros(NumberOfDrifterControlTimes,NumberOfControlDrifters);

        % Initialize array of trajectories
        InquiryTrajectoriesAtControlTimesT(1)   = ControlDriftersTimes(1);
        InquiryTrajectoriesAtControlTimesX(1,:) = TrajectoriesStartX;
        InquiryTrajectoriesAtControlTimesY(1,:) = TrajectoriesStartY;

        ControlTrajectoriesAtControlTimesX(1,:) = ControlDriftersX(1,:);
        ControlTrajectoriesAtControlTimesY(1,:) = ControlDriftersY(1,:);
        
        % Observation matrix
        H_ux = zeros(NumberOfControlDrifters,NumberOfAllTrajectories);
        for i = 1:NumberOfControlDrifters

            % Leave out inquiry points and only retain the drifter control points for observation
            H_ux(i,i+NumberOfInquiryPoints) = 1;
        end
        H_vy = H_ux;
        H_uy = zeros(size(H_ux));
        H_vx = zeros(size(H_ux));
        H = [H_ux,H_uy
             H_vx,H_vy];

        % Initialize covariance
        % P_prediction = zeros(StateVectorSize,StateVectorSize);
        P_current = zeros(StateVectorSize,StateVectorSize);
 
        % Ietrate over drifter control times
        for IntervalIterator = 1:NumberOfDrifterControlTimes - 1

            % Start and End Times
            TrajectoryStartTime = ControlDriftersTimes(IntervalIterator);
            TrajectoryEndTime = ControlDriftersTimes(IntervalIterator+1);

            % Append control trajectories to the list of trajectories
            AllTrajectoriesStartX = [InquiryTrajectoriesAtControlTimesX(IntervalIterator,:),ControlTrajectoriesAtControlTimesX(IntervalIterator,:)];
            AllTrajectoriesStartY = [InquiryTrajectoriesAtControlTimesY(IntervalIterator,:),ControlTrajectoriesAtControlTimesY(IntervalIterator,:)];

            % Covariance of observation (matrix R)
            R = Config.KalmanFilter.ObservationStandardDeviation^2 * ones(Dimension*NumberOfControlDrifters,Dimension*NumberOfControlDrifters);

            % Prediction for state (x and y of trajectories)
            [AllTrajectoriesX_prediction,AllTrajectoriesY_prediction] = KalmanFilter.PredictState( ...
                TrajectoryStartTime,TrajectoryEndTime, ...
                AllTrajectoriesStartX,AllTrajectoriesStartY, ...
                OceanVelocity);

            % Prediction for state covariance (covariance P between x and y of all trajectories)
            P_prediction = KalmanFilter.PredictStateCovariance( ...
                Config, ...
                TrajectoryStartTime,TrajectoryEndTime, ...
                AllTrajectoriesStartX,AllTrajectoriesStartY,P_current, ...
                BasisPolynomialOrder,ScaleLength,a_u,a_v, ...
                VelocityStandardDeviationAverage,SpatialDecorrelationScale,SpatialDecorrelationCutoff);

            % Choose whether update the prediction based on observations or not
            if Config.KalmanFilter.DoTheUpdateStep == false

                % Store output (from predicton) COMMNET these lines if you want to update the prediction from observations
                InquiryTrajectoriesAtControlTimesT(IntervalIterator+1) = TrajectoryEndTime;
                InquiryTrajectoriesAtControlTimesX(IntervalIterator+1,:) = AllTrajectoriesX_prediction(1:NumberOfInquiryPoints);
                InquiryTrajectoriesAtControlTimesY(IntervalIterator+1,:) = AllTrajectoriesY_prediction(1:NumberOfInquiryPoints);

                ControlTrajectoriesAtControlTimesX(IntervalIterator+1,:) = AllTrajectoriesX_prediction(NumberOfInquiryPoints+1:end);
                ControlTrajectoriesAtControlTimesY(IntervalIterator+1,:) = AllTrajectoriesY_prediction(NumberOfInquiryPoints+1:end);

                InquiryTrajectoriesXStandardDeviationAtControlTimes(IntervalIterator+1,:) = sqrt(diag(P_prediction(1:M,1:M)));
                InquiryTrajectoriesYStandardDeviationAtControlTimes(IntervalIterator+1,:) = sqrt(diag(P_prediction(N+1:N+M,N+1:N+M)));

                % Update covariance
                P_current = P_prediction;

            else

                % Observations at end of time interval
                ControlDriftersXAtUpdateTime = ControlDriftersX(IntervalIterator+1,:);
                ControlDriftersYAtUpdateTime = ControlDriftersY(IntervalIterator+1,:);

                % Update state (trajecotries) and state covariance (matrix P)
                [AllTrajectoriesX_update,AllTrajectoriesY_update,P_update] = KalmanFilter.UpdateStateAndStateCovariance( ...
                    AllTrajectoriesX_prediction,AllTrajectoriesY_prediction, ...
                    ControlDriftersXAtUpdateTime,ControlDriftersYAtUpdateTime, ...
                    H,P_prediction,R);
                
                % Store output (from updates)
                InquiryTrajectoriesAtControlTimesT(IntervalIterator+1) = TrajectoryEndTime;
                InquiryTrajectoriesAtControlTimesX(IntervalIterator+1,:) = AllTrajectoriesX_update(1:NumberOfInquiryPoints);
                InquiryTrajectoriesAtControlTimesY(IntervalIterator+1,:) = AllTrajectoriesY_update(1:NumberOfInquiryPoints);

                ControlTrajectoriesAtControlTimesX(IntervalIterator+1,:) = AllTrajectoriesX_update(NumberOfInquiryPoints+1:end);
                ControlTrajectoriesAtControlTimesY(IntervalIterator+1,:) = AllTrajectoriesY_update(NumberOfInquiryPoints+1:end);

                InquiryTrajectoriesXStandardDeviationAtControlTimes(IntervalIterator+1,:) = sqrt(diag(P_update(1:M,1:M)));
                InquiryTrajectoriesYStandardDeviationAtControlTimes(IntervalIterator+1,:) = sqrt(diag(P_update(N+1:N+M,N+1:N+M)));

                % Update covariance
                P_current = P_update;

            end

        end

    end

    % =============
    % Predict State
    % =============

    function [AllTrajectoriesX_prediction,AllTrajectoriesY_prediction] = PredictState( ...
        TrajectoryStartTime,TrajectoryEndTime, ...
        AllTrajectoriesStartX,AllTrajectoriesStartY, ...
        OceanVelocity)

        % The ODE for positions (state vector) is integrated using Runge-Kutta method.

        % Number of trajectories
        AllNumberOfTrajectories = length(AllTrajectoriesStartX);

        % Trace trajectory
        opts = odeset('RelTol',1e-5);  % SETTING
        [TrajectoryTimes,TrajectoryCoordinates] = ...
            ode45(OceanVelocity,[TrajectoryStartTime,TrajectoryEndTime],[AllTrajectoriesStartX,AllTrajectoriesStartY],opts);

        % Output
        AllTrajectoriesX_prediction = TrajectoryCoordinates(end,1:AllNumberOfTrajectories);
        AllTrajectoriesY_prediction = TrajectoryCoordinates(end,AllNumberOfTrajectories+1:end);

    end

    % ========================
    % Predict State Covariance
    % ========================

    function P_prediction = PredictStateCovariance( ...
        Config, ...
        TrajectoryStartTime,TrajectoryEndTime, ...
        AllTrajectoriesStartX,AllTrajectoriesStartY,P_current, ...
        BasisPolynomialOrder,ScaleLength,a_u,a_v, ...
        VelocityStandardDeviationAverage,SpatialDecorrelationScale,SpatialDecorrelationCutoff)

        % The covariance ODE is integrated based on Euler backward method.

        % Concatenate X and Y coordinates in one array
        X_InquiryPoints = [AllTrajectoriesStartX',AllTrajectoriesStartY'];

        % Jacobian
        J = BasisFunctions.JacobianOfBasisFunctions( ...
            BasisPolynomialOrder, ...
            X_InquiryPoints, ...
            ScaleLength, ...
            a_u,a_v);

        % Covariance either east or north velocity between all trajectory points (inquiry points and control drifter points)
        Q_uu = Kriging.CovarianceOfInquiryPointsWithItself( ...
            Config, ...
            X_InquiryPoints, ...
            VelocityStandardDeviationAverage, ...
            SpatialDecorrelationScale,SpatialDecorrelationCutoff);

        % We assumed covariance of east velocity and covariance of north velocty are the same
        Q_vv = Q_uu;

        % Covariance between u and v is assumed to be zero
        Q_uv = zeros(size(Q_uu));

        % Concatenate
        Q = [Q_uu,Q_uv;
             Q_uv,Q_vv];

        % Forward covariance error using forward Euler
        Delta_t = TrajectoryEndTime - TrajectoryStartTime;
        P_prediction = P_current + Delta_t * (J*P_current + P_current*J' + Q);

    end

    % ============
    % Update State
    % ============

    function [AllTrajectoriesX_update,AllTrajectoriesY_update,P_update] = UpdateStateAndStateCovariance( ...
        AllTrajectoriesX_prediction,AllTrajectoriesY_prediction, ...
        ControlDriftersXAtUpdateTime,ControlDriftersYAtUpdateTime, ...
        H,P_prediction,R)

        % This is the second stage of kalman filter, where the prediction from the first stage is 
        % updated based on the "innovtion", or misfit of the projection of state vector with the
        % observed drifter data.

        % State vector prediction
        State_prediction = [AllTrajectoriesX_prediction,AllTrajectoriesY_prediction]';

        % Observation vector
        Observation = [ControlDriftersXAtUpdateTime,ControlDriftersYAtUpdateTime]';

        % Innovation
        Innovation = Observation - H*State_prediction;

        % Observation covariance
        S = H * P_prediction * H' + R;

        % Kalman gain
        K = P_prediction * H' * inv(S);

        % Update state vector
        State_update = State_prediction + K * Innovation;

        % Output
        NumberOfAllTrajectories = length(AllTrajectoriesX_prediction);
        AllTrajectoriesX_update = State_update(1:NumberOfAllTrajectories);
        AllTrajectoriesY_update = State_update(NumberOfAllTrajectories+1:end);

        % Update state covariance
        P_update = (eye(size(P_prediction)) - K*H) * P_prediction;

    end

    % ---------------------
    % End of Static methods

    end
end
