% =======
% Manager
% =======

function Manager(Config)

    % Read all data
    if (exist(Config.InputLogFilename,'file') == 2) && (Config.UseInputLogFileFromPreviousRun == true)

        % Input data has already read into mat file.
        InputLog = load(Config.InputLogFilename);

        if isfield(InputLog,'Data')

            fprintf('Reading input data from: %s\n',Config.InputLogFilename);
            Data = InputLog.Data;

        else
            error('InputLog.mat does not contain Data struct.')
        end
    else
        % Input data are not exist in log file. Read all files from ADCP raw data.
        disp('Reading data from raw ADCP files. May take moments ...');
        Data = ReadData.ConcatenateData(Config);

        % Save input data into input log file, which makes the future calls faster.
        if Config.WriteInputDataToInputLogFile == true
            save(Config.InputLogFilename,'Data')
            fprintf('Wrote input data to: %s.\n',Config.InputLogFilename)
        end
    end

    % Adjust time shift (from original to UTC or anything that makes this array usable)
    TimeShiftInSeconds = Config.TimeShift * (3600.0);    % Config.TimeShift is assumed to be in hours
    Data.Time = Data.Time + TimeShiftInSeconds;

    % Save Ship Tracks
    % ShipTrack.Time = Data.Time / (24 * 3600);
    % ShipTrack.TimeUnit = 'Days since Jan 0, 0000 at 00:00:00 (UTC)';
    % ShipTrack.Longitude = Data.Ship.Longitude;
    % ShipTrack.Latitude = Data.Ship.Latitude;
    % save('ShipTrack.mat','ShipTrack')
    % fprintf('Wrote to: ShipTrack.mat\n')

    % Cut data in time
    Data = CutDataInTime(Data,Config);

    % Print times
    fprintf('------------------------------------------------\n')
    fprintf('Times summary:\n')
    fprintf('First data time:            %s\n',datestr(Data.Time(1) / (24*3600)))
    fprintf('Last  data time:            %s\n',datestr(Data.Time(end) / (24*3600)))
    fprintf('Start plume release time:   %s\n',datestr(Config.PlumeRelease.StartReleaseTime))
    fprintf('End   plume release time:   %s\n',datestr(Config.PlumeRelease.EndReleaseTime))
    fprintf('Inquiry time:               %s\n',datestr(Config.PlumeInquiry.InquiryTime))
    fprintf('-------------------------------------------------\n\n')

    % Plot Data
    if Config.Plots.ExtraPlots == true
        PlotData.PlotMeasurementData(Data,Config);
        PlotData.PlotVelocityColumnAtGivenTime(Data,Config,Config.PlumeInquiry.InquiryTime);
        PlotData.PlotShipAndAveragedOceanVelocities(Data,Config)
    end

    % Read wind data
    if Config.Wind.IncludeWindData == true
        Data.Wind.Field = ReadWind.ReadWindFile(Config);
        Data.Wind.OnShipTrack = ReadWind.InterpolateWindOnShipTracks(Config,Data);
    end

    % Load Log
    ComputeStreakLineFlag = true;
    StoreStreakLineFlag = true;
    % if (exist(Config.OutputLogFilename,'file') == 2) && (Config.UseOutputLogFileFromPreviousRun == true)
    %     
    %     % Already exists. Load form previous log
    %     Log = load(Config.OutputLogFilename);
    % 
    %     if isfield(Log,'Output')
    %         
    %         disp('Streakline already computed. Loading form log file.');
    %         Output = Log.Output;
    % 
    %         if isfield(Output,'PlumeInitialStreakLine')
    % 
    %         % Check streakline stored before
    %             if Output.PlumeInitialStreakLine.PlumeReleaseEnded == true
    % 
    %                 % do not sore current steakline computation if it already computed and reached to its end
    %                 StoreStreakLineFlag = false;
    % 
    %             end
    % 
    %             % Check if previous streakline computation was complete
    %             if Output.PlumeInitialStreakLine.FinalStreakLineTime >= Output.PlumeInitialStreakLine.EndReleaseTime - 10*eps
    % 
    %                 % Streakline is traced till end of plume release. No recomputation needed.
    %                 ComputeStreakLineFlag = false;
    %             end
    %         end
    %     end
    % end

    if ComputeStreakLineFlag == true
         
        % Trace Initial Plume StreakLine
        Output = StreakLine.TracePlumeInitialStreakLine(Data,Config);

        % Plot Plume Initial StreakLine
        if Config.Plots.ExtraPlots == true
            StreakLine.PlotPlumeInitialStreakLine(Output.PlumeInitialStreakLine,Data,Config);
        end
        
    end

    % Find Plume At Inquiry Time
    if Output.PlumeInitialStreakLine.PlumeReleaseEnded == true

        % Use stochastic or deterministic model
        if Config.StochasticModel.UseStochasticModel == true

            % Stochastic model
            [Output,ErrorCode] = StochasticModel.StochasticPlumeTrajectoryAtInquiryTime(Data,Config,Output);

            % Plot trajectory
            if ErrorCode == 0
                StochasticModel.PlotPlumeTrajectory(Data,Config,Output);
            else
                return
            end

        else

            % Deterministic model
            Output = DeterministicModel.FindPlumeAtInquiryTime(Data,Config,Output);
         
            % Plot trajectory
            DeterministicModel.PlotPlumeTrajectory(Data,Config,Output);
        end
    end

    % Save results to log
    if Config.WriteOutputDataToOutputLogFile == true
        save(Config.OutputLogFilename,'Output');
        fprintf('Wrote output data to: %s.\n',Config.InputLogFilename)
    end

    % Save all data for python plots (external codes)
    if Config.SaveForPythonPlots == true
        save(fullfile(Config.SaveForPythonDirectory,'/Data'),'Data')
        save(fullfile(Config.SaveForPythonDirectory,'/Config'),'Config')
        save(fullfile(Config.SaveForPythonDirectory,'/Output'),'Output')

        % Save drifters iteratively. It is important to save each drifter as a separate "scalar" struct, not a vector struct, since they cannot be opened in python while keeping the struture properly
        DriftersFileData = load(Config.KalmanFilter.DriftersFileData);
        Drifters = DriftersFileData.Drifters;
        for i = 1:size(Drifters,2)
            eval(sprintf('Drifter%d = Drifters(%d);',i,i));
            save(fullfile(Config.SaveForPythonDirectory,'drifters',sprintf('/Drifters%d.mat',i)),sprintf('Drifter%d',i))
        end
    end

end

% ================
% Cut Data In Time
% ================

function Data = CutDataInTime(Data,Config)

    % Cuts all arrays in Data struct that have time dimension. The time dimension
    % is cut before initial plume release and after inquiry time of plume to shorten
    % the long array of data. This enhances the interpolation when several (or all)
    % adcp files are loaded.

    DayToSecond = 24.0 * 60.0 * 60.0;

    % Find earliest time used in computation
    if Config.PlumeRelease.ComputeStreakLine == true
        % Start from plume release
        FirstTimeInUse = datenum(Config.PlumeRelease.StartReleaseTime) * DayToSecond;
    else
        % Start from tracking plume (plume already released)
        if Config.PlumeLastSeen.UseUpdatedPosition == true
            % Use user-defined time to track plume
            FirstTimeInUse = datenum(Config.PlumeLastSeen.LastObservedTime) * DayToSecond;
        else
            % Use end time of plume release as start of plume track
            FirstTimeInUse = datenum(Config.PlumeRelease.EndReleaseTime) * DayToSecond;
        end
    end

    % Find last time used in computation
    if Config.PlumeInquiry.UseLastAvailableTime == true
        % Use the end of whole data as the last time
        LastTimeInUse = Data.Time(end); % no day to second conversion is needed
    else
        % Use user-defined inquiry time
        LastTimeInUse = datenum(Config.PlumeInquiry.InquiryTime) * DayToSecond;
    end

    % Check if time is monotonically increasing
    TimeIndexNotIncreasing = find(diff(Data.Time) <= 0);
    if any(TimeIndexNotIncreasing)

        disp('Time indices that Data.Time is not monotonically increasing are:')
        disp(TimeIndexNotIncreasing')
        disp('The size of Data.Time array:')
        disp(size(Data.Time))
        error('Data.Time is not monotonically increasing. Program terminated.')
        return
    end

    % Find index of the first time in use
    FirstTimeIndexInUse = find(Data.Time <= FirstTimeInUse,1,'last') - Config.ExtraTimeCutMarginIndex;
    LastTimeIndexInUse = find(Data.Time >= LastTimeInUse,1,'first') + Config.ExtraTimeCutMarginIndex;

    % Check if indices are within the range indices of data
    if FirstTimeIndexInUse < 1
        FirstTimeIndexInUse = 1;
    end
    if LastTimeIndexInUse > size(Data.Time,2)
        LastTimeIndexInUse = size(Data.Time,2);
    end

    % Check if first and last cut indices are not overlapping
    if LastTimeIndexInUse <= FirstTimeIndexInUse
        error('First and last time indices in use are overlaping. First index: %d, last index: %d',FirstTimeIndexInUse,LastTimeIndexInUse)
        return
    end

    % Cut Data in time
    Data.Time                = Data.Time(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Ship.Longitude      = Data.Ship.Longitude(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Ship.Latitude       = Data.Ship.Latitude(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Ship.Heading        = Data.Ship.Heading(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Ship.Pitch          = Data.Ship.Pitch(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Ship.Roll           = Data.Ship.Roll(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Ship.Madegd         = Data.Ship.Madegd(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Ship.Vel_u          = Data.Ship.Vel_u(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Ship.Vel_v          = Data.Ship.Vel_v(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Ocean.AbsoluteVel_u = Data.Ocean.AbsoluteVel_u(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Ocean.AbsoluteVel_v = Data.Ocean.AbsoluteVel_v(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Ocean.AbsoluteVel_w = Data.Ocean.AbsoluteVel_w(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Beam.Intensity      = Data.Beam.Intensity(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Beam.Range          = Data.Beam.Range(:,FirstTimeIndexInUse:LastTimeIndexInUse);
    Data.Beam.Correlation    = Data.Beam.Correlation(:,:,FirstTimeIndexInUse:LastTimeIndexInUse);

end
