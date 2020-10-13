classdef ReadData

    methods(Static)

    % ================
    % Concatenate Data
    % ================

    function Data = ConcatenateData(Config)

        % Reads all files and returns a single Data struct including all data

        % FullPath filename
        FullPathFilenames = fullfile(Config.DataDirectory,Config.Filename);

        ListOfFilesStruct = dir(FullPathFilenames);

        if isempty(ListOfFilesStruct)
            error('List of files is empty. Filenames: %s',FullPathFilenames)
        end

        ListOfDataCell = cell(1,size(ListOfFilesStruct,1));
        FileBeginTimes = zeros(1,size(ListOfFilesStruct,1));
        FileMedianTimes = zeros(1,size(ListOfFilesStruct,1));
        FileEndTimes = zeros(1,size(ListOfFilesStruct,1));

        % Store data of each file in a separate struct.
        for FileId = 1:size(ListOfFilesStruct,1)

            % Get Filename
            FullPathFilename = fullfile(Config.DataDirectory,ListOfFilesStruct(FileId).name);

            % Read variables from ADCP file (Choose one of the two functions below)
            if Config.ReadFromMatFiles == true
                RawData_File = ReadData.GetRawDataFromMATFile(FullPathFilename,Config);      % Uses MAT files. Note: velocities are in earth coordinates
            else
                RawData_File = ReadData.GetRawDataFromLTAFile(FullPathFilename,Config);      % Uses LTA files. Note: velocities are relative to ship. Transformation needed.
            end

            % Remove corrupted time indices from all arrays
            RawData_File = ReadData.RemoveCorruptedTimeIndicesFromRawData(RawData_File);

            % Fill missing
            RawData_File = ReadData.FillMissingInTime(RawData_File,Config);

            % Process raw data (that that we can use)
            Data_File = ReadData.ProcessRawData(RawData_File,Config);

            ListOfDataCell{FileId} = Data_File;

            FileBeginTimes(FileId) = Data_File.Time(1);
            FileMedianTimes(FileId) = median(Data_File.Time);
            FileEndTimes(FileId) = Data_File.Time(end);

        end

        % Sort times
        [~,SortingIndex] = sort(FileMedianTimes);

        % Add the first file to the Data struct
        Data = ListOfDataCell{SortingIndex(1)};

        % Display times
        disp(' ')
        disp('------------------------------------------------')
        disp('File concatenation summary:')
        disp(sprintf('File %03d: start: %s, end: %s',1,datestr(FileBeginTimes(SortingIndex(1))/(24*3600)),datestr(FileEndTimes(SortingIndex(1))/(24*3600))))

        % Concatenate only variables that have time argument. Other variables exist in Data struct from the first file.
        if size(SortingIndex,2) > 1
            for i = 2:size(SortingIndex,2)

                disp(sprintf('File %03d: start: %s, end: %s',i,datestr(FileBeginTimes(SortingIndex(i))/(24*3600)),datestr(FileEndTimes(SortingIndex(i))/(24*3600))))
                FileId = SortingIndex(i);
                Data_Next = ListOfDataCell{FileId};

                % Concatenate
                Data.Time = cat(2,Data.Time,Data_Next.Time);
                Data.Ship.Longitude = cat(2,Data.Ship.Longitude,Data_Next.Ship.Longitude);                   % Concatenate second (time) dimension
                Data.Ship.Latitude = cat(2,Data.Ship.Latitude,Data_Next.Ship.Latitude);                      % Concatenate second (time) dimension
                Data.Ship.Heading = cat(2,Data.Ship.Heading,Data_Next.Ship.Heading);                         % Concatenate second (time) dimension 
                Data.Ship.Pitch = cat(2,Data.Ship.Pitch,Data_Next.Ship.Pitch);                               % Concatenate second (time) dimension 
                Data.Ship.Roll = cat(2,Data.Ship.Roll,Data_Next.Ship.Roll);                                  % Concatenate second (time) dimension 
                Data.Ship.Madegd = cat(2,Data.Ship.Madegd,Data_Next.Ship.Madegd);                            % Concatenate second (time) dimension 
                Data.Ship.Vel_u = cat(2,Data.Ship.Vel_u,Data_Next.Ship.Vel_u);                               % Concatenate second (time) dimension 
                Data.Ship.Vel_v = cat(2,Data.Ship.Vel_v,Data_Next.Ship.Vel_v);                               % Concatenate second (time) dimension
                Data.Ocean.AbsoluteVel_u = cat(2,Data.Ocean.AbsoluteVel_u,Data_Next.Ocean.AbsoluteVel_u);    % Concatenate second (time) dimension 
                Data.Ocean.AbsoluteVel_v = cat(2,Data.Ocean.AbsoluteVel_v,Data_Next.Ocean.AbsoluteVel_v);    % Concatenate second (time) dimension 
                Data.Ocean.AbsoluteVel_w = cat(2,Data.Ocean.AbsoluteVel_w,Data_Next.Ocean.AbsoluteVel_w);    % Concatenate second (time) dimension 
                Data.Beam.Intensity = cat(2,Data.Beam.Intensity,Data_Next.Beam.Intensity);                   % Concatenate second (time) dimension  
                Data.Beam.Range = cat(2,Data.Beam.Range,Data_Next.Beam.Range);                               % Concatenate second (time) dimension  
                Data.Beam.Correlation = cat(3,Data.Beam.Correlation,Data_Next.Beam.Correlation);             % Concatenate third  (time) dimension 

            end
        end

        % --------------------------------
        % Check Concatenated Variable Size (nested function)
        % --------------------------------

        function CheckConcatenatedVariableSize(Variable,VariableName,ExpectedVariableSize)

            fprintf('Checking concatenated %s',VariableName)

            if any(size(Variable) ~= ExpectedVariableSize)
                ErrorMessage = sprintf('The concatenated variable %s has size [%s] but expected size is [%s].', ...
                VariableName,num2str(size(Variable)),num2str(ExpectedVariableSize));

                error(ErrorMessage)
            else
                fprintf(' ... OK.\n')
            end

        end

        % --------------------

        % Check dimensions
        D = size(Data.Ocean.AbsoluteVel_u,1);   % Number of depth bins
        T = size(Data.Ocean.AbsoluteVel_u,2);   % Number of times

        CheckConcatenatedVariableSize(Data.Time,'Time',[1,T])
        CheckConcatenatedVariableSize(Data.Ship.Longitude,'Longitude',[1,T])
        CheckConcatenatedVariableSize(Data.Ship.Latitude,'Latitude',[1,T])
        CheckConcatenatedVariableSize(Data.Ship.Heading,'Heading',[1,T])
        CheckConcatenatedVariableSize(Data.Ship.Pitch,'Pitch',[1,T])
        CheckConcatenatedVariableSize(Data.Ship.Roll,'Roll',[1,T])
        CheckConcatenatedVariableSize(Data.Ship.Madegd,'Madegd',[1,T])
        CheckConcatenatedVariableSize(Data.Ship.Vel_u,'Vel_u',[1,T])
        CheckConcatenatedVariableSize(Data.Ship.Vel_v,'Vel_v',[1,T])
        CheckConcatenatedVariableSize(Data.Ocean.AbsoluteVel_u,'AbsoluteVel_u',[D,T])
        CheckConcatenatedVariableSize(Data.Ocean.AbsoluteVel_v,'AbsoluteVel_v',[D,T])
        CheckConcatenatedVariableSize(Data.Ocean.AbsoluteVel_w,'AbsoluteVel_w',[D,T])
        CheckConcatenatedVariableSize(Data.Beam.Intensity,'Intensity',[D,T])
        CheckConcatenatedVariableSize(Data.Beam.Range,'BeamRange',[1,T])
        CheckConcatenatedVariableSize(Data.Beam.Correlation,'Correlation',[D,4,T])

        disp('')

    end

    % ==========================
    % Get Raw Data From LTA File
    % ==========================

    function RawData = GetRawDataFromLTAFile(FullPathFilename,Config)

        % There are two functions that may be used to read raw data.
        %
        % 1. ReadData.GetRawDataFromLTAFile (this function)
        % 2. ReadData.GetRawDataFromMatFile 
        %
        % This function reads variables from LTA or STA file using external code (rdradcp.m)
        % The rdradcp.m has some variables to set. For instance, you may set num_av, the number of
        % ensembles to block together and average. Default in that code is num_av = 5, which shrinks
        % the time length by the factor of 5. For instance, an LTA file of time length 726 data points
        % becomes 145 length. 
        %
        % The file is read into the adcp struct that contains all variables.%
        % This function a struct RawData containing some of the variables that
        % are in the adcp struct.

        % Velocities in LTA files are relative to ship coordinates.
        if Config.OceanVelocityIsRelativeToShip == false;
            error('When using *.LTA files, enable velocity coordinate transformation since they are in ship coordinate.')
        end

        % Read ADCP struct (external code)
        adcp = rdradcp(FullPathFilename);

        % Ship coordinates (degrees)
        RawData.Ship.StartLongitude = adcp.nav_slongitude;    % Degrees
        RawData.Ship.StartLatitude  = adcp.nav_slatitude;     % Degrees
        RawData.Ship.EndLongitude   = adcp.nav_elongitude;    % Degrees
        RawData.Ship.EndLatitude    = adcp.nav_elatitude;     % Degrees

        % Ship Speed and Directions, navigation angles
        RawData.Ship.Speed   = adcp.nav_speedmadegd;          % m/s. Speed made good based on GPS.
        RawData.Ship.Madegd  = adcp.nav_direcmadegd;          % Degrees, zero is north, clockwise. Direction made good, based on GPS.
        RawData.Ship.Heading = adcp.heading;                  % Degrees, zero is north, clockwise
        RawData.Ship.Pitch   = adcp.pitch;                    % Degrees
        RawData.Ship.Roll    = adcp.roll;                     % Degrees

        % Data times (days from 0000 January 01) in unit of "days". From time variabels below, choose adcp.mtime
        RawData.Time = adcp.mtime;           % Shows 3:53' hours behind EST clock. Rate is the same as time. Use this time with +4 shift.
        % RawData.Time = adcp.nav_smtime;        % Shows 4:00' hours after  EST clock, Rate is the same as time. Use this time with 0 shift. (Not monotonocally increasing)
        % RawData.Time = adcp.nav_emtime;      % Shows 4:00' hours after  EST clock, Rate is the same as time. 
        % RawData.Time = adcp.nav_mtime;       % This seems to be wrong. Rate is 1/100 of the time rate. Do NOT use this time.

        % Ocean data depths at center of cell bins (meters, positive downward)
        RawData.Ocean.Range = adcp.config.ranges;

        % Ocean measured velocity (m/s, relative to ship)
        RawData.Ocean.MeasuredVel_u = adcp.east_vel;
        RawData.Ocean.MeasuredVel_v = adcp.north_vel;
        RawData.Ocean.MeasuredVel_w = adcp.vert_vel;

        % BottomTrack (bt) velocity of ground relative to ship. bt_vel is negative of ship velocity
        MillimeterToMeter = 1000;
        RawData.Ocean.BottomTrack_Vel_u = adcp.bt_vel(1,:) / MillimeterToMeter;   % in m/s unit
        RawData.Ocean.BottomTrack_Vel_v = adcp.bt_vel(2,:) / MillimeterToMeter;   % in m/s unit

        % BottomTrack (bt) data: Ocean range for each beam
        RawData.Beam.Range = adcp.bt_range;      % meters

        % Beam Correlations
        RawData.Beam.TotalCorrelationCounts = adcp.bt_corr;   % Total correlation count, which is around 255
        RawData.Beam.CorrelationCounts = adcp.corr;           % Correlation count for each beam at each cell (between adcp.config.corr_threshod=64 and 255)

        % Beam intensity
        RawData.Beam.Intensity = adcp.intens;

        % Pings per ensemble
        RawData.Beam.PingsPerEnsemble = adcp.config.pings_per_ensemble;

        % Beam configuration
        RawData.Beam.Angle = adcp.config.beam_angle;          % Degrees
        RawData.Beam.Frequency = adcp.config.beam_freq;       % Hertz
        RawData.Beam.CellSize = adcp.config.cell_size;        % meters
        RawData.Beam.ProfileMode = adcp.config.prof_mode;     % integer
        RawData.Beam.LagLength = adcp.config.laglength;       % Number of bits in each period

        % -------------------
        % Check Variable Size (nested function)
        % -------------------

        function CheckVariableSize(Filename,Variable,VariableName,ExpectedVariableSize)

            fprintf('Checking %s',VariableName)

            if any(size(Variable) ~= ExpectedVariableSize)
                ErrorMessage = sprintf('In file %s, the variable %s has size [%s] but expected size is [%s].', ...
                Filename,VariableName,num2str(size(Variable)),num2str(ExpectedVariableSize));

                error(ErrorMessage)
            else
                fprintf(' ... OK.\n')
            end

        end

        % --------------------

        % Use the size of east velocity array to check with other variables sizes
        D = size(RawData.Ocean.MeasuredVel_u,1);  % Number of depth bins
        T = size(RawData.Ocean.MeasuredVel_u,2);  % Number of times

        % Check consistency of variable sizes
        CheckVariableSize(FullPathFilename,adcp.nav_slongitude,'nav_slongitude',[1,T])
        CheckVariableSize(FullPathFilename,adcp.nav_elongitude,'nav_elongitude',[1,T])
        CheckVariableSize(FullPathFilename,adcp.nav_slatitude,'nav_slatitude',[1,T])
        CheckVariableSize(FullPathFilename,adcp.nav_elatitude,'nav_elatitude',[1,T])
        CheckVariableSize(FullPathFilename,adcp.nav_speedmadegd,'nav_speedmadegd',[1,T])
        CheckVariableSize(FullPathFilename,adcp.nav_direcmadegd,'nav_direcmadegd',[1,T])
        CheckVariableSize(FullPathFilename,adcp.heading,'heading',[1,T])
        CheckVariableSize(FullPathFilename,adcp.pitch,'pitch',[1,T])
        CheckVariableSize(FullPathFilename,adcp.roll,'roll',[1,T])
        CheckVariableSize(FullPathFilename,adcp.mtime,'mtime',[1,T])
        CheckVariableSize(FullPathFilename,adcp.config.ranges,'ranges',[D,1])
        CheckVariableSize(FullPathFilename,adcp.east_vel,'east_vel',[D,T])
        CheckVariableSize(FullPathFilename,adcp.north_vel,'north_vel',[D,T])
        CheckVariableSize(FullPathFilename,adcp.vert_vel,'vert_vel',[D,T])
        CheckVariableSize(FullPathFilename,adcp.intens,'intens',[D,4,T])
        CheckVariableSize(FullPathFilename,adcp.bt_range,'bt_range',[4,T])
        CheckVariableSize(FullPathFilename,adcp.bt_vel,'bt_vel',[4,T])
        CheckVariableSize(FullPathFilename,adcp.bt_corr,'bt_corr',[4,T])
        CheckVariableSize(FullPathFilename,adcp.corr,'corr',[D,4,T])

    end

    % ==========================
    % Get Raw Data From MAT File
    % ==========================

    function RawData = GetRawDataFromMATFile(FullPathFilename,Config)

        % There are two functions that may be used to read raw data.
        %
        % 1. ReadData.GetRawDataFromLTAFile
        % 2. ReadData.GetRawDataFromMatFile (this function)
        %
        % This function reads variables from matlab's mat file.
        % The mat files are produced from the Teledyne's company software, called "Velocity Software".
        % The software exports LTA files to mat files.
        %
        % The file is read into the adcp struct that contains all variables.
        % This function outputs a struct, RawData, containing some of the variables that
        % are in the adcp struct.

        % Velocities in MAT files are relative to earth coordinates.
        if Config.OceanVelocityIsRelativeToShip == true;
            error('When using *.mat files, disable velocity coordinate transformation since they are already in earth coordinate.')
        end

        % Get the LTA filename without directory and extention
        FindSlash = strfind(FullPathFilename,'/');
        LTAFilename = FullPathFilename(FindSlash(end)+1:end);
        FindDot = strfind(LTAFilename,'.');
        Filename = LTAFilename(1:FindDot(1)-1);

        % Construct Mat filename
        MatFileDirectory = Config.MATFileDirectory;
        MatFileExtension = '.mat';
        MatFilename = fullfile(Filename,MatFileExtension);
        FullPathMatFilename = fullfilw(MatFileDirectory,MatFilename);

        % Read ADCP struct (external code)
        NumAv = 1;                                            % This must be set to 1 in rdradcp.m so that time indices are compatible with mat files
        adcp = rdradcp(FullPathFilename,NumAv);               % Loads from LTA files
        adcp_mat = load(FullPathMatFilename);                 % Loads from MAT files

        % Ship coordinates (degrees)
        RawData.Ship.StartLongitude = adcp.nav_slongitude;    % Degrees
        RawData.Ship.StartLatitude  = adcp.nav_slatitude;     % Degrees
        RawData.Ship.EndLongitude   = adcp.nav_elongitude;    % Degrees
        RawData.Ship.EndLatitude    = adcp.nav_elatitude;     % Degrees

        % Ship Speed and Directions, navigation angles
        RawData.Ship.Speed   = adcp.nav_speedmadegd;          % m/s. Speed made good based on GPS.
        RawData.Ship.Madegd  = adcp.nav_direcmadegd;          % Degrees, zero is north, clockwise. Direction made good, based on GPS.
        RawData.Ship.Heading = adcp.heading;                  % Degrees, zero is north, clockwise
        RawData.Ship.Pitch   = adcp.pitch;                    % Degrees
        RawData.Ship.Roll    = adcp.roll;                     % Degrees

        % Data times (days from 0000 January 01) in unit of "days". From time variabels below, choose adcp.mtime
        % RawData.Time = adcp.mtime;               % Shows 3:53' hours behind EST clock. Rate is the same as time. Use this time with +4 shift.
        % RawData.Time = adcp.nav_smtime;          % Shows 4:00' hours after  EST clock, Rate is the same as time. Use this time with 0 shift. (Not monotonocally increasing)
        % RawData.Time = adcp.nav_emtime;          % Shows 4:00' hours after  EST clock, Rate is the same as time. 
        % RawData.Time = adcp.nav_mtime;           % This seems to be wrong. Rate is 1/100 of the time rate. Do NOT use this time.
        TimeOffset = datenum([1970 1 1 0 0 0]);                  % Time in mat files start in Jan 01, 1970 at 00:00:00
        DayToSecond = 24.0 * 3600.0;
        FileTimes = adcp_mat.sens.time';                         % In seconds since Jan 01, 1970, 00:00:00
        RawData.Time = TimeOffset + (FileTimes / DayToSecond);   % In days since Jan 01, 0000, 00:00:00

        % Ocean data depths at center of cell bins (meters, positive downward)
        RawData.Ocean.Range = adcp.config.ranges;

        % Ocean measured velocity (m/s, relative to ship)
        % RawData.Ocean.MeasuredVel_u = adcp.east_vel;
        % RawData.Ocean.MeasuredVel_v = adcp.north_vel;
        % RawData.Ocean.MeasuredVel_w = adcp.vert_vel;
        RawData.Ocean.MeasuredVel_u = adcp_mat.wt.vel(:,:,1)';    % Using mat file. in m/s. In mat file velocity array size is [T,D,4]. We transpose it to be [D,T]
        RawData.Ocean.MeasuredVel_v = adcp_mat.wt.vel(:,:,2)';
        RawData.Ocean.MeasuredVel_w = adcp_mat.wt.vel(:,:,3)';

        % BottomTrack (bt) velocity of ground relative to ship. bt_vel is negative of ship velocity
        MillimeterToMeter = 1000;
        RawData.Ocean.BottomTrack_Vel_u = adcp.bt_vel(1,:) / MillimeterToMeter;   % in m/s unit
        RawData.Ocean.BottomTrack_Vel_v = adcp.bt_vel(2,:) / MillimeterToMeter;   % in m/s unit

        % BottomTrack (bt) data: Ocean range for each beam
        RawData.Beam.Range = adcp.bt_range;      % meters

        % Beam Correlations
        RawData.Beam.TotalCorrelationCounts = adcp.bt_corr;   % Total correlation count, which is around 255
        RawData.Beam.CorrelationCounts = adcp.corr;           % Correlation count for each beam at each cell (between adcp.config.corr_threshod=64 and 255)

        % Beam intensity
        RawData.Beam.Intensity = adcp.intens;

        % Pings per ensemble
        RawData.Beam.PingsPerEnsemble = adcp.config.pings_per_ensemble;

        % Beam configuration
        RawData.Beam.Angle = adcp.config.beam_angle;          % Degrees
        RawData.Beam.Frequency = adcp.config.beam_freq;       % Hertz
        RawData.Beam.CellSize = adcp.config.cell_size;        % meters
        RawData.Beam.ProfileMode = adcp.config.prof_mode;     % integer
        RawData.Beam.LagLength = adcp.config.laglength;       % Number of bits in each period

        % -------------------
        % Check Variable Size (nested function)
        % -------------------

        function CheckVariableSize(Filename,Variable,VariableName,ExpectedVariableSize)

            fprintf('Checking %s',VariableName)

            if any(size(Variable) ~= ExpectedVariableSize)
                ErrorMessage = sprintf('In file %s, the variable %s has size [%s] but expected size is [%s].', ...
                Filename,VariableName,num2str(size(Variable)),num2str(ExpectedVariableSize));

                error(ErrorMessage)
            else
                fprintf(' ... OK.\n')
            end

        end

        % --------------------

        % Use the size of east velocity array to check with other variables sizes
        D = size(RawData.Ocean.MeasuredVel_u,1);  % Number of depth bins
        T = size(RawData.Ocean.MeasuredVel_u,2);  % Number of times

        % Check consistency of variable sizes
        CheckVariableSize(FullPathFilename,adcp.nav_slongitude,'nav_slongitude',[1,T])
        CheckVariableSize(FullPathFilename,adcp.nav_elongitude,'nav_elongitude',[1,T])
        CheckVariableSize(FullPathFilename,adcp.nav_slatitude,'nav_slatitude',[1,T])
        CheckVariableSize(FullPathFilename,adcp.nav_elatitude,'nav_elatitude',[1,T])
        CheckVariableSize(FullPathFilename,adcp.nav_speedmadegd,'nav_speedmadegd',[1,T])
        CheckVariableSize(FullPathFilename,adcp.nav_direcmadegd,'nav_direcmadegd',[1,T])
        CheckVariableSize(FullPathFilename,adcp.heading,'heading',[1,T])
        CheckVariableSize(FullPathFilename,adcp.pitch,'pitch',[1,T])
        CheckVariableSize(FullPathFilename,adcp.roll,'roll',[1,T])
        CheckVariableSize(FullPathFilename,adcp.mtime,'mtime',[1,T])
        CheckVariableSize(FullPathFilename,adcp.config.ranges,'ranges',[D,1])
        CheckVariableSize(FullPathFilename,adcp.east_vel,'east_vel',[D,T])
        CheckVariableSize(FullPathFilename,adcp.north_vel,'north_vel',[D,T])
        CheckVariableSize(FullPathFilename,adcp.vert_vel,'vert_vel',[D,T])
        CheckVariableSize(FullPathFilename,adcp.intens,'intens',[D,4,T])
        CheckVariableSize(FullPathFilename,adcp.bt_range,'bt_range',[4,T])
        CheckVariableSize(FullPathFilename,adcp.bt_vel,'bt_vel',[4,T])
        CheckVariableSize(FullPathFilename,adcp.bt_corr,'bt_corr',[4,T])
        CheckVariableSize(FullPathFilename,adcp.corr,'corr',[D,4,T])

    end

    % ===========================================
    % Remove Corrupted Time Indices From Raw Data
    % ===========================================

    function RawData = RemoveCorruptedTimeIndicesFromRawData(RawData)

        % Often, when the rdradcp.m is called to read adcp data, if we set the num_av (number of averagings inthe code) to 1,2,3, or 4, 
        % some of the time indices (especially the last time index in each adcp array) is not set, and the values are zero.
        % This function removes the columns from arrays that their values are zero. We find which columns are corrupted
        % by looking at the intensity array (adcp.intens).
        %
        % adcp.intens (here, RawData.Beam.Intensity) is array of size [D,4,T], where D is depth, and T is time.

        CorruptedTimeLogical = squeeze(all(all(RawData.Beam.Intensity,1),2));
        CorruptedTimeIndices = find(~CorruptedTimeLogical);

        if ~isempty(CorruptedTimeIndices)

            % Remove corrupted columns from arrays
            RawData.Ship.StartLongitude(CorruptedTimeIndices) = [];
            RawData.Ship.StartLatitude(CorruptedTimeIndices) = [];
            RawData.Ship.EndLongitude(CorruptedTimeIndices) = [];
            RawData.Ship.EndLatitude(CorruptedTimeIndices) = [];

            RawData.Ship.Speed(CorruptedTimeIndices) = [];
            RawData.Ship.Madegd(CorruptedTimeIndices) = [];
            RawData.Ship.Heading(CorruptedTimeIndices) = [];
            RawData.Ship.Pitch(CorruptedTimeIndices) = [];
            RawData.Ship.Roll(CorruptedTimeIndices) = [];

            RawData.Time(CorruptedTimeIndices) = [];

            RawData.Ocean.MeasuredVel_u(:,CorruptedTimeIndices) = [];
            RawData.Ocean.MeasuredVel_v(:,CorruptedTimeIndices) = [];
            RawData.Ocean.MeasuredVel_w(:,CorruptedTimeIndices) = [];

            RawData.Ocean.BottomTrack_Vel_u(CorruptedTimeIndices) = [];
            RawData.Ocean.BottomTrack_Vel_v(CorruptedTimeIndices) = [];

            RawData.Beam.Range(:,CorruptedTimeIndices) = [];
            RawData.Beam.TotalCorrelationCounts(:,CorruptedTimeIndices) = [];
            RawData.Beam.CorrelationCounts(:,:,CorruptedTimeIndices) = [];
            RawData.Beam.Intensity(:,:,CorruptedTimeIndices) = [];

        end

    end

    % ====================
    % Fill Missing In Time
    % ====================

    function RawData = FillMissingInTime(RawData,Config)

        % Fills missing values (volumns in a time frame) using nearby columns by moving average

        MovingAverageWindow = Config.FillMissingAlongTimeMovingAverageWindow;

        if MovingAverageWindow > 0

            % Ocean velocity. Size is [D,T]
            RawData.Ocean.MeasuredVel_u = fillmissing(RawData.Ocean.MeasuredVel_u,'movmean',MovingAverageWindow,2);
            RawData.Ocean.MeasuredVel_v = fillmissing(RawData.Ocean.MeasuredVel_v,'movmean',MovingAverageWindow,2);
            RawData.Ocean.MeasuredVel_w = fillmissing(RawData.Ocean.MeasuredVel_w,'movmean',MovingAverageWindow,2);

            % Bottom Track. Size is [T]
            RawData.Ocean.BottomTrack_Vel_u = fillmissing(RawData.Ocean.BottomTrack_Vel_u,'movmean',MovingAverageWindow);
            RawData.Ocean.BottomTrack_Vel_v = fillmissing(RawData.Ocean.BottomTrack_Vel_v,'movmean',MovingAverageWindow);

        end

    end

    % ================
    % Process Raw Data
    % ================

    function Data = ProcessRawData(RawData,Config)

        % Data depths (meters, positive downward)
        Data.Ocean.Range = RawData.Ocean.Range;

        % Ocean range for each beam
        Data.Beam.Range = mean(RawData.Beam.Range,1);                  % Average ranges for 4 beams.
        Data.Beam.Range = fillmissing(Data.Beam.Range,'movmean',5);    % replace NaNs with nearby values by interpolation

        % Beam Correlation
        TempTotalCorrelationCounts = repmat(RawData.Beam.TotalCorrelationCounts,[1,1,size(RawData.Beam.CorrelationCounts,1)]); % add new dimension at end
        TotalCorrelationCounts = permute(TempTotalCorrelationCounts,[3,1,2]);  % Move last dimension to the begining
        Data.Beam.Correlation = RawData.Beam.CorrelationCounts ./ TotalCorrelationCounts;   % Normalize with total counts in each cell

        % Beam configuration
        DegreeToRadian = pi/180;
        Data.Beam.Angle = RawData.Beam.Angle * DegreeToRadian;               % Radian
        Data.Beam.Frequency = RawData.Beam.Frequency;                        % Hertz
        Data.Beam.CellSize = RawData.Beam.CellSize;                          % meters
        Data.Beam.ProfileMode = RawData.Beam.ProfileMode;                    % integer
        Data.Beam.LagLength = RawData.Beam.LagLength;                        % Number of bits in each period
        Data.Beam.PingsPerEnsemble = RawData.Beam.PingsPerEnsemble;          % 1 or 3

        % Data times in seconds from 0000 Jan 00.
        DayToSecond = 24.0 * 60.0 * 60.0;
        Data.Time = RawData.Time * DayToSecond;                              % In second

        % Average start and end longitude and latitude
        Data.Ship.Longitude= 0.5 * (RawData.Ship.StartLongitude + RawData.Ship.EndLongitude);
        Data.Ship.Latitude = 0.5 * (RawData.Ship.StartLatitude + RawData.Ship.EndLatitude);

        % Ship heading and madeg angles (counter-clockwise, zero eastward, radian)
        Data.Ship.Pitch = RawData.Ship.Pitch * DegreeToRadian;
        Data.Ship.Roll = RawData.Ship.Roll * DegreeToRadian;
        Data.Ship.Heading = RawData.Ship.Heading * DegreeToRadian;          % Zero is north, clockwise (same as compass)
        Data.Ship.Madegd = (90 - RawData.Ship.Madegd) * DegreeToRadian;     % Zero is shifted from north cw to east ccw.

        % Ship velocity (eastward and northward directions)
        Data.Ship.Vel_u = RawData.Ship.Speed .* cos(Data.Ship.Madegd) * Config.ShipVelocityUnitToMeterPerSecond;
        Data.Ship.Vel_v = RawData.Ship.Speed .* sin(Data.Ship.Madegd) * Config.ShipVelocityUnitToMeterPerSecond;

        % Adjust unit of velocity
        Ocean_MeasuredVel_u = RawData.Ocean.MeasuredVel_u * Config.OceanVelocityUnitToMeterPerSecond;
        Ocean_MeasuredVel_v = RawData.Ocean.MeasuredVel_v * Config.OceanVelocityUnitToMeterPerSecond;
        Ocean_MeasuredVel_w = RawData.Ocean.MeasuredVel_w * Config.OceanVelocityUnitToMeterPerSecond;

        % Ocean velocity: fill gaps vertically (replace missing row values in each columns from nearby row data)
        % Ocean_MeasuredVel_u = fillgapsOcean_MeasuredVel_u,20,8);
        % Ocean_MeasuredVel_v = fillgapsOcean_MeasuredVel_v,20,8);
        % Ocean_MeasuredVel_w = fillgapsOcean_MeasuredVel_w,20,8);

        % Ocean velocity relative to the ship (smooth along ranges)
        if Config.SmoothingWindowAlongVerticalRange > 1
            Ocean_MeasuredVel_u = smoothdata(Ocean_MeasuredVel_u,'gaussian',Config.SmoothingWindowAlongVerticalRange);
            Ocean_MeasuredVel_v = smoothdata(Ocean_MeasuredVel_v,'gaussian',Config.SmoothingWindowAlongVerticalRange);
            Ocean_MeasuredVel_w = smoothdata(Ocean_MeasuredVel_w,'gaussian',Config.SmoothingWindowAlongVerticalRange);
            % SmoothingKernel = ones(Config.SmoothingWindowAlongVerticalRange) / Config.SmoothingWindowAlongVerticalRange;
            % Ocean_MeasuredVel_u = conv2(Ocean_MeasuredVel_u,SmoothingKernel,'same');
            % Ocean_MeasuredVel_v = conv2(Ocean_MeasuredVel_v,SmoothingKernel,'same');
            % Ocean_MeasuredVel_w = conv2(Ocean_MeasuredVel_w,SmoothingKernel,'same');
        end

        % ---------------------------------------------------------

        % Ocean velocity on earth frame (only u and v, but not w)
        if Config.OceanVelocityIsRelativeToShip == true

            % Option 1: Use ship velocity that is obtained from madegood data (gps) (add to measured ocean velocity)
            % Data.Ocean.AbsoluteVel_u = Ocean_MeasuredVel_u + repmat(Data.Ship.Vel_u,size(Ocean_MeasuredVel_u,1),1) * Config.ShipVelocityScale;
            % Data.Ocean.AbsoluteVel_v = Ocean_MeasuredVel_v + repmat(Data.Ship.Vel_v,size(Ocean_MeasuredVel_u,1),1) * Config.ShipVelocityScale;

            % Option 2: Use bottom track velocity (subtract from measured ocean velocity)
            Data.Ocean.AbsoluteVel_u = Ocean_MeasuredVel_u - repmat(RawData.Ocean.BottomTrack_Vel_u,size(Ocean_MeasuredVel_u,1),1) * Config.ShipVelocityScale;
            Data.Ocean.AbsoluteVel_v = Ocean_MeasuredVel_v - repmat(RawData.Ocean.BottomTrack_Vel_v,size(Ocean_MeasuredVel_u,1),1) * Config.ShipVelocityScale;
        else
            Data.Ocean.AbsoluteVel_u = Ocean_MeasuredVel_u;
            Data.Ocean.AbsoluteVel_v = Ocean_MeasuredVel_v;
        end

        % ---------------------------------------------------------
        
        % % Use projection to ship direction and normal to ship direction
        % OceanVel_ParallelToShip = zeros(size(Data.Time));
        % OceanVel_TransversalToShip = zeros(size(Data.Time));

        % Data.Ocean.AbsoluteVel_u = zeros(size(Ocean_MeasuredVel_u));
        % Data.Ocean.AbsoluteVel_v = zeros(size(Ocean_MeasuredVel_v));

        % for t = 1:size(Data.Time,2)

        %     ShipVelocityMagnitude = hypot(Data.Ship.Vel_u(t),Data.Ship.Vel_v(t));

        %     % Ship 
        %     ShipVel = [Data.Ship.Vel_u(t),Data.Ship.Vel_v(t)];
        %     ShipUnitDirection = ShipVel / hypot(ShipVel(1),ShipVel(2));
        %     ShipTransverseUnitDirection = [-ShipUnitDirection(2),ShipUnitDirection(1)];

        %     % Ocean velocity along ship direction
        %     OceanVels_ParallelToShip = Ocean_MeasuredVel_u(:,t)*ShipUnitDirection(1) + Ocean_MeasuredVel_v(:,t)*ShipUnitDirection(2);
        %     OceanVels_TransversalToShip = Ocean_MeasuredVel_u(:,t)*ShipTransverseUnitDirection(1) + Ocean_MeasuredVel_v(:,t)*ShipTransverseUnitDirection(2);

        %     % Add ship refrerence frame
        %     OceanVels_ParallelToShip = OceanVels_ParallelToShip + ShipVelocityMagnitude * Config.ShipVelocityScale;

        %     % Return back to x-y coordinates
        %     Data.Ocean.AbsoluteVel_u(:,t) = OceanVels_ParallelToShip * ShipUnitDirection(1) + OceanVels_TransversalToShip * ShipTransverseUnitDirection(1);
        %     Data.Ocean.AbsoluteVel_v(:,t) = OceanVels_ParallelToShip * ShipUnitDirection(2) + OceanVels_TransversalToShip * ShipTransverseUnitDirection(2);

        % end

        % -----------------------------------------------------------

        % w velocity is already in earth coordinate
        Data.Ocean.AbsoluteVel_w = Ocean_MeasuredVel_w;

        % Intensity
        Data.Beam.Intensity = squeeze(mean(RawData.Beam.Intensity,2));  % average over all beam measurements
        % Data.Beam.Intensity = fillgaps(Data.Beam.Intensity,20,8);
        Data.Beam.Intensity = smoothdata(Data.Beam.Intensity,'gaussian',Config.SmoothingWindowAlongVerticalRange);
        % Data.Beam.Intensity = conv2(Data.Beam.Intensity,SmoothingKernel,'same');

        % Estimate Bathymetry from intensity
        % % Method 1: Find Bathymetry from max of intensity
        % % IgnoreDepthIndex = 1;    % Exclude top depths from searching for spike of intensity
        % % [~,MaxIntensityIndex] = max(Data.Beam.Intensity(IgnoreDepthIndex:end,:));
        % % Data.Ocean.BathymetryIndex = MaxIntensityIndex + IgnoreDepthIndex;
        % 
        % % Method 2: Find Bathymetry from max of prominence-to-width ratio of spike peaks of intensity
        % Data.Ocean.BathymetryIndex = zeros(size(Data.Beam.Intensity,2),1);
        % for TimeIndex = 1:size(Data.Beam.Intensity,2)
        %     [PeakHeight,PeakLocation,PeakWidth,PeakProminence] = findpeaks(Data.Beam.Intensity(:,TimeIndex),'MinPeakProminence',10,'MinPeakDistance',5);
        %     PeakShape = PeakProminence ./ PeakWidth;   % larger ratio is sharper peak
        % 
        %     % Method 2: Select peak from the height
        %     % [~,SelectedPeakIndex] = max(PeakHeight);
        % 
        %     % Method 3: Select peak by max PeakShape
        %     % [~,SelectedPeakIndex] = max(PeakShape);
        % 
        %     % Method 4: Select peak by height, from the first two max peak shape
        %     if length(PeakShape) == 1
        %         SelectedPeakIndex = 1;
        %     else
        %         % Choose the peak with highest height from the first two Peakshapes
        %         [SortedPeakShape,SortingIndex] = sort(PeakShape,'descend');
        %         if PeakHeight(SortingIndex(1)) > PeakHeight(SortingIndex(2))
        %         % if PeakLocation(SortingIndex(1)) < PeakLocation(SortingIndex(2))
        %             SelectedPeakIndex = SortingIndex(1);
        %         else
        %             SelectedPeakIndex = SortingIndex(2);
        %         end
        %     end
        % 
        %     Data.Ocean.BathymetryIndex(TimeIndex) = PeakLocation(SelectedPeakIndex);
        % end

        % Find Bathymetry index in the range
        Data.Ocean.BathymetryIndex = zeros(size(Data.Beam.Intensity,2),1);
        for TimeIndex = 1:size(Data.Beam.Intensity,2)

            % Compare BeamRange with cell range intervals to find index
            for RangeIndex = 1:length(Data.Ocean.Range)-1
                if (Data.Beam.Range(TimeIndex) >= Data.Ocean.Range(RangeIndex)) & (Data.Beam.Range(TimeIndex) < Data.Ocean.Range(RangeIndex+1))
                    Data.Ocean.BathymetryIndex(TimeIndex) = RangeIndex;
                    break
                end
            end

            % Check
            if Data.Ocean.BathymetryIndex(TimeIndex) == 0
                error('BathymetryIndex is zero.')
            end
        end

        % For all variables with depth argument, make below bathymetry to be NaN.
        if Config.DiscardOceanVelocitiesBelowIntensitySpike == true
            % Find where the intensity has a spike (to find ground)
            for TimeIndex = 1:size(Data.Beam.Intensity,2)
                
                % Reomve below bathymetry
                Data.Ocean.AbsoluteVel_u(Data.Ocean.BathymetryIndex(TimeIndex):end,TimeIndex) = NaN;
                Data.Ocean.AbsoluteVel_v(Data.Ocean.BathymetryIndex(TimeIndex):end,TimeIndex) = NaN;
                Data.Beam.Correlation(Data.Ocean.BathymetryIndex(TimeIndex):end,:,TimeIndex) = NaN;
            end
        end

    end

    % ---------------------
    % End of Static methods

    end
end
