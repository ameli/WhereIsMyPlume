% ----------------------------------------------------------------
%          filename: Drifters.m
%           purpose: Extracts drifter data from raw xlsx files.
%            author: Siavash Ameli
%              date: 2019/02/27
%           project: NSF ALPHA, full study experiment, August 2018
%   acknowledgement: This work is supported by NSF No. 1520825
% -----------------------------------------------------------------

% Find the project directry. We assume this script is in ProjectDirectory/utilities/
[CurrentDirectory,~,~,] = fileparts(mfilename('fullpath'));
[ProjectDirectory,~,~] = fileparts(CurrentDirectory);

% Directory of all drifter xlsx files
DataDirectory = fullfile(ProjectDirectory,'data','drifters','xlsx');
Filenames = '*.xlsx';

% Drifter Ids
DrifterIdsFilename = fullfile(ProjectDirectory,'data','drifters','DrifterIds.txt');
DrifterIdsTable = readtable(DrifterIdsFilename);

% Get list of all files
FullPathFilenames = fullfile(DataDirectory,Filenames);
ListOfFilesStruct = dir(FullPathFilenames);

% On the xlsx document, the first row where data start
FirstValidRow = 7;

% Initialize output struct
Drifters = [];

% Iterate over each file, read data into array of structs, 'Drifters'.
for FileId = 1:size(ListOfFilesStruct,1)

    fprintf('%2d/%d Reading: %s ',FileId,size(ListOfFilesStruct,1),ListOfFilesStruct(FileId).name);
    FullPathFilename = fullfile(DataDirectory,ListOfFilesStruct(FileId).name);

    % Read xlsx file
    [Numerics,Texts,RawData] = xlsread(FullPathFilename);

    % Initialize arrays
    DrifterName = '';
    ValidRows = boolean(ones(size(RawData,1),1));   % zero is invalid, one is valid
    Longitude = zeros(size(RawData,1),1);
    Latitude = zeros(size(RawData,1),1);
    Time = zeros(size(RawData,1),1);
    Date = cell(size(RawData,1),1);

    % Columns in xlsx file corresponding to each variable
    NameColumn = 1;
    LatitudeColumn = 2;
    LongitudeColumn = 3;
    DateColumn = 4;

    % Iterate over rows
    for Row = 1:size(RawData,1)

        if mod(Row,floor(size(RawData,1)/10)) == 0
            fprintf('.')
        end
        
        % Read Latitude
        if isnumeric(RawData{Row,LatitudeColumn}) & ~isnan(RawData{Row,LatitudeColumn})
            if RawData{Row,LatitudeColumn} >= -90.0 & RawData{Row,LatitudeColumn} <= 90.0
                Latitude(Row) = RawData{Row,LatitudeColumn};
                DrifterName = RawData{Row,NameColumn};
            else
                ValidRows(Row) = false;
                continue
            end
        else
            ValidRows(Row) = false;
            continue
        end

        % Read Longitude
        if isnumeric(RawData{Row,LongitudeColumn}) & ~isnan(RawData{Row,LongitudeColumn})
            if RawData{Row,LongitudeColumn} >= -360.0 & RawData{Row,LongitudeColumn} <= 360.0
                Longitude(Row) = RawData{Row,LongitudeColumn};
            else
                ValidRows(Row) = false;
                continue
            end
        else
            ValidRows(Row) = false;
            continue
        end

        % Read Date
        if ischar(RawData{Row,DateColumn})
            Date{Row} = RawData{Row,DateColumn};
            Time(Row) = datenum(RawData{Row,DateColumn},'dd-mm-yyyy HH:MM:SS');
        else
            ValidRows(Row) = false;
            continue
        end
    end

    % Find drifter Id (integrers, 1,2,...,48) from drifter name (string, like '0-3103621')
    DrifterId = NaN;
    for TableRow = 1:size(DrifterIdsTable,1)
        if strcmp(DrifterName,DrifterIdsTable.Var2(TableRow))
            DrifterId = DrifterIdsTable.Var1(TableRow);
            break
        end
    end

    % Filter arrays by valid entries
    Latitude = Latitude(ValidRows);
    Longitude = Longitude(ValidRows);
    Time = Time(ValidRows);
    Date = Date(ValidRows);

    % Remove time redundancy
    [Time,UniqueIndex] = unique(Time);
    Latitude = Latitude(UniqueIndex);
    Longitude = Longitude(UniqueIndex);
    Date = Date(UniqueIndex);

    % Reverse order since they are in descending time order
    Time = Time(end:-1:1);
    Longitude = Longitude(end:-1:1);
    Latitude = Latitude(end:-1:1);
    Date = Date(end:-1:1);

    % Sort again time in chronologically increasing order
    [Time,SortingIndex] = sort(Time,'ascend');
    Longitude = Longitude(SortingIndex);
    Latitude = Latitude(SortingIndex);
    Date = Date(SortingIndex);

    % Remove outliers, sudden jumpts on the drifter trajectory, SETTINGS
    RemoveOutlier = true;                         % perform removing (or removing and filling) outliers
    FillOutlier = true;                           % instead of removing, replace it with an interpolation value
    Window = 11;                                  % move mean window
    ThresholdFactor = 1;                          % factor multiplied by standard deviation in the local window
    DetectOutlierMethod = 'movmedian';            % 'movmean' or 'movmedian'
    FillOutlierMethod = 'pchip';                  % 'linear', 'spline', 'pchip', etc     

    if RemoveOutlier == true

        % Remove outler
        X = [Longitude(:),Latitude(:)];
        [X_new,RemovedIndices] = rmoutliers(X,DetectOutlierMethod,Window,'ThresholdFactor',ThresholdFactor);

        if FillOutlier == true

            % Fill the removed outliers. In this case, the number of original data (X) and output data are the same since the removed data are filled
            X_copy = X;
            for Dim = 1:2
                X_copy(RemovedIndices,Dim) = interp1(Time(~RemovedIndices),X(~RemovedIndices,Dim),Time(RemovedIndices));
            end

            Longitude = X_copy(:,1);
            Latitude = X_copy(:,2);


        else

            % Remove outliers (do not fill them). In this case the number of output data are less since some data are removed
            Longitude = X_new(:,1);
            Latitude = X_new(:,2);
            Time(RemovedIndices) = [];
            Date(RemovedIndices) = [];
        end


        % Check lengths for consistency
        if length(Longitude) ~= length(Date)
            error('Lengths are not the same')
        end
    end

    % Save in Drifters struct in reverse order to be chronologically increasing
    Drifters(FileId).Name = DrifterName;
    Drifters(FileId).Id = DrifterId;
    Drifters(FileId).Longitude = Longitude;
    Drifters(FileId).Latitude = Latitude;
    Drifters(FileId).Time = Time;
    Drifters(FileId).Date = Date;
    Drifters(FileId).LongitudeUnit = 'Degrees, positive eastward';
    Drifters(FileId).LatitudeUnit = 'Degrees, positive northward';
    Drifters(FileId).TimeUnit = 'Days since Jan 0, 0000 at 00:00:00 (UTC)';
    Drifters(FileId).DateUnit = 'UTC';

    fprintf(' Done.\n')
end

% Save to mat file
OutputDirectory = fullfile(ProjectDirectory,'data','drifters');
OutputFilename = 'Drifters.mat';
FullPathOutputFilename = fullfile(OutputDirectory,OutputFilename);
save(FullPathOutputFilename,'Drifters')
fprintf('Wrote to: %s.\n',FullPathOutputFilename)

% Plots
figure()
for DrifterId = 1:size(Drifters,2)
    plot(Drifters(DrifterId).Longitude,Drifters(DrifterId).Latitude)
    hold on
end

figure()
for DrifterId = 1:size(Drifters,2)
    plot(Drifters(DrifterId).Time)
    hold on
end

