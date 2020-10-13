classdef ReadWind

    methods(Static)

    % ==============
    % Read Wind File
    % ==============

    function WindField = ReadWindFile(Config)

    % Note: Wind data time is in seconds since 1970. We convert then to be since year 0000.
    % Also longitudes are from 0 to 360. We convert them to be from -180 to 180.

    % Open file object
    ncid = netcdf.open(Config.Wind.WindFilename,'NC_NOWRITE');
    
    % Read variables
    LongitudesGrid = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'longitude'));
    LatitudesGrid = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latitude'));
    Time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
    EastWind = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'UGRD_10maboveground'));
    NorthWind = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'VGRD_10maboveground'));

    % Adjust longitude to be from -180 to +180
    LongitudesGrid = LongitudesGrid - 360;

    % Adjust time
    OriginTimeInDays = datenum(Config.Wind.WindDataOriginTime);
    OriginTimeInSeconds = OriginTimeInDays * (24*3600);
    Time = Time + OriginTimeInSeconds;   % Still in seconds unit.

    % Store in struct
    WindField.LongitudesGrid = LongitudesGrid;
    WindField.LatitudesGrid = LatitudesGrid;
    WindField.Time = Time;
    WindField.EastWind = EastWind;
    WindField.NorthWind = NorthWind;

    % Close file object
    netcdf.close(ncid)

    % Plot
    if Config.Plots.ExtraPlots == true
        ReadWind.PlotWindData(WindField,Config)
    end

    end

    % ===============
    % Plot Wind Field
    % ===============

    function PlotWindData(WindField,Config)

    figure()

    % Set the TimeIndex below to plot
    % Each time index correspond to one hour in nam wind data. For instance, in nam.conusnest.hiresf.20180816.t00z, 
    % the data starts at 2018/08/16 at 00 UTC time. To plot wind at 2018/08/16 at 17:00 UTC, set TimeIndex to 18.
    TimeIndex = 18;

    quiver(WindField.LongitudesGrid,WindField.LatitudesGrid,WindField.EastWind(:,:,TimeIndex),WindField.NorthWind(:,:,TimeIndex),'color','black')
    xlabel('Longitude (degrees)')
    ylabel('Latitude (degrees)')
    title(sprintf('Wind field on %s',datestr(WindField.Time(TimeIndex)/(24*3600))))

    % Font
    set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

    end

    % ===============================
    % Interpolate Wind On Ship Tracks
    % ===============================

    function WindOnShipTrack = InterpolateWindOnShipTracks(Config,Data)
    
    % Check wind data time is within the velocity data time
    if (Data.Time(1) < Data.Wind.Field.Time(1)) || (Data.Time(end) > Data.Wind.Field.Time(end))
    SecondToDay = 1/(24*3600);
        fprintf('Velocity data time from %s to %s.\n',datestr(Data.Time(1)*SecondToDay),datestr(Data.Time(end)*SecondToDay))
        fprintf('Wind data time from %s to %s.\n',datestr(Data.Wind.Field.Time(1)*SecondToDay),datestr(Data.Wind.Field.Time(end)*SecondToDay))
        error('Velocity time is not within the wind data time.')
    end

    % Wind grid points
    LongitudePoints = Data.Wind.Field.LongitudesGrid(:);
    LatitudePoints = Data.Wind.Field.LatitudesGrid(:);

    % Initialize arrays
    NumTimes = length(Data.Ship.Longitude);
    WindOnShipTrack.Time = zeros(1,NumTimes);
    WindOnShipTrack.Longitude = zeros(1,NumTimes);
    WindOnShipTrack.Latitude = zeros(1,NumTimes);
    WindOnShipTrack.EastWind = zeros(1,NumTimes);
    WindOnShipTrack.NorthWind = zeros(1,NumTimes);

    % Iterate over points
    for i = 1:NumTimes

        % Get time index in wind time data
        TimeIndexInWindData = find(Data.Wind.Field.Time < Data.Time(i),1,'last');
        if isempty(TimeIndexInWindData)
            error('Cannot find a time index from wind time data.')
        end

        % Time interpolating coefficient
        TimeInterpolatingCoefficient = (Data.Time(i) - Data.Wind.Field.Time(TimeIndexInWindData)) / (Data.Wind.Field.Time(TimeIndexInWindData+1) - Data.Wind.Field.Time(TimeIndexInWindData));

        % Inquiry points (on the ship)
        InquiryLongitude = Data.Ship.Longitude(i);
        InquiryLatitude = Data.Ship.Latitude(i);

        % Data fields
        EastWindFieldAtLeftTime = double(Data.Wind.Field.EastWind(:,:,TimeIndexInWindData));
        NorthWindFieldAtLeftTime = double(Data.Wind.Field.NorthWind(:,:,TimeIndexInWindData));
        EastWindFieldAtRightTime = double(Data.Wind.Field.EastWind(:,:,TimeIndexInWindData+1));
        NorthWindFieldAtRightTime = double(Data.Wind.Field.NorthWind(:,:,TimeIndexInWindData+1));

        % Interpolate in space at left time
        EastWindAtLeftTime = griddata(LongitudePoints,LatitudePoints,EastWindFieldAtLeftTime(:),InquiryLongitude,InquiryLatitude);
        NorthWindAtLeftTime = griddata(LongitudePoints,LatitudePoints,NorthWindFieldAtLeftTime(:),InquiryLongitude,InquiryLatitude);

        % Interpolate in space in right time
        EastWindAtRightTime = griddata(LongitudePoints,LatitudePoints,EastWindFieldAtRightTime(:),InquiryLongitude,InquiryLatitude);
        NorthWindAtRightTime = griddata(LongitudePoints,LatitudePoints,NorthWindFieldAtRightTime(:),InquiryLongitude,InquiryLatitude);

        % Interpolate in time
        WindOnShipTrack.Time(i) = Data.Time(i);
        WindOnShipTrack.Longitude(i) = InquiryLongitude;
        WindOnShipTrack.Latitude(i) = InquiryLatitude;
        WindOnShipTrack.EastWind(i) = (1-TimeInterpolatingCoefficient) * EastWindAtLeftTime + TimeInterpolatingCoefficient * EastWindAtRightTime;
        WindOnShipTrack.NorthWind(i) = (1-TimeInterpolatingCoefficient) * NorthWindAtLeftTime + TimeInterpolatingCoefficient * NorthWindAtRightTime;
    end

    % Plot
    if Config.Plots.ExtraPlots == true
        ReadWind.PlotWindDataAlongShipTracks(WindOnShipTrack,Config)
    end

    end

    % =================================
    % Plot Wind Field Along Ship Tracks
    % =================================

    function PlotWindDataAlongShipTracks(WindOnShipTrack,Config)

    figure()

    TimeIndex = 18;
    quiver(WindOnShipTrack.Longitude,WindOnShipTrack.Latitude,WindOnShipTrack.EastWind,WindOnShipTrack.NorthWind,'color','black','DisplayName','Wind')
    hold on
    plot(WindOnShipTrack.Longitude,WindOnShipTrack.Latitude,'-o','color',[0.5,0.5,0.5],'DisplayName','Ship track')
    xlabel('Longitude (degrees)')
    ylabel('Latitude (degrees)')
    legend()
    title(sprintf('Wind field on %s to %s',datestr(WindOnShipTrack.Time(1)/(24*3600)),datestr(WindOnShipTrack.Time(end)/(24*3600))))

    % Font
    set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

    end

    % ---------------------
    % End of static methods

    end
end
