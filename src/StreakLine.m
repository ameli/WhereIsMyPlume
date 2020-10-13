classdef StreakLine

    methods(Static)

    % ==============================
    % Trace Initial Plume StreakLine
    % ==============================

    function Output = TracePlumeInitialStreakLine(Data,Config)

    disp('Computing streakline ...');

    % Convert datetimes to times (second from 0000 Jan 00)
    DayToSecond = 24.0 * 60.0 * 60.0;
    PlumeStartReleaseTime = datenum(Config.PlumeRelease.StartReleaseTime) * DayToSecond;
    PlumeEndReleaseTime = datenum(Config.PlumeRelease.EndReleaseTime) * DayToSecond;

    % Inquiry time
    if Config.PlumeInquiry.UseLastAvailableTime == true
        InquiryTime = datenum(Data.Time(end)) * DayToSecond;
    else
        InquiryTime = datenum(Config.PlumeInquiry.InquiryTime) * DayToSecond;
    end


    % Check release times with data times
    if PlumeStartReleaseTime < Data.Time(1)
        error('ERROR: File time starts at: %s, but plume start release time is %s', ...
            datestr(Data.Time(1)/DayToSecond),datestr(PlumeStartReleaseTime/DayToSecond));
    end

    % Check release times with data times
    if PlumeEndReleaseTime > Data.Time(end)
        error('ERROR: File time ends at: %s, but plume end release time is %s', ...
            datestr(Data.Time(end)/DayToSecond),datestr(PlumeEndReleaseTime/DayToSecond));
    end

    % Check Inquiry time with Start Release Time
    if InquiryTime < PlumeStartReleaseTime
        error('Inquiry time: %s cannot be earlier than plume start release time: %s', ...
        datestr(InquiryTime/DayToSecond),datestr(PlumeStartReleaseTime/DayToSecond));
    end

    % Terminate streakline earlier of inquity time is earlier than end release time
    if InquiryTime < PlumeEndReleaseTime

        sprintf('NOTE: Streakline will be terminated at %s, which is earlier than end of plume release time at %s.', ...
        datestr(InquiryTime/DayToSecond),datestr(PlumeEndReleaseTime/DayToSecond))
        FinalStreakLineTime = InquiryTime;
        PlumeReleaseEndedAtInquiryTime = false;

    else

        FinalStreakLineTime = PlumeEndReleaseTime;
        PlumeReleaseEndedAtInquiryTime = true;

    end

    % Check Depth
    % if Config.PlumeRelease.ReleaseDepth < Data.Ocean.Range(1) || Config.PlumeRelease.ReleaseDepth > Data.Ocean.Range(end)
    %     error('Plume depth: %f is not in the data range from %f(m) to %f(m).', ...
    %     Config.PlumeRelease.ReleaseDepth,Data.Ocean.Range(1),Data.Ocean.Range(end))
    % end

    % Convert Ship Lat/Lon to XY, origin at initial position of ship
    OriginLongitude = Data.Ship.Longitude(1);
    OriginLatitude = Data.Ship.Latitude(1);

    % Interpolate Velocity At Depth
    Ocean_VelAtDepth_u = interp1(Data.Ocean.Range,Data.Ocean.AbsoluteVel_u,Config.PlumeRelease.ReleaseDepth,'pchip','extrap');
    Ocean_VelAtDepth_v = interp1(Data.Ocean.Range,Data.Ocean.AbsoluteVel_v,Config.PlumeRelease.ReleaseDepth,'pchip','extrap');

    % Include wind to the ocean surface data
    if Config.Wind.IncludeWindData == true
        Ocean_VelAtDepth_u = Ocean_VelAtDepth_u + Config.Wind.WindageCoefficient * Data.Wind.OnShipTrack.EastWind;
        Ocean_VelAtDepth_v = Ocean_VelAtDepth_v + Config.Wind.WindageCoefficient * Data.Wind.OnShipTrack.NorthWind;
    end

    % Function handle to interpolate velocity in time
    OceanVelocity =@(t,x) [interp1(Data.Time,Ocean_VelAtDepth_u,t);interp1(Data.Time,Ocean_VelAtDepth_v,t)];

    NumberOfStreakLineParticles = Config.NumberOfStreakLineParticles;
    PlumeInjectionTime = linspace(PlumeStartReleaseTime,FinalStreakLineTime,NumberOfStreakLineParticles);
    PlumeInjectionLongitude = zeros(1,NumberOfStreakLineParticles);
    PlumeInjectionLatitude = zeros(1,NumberOfStreakLineParticles);
    PlumeStreakLineX = zeros(1,NumberOfStreakLineParticles);
    PlumeStreakLineY = zeros(1,NumberOfStreakLineParticles);

    for i = 1:size(PlumeInjectionTime,2)
        
        % Initial position of trajectory (lat/lon)
        PlumeInjectionLongitude(i) = interp1(Data.Time,Data.Ship.Longitude,PlumeInjectionTime(i),'pchip');
        PlumeInjectionLatitude(i) = interp1(Data.Time,Data.Ship.Latitude,PlumeInjectionTime(i),'pchip');
        
        % Initial Position of traejctoty (X,Y)
        [PlumeInjectionX,PlumeInjectionY] = Utilities.ConvertLonLatToXY(OriginLongitude,OriginLatitude,PlumeInjectionLongitude(i),PlumeInjectionLatitude(i));
            
        if PlumeInjectionTime(i) < FinalStreakLineTime
            
            % Integrate trajectory
            opts = odeset('RelTol',1e-7);
            [~,Trajectory_Coordinates] = ...
                ode45(OceanVelocity,[PlumeInjectionTime(i),FinalStreakLineTime],[PlumeInjectionX,PlumeInjectionY],opts);

            PlumeStreakLineX(i) = Trajectory_Coordinates(end,1);
            PlumeStreakLineY(i) = Trajectory_Coordinates(end,2);
            
        elseif abs(PlumeInjectionTime(i) - FinalStreakLineTime) < 10*eps
            
            % Plume is just injected. Nothing to integrate.
            PlumeStreakLineX(i) = PlumeInjectionX;
            PlumeStreakLineY(i) = PlumeInjectionY;
            
        else
            error('ERORR: PlumeInjectionTime is after PlumeReleaseTime')
        end
    end

    % Convert back XY to LonLat
    [PlumeStreakLineLongitude,PlumeStreakLineLatitude] = ...
        Utilities.ConvertXYToLonLat(OriginLongitude,OriginLatitude,PlumeStreakLineX,PlumeStreakLineY);
        
    % Outputs
    Output.PlumeInitialStreakLine.LongitudeAtStartReleaseTime = PlumeInjectionLongitude;
    Output.PlumeInitialStreakLine.LatitudeAtStartReleaseTime = PlumeInjectionLatitude;

    Output.PlumeInitialStreakLine.LongitudeAtEndReleaseTime = PlumeStreakLineLongitude;
    Output.PlumeInitialStreakLine.LatitudeAtEndReleaseTime = PlumeStreakLineLatitude;

    % Output.PlumeInitialStreakLine.CenterLongitudeAtEndReleaseTime = PlumeStreakLineLongitude(floor(size(PlumeStreakLineLongitude,2)/2));    % Use center point on streakline
    % Output.PlumeInitialStreakLine.CenterLatitudeAtEndReleaseTime = PlumeStreakLineLatitude(floor(size(PlumeStreakLineLatitude,2)/2));
    Output.PlumeInitialStreakLine.CenterLongitudeAtEndReleaseTime = mean(PlumeStreakLineLongitude);                                           % Use geometric center of streakline curve
    Output.PlumeInitialStreakLine.CenterLatitudeAtEndReleaseTime = mean(PlumeStreakLineLatitude);

    Output.PlumeInitialStreakLine.StartReleaseTime = PlumeStartReleaseTime;     % Time plume release starts
    Output.PlumeInitialStreakLine.EndReleaseTime = PlumeEndReleaseTime;         % Time plume release ends
    Output.PlumeInitialStreakLine.FinalStreakLineTime = FinalStreakLineTime;    % Time we ended integration
    Output.PlumeInitialStreakLine.PlumeReleaseEnded = PlumeReleaseEndedAtInquiryTime;

    end

    % =============================
    % Plot Plume Initial StreakLine
    % =============================

    function PlotPlumeInitialStreakLine(PlumeInitialStreakLine,Data,Config)

    figure()
    oceanColor = [.6 .8 .9];
    ax = axesm('MapProjection','mercator');
    setm(ax,'FFaceColor',oceanColor);

    % Ship total trajectory
    h1 = plotm(Data.Ship.Latitude,Data.Ship.Longitude,'-o','color','black','DisplayName','Ship trajectory');
    hold on

    % Ship trajectory during plume streakline
    h2 = plotm(PlumeInitialStreakLine.LatitudeAtStartReleaseTime,PlumeInitialStreakLine.LongitudeAtStartReleaseTime,'color','green','linewidth',2,'DisplayName','Ship trajectory during plume release');
    plotm(PlumeInitialStreakLine.LatitudeAtStartReleaseTime(1),PlumeInitialStreakLine.LongitudeAtStartReleaseTime(1),'o','MarkerEdgeColor','blue','MarkerFaceColor','green','MarkerSize',8)
    plotm(PlumeInitialStreakLine.LatitudeAtStartReleaseTime(end),PlumeInitialStreakLine.LongitudeAtStartReleaseTime(end),'o','MarkerEdgeColor','blue','MarkerFaceColor','green','MarkerSize',8)
    textm(PlumeInitialStreakLine.LatitudeAtStartReleaseTime(1),PlumeInitialStreakLine.LongitudeAtStartReleaseTime(1),datestr(PlumeInitialStreakLine.StartReleaseTime/(24*3600)))
    textm(PlumeInitialStreakLine.LatitudeAtStartReleaseTime(end),PlumeInitialStreakLine.LongitudeAtStartReleaseTime(end),datestr(PlumeInitialStreakLine.FinalStreakLineTime/(24*3600)))

    % Plume Streakline
    h3 = plotm(PlumeInitialStreakLine.LatitudeAtEndReleaseTime,PlumeInitialStreakLine.LongitudeAtEndReleaseTime,'color','red','linewidth',2,'DisplayName','Plume streakline during release');
    plotm(PlumeInitialStreakLine.LatitudeAtEndReleaseTime(1),PlumeInitialStreakLine.LongitudeAtEndReleaseTime(1),'o','MarkerEdgeColor','red','MarkerFaceColor','red','MarkerSize',8)
    plotm(PlumeInitialStreakLine.LatitudeAtEndReleaseTime(end),PlumeInitialStreakLine.LongitudeAtEndReleaseTime(end),'o','MarkerEdgeColor','red','MarkerFaceColor','red','MarkerSize',8)
    textm(PlumeInitialStreakLine.LatitudeAtEndReleaseTime(1),PlumeInitialStreakLine.LongitudeAtEndReleaseTime(1),datestr(PlumeInitialStreakLine.StartReleaseTime/(24*3600)));

    MinLon = min([Data.Ship.Longitude(:);PlumeInitialStreakLine.LongitudeAtEndReleaseTime(:);PlumeInitialStreakLine.LongitudeAtStartReleaseTime(:)]);
    MaxLon = max([Data.Ship.Longitude(:);PlumeInitialStreakLine.LongitudeAtEndReleaseTime(:);PlumeInitialStreakLine.LongitudeAtStartReleaseTime(:)]);
    MinLat = min([Data.Ship.Latitude(:);PlumeInitialStreakLine.LatitudeAtEndReleaseTime(:);PlumeInitialStreakLine.LatitudeAtStartReleaseTime(:)]);
    MaxLat = max([Data.Ship.Latitude(:);PlumeInitialStreakLine.LatitudeAtEndReleaseTime(:);PlumeInitialStreakLine.LatitudeAtStartReleaseTime(:)]);

    LatExt = MaxLat - MinLat;
    LonExt = MaxLon - MinLon;
    Ext = max(LatExt,LonExt);
    Ratio = 0.1;

    LatLim = [MinLat-Ratio*Ext,MaxLat+Ratio*Ext];
    LonLim = [MinLon-Ratio*Ext,MaxLon+Ratio*Ext];

    setm(gca,'MapLatLimit',LatLim,'MapLonLimit',LonLim)

    setm(gca,'MLineLocation',0.005,'MLabelLocation',0.005,'MLabelRound',-3,'MeridianLabel','on','MLabelParallel','south');
    setm(gca,'PLineLocation',0.005,'PLabelLocation',0.005,'PLabelRound',-3,'ParallelLabel','on','PLabelMeridian','west');

    tightmap; %Tighten up the map area shown
    framem on; %Turn on the black frame
    gridm on; %Turn on grid lines
    mlabel on;
    plabel on;
    showaxes;
    grid off;
    axis off;

    legend([h1,h2,h3])

    % xlim([MinLon,MaxLon])
    % ylim([MinLat,MaxLat])

    xlabel('Longitude(deg)')
    ylabel('Latitude (deg)')
    title('Plume initial streakline')

    % Font
    set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

    end

    % ---------------------
    % End of static methods

    end
end
