classdef PlotData

    methods(Static)

    % ====================
    % Plot Measurment Data
    % ====================

    function PlotMeasurementData(Data,Config)

    % Plots all concatenated data from the begining to the end of times.

    % Find max bathymetry
    NonZero = find(all(isnan(Data.Ocean.AbsoluteVel_u),2));
    if isempty(NonZero)
        MaxBathymetry = size(Data.Ocean.AbsoluteVel_u,1);
    else
        MaxBathymetry = NonZero(1);
    end

    SecondToDay = 1.0 / (24.0 * 3600.0);
    Datetime = Data.Time * SecondToDay;

    figure()

    % Plot east velocity
    subplot(2,3,1)
    pcolor(Datetime,Data.Ocean.Range(1:MaxBathymetry,1),Data.Ocean.AbsoluteVel_u(1:MaxBathymetry,:))
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    shading interp;
    h = colorbar;
    set(gca,'Ydir','reverse');
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Ranges (m)');
    ylabel(h,'(m/s)');
    title('East velocity');
    set(gca,'XGrid','on')
    set(gca,'color',[0.8,0.52,0.24])

    % Plot north velocity
    subplot(2,3,4)
    pcolor(Datetime,Data.Ocean.Range(1:MaxBathymetry,1),Data.Ocean.AbsoluteVel_v(1:MaxBathymetry,:))
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    shading interp;
    h = colorbar;
    set(gca,'Ydir','reverse');
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Ranges (m)');
    ylabel(h,'(m/s)');
    title('North velocity');
    set(gca,'XGrid','on')
    set(gca,'color',[0.8,0.52,0.24])

    % Plot Velocity magnitude
    subplot(2,3,2)
    ComplexVelocity = Data.Ocean.AbsoluteVel_u + i * Data.Ocean.AbsoluteVel_v;
    VelocityMagnitude = abs(ComplexVelocity);
    pcolor(Datetime,Data.Ocean.Range(1:MaxBathymetry,1),VelocityMagnitude(1:MaxBathymetry,:))
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    shading interp;
    h = colorbar;
    set(gca,'Ydir','reverse');
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Ranges (m)');
    ylabel(h,'(m/s)');
    title('Velocity magnitude');
    set(gca,'XGrid','on')
    set(gca,'color',[0.8,0.52,0.24])

    % Velocity direction
    subplot(2,3,5)
    VelocityDirection = mod(angle(ComplexVelocity),2*pi) * 180.0 / pi;
    pcolor(Datetime,Data.Ocean.Range(1:MaxBathymetry,1),VelocityDirection(1:MaxBathymetry,:))
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    shading interp;
    h = colorbar;
    h.Ticks = linspace(0,360,5);
    set(gca,'Ydir','reverse');
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Ranges (m)');
    ylabel(h,'(degrees)');
    title('Velocity direction');
    set(gca,'XGrid','on')
    set(gca,'color',[0.8,0.52,0.24])

    % Plot intensity
    subplot(2,3,3)
    pcolor(Datetime,Data.Ocean.Range,Data.Beam.Intensity)
    hold on
    plot(Datetime,Data.Beam.Range,'color','black','linewidth',1,'DisplayName','BeamRange')
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    shading interp;
    h = colorbar;
    set(gca,'Ydir','reverse');
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Ranges (m)');
    title('Intensity');
    set(gca,'XGrid','on')

    % Plot ship trajectory
    subplot(2,3,6)
    N = size(Data.Ship.Longitude,2)-1;
    Speed = hypot(Data.Ship.Vel_u,Data.Ship.Vel_v);;
    ScaledSpeed = floor(N * Speed / max(Speed));
    ScaledSpeed(ScaledSpeed < 1) = 1;
    ScaledSpeed(ScaledSpeed > N) = N;
    cmap = jet(N);
    hold on
    for n = 1:N
        plot(Data.Ship.Longitude(n:n+1),Data.Ship.Latitude(n:n+1),'-o','color',cmap(ScaledSpeed(n),:));
    end
    text(Data.Ship.Longitude(1),Data.Ship.Latitude(1),datestr(Datetime(1),' mm/dd HH:MM'));
    text(Data.Ship.Longitude(end),Data.Ship.Latitude(end),datestr(Datetime(end),' mm/dd HH:MM'));
    set(gca,'CLim',[0,max(Speed)]);
    colormap(gca,cmap);
    h = colorbar;
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');
    ylabel(h,'Ship speed (m/s)');
    title('Ship trajectory');
    set(gca,'XGrid','on')
    set(gca,'YGrid','on')

    % Font
    set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

    % Save
    saveas(gcf,fullfile(Config.Plots.FiguresDirectory,strcat('MeasurementData',Config.Plots.FiguresFormat)))
    
    end

    % ==================================
    % Plot Velocity Column At Given Time
    % ==================================

    function PlotVelocityColumnAtGivenTime(Data,Config,GivenDatetime)

    % Plots velocity columns at a given date time.
    % The GivenDatetime should be for the format [yyyy,mm,dd,HH,MM,SS].

    DataMinTime = Data.Time(1);                           % In seconds
    DataMaxTime = Data.Time(end);                         % In seconds
    GivenTime = datenum(GivenDatetime) * (24 * 3600);     % In seconds

    if (GivenTime > DataMaxTime) || (GivenTime < DataMinTime)
        error('GivenTime: %s is not between data times: %s and %s', ...
        datestr(GivenTime / (24*3600)), ...
        datestr(DataMinTime / (24*3600)), ...
        datestr(DataMaxTime / (24*3600)))
    end

    t_diff = GivenTime - Data.Time';
    [~,Index] = min(abs(t_diff));
    TimeString = sprintf('%s',datestr(Data.Time(Index) / (24 * 3600)));

    % Velocity column at lookup time
    u_at_t = Data.Ocean.AbsoluteVel_u(:,Index);
    v_at_t = Data.Ocean.AbsoluteVel_v(:,Index);
    w_at_t = Data.Ocean.AbsoluteVel_w(:,Index);

    % Bathymetry
    NonZero = find(isnan(u_at_t));
    BathymetryIndex = NonZero(1);

    u_at_t = u_at_t(1:BathymetryIndex);
    v_at_t = v_at_t(1:BathymetryIndex);
    w_at_t = w_at_t(1:BathymetryIndex);
    range = Data.Ocean.Range(1:BathymetryIndex);

    % norm and angle for planar velocity (on x-y plane)
    vel_plane = u_at_t + i * v_at_t;
    vel_plane_mag = abs(vel_plane);
    vel_plane_direction = mod(angle(vel_plane),2*pi) * 180.0 / pi;

    % Full velocity in (x,y,z) space
    vel_3d_mag = sqrt(u_at_t.^2 + v_at_t.^2 + w_at_t.^2);
    vel_vertical_angle = atan2(w_at_t,vel_plane_mag) * 180.0 / pi;

    % Font
    set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

    % Figure
    figure()

    % Plot u,v velocities
    subplot(1,3,1)
    plot(u_at_t,range,'color','red','linewidth',1.5,'DisplayName','East')
    hold on
    plot(v_at_t,range,'color','blue','linewidth',1.5,'DisplayName','North')
    plot(w_at_t,range,'color',[0.6,0.8,0.9],'linewidth',1.5,'DisplayName','Vertical')
    legend()
    title('Velocities')
    ylabel('Depth (m)')
    xlabel(sprintf('u and v (m/s)\n%s',TimeString))
    set(gca,'Ydir','reverse')
    set(gca,'YGrid','on')
    set(gca,'XGrid','on')

    % Plot velocity magnitude
    subplot(1,3,2)
    plot(vel_plane_mag,range,'color','black','linewidth',1.5,'DisplayName','Planar velocity')
    hold on
    plot(vel_3d_mag,range,'color',[0.6,0.8,0.9],'linewidth',1.5,'DisplayName','3D velocity')
    ylabel('Depth (m)')
    xlabel(sprintf('Magnitude (m/s)\n%s',TimeString))
    title('Velocity Magnitude')
    set(gca,'Ydir','reverse')
    set(gca,'YGrid','on')
    legend()

    % Plot velocity direction
    subplot(1,3,3)
    plot(vel_plane_direction,range,'color','black','linewidth',1.5,'DisplayName','Planar direction')
    hold on
    plot(vel_vertical_angle,range,'color',[0.6,0.8,0.9],'linewidth',1.5,'DisplayName','Vertical angle')
    xlim([-90,360])
    ylabel('Depth (m)')
    xlabel(sprintf('Direction (degrees)\n%s',TimeString))
    title('Velocity Direction')
    set(gca,'Ydir','reverse')
    set(gca,'xtick',-90:90:360)
    set(gca,'YGrid','on');
    set(gca,'XGrid','on');
    legend()

    % Font
    set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

    end

    % =======================================
    % Plot Ship and Averaged Ocean Velocities
    % =======================================

    function PlotShipAndAveragedOceanVelocities(Data,Config)

    % Find max bathymetry
    % NonZero = find(all(isnan(Data.Ocean.AbsoluteVel_u),2));
    % MaxBathymetry = NonZero(1);

    SecondToDay = 1.0 / (24.0 * 3600.0);
    Datetime = Data.Time * SecondToDay;

    ShipVelocityMagnitude = hypot(Data.Ship.Vel_u,Data.Ship.Vel_v);

    U_averaged = zeros(size(Data.Time));
    V_averaged = zeros(size(Data.Time));
    OceanVel_ParallelToShip = zeros(size(Data.Time));
    OceanVel_TransversalToShip = zeros(size(Data.Time));

    for t = 1:size(Data.Time,2)

        NonZero = find(isnan(Data.Ocean.AbsoluteVel_u(:,t)));
        BathymetryIndex = NonZero(1) - 1;

        U = Data.Ocean.AbsoluteVel_u(1:BathymetryIndex,t);
        V = Data.Ocean.AbsoluteVel_v(1:BathymetryIndex,t);

        U_averaged(t) = mean(U);
        V_averaged(t) = mean(V);

        % Ship 
        ShipVel = [Data.Ship.Vel_u(t),Data.Ship.Vel_v(t)];
        ShipUnitDirection = ShipVel / hypot(ShipVel(1),ShipVel(2));
        ShipTransverseUnitDirection = [-ShipUnitDirection(2),ShipUnitDirection(1)];

        % Ocean velocity along ship direction
        OceanVels_ParallelToShip = U*ShipUnitDirection(1) + V*ShipUnitDirection(2);
        OceanVels_TransversalToShip = U*ShipTransverseUnitDirection(1) + V*ShipTransverseUnitDirection(2);

        % Averages
        OceanVel_ParallelToShip(t) = mean(OceanVels_ParallelToShip);
        OceanVel_TransverseToShip(t) = mean(OceanVels_TransversalToShip);

    end

    % Font
    set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

    figure()

    subplot(1,3,1)
    hold on
    plot(Datetime,Data.Ship.Vel_u,'linewidth',1.5,'color','red')
    plot(Datetime,U_averaged,'linewidth',1.5,'color','blue')
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Velocity (m)');
    xlim([Datetime(1),Datetime(end)])
    title('East velocity');
    set(gca,'XGrid','on')
    legend('Ship','Ocean (averaged over column)')

    subplot(1,3,2)
    hold on
    plot(Datetime,Data.Ship.Vel_v,'linewidth',1.5,'color','red')
    plot(Datetime,V_averaged,'linewidth',1.5,'color','blue');
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Velocity (m)');
    xlim([Datetime(1),Datetime(end)])
    title('North velocity');
    set(gca,'XGrid','on')
    legend('Ship','Ocean (averaged over column)')

    subplot(1,3,3)
    hold on
    plot(Datetime,ShipVelocityMagnitude,'linewidth',1.5,'color','red')
    plot(Datetime,OceanVel_ParallelToShip,'linewidth',1.5,'color','blue')
    plot(Datetime,OceanVel_TransverseToShip,'linewidth',1.5,'color',[0.9290, 0.6940, 0.1250])
    xticks(linspace(Datetime(1),Datetime(end),4))
    datetick('x','mm/dd HH:MM','keepticks','keeplimits')
    xtickangle(30)
    xlabel('Time (mm/dd HH:MM)');
    ylabel('Velocity (m)');
    xlim([Datetime(1),Datetime(end)])
    title('Velocity projection to ship course');
    set(gca,'XGrid','on')
    legend('Ship vel magnitude','Ocean vel parallel to ship','Ocean vel transverse to ship')

    % Font
    set(findall(gcf,'-property','FontName'),'FontName',Config.Plots.FontName)

    end

    % ---------------------
    % End of static methods

    end
end
