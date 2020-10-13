% 
% Plot all drifters between a give time window
% 

% ---------------- Settings ---------------------

% Filename, start and end of time to plot
DrifterFilename   = '/home/sia/work/alpha/WhereIsMyPlume/Drifters.mat';
ShipTrackFilename = '/home/sia/work/alpha/WhereIsMyPlume/ShipTrack.mat';

% First dye relase
% StartPlotDate = '14-AUG-2018 12:00:00';        % In UTC
% EndPlotDate = '14-AUG-2018 17:00:00';          % In UTC

% Second dye release
% StartPlotDate = '16-AUG-2018 11:00:00';      % In UTC
% EndPlotDate = '16-AUG-2018 17:00:00';        % In UTC

StartPlotDate = '16-AUG-2018 08:50:00';      % In UTC
EndPlotDate = '16-AUG-2018 18:00:00';        % In UTC

% -----------------------------------------------

% Plot Drifters
DriftersFileData = load(DrifterFilename);
Drifters = DriftersFileData.Drifters;

% Convert dates to times
StartPlotTime = datenum(StartPlotDate);
EndPlotTime = datenum(EndPlotDate);

PlotHandles = [];
Counter = 0;
NumberOfDrifters = size(Drifters,2);
cmap = jet(NumberOfDrifters);
DrifterColor = [0.8,0.2,0.2];

% Plot each drifter
for DrifterId = 1:NumberOfDrifters

    DrifterTime = Drifters(DrifterId).Time;
    
    % Find start and end of time indices of plots
    StartPlotTimeIndex = find(DrifterTime >= StartPlotTime,1,'first');
    EndPlotTimeIndex = find(DrifterTime <= EndPlotTime,1,'last');

    if (~isempty(StartPlotTimeIndex)) && (~isempty(EndPlotTimeIndex))
        
        % Get Drifter info
        DrifterLongitude = Drifters(DrifterId).Longitude(StartPlotTimeIndex:EndPlotTimeIndex);
        DrifterLatitude = Drifters(DrifterId).Latitude(StartPlotTimeIndex:EndPlotTimeIndex);
        
        if ~isempty(DrifterLongitude) && ~isempty(DrifterLatitude)
            Counter = Counter + 1;

            fprintf('Plotting drifter: %d\n',DrifterId)

            PlotHandles(Counter) = plot(DrifterLongitude,DrifterLatitude,'color',cmap(DrifterId,:),'DisplayName',strcat('Drifter ',num2str(DrifterId)));
            hold on
            plot(DrifterLongitude(1),DrifterLatitude(1),'o','MarkerEdgeColor','black','MarkerFaceColor',DrifterColor,'MarkerSize',4)
            plot(DrifterLongitude(end),DrifterLatitude(end),'v','MarkerEdgeColor','black','MarkerFaceColor',DrifterColor,'MarkerSize',4)
            text(DrifterLongitude(floor(end/2+0.5)),DrifterLatitude(floor(end/2+0.5)),num2str(DrifterId),'HorizontalAlignment','center','FontSize',10)
        
        end
    end
end

% Plot ship track
ShipTrackFile = load(ShipTrackFilename);
ShipTrack = ShipTrackFile.ShipTrack;
StartPlotTimeIndex = find(ShipTrack.Time >= StartPlotTime,1,'first');
EndPlotTimeIndex = find(ShipTrack.Time <= EndPlotTime,1,'last');
ShipTrackLongitude = ShipTrack.Longitude(StartPlotTimeIndex:EndPlotTimeIndex);
ShipTrackLatitude = ShipTrack.Latitude(StartPlotTimeIndex:EndPlotTimeIndex);
PlotHandles(Counter+1) = plot(ShipTrackLongitude,ShipTrackLatitude,'--','linewidth',2,'color','black','DisplayName','Ship track');

text(0.05,0.15,sprintf('Start time: %s',StartPlotDate),'units','normalized','Backgroundcolor',[1,1,0])
text(0.05,0.05,sprintf('End time: %s',EndPlotDate),'units','normalized','Backgroundcolor',[1,1,0])

daspect([1,cos(ShipTrack.Latitude(1)*pi/180),1])
title('Drifter and ship trajectories')
xlabel('Longitude')
ylabel('Latitude')
legend(PlotHandles)

if Counter == 0
    fprintf('No drifter was found within %s to %s.\n',StartPlotDate,EndPlotDate);
end
