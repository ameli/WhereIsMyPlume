% ---------------------------------------
%
% This class is used in StochasticModel.m
%
% ---------------------------------------

classdef ProbabilityDensityFunctions

    methods(Static)

    % ===================================
    % Estimate 2D Kernel Density Function
    % ===================================

    function [PDF_ValuesOnGrid,PDF_LongitudeGrid,PDF_LatitudeGrid,X_Grid,Y_Grid] = Estimate2DKernelDensityFunction(StochasticTrajectoriesLongitudes,StochasticTrajectoriesLatitudes)

        % Get the location of end points
        Lons = StochasticTrajectoriesLongitudes(end,:);
        Lats = StochasticTrajectoriesLatitudes(end,:);

        % Center and size of grid
        Lons_Center = mean(Lons);
        Lats_Center = mean(Lats);
        Lons_Span = max(Lons) - min(Lons);
        Lats_Span = max(Lats) - min(Lats);

        % Min and max of grid
        Scale = 2.5;
        MinLon = Lons_Center - Scale * Lons_Span/2.0;
        MaxLon = Lons_Center + Scale * Lons_Span/2.0;
        MinLat = Lats_Center - Scale * Lats_Span/2.0;
        MaxLat = Lats_Center + Scale * Lats_Span/2.0;

        % Create a grid
        Resolution = 500;
        Lon_Axis = linspace(MinLon,MaxLon,Resolution);
        Lat_Axis = linspace(MinLat,MaxLat,Resolution);
        [Lon_Grid,Lat_Grid] = meshgrid(Lon_Axis,Lat_Axis);

        % Grid in XY coordinate (needed for correct estimation of PDF based on actual lenght distances)
        OriginLongitude = Lons_Center;
        OriginLatitude = Lats_Center;
        [X_Data,Y_Data] = Utilities.ConvertLonLatToXY(OriginLongitude,OriginLatitude,Lons,Lats);
        [X_Grid,Y_Grid] = Utilities.ConvertLonLatToXY(OriginLongitude,OriginLatitude,Lon_Grid,Lat_Grid);

        % Scatter point of the grid
        GridPoints = [X_Grid(:), Y_Grid(:)];

        % Kernel density estimation
        DataPoints = [X_Data(:),Y_Data(:)];
        rng('default');
        [PDF,PDF_Points_XY] = ksdensity(DataPoints,GridPoints,'Function','pdf');

        % Reshape to grid
        PDF_ValuesOnGrid = reshape(PDF,size(Lon_Grid));
        PDF_LongitudeGrid = Lon_Grid;
        PDF_LatitudeGrid = Lat_Grid;
    end

    % ===========
    % Compute CDF
    % ===========

    function CDF_OnLevelSet = ComputeCDF(X_Grid,Y_Grid,PDF,LevelSetValue)
        
        % Make points below levelset to be zero
        PDF_InsideLevelSet = PDF(:,:);
        PDF_InsideLevelSet(PDF_InsideLevelSet <= LevelSetValue) = 0.0;

        % Integrate
        dx = X_Grid(1,2) - X_Grid(1,1);
        dy = Y_Grid(2,1) - Y_Grid(1,1);
        CDF_OnLevelSet = trapz(trapz(PDF_InsideLevelSet)) * dx * dy;
    end

    % ==============================================
    % Complementary Cummulative Distribution Functon
    % ==============================================

    function CCDF_ValuesOnGrid = ComplementaryCumulativeDistributionFunction(X_Grid,Y_Grid,PDF)

    % Level sets, from zero to the maximum of PDF
    LevelSets = linspace(min(PDF(:)),max(PDF(:)),500);

    % CDF, this is an array that each element is the CDF integral inside the levele sets of PDF.
    CDF_OnLevelSet = zeros(1,length(LevelSets));
    for i = 1:length(CDF_OnLevelSet)
        CDF_OnLevelSet(i) = ProbabilityDensityFunctions.ComputeCDF(X_Grid,Y_Grid,PDF,LevelSets(i));
    end

    % Initialize
    AllLevelSetPoints_X = [];
    AllLevelSetPoints_Y = [];
    AllLevelSet_CDF_Values = [];

    % Get points on the level sets
    for i = 1:length(LevelSets)

        % Get xy points of the level set
        LevelSetPoints = contourc(X_Grid(1,:),Y_Grid(:,1),PDF,[LevelSets(i),LevelSets(i)]);

        % Extract points
        LevelSetPoints_X = LevelSetPoints(1,2:end);
        LevelSetPoints_Y = LevelSetPoints(2,2:end);

        % CDF values for this levelset
        LevelSet_CDF_Values = ones(1,length(LevelSetPoints_X)) * CDF_OnLevelSet(i);

        % Append
        AllLevelSetPoints_X = [AllLevelSetPoints_X,LevelSetPoints_X(:)'];
        AllLevelSetPoints_Y = [AllLevelSetPoints_Y,LevelSetPoints_Y(:)'];
        AllLevelSet_CDF_Values = [AllLevelSet_CDF_Values,LevelSet_CDF_Values(:)'];

    end

    % Include lower boundary points as ones
    AllLevelSetPoints_X = [AllLevelSetPoints_X,X_Grid(1,:)];
    AllLevelSetPoints_Y = [AllLevelSetPoints_Y,Y_Grid(1,:)];
    AllLevelSet_CDF_Values = [AllLevelSet_CDF_Values,ones(1,size(X_Grid,2))];

    % Include uper boundary points as ones
    AllLevelSetPoints_X = [AllLevelSetPoints_X,X_Grid(end,:)];
    AllLevelSetPoints_Y = [AllLevelSetPoints_Y,Y_Grid(end,:)];
    AllLevelSet_CDF_Values = [AllLevelSet_CDF_Values,ones(1,size(X_Grid,2))];

    % Include left boundary points as ones
    AllLevelSetPoints_X = [AllLevelSetPoints_X,X_Grid(:,1)'];
    AllLevelSetPoints_Y = [AllLevelSetPoints_Y,Y_Grid(:,1)'];
    AllLevelSet_CDF_Values = [AllLevelSet_CDF_Values,ones(1,size(X_Grid,1))];

    % Include right boundary points as ones
    AllLevelSetPoints_X = [AllLevelSetPoints_X,X_Grid(:,end)'];
    AllLevelSetPoints_Y = [AllLevelSetPoints_Y,Y_Grid(:,end)'];
    AllLevelSet_CDF_Values = [AllLevelSet_CDF_Values,ones(1,size(X_Grid,1))];

    % Regrid the scattered collection of level set points onto the original X_Grid and Y_Grid
    CDF_ValuesOnGrid = griddata(AllLevelSetPoints_X,AllLevelSetPoints_Y,AllLevelSet_CDF_Values,X_Grid,Y_Grid);

    % CCDF on grid
    CCDF_ValuesOnGrid = 1 - CDF_ValuesOnGrid;

    end

    % ---------------------
    % End of static methods

    end
end
