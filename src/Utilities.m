classdef Utilities

    methods(Static)

    % ====================
    % Convert LonLat To XY
    % ====================

    function [X,Y] = ConvertLonLatToXY(OriginLongitude,OriginLatitude,LongitudesArray,LatitudesArray)

    % LongitudeArray and LatitudeArray are N*M arrays. 
    % Outputs X and Y are also N*M arrays.

    % Average Earth radius (meters)
    EarthRadius = 6371008.8;

    % Degree To Radian
    DegreeToRadian = pi / 180.0;

    % Latitude to Y
    Y = (LatitudesArray - OriginLatitude) * DegreeToRadian * EarthRadius;

    % Longitude to X
    EarthAxisRadiusAtOriginLatitude = EarthRadius * cos(OriginLatitude * DegreeToRadian);
    X = (LongitudesArray - OriginLongitude) * DegreeToRadian * EarthAxisRadiusAtOriginLatitude;

    end

    % ====================
    % Convert XY To LonLat
    % ====================

    function [LongitudesArray,LatitudesArray] = ConvertXYToLonLat(OriginLongitude,OriginLatitude,X,Y)

    % X and Y are N*M array.
    % Outputs LongitudesArray and LatitudesArray are also N*M matrices.

    % Average Eath radius (meters)
    EarthRadius = 6371008.8;

    % Degree To Radian
    RadianToDegree = 180.0 / pi;

    % Y to Latitude
    LatitudesArray = OriginLatitude + (Y / EarthRadius) * RadianToDegree;

    % X to Longitude
    EarthAxisRadiusAtOriginLatitude = EarthRadius * cos(OriginLatitude / RadianToDegree);
    LongitudesArray = OriginLongitude + (X / EarthAxisRadiusAtOriginLatitude) * RadianToDegree;

    end 

    % ---------------------
    % End of static methods

    end
end
