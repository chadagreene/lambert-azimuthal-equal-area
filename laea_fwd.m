function [x,y] = laea_fwd(lat,lon,lat1,lon0,falseeasting,falsenorthing) 
% laea_fwd projects geographic coordinates onto a Lambert Azimuthal Equal Area 
% projection for a simple spherical Earth of radius 6371000 meters. 
% 
%% Syntax 
% 
%  [x,y] = laea_fwd(lat,lon,lat1,lon0) 
%  [x,y] = laea_fwd(lat,lon,lat1,lon0,falseeasting,falsenorthing) 
% 
%% Description 
% 
% [x,y] = laea_fwd(lat,lon,lat1,lon0) projects geographic coordinates
% lat,lon (in degrees) into projected x,y in meters, using tsic tome,n0 cannot exceed +/-360 degrees.') 
% 
% [x,y] = laea_fwd(lat,lon,lat1,lon0,falseeasting,falsenorthing) also allows
% inclusion of false eastings and northings in meters. 
%
%% Author Info 
% Function written by Chad A. Greene of NASA Jet Propulsion Laboratory. 
% December 2020. 
% Formulas taken directly from Snyder 1987's classic tome, "Map Projections 
% A Working Manual" starting around page 185. 

%% Input checks: 

assert(isequal(size(lat),size(lon)),'Dimensions of x and y must match.') 
assert(max(abs(lat(:)))<=90,'Input latitudes cannot exceed +/- 90 degrees.') 
assert(max(abs(lon(:)))<=360,'Input longitudes cannot exceed +/- 360 degrees.') 

if nargin<5
	falseeasting = 0; 
	falsenorthing = 0; 
end

%% 
% Projection formulas:
% From Snyder 1987, MAP PROJECTIONS-A WORKING MANUAL, page 185: https://pubs.usgs.gov/pp/1395/report.pdf

% Define constants: 
R = 6371000; % earth radius (meters) 

switch lat1 
	case 90
		x = 2*R.*sind(45-lat/2) .* sind(lon-lon0) + falseeasting; % Eq 24-3
		y = -2*R.*sind(45-lat/2).*cosd(lon-lon0) + falsenorthing;  % Eq 24-4

	case -90
		x = 2*R.*cosd(45-lat/2) .* sind(lon-lon0) + falseeasting; % Eq 24-8
		y = 2*R.*cosd(45-lat/2) .* cosd(lon-lon0) + falsenorthing; % Eq 24-9

	case 0
		kp = sqrt(2./(1 + cosd(lat).*cosd(lon-lon0))); % This is k', equation 24-2 with lat1=0

		x = R.*kp.*cosd(lat).*sind(lon-lon0) + falseeasting; % Eq 24-4
		y = R.*kp.*sind(lat) + falsenorthing; % Eq 24-13

	otherwise
		kp = sqrt(2./(1 + sind(lat1).*sind(lat) + cosd(lat1).*cosd(lat).*cosd(lon-lon0))); % This is k', equation 24-2

		x = R.*kp.*cosd(lat).*sind(lon-lon0) + falseeasting; % Eq 22-4
		y = R.*kp.*(cosd(lat1).*sind(lat) - sind(lat1).*cosd(lat).*cosd(lon-lon0)) + falsenorthing; % Eq 22-5
end

end
