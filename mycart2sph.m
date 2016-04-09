function [az,elev,r] = mycart2sph(x,y,z)

%   3D Cartesian to Spherical Coordinates
%
%   [az,elev,r] = mycart2sph(x,y,z)
%
%   Inputs:
%       x       x-coordinate
%       y       y-coordinate
%       z       z-coordinate
%
%   Outputs:
%       az      Azimuth angle in radians, measured counterclockwise
%               from the positive x axis (otherwise referred to as phi)
%       elev    Elevation angle in radians, from xy plane (otherwise 
%               referred to as theta) 
%       r       Radius
%
%   Notes:
%       The MATLAB function cart2sph reverses phi and theta.
%
%**************************************************************************
% Author:           E. A. P. Habets and M. R. P. Thomas
% Date:             27 July 2010
%**************************************************************************

hypotxy = hypot(x,y);
r = hypot(hypotxy,z);
elev = acos(z./r);
az = atan2(y,x);
