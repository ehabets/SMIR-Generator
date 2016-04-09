function [x,y,z] = mysph2cart(az,inc,r)

%   Spherical Coordinates to 3D Cartesian
%
%   [x,y,z] = mysph2cart(az,inc,r)
%
%   Inputs:
%       az      Azimuth angle in radians, measured counterclockwise
%               from the positive x axis (otherwise referred to as phi)
%       inc     Inclination angle in radians, from positive z axis 
%               (otherwise referred to as theta) 
%       r       Radius
%
%   Outputs:
%       x       x-coordinate
%       y       y-coordinate
%       z       z-coordinate
%
%   Notes:
%       The MATLAB function cart2sph reverses phi and theta.
%
%**************************************************************************
% Author:           E. A. P. Habets, M. R. P. Thomas and D. P. Jarrett
% Date:             24 August 2011
%**************************************************************************

z = r .* cos(inc);
rcosinc = r .* sin(inc);
x = rcosinc .* cos(az);
y = rcosinc .* sin(az);