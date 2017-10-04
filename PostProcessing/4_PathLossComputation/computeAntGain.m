function [ gain ] = computeAntGain(lat0, lon0, alt0, ...
    az0, el0, ...
    lat, lon, alt, ...
    antPatAz, antPatEl)
%COMPUTEANTGAIN Compute the antenna gain for the antenna located at (lat0,
%lon0, alt0) facing at (az0, el0), with its pattern specified by the
%Azimuth and Elevation sweep results (antPatAz, antPatEl), from a
%transmitter / to a receiver located at (lat, lon, alt).
%
%   - antPatAz, antPatEl
%     The antenna patterns, for the Azimuth and the Elevation sweeps; Each
%     of them is a struct containing fields:
%       - azs
%         The azimuth angles in degree.
%       - els
%         The elevation angles in degree.
%       - amps
%         The linear amplitudes of the samples.
%       - phases
%         The phases of the samples.
%
% Yaguang Zhang, Purdue, 10/04/2017

% Convert GPS to UTM so that the location coordinates will all in the unit
% of meter.
[x0, y0, zone0] = deg2utm(lat0, lon0);
[x, y, zone] = deg2utm(lat, lon);

if ~strcmp(zone0, zone)
    error('computeAntGain: The two locations specified are not in the same UTM zone!');
end

% Compute the direction vector from (x0,y0,alt0) to (x,y,alt).
xDiff = x-x0;
yDiff = y-y0;
zDiff = alt-alt0;
% Compute the (az,el,r) for that direction. Note that the resulted angles
% are in radians.
[az1, el1, r1] = cart2sph(xDiff, yDiff, zDiff);

% Get the spherical coordinates for that direction with respect of the
% origin antenna's point of view. One can think of the operation as:
%   1. Plot the vector for the direction from the origin antenna to the
%   destination antenna, which is (xDiff, yDiff, zDiff) here;
%  2.  Plot in the same Cartesian coordinate system a vector for the
%  direction to which origin antenna is facing; We essentially use
%  sph2cart(az0, el0, 1) here.
%   3. Rotate these two vectors at the same time, first along Azimuth then
%   along Elevation (of vector in (2)), until the vector in (2) align perfectly with (1,0,0),
%   i.e. the position x axis.
%  4.  The resulted vector in (1) will give us the (az, el) needed.

% Rotate along Azimuth.
az1 = az1-deg2rad(az0);
[x1, y1, z1] = sph2cart(az1, el1, r1);

% Rotate along Elevation of vector in (2). Essentially, we keep y1
% untouched, and rotate (x1, 0, z1) to get the new x1 and z1.
[az2, el2, r2] = cart2sph(x1, 0, z1);
[x1, ~, z1] = sph2cart(az2, el2-deg2rad(el0), r2);

% Get the final result.
[az, el, ~] = cart2sph(x1, y1, z1);

% Now get the antenna gain from the interpolated patten. We will
% interpolate directly using the linear amplitudes of the antenna pattern
% with the "WeightedSum" method.
amps = antPatInter(antPatAz, antPatEl, az, el, 'WeightedSum');

gain = 10.*log10(amps);

end

%EOF