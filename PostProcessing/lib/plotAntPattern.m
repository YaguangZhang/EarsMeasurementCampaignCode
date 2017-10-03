function [ hPat2DRef, hPat3DRef, ...
    hInterPat3DOnLineLinear, hInterPat3DOnLineDb, ...
    hInterPat3DWeightedSumLinear, hInterPat3DWeightedSumDb ] ...
    = plotAntPattern( patAz, patEl )
%PLOTANTPATTERN Plot the antenna pattern specified by the inputs patAz and
%patEl.
%
% Inputs:
%   - patAz, patEl
%     The antenna patterns, for the Azimuth and Elevation sweeps,
%     respectively; Each of which is a struct containing fields:
%       - azs
%         The azimuth angles in degree from set [0, 360).
%       - els
%         The elevation angles in degree from set [0, 360).
%       - amps
%         The linear amplitudes of the samples.
%       - phases
%         The phases of the samples.
%     All of these fields contains a column vector with each row
%     corresponding to a sweep sample.
%
% For the reference input sweep data, both a 2D illsutration and a 3D
% illustrations will be generated. For the interpolated data, a 3D
% illustration will be generated for each interpolation method.
%
% Ref:
% http://antennatutorials.blogspot.com/2013/05/radiation-pattern-of-half-wave-dipole.html
%
% Yaguang Zhang, Purdue, 10/02/2017

%% 2D Plot
hPat2DRef = figure('units','normalized', ...
    'outerposition',[0.1 0.05 0.8 0.9], 'Name','hPat2DRef');
% Azimuth.
polarAx = subplot(2,2,1,polaraxes);
plotAntPlanePat(polarAx, patAz.azs, 10.*log10(patAz.amps));
title('Azimuth Plane Pattern (Relative to the Minimum Amplitude)');

subplot(2,2,3);
angles = patAz.azs;
plot(angles, 10.*log10(patAz.amps));
title('Azimuth Sweep Data');
curAxis = axis; axis([min(angles), max(angles), curAxis(3:4)]);
xlabel('Azimuth'); ylabel('Normalized Amplitude (dB)');

% Elevation.
polarAx = subplot(2,2,2,polaraxes);
plotAntPlanePat(polarAx, patEl.els, 10.*log10(patEl.amps));
title('Elevation Plane Pattern (Relative to the Minimum Amplitude)');

subplot(2,2,4);
angles = patEl.els;
plot(angles, 10.*log10(patEl.amps));
title('Elevation Sweep Data');
curAxis = axis; axis([min(angles), max(angles), curAxis(3:4)]);
xlabel('Elevation'); ylabel('Normalized Amplitude (dB)');

%% 3D Plot

% Plot the interpolation results.
numPtsPerDim = 1000;
hInterPat3DOnLineLinear = plotInterPat3D( patAz, patEl, ...
    'OnLine', false, numPtsPerDim);
set(hInterPat3DOnLineLinear, 'Name', 'interPat3DOnLineLinear');
hInterPat3DOnLineDb = plotInterPat3D( patAz, patEl, ...
    'OnLine', true, numPtsPerDim);
set(hInterPat3DOnLineDb, 'Name', 'interPat3DOnLineDb');
hInterPat3DWeightedSumLinear = plotInterPat3D( patAz, patEl, ...
    'WeightedSum', false, numPtsPerDim);
set(hInterPat3DWeightedSumLinear, 'Name', 'interPat3DWeightedSumLinear');
hInterPat3DWeightedSumDb = plotInterPat3D( patAz, patEl, ...
    'WeightedSum', true, numPtsPerDim);
set(hInterPat3DWeightedSumDb, 'Name', 'interPat3DWeightedSumDb');

% We will also plot the sweep data, just for reference.
hPat3DRef = plotRefPat3D( patAz, patEl);

end
% EOF