% PLOTCALIBRATIONLINES Linearly fill calibrationlines needed by our data
% set and plot all of them in one figure.
%
% Yaguang Zhang, Purdue, 09/14/2017

clear; clc; close all;

%% Configurations

% The absolute path to the shared Google Drive folder "Annapolis
% Measurement Campaign". Please make sure it is correct for the machine
% which will run this script.
%  On Mac Lemma:
%    '/Users/zhan1472/Google Drive/Annapolis Measurement Campaign'
%  On Windows Dell:
%    'C:\Users\Zyglabs\Documents\MEGAsync\EARS'
ABS_PATH_TO_EARS_SHARED_FOLDER = ...
    'C:\Users\Zyglabs\Documents\MEGAsync\EARS';

ABS_PATH_TO_CALIBRATION_REF_POLYGONS = fullfile(...
    ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'Calibration');

ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation');

% Add libs to current path.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

%% Before Processing the Data

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

%% Get All the Gains Covered by the Measurement Data

gpsLogs = rdir(fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, 'Data', '**', ...
    '*_GPS.log'));
% Load the GPS samples.
gpsLogsParsed = arrayfun(...
    @(log) parseGpsLog(log.name), gpsLogs);
gains_needed = cell2mat(arrayfun( ...
    @(logP) str2num(logP.rxChannelGain), gpsLogsParsed, ...
    'UniformOutput', false));
gains_needed = unique(gains_needed);

%% Read the Polynomials for the Reference Calibration
% WE will get lsLinesPolys, lsLinesPolysInv, fittedMeaPs, fittedCalPs, and
% rxGains.

fileLsLinesPolys = fullfile(ABS_PATH_TO_CALIBRATION_REF_POLYGONS, ...
    'lsLinesPolys.mat');
if (exist(fileLsLinesPolys, 'file') ==2)
    load(fileLsLinesPolys);
else
    error('plotCalibrationLines:fileLsLinesPolysNotAvailable', ...
        ['Calibration results not available! \n', ...
        'Please run first: \n', ...
        '    EarsMeasurementCampaignCode/PostProcessing/2_0_Calibration/clibrateRx.m \n', ...
        'first and retry.']);
end

%% Compute and Plot Calibration Lines

% Because all the lines are linear, two points are enough for plotting.
[xRangeToShow, yRangeToShow] = deal([-100, 0]);
powerShiftsForCali = genCalibrationFct(lsLinesPolysInv, rxGains, ...
    gains_needed);

% Fix the colors to use.
rng(1);
% Plot.
hCalibrationLines = figure; hold on;
% Plot the calibration points used for coming up with the first two
% calibration lines.
hCaliPts = scatter(vertcat(fittedCalPs{1:end}), ...
    vertcat(fittedMeaPs{1:end}), '*', ...
    'MarkerEdgeColor', 'blue', ...
    'LineWidth',1.5);
% Plot these two calibrations lines, too.
for idxRefCaliLine = 1:length(lsLinesPolysInv)
    hLastRefCaliLine = ...
        plot(xRangeToShow, ...
        polyval(lsLinesPolysInv{idxRefCaliLine}, xRangeToShow), ...
        'Color', 'black','LineStyle', '--');
end
% Plot newly computed calibration lines.
for idxGain = 1:length(gains_needed)
    hLastNewCaliLine = plot(xRangeToShow, ...
        polyval([1 powerShiftsForCali(idxGain)],xRangeToShow), ...
        'LineStyle', '-');
end
axis([xRangeToShow, yRangeToShow]);
hold off;
title('Calibration Lines');
xlabel('Calculated Power (dB)'); ylabel('Measured Power (dB)');
legend([hCaliPts, hLastRefCaliLine, hLastNewCaliLine], ...
    'Calibration Pts', 'Ref Cali Lines', 'Cali Lines Needed');
transparentizeCurLegends;
grid minor; axis square;

% Save the plot.
pathCaliLinesFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'caliLinesNeeded');
saveas(hCalibrationLines, [pathCaliLinesFileToSave, '.png']);
saveas(hCalibrationLines, [pathCaliLinesFileToSave, '.fig']);

% EOF