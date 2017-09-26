% COMPUTEPATHLOSSES Compute the path losses in dB for all the locations
% covered by our measurement data set.
%
% We will consider both the TX calibration and the antenna normalization.
% The results from
%   - PostProcessing/1_SummaryReport/genPlots.m and
%     Output file plotInfo.mat contains the information for all the
%     measurement data files found locally on the machine. Note that only
%     the information for _LargeScale, _SIMO and _Conti folders was saved.
%   - PostProcessing/2_0_Calibration/calibrateRx.m
%     Output file plotInfo.mat contains the information for all the
%     measurement data files found locally on the machine.
% will be reused.
%
% Yaguang Zhang, Purdue, 09/26/2017

clear; clc; close all;

%% Configurations

% Add libs to current path and set ABS_PATH_TO_EARS_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% We also need the function thresholdWaveform.m for noise elimination.
addpath(fullfile(pwd, '2_0_Calibration'));

% Configure other paths accordingly.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation');

% Reuse results from plotInfo.m and calibrateRx.m.
ABS_PATH_TO_DATA_INFO_FILE = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'SummaryReport', 'plots', 'plotInfo.mat');
ABS_PATH_TO_CALI_LINES_FILE = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'Calibration', 'lsLinesPolys.mat');

% Sample rate used for GnuRadio.
Fs = 1.04 * 10^6;

% For setting the threshold during the noise elimination.
NUM_SIGMA_FOR_THRESHOLD = 3.5;

%% Before Processing the Data

disp(' ----------------------- ')
disp('  computePathLosses')
disp(' ----------------------- ')

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

%% Get Info for Measurement Data Files and Calibration Polynomials

disp(' ')
disp('    Loading results from plotInfo.m and calibrateRx.m ...')

assert(exist(ABS_PATH_TO_DATA_INFO_FILE, 'file')==2, ...
    'Couldn''t find plotInfo.mat! Please run PostProcessing/1_SummaryReport/genPlots.m first.');
assert(exist(ABS_PATH_TO_CALI_LINES_FILE, 'file')==2, ...
    'Couldn''t find lsLinesPolys.mat! Please run PostProcessing/2_0_Calibration/calibrateRx.m  first.');

% The data have been processed before and the result files have been found.
disp('    Found plotInfo.mat and lsLinesPolys.mat');
disp('        Loading the results...')
% Get 'allSeriesParentDirs' and 'allSeriesDirs'.
load(ABS_PATH_TO_DATA_INFO_FILE);
% Get 'lsLinesPolys', 'lsLinesPolysInv', 'fittedMeaPs', 'fittedCalPs', and
% 'rxGains'.
load(ABS_PATH_TO_CALI_LINES_FILE);

disp('    Done!')

%% Search for the Measurement Data Files

disp(' ')
disp('    Searching for measurement data files ...')

% Here we don't care about when the data were collected so let's rearrange
% all the dir struct into one array.
allSeriesDirsArray = vertcat(allSeriesDirs{1:end});
numSeries = length(allSeriesDirsArray);

% Scan the series folder for Gnu Radio out files, as well as the
% corresponding GPS log files.
[allOutFilesDirs, allGpsFilesDirs] = deal(cell(numSeries,1));

for idxSeries = 1:numSeries
    disp(['        Scanning series folder ', num2str(idxSeries), '/', ...
        num2str(numSeries), '...']);
    
    [allOutFilesDirs{idxSeries}, allGpsFilesDirs{idxSeries}] = ...
        loadFilesFromSeriesDir(allSeriesDirsArray(idxSeries));
end

disp('    Done!')

%% Compute Path Losses

disp(' ')
disp('    Computing path losses...')

% Compute the path losses and save them into a matrix together with the GPS
% info.
numOutFiles = sum(cellfun(@(d) length(d), allOutFilesDirs));
% More specifically, each row is a [path loss (dB), lat, lon] array.
pathLossesWithGpsInfo = nan(numOutFiles, 3);
pathLossCounter = 1;
for idxSeries = 1:numSeries
    disp(['        Processing series ', num2str(idxSeries), '/', ...
        num2str(numSeries), '...']);
    for idxOutFile = 1:numOutFiles
        disp(['            Outfile ', num2str(idxOutFile), '/', ...
            num2str(numOutFiles), '...']);
        
        curOutFileDir = allOutFilesDirs{idxSeries}(idxOutFile);
        [lat, lon, gpsLog] = fetchGpsForOutFileDir(curOutFileDir);
        
        % Compute b for the calibration line corresponding to the RX gain.
        usrpGain = str2num(gpsLog.rxChannelGain);
        powerShiftsForCali = genCalibrationFct( lsLinesPolysInv, ...
            rxGains, usrpGain);

        % Compute path loss. We will use the amplitude version of
        % thresholdWaveform.m without plots for debugging as the noise
        % eliminiation function.
        noiseEliminationFct = @(waveform) thresholdWaveform(abs(waveform));
        pathLossInDb = computePathLossForOutFileDir(curOutFileDir, ...
            usrpGain, noiseEliminationFct, powerShiftsForCali);

        % Store the results.
        pathLossesWithGpsInfo(pathLossCounter,:) = [pathLossInDb, lat, lon];
        pathLossCounter = pathLossCounter+1;
    end
end

disp('    Done!')

%% Plot

disp(' ')
disp('    Plotting...')

% Fix the colors to use.
rng(1);
% Plot.
hPathLossesOnMap = figure; hold on;
plot(pathLossesWithGpsInfo(:,3), pathLossesWithGpsInfo(:,2), 'w.');
plot_google_map;
plot3k([pathLossesWithGpsInfo(:,3), pathLossesWithGpsInfo(:,2), ...
    pathLossesWithGpsInfo(:,1)]);

% Save the plot.
pathCaliLinesFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'pathLossesOnMap');
saveas(hPathLossesOnMap, [pathCaliLinesFileToSave, '.png']);
saveas(hPathLossesOnMap, [pathCaliLinesFileToSave, '.fig']);

disp('    Done!')

% EOF