% COMPUTEPATHLOSSES Compute the path losses in dB for all the locations
% covered by our measurement data set (excluding the Conti case).
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
% Note that we only load and process the LargeScale and SIMO measurements
% in this script.
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

% Transmitter location.
TX_LAT = 38.983899;
TX_LON = -76.486682;

% TX power into upconverter in dBm.
txPower  = -8;

%% Before Processing the Data

disp(' ------------------- ')
disp('  computePathLosses ')
disp(' ------------------- ')

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
    
    % Here it doesn't make much sense to load the conti measurements and
    % come up with only 1 path loss value for each long sequence. We will
    % deal with them separately with another script.
    regexPattern = '\d+_(LargeScale|SIMO)';
    [allOutFilesDirs{idxSeries}, allGpsFilesDirs{idxSeries}] = ...
        loadFilesFromSeriesDir(allSeriesDirsArray(idxSeries), regexPattern);
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
% Also save the meta info needed to map the path loss back to the
% measurements. We choose to save the full file path to the .out file for
% convenience.
absPathsOutFiles = cell(numOutFiles, 1);
for idxSeries = 1:numSeries
    disp(['        Processing series ', num2str(idxSeries), '/', ...
        num2str(numSeries), '...']);
    
    numOutFileCurSeries = length(allOutFilesDirs{idxSeries});
    for idxOutFile = 1:numOutFileCurSeries
        disp(['            Outfile ', num2str(idxOutFile), '/', ...
            num2str(numOutFileCurSeries), '...']);
        
        curOutFileDir = allOutFilesDirs{idxSeries}(idxOutFile);
        [lat, lon, gpsLog] = fetchGpsForOutFileDir(curOutFileDir);
        
        % Compute b for the calibration line corresponding to the RX gain.
        usrpGain = str2double(gpsLog.rxChannelGain);
        powerShiftsForCali = genCalibrationFct( lsLinesPolysInv, ...
            rxGains, usrpGain);
        
        % Compute path loss. We will use the amplitude version of
        % thresholdWaveform.m without plots for debugging as the noise
        % eliminiation function.
        noiseEliminationFct = @(waveform) thresholdWaveform(abs(waveform));
        [ pathLossInDb, absPathOutFile] ...
            = computePathLossForOutFileDir(curOutFileDir, ...
            usrpGain, noiseEliminationFct, powerShiftsForCali);
        
        % Store the results.
        pathLossesWithGpsInfo(pathLossCounter,:) = [pathLossInDb, lat, lon];
        absPathsOutFiles{pathLossCounter} = absPathOutFile;
        pathLossCounter = pathLossCounter+1;
    end
end
assert(all(~isnan(pathLossesWithGpsInfo(1:end))));

disp('    Saving the results...')
% For absPathsOutFiles, convert it to relative paths under the data folder,
% which will already contain enough information we need.
relPathsOutFilesUnderDataFolder = ...
    cellfun(@(p) regexp(p, 'Data[\/\\]([a-zA-Z\d\/\\_]+.out)$', ...
    'tokens'), absPathsOutFiles);
pathPathLossFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'pathLossesWithGpsInfo.mat');
save(pathPathLossFileToSave, ...
    'pathLossesWithGpsInfo', 'relPathsOutFilesUnderDataFolder');

disp('    Done!')

%% Plot

disp(' ')
disp('    Plotting...')

boolsInvalidCoor = pathLossesWithGpsInfo(:,2)==0 ...
    & pathLossesWithGpsInfo(:,3)==0;
if any(boolsInvalidCoor)
    warning([num2str(sum(boolsInvalidCoor)), ...
        ' invalid (lat, lon) pairs detected (both are 0).', ...
        ' We will ignore these points together with their path losses.']);
end
pathLossesWithValidGps = pathLossesWithGpsInfo(~boolsInvalidCoor,:);

boolsInfPathloss = isinf(pathLossesWithValidGps(:,1));
if any(boolsInfPathloss)
    warning([num2str(sum(boolsInfPathloss)), ...
        ' inf path loss detected.', ...
        ' We will show these points at z=0 with different markers.']);
end
validPathLossesWithValidGps = pathLossesWithValidGps(~boolsInfPathloss,:);

% Plot path losses on map.
hPathLossesOnMap = figure; hold on; colormap jet;
plot(validPathLossesWithValidGps(:,3), validPathLossesWithValidGps(:,2), 'w.');
plot(pathLossesWithValidGps(boolsInfPathloss,3), ...
    pathLossesWithValidGps(boolsInfPathloss,2), 'kx');
hTx = plot(TX_LON, TX_LAT, '^b');
plot_google_map('MapType','satellite');
plot3k([validPathLossesWithValidGps(:,3), validPathLossesWithValidGps(:,2), ...
    validPathLossesWithValidGps(:,1)], 'Marker', {'.', 12});
% The command plot_google_map messes up the color legend of plot3k, so we
% will have to fix it here.
hCb = findall( allchild(hPathLossesOnMap), 'type', 'colorbar');
hCb.Ticks = linspace(1,length(colormap),length(hCb.TickLabels));
hold off; grid on; view(0, 90); legend(hTx, 'TX');
title('Path Losses on Map (Large Scale & SIMO)');
xlabel('Lon'); ylabel('Lat'); zlabel('Path Loss (dB)');

% Plot path losses over distance from Tx.
validPLWithValidGPSCell = num2cell(validPathLossesWithValidGps, 2);
distsFromTx = cellfun(@(s) 1000.*lldistkm([s(2) s(3)],[TX_LAT,TX_LON]), ...
    validPLWithValidGPSCell);

hPathLossesOverDist = figure; hold on; colormap jet;
plot3k([distsFromTx, zeros(length(distsFromTx),1), ...
    validPathLossesWithValidGps(:,1)], 'Marker', {'.', 6});
curAxis = axis;
axis([min(distsFromTx)-10, max(distsFromTx)+100, curAxis(3:6)]);
view(0, 0); set(gca, 'XScale', 'log'); grid on;
newXTicks = [10,100,200,500,1000];
set(gca, 'XTickLabels',newXTicks);
set(gca, 'XTick',newXTicks);
title('Path Losses over Distance (Large Scale & SIMO)');
xlabel('Distance to Tx (m)'); ylabel(''); zlabel('Path Loss (dB)');

% Save the plots.
pathPathossesOnMapFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'pathLossesOnMap');
saveas(hPathLossesOnMap, [pathPathossesOnMapFileToSave, '.png']);
saveas(hPathLossesOnMap, [pathPathossesOnMapFileToSave, '.fig']);
% Save the plot.
pathPathLossesOverDistFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'pathLossesOverDist');
saveas(hPathLossesOverDist, [pathPathLossesOverDistFileToSave, '.png']);
saveas(hPathLossesOverDist, [pathPathLossesOverDistFileToSave, '.fig']);

disp('    Done!')

% EOF