% EVALPATHLOSSES Evaluate the path losses in dB for all the locations
% covered by our measurement data set (excluding the Conti case).
%
% We will consider both the TX calibration and the antenna normalization,
% (i.e. the resulted path losses are antenna-independent, which are known
% as the Basic Transmission Losses).
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

clear; clc; close all; dbstop if error;

%% Configurations

warning('on');

% Add libs to current path and set ABS_PATH_TO_EARS_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% We will need the function thresholdWaveform.m for noise elimination.
addpath(fullfile(pwd, '2_0_Calibration'));
% We also need the function antPatInter.m for antenna gain calculation.
addpath(fullfile(pwd, '3_AntennaPattern'));

% Configure other paths accordingly.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation');

% Reuse results from plotInfo.m, calibrateRx.m, fetchAntennaPattern.m, and
% loadMeasCampaignInfo.m.
ABS_PATH_TO_DATA_INFO_FILE = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'SummaryReport', 'plots', 'plotInfo.mat');
ABS_PATH_TO_CALI_LINES_FILE = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'Calibration', 'lsLinesPolys.mat');
ABS_PATH_TO_ANT_PAT_FILE = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'AntennaPattern', 'antennaPattern.mat');
ABS_PATH_TO_TX_INFO_LOGS_FILE= fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation', 'txInfoLogs.mat');

% For setting the threshold during the noise elimination.
NUM_SIGMA_FOR_THRESHOLD = 3.5;

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
disp('    Loading results from: ')
disp('      - plotInfo.m')
disp('      - calibrateRx.m')
disp('      - antennaPattern.m')
disp('      - loadMeasCampaignInfo.m')

assert(exist(ABS_PATH_TO_DATA_INFO_FILE, 'file')==2, ...
    'Couldn''t find plotInfo.mat! Please run PostProcessing/1_SummaryReport/genPlots.m first.');
assert(exist(ABS_PATH_TO_CALI_LINES_FILE, 'file')==2, ...
    'Couldn''t find lsLinesPolys.mat! Please run PostProcessing/2_0_Calibration/calibrateRx.m first.');
assert(exist(ABS_PATH_TO_ANT_PAT_FILE, 'file')==2, ...
    'Couldn''t find antennaPattern.mat! Please run PostProcessing/3_AntennaPattern/fetchAntennaPattern.m first.');
assert(exist(ABS_PATH_TO_TX_INFO_LOGS_FILE, 'file')==2, ...
    'Couldn''t find txInfoLogs.mat! Please run PostProcessing/4_0_PathLossComputation/loadMeasCampaignInfo.m first.');

% The data have been processed before and the result files have been found.
disp('    Found all .mat files required.');
disp('        Loading the results...')
% Get 'allSeriesParentDirs' and 'allSeriesDirs'.
load(ABS_PATH_TO_DATA_INFO_FILE);
% Get 'lsLinesPolys', 'lsLinesPolysInv', 'fittedMeaPs', 'fittedCalPs', and
% 'rxGains'.
load(ABS_PATH_TO_CALI_LINES_FILE);
% Get 'pat28GAzNorm', and 'pat28GElNorm'.
load(ABS_PATH_TO_ANT_PAT_FILE);

% Get records of the TxInfo.txt files (among other contant parameters for
% the measurement campaign, e.g. F_S, TX_LAT, TX_LON, and TX_POWER_DBM):
% 'TX_INFO_LOGS' and 'TX_INFO_LOGS_ABS_PAR_DIRS'.
load(ABS_PATH_TO_TX_INFO_LOGS_FILE);
% Sample rate used for GnuRadio. Needed by computePathLossForOutFileDir.m.
Fs = F_S;
% TX power into upconverter in dBm.
txPower  = TX_POWER_DBM;

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
% More specifically, each row is a [path loss (dB), lat, lon, latM, lonM]
% array, where (lat, lon) is the GPS coordinates for the individule .out
% file, while (latM, latM) is the average (via median) coordinates for all
% the locked GPS samples on that site.
pathLossesWithGpsInfo = nan(numOutFiles, 5);
pathLossCounter = 1;
% Also save the meta info needed to map the path loss back to the
% measurements. We choose to save the full file path to the .out file for
% convenience.
absPathsOutFiles = cell(numOutFiles, 1);
for idxSeries = 1:numSeries
    disp(['        Processing series ', num2str(idxSeries), '/', ...
        num2str(numSeries), '...']);
    disp(['            Folder: ', allSeriesDirsArray(idxSeries).name]);
    
    numOutFileCurSeries = length(allOutFilesDirs{idxSeries});
    for idxOutFile = 1:numOutFileCurSeries
        disp(['            Outfile ', num2str(idxOutFile), '/', ...
            num2str(numOutFileCurSeries), '...']);
        
        curOutFileDir = allOutFilesDirs{idxSeries}(idxOutFile);
        [lat, lon, alt, gpsLog] = fetchGpsForOutFileDir(curOutFileDir);
        
        % Get the median RX (lat, lon, alt) for all the GPS samples in this
        % series, which will be needed for the antenna gain calculation.
        [latM, lonM, altM] = fetchMedianGpsForSeriesDir(...
            curOutFileDir.folder);
        
        % Compute b for the calibration line corresponding to the RX gain.
        usrpGain = str2double(gpsLog.rxChannelGain);
        powerShiftsForCali = genCalibrationFct( lsLinesPolysInv, ...
            rxGains, usrpGain);
        
        % Compute path loss (without considering the antenna gain). We will
        % use the amplitude version of thresholdWaveform.m without plots
        % for debugging as the noise eliminiation function.
        noiseEliminationFct = @(waveform) thresholdWaveform(abs(waveform));
        [ pathLossInDb, absPathOutFile] ...
            = computePathLossForOutFileDir(curOutFileDir, ...
            usrpGain, noiseEliminationFct, powerShiftsForCali);
        
        % Fetch the measurement campaign meta records.
        [absCurParDir, curSeries] = fileparts(curOutFileDir.folder);
        idxParDir = find(cellfun(@(d) strcmp(d, absCurParDir), ...
            TX_INFO_LOGS_ABS_PAR_DIRS));
        idxCurSer = str2double(curSeries((length('Series_')+1):end));
        assert(length(idxParDir)==1, ...
            'Error: More than 1 parent folder match with the meta information of the current Series data!');
        curTxInfoLog = TX_INFO_LOGS{idxParDir}(idxCurSer);
        assert(curTxInfoLog.seriesNum==idxCurSer, ...
            'The series index number in the matched Tx info log does not agree with idxCurSer.');
        
        % Compute the antenna gains.
        %    function [ gain ] = computeAntGain(lat0, lon0, alt0, ...
        %        az0, el0, ...
        %         lat, lon, alt, ...
        %        antPatAz, antPatEl)
        txGain = computeAntGain(TX_LAT, TX_LON, TX_HEIGHT_M, ...
            curTxInfoLog.txAz, curTxInfoLog.txEl, ...
            latM, lonM, altM, ...
            pat28GAzNorm, pat28GElNorm);
        rxGain = computeAntGain(latM, lonM, altM, ...
            curTxInfoLog.rxAz, curTxInfoLog.rxEl, ...
            TX_LAT, TX_LON, TX_HEIGHT_M, ...
            pat28GAzNorm, pat28GElNorm);
        
        % Store the results, considering the antenna gains.
        pathLossesWithGpsInfo(pathLossCounter,:) ...
            = [pathLossInDb + txGain + rxGain, lat, lon, latM, lonM];
        absPathsOutFiles{pathLossCounter} = absPathOutFile;
        pathLossCounter = pathLossCounter+1;
    end
end
assert(all(~isnan(pathLossesWithGpsInfo(1:end))));

%% Just for Fun: Get the Maximum Measurable Path Loss

maxMeasurablePathLossInfo = struct(...
    'boltzmannsConst', physconst('Boltzmann'), ...
    'rxTemperatureInKel', 290, ...
    'noiseFigureInDb', 2.4, ...
    'ifBandwidth', 200*10^6, ...
    'snrRequiredInDb', 3, ...
    'txPowerInDbm', 23, ...
    'txGain', 22, ...
    'rxGain', 22, ...
    'maxMeasurablePathLossInDb', nan ...
    );
maxMeasurablePathLossInfo.maxMeasurablePathLossInDb = ...
    - (10*log10(maxMeasurablePathLossInfo.boltzmannsConst ...
    *maxMeasurablePathLossInfo.rxTemperatureInKel...
    *1000) + maxMeasurablePathLossInfo.noiseFigureInDb ...
    +10*log10(maxMeasurablePathLossInfo.ifBandwidth) ...
    + maxMeasurablePathLossInfo.snrRequiredInDb) ...
    + maxMeasurablePathLossInfo.txPowerInDbm ...
    + maxMeasurablePathLossInfo.txGain ...
    + maxMeasurablePathLossInfo.rxGain;

%% Save the Results
disp('    Saving the results...')
% For absPathsOutFiles, convert it to relative paths under the data folder,
% which will already contain enough information we need.
relPathsOutFilesUnderDataFolder = ...
    cellfun(@(p) regexp(p, 'Data[\/\\]([a-zA-Z\d\/\\_]+.out)$', ...
    'tokens'), absPathsOutFiles);
pathPathLossFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'pathLossesWithGpsInfo.mat');
save(pathPathLossFileToSave, ...
    'pathLossesWithGpsInfo', 'relPathsOutFilesUnderDataFolder',...
    'maxMeasurablePathLossInfo');

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

% Plot path losses on map with individual GPS coordinates.
hPathLossesOnMapIndi = figure; hold on; colormap jet;
plot(validPathLossesWithValidGps(:,3), validPathLossesWithValidGps(:,2), 'w.');
plot(pathLossesWithValidGps(boolsInfPathloss,3), ...
    pathLossesWithValidGps(boolsInfPathloss,2), 'kx');
hTx = plot(TX_LON, TX_LAT, '^w', 'MarkerFaceColor', 'b');
plot_google_map('MapType','satellite');
plot3k([validPathLossesWithValidGps(:,3), validPathLossesWithValidGps(:,2), ...
    validPathLossesWithValidGps(:,1)], 'Marker', {'.', 12});
% The command plot_google_map messes up the color legend of plot3k, so we
% will have to fix it here.
hCb = findall( allchild(hPathLossesOnMapIndi), 'type', 'colorbar');
hCb.Ticks = linspace(1,length(colormap),length(hCb.TickLabels));
hold off; grid on; view(0, 90); legend(hTx, 'TX');
title('Path Losses on Map (Large Scale & SIMO)');
xlabel('Lon'); ylabel('Lat'); zlabel('Path Loss (dB)');

% Plot path losses on map with average GPS coordinates.
hPathLossesOnMap = figure; hold on; colormap jet;
plot(validPathLossesWithValidGps(:,5), validPathLossesWithValidGps(:,4), 'w.');
plot(pathLossesWithValidGps(boolsInfPathloss,3), ...
    pathLossesWithValidGps(boolsInfPathloss,2), 'kx');
hTx = plot(TX_LON, TX_LAT, '^w', 'MarkerFaceColor', 'b');
plot_google_map('MapType','satellite');
plot3k([validPathLossesWithValidGps(:,5), validPathLossesWithValidGps(:,4), ...
    validPathLossesWithValidGps(:,1)], 'Marker', {'.', 12});
% The command plot_google_map messes up the color legend of plot3k, so we
% will have to fix it here.
hCb = findall( allchild(hPathLossesOnMap), 'type', 'colorbar');
hCb.Ticks = linspace(1,length(colormap),length(hCb.TickLabels));
hold off; grid on; view(45, 45); legend(hTx, 'TX');
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
pathPathossesOnMapIndiFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'pathLossesOnMapIndi');
saveas(hPathLossesOnMapIndi, [pathPathossesOnMapIndiFileToSave, '.png']);
saveas(hPathLossesOnMapIndi, [pathPathossesOnMapIndiFileToSave, '.fig']);
pathPathossesOnMapFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'pathLossesOnMap');
saveas(hPathLossesOnMap, [pathPathossesOnMapFileToSave, '.png']);
saveas(hPathLossesOnMap, [pathPathossesOnMapFileToSave, '.fig']);
pathPathLossesOverDistFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'pathLossesOverDist');
saveas(hPathLossesOverDist, [pathPathLossesOverDistFileToSave, '.png']);
saveas(hPathLossesOverDist, [pathPathLossesOverDistFileToSave, '.fig']);

disp('    Done!')

% EOF