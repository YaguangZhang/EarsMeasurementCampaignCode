% PATHLOSSESFORCONTITRACKS Load the Conti track measurements and their
% cooresponding GPS logs, and generate for each track a path loss over
% track (on map) plot.
%
% Yaguang Zhang, Purdue, 09/27/2017

clear; clc; close all;

%% Configurations

warning('on');

% Add libs to current path and set ABS_PATH_TO_EARS_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% We also need the functions thresholdWaveform.m for noise elimination, as
% well as genCalibrationFct.m for Rx calibration.
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

disp(' -------------------------- ')
disp('  pathLossesForContiTracks')
disp(' -------------------------- ')

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

%% Search for the Conti Measurement Data Files

disp(' ')
disp('    Searching for measurement data files ...')

% Here we don't care about when the data were collected so let's rearrange
% all the dir struct into one array.
allSeriesDirsArray = vertcat(allSeriesDirs{1:end});
numSeries = length(allSeriesDirsArray);

% Scan the series folder for Gnu Radio out files, as well as the
% corresponding GPS log files.
[contiOutFilesDirs, contiGpsFilesDirs] = deal(cell(numSeries,1));

for idxSeries = 1:numSeries
    disp(['        Scanning series folder ', num2str(idxSeries), ...
        '/', num2str(numSeries), '...']);
    
    % Here we only load the _Conti .out files.
    regexPattern = '\d+_Conti';
    [contiOutFilesDirs{idxSeries}, contiGpsFilesDirs{idxSeries}] = ...
        loadFilesFromSeriesDir(allSeriesDirsArray(idxSeries), regexPattern);
end

disp('    Done!')

%% Compute Path Losses

disp(' ')
disp('    Computing path losses...')

% Compute the path losses and save them into a matrix together with the GPS
% info.
numContiOutFiles = sum(cellfun(@(d) length(d), contiOutFilesDirs));
% More specifically, we will create one cell for each conti track to store
% a matrix, where each row is a [path loss (dB), lat, lon] array for a GPS
% position.
contiPathLossesWithGpsInfo = cell(numContiOutFiles, 1);
% Also save the meta info needed to map the path loss back to the
% measurements. We choose to save the full file path to the .out file for
% convenience.
absPathsContiOutFiles = cell(numContiOutFiles, 1);
% For indexing and storing the results into cells.
contiOutFileCounter = 1;
for idxSeries = 1:numSeries
    disp(['        Processing series ', num2str(idxSeries), '/', ...
        num2str(numSeries), '...']);
    
    numOutFileCurSeries = length(contiOutFilesDirs{idxSeries});
    for idxOutFile = 1:numOutFileCurSeries
        disp(['            Outfile ', num2str(idxOutFile), '/', ...
            num2str(numOutFileCurSeries), '...']);

        curOutFileDir = contiOutFilesDirs{idxSeries}(idxOutFile);
        contiGpsFilesDirs = contiGpsFilesDirs{idxSeries};
        
        % Compute path losses for this track. We will use the amplitude
        % version of thresholdWaveform.m without plots for debugging as the
        % noise eliminiation function.
        noiseEliminationFct = @(waveform) thresholdWaveform(abs(waveform));
        [ curContiPathLossesWithGpsInfo, absPathOutFile] ...
            = computePathLossesForContiOutFileDir(...
            curOutFileDir, contiGpsFilesDirs, ...
            noiseEliminationFct);
        
        % Store the results.
        contiPathLossesWithGpsInfo{contiOutFileCounter} ...
            = curContiPathLossesWithGpsInfo;
        absPathsContiOutFiles{contiOutFileCounter} = absPathOutFile;
        contiOutFileCounter = contiOutFileCounter+1;
    end
end
assert(all(~isempty(contiPathLossesWithGpsInfo{1:end})));

disp('    Saving the results...')
% For absPathsOutFiles, convert it to relative paths under the data folder,
% which will already contain enough information we need.
relPathsContiOutFilesUnderDataFolder = ...
    cellfun(@(p) regexp(p, 'Data[\/\\]([a-zA-Z\d\/\\_]+.out)$', ...
    'tokens'), absPathsContiOutFiles);
pathContiPathLossesFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'contiPathLossesWithGpsInfo.mat');
save(pathContiPathLossesFileToSave, ...
    'contiPathLossesWithGpsInfo', 'relPathsContiOutFilesUnderDataFolder');

disp('    Done!')

%% Plot

disp(' ')
disp('    Plotting...')

numTracks = length(contiPathLossesWithGpsInfo);
for idxTrack = 1:numTracks
    
    curRelPath = cell2mat(relPathsContiOutFilesUnderDataFolder{idxTrack});
    disp(['        Track ', num2str(idxTrack), '/', ...
        num2str(numSeries), ' - ', ...
        curRelPath, '...']);
    % Get the date, series number and timestamp for this track, and use
    % them to construct the filename.
    curFileName = strrep(strrep( ...
        strrep(strrep( curRelPath, '\', '_' ),'/','_'), ...
        'measureSignalCont','Epoch'), '.out', '');    
    
    curPathLossesWithGpsInfo = contiPathLossesWithGpsInfo{idxTrack};
    
    boolsInvalidCoor = curPathLossesWithGpsInfo(:,2)==0 ...
        & curPathLossesWithGpsInfo(:,3)==0;
    if any(boolsInvalidCoor)
        warning([num2str(sum(boolsInvalidCoor)), ...
            ' invalid (lat, lon) pairs detected (both are 0).', ...
            ' We will ignore these points together with their path losses.']);
    end
    pathLossesWithValidGps = curPathLossesWithGpsInfo(~boolsInvalidCoor,:);

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
    hTx = plot(TX_LON, TX_LAT, '^w', 'MarkerFaceColor', 'b');
    plot_google_map('MapType','satellite');
    plot3k([validPathLossesWithValidGps(:,3), validPathLossesWithValidGps(:,2), ...
        validPathLossesWithValidGps(:,1)], 'Marker', {'.', 12});
    % The command plot_google_map messes up the color legend of plot3k, so we
    % will have to fix it here.
    hCb = findall( allchild(hPathLossesOnMap), 'type', 'colorbar');
    hCb.Ticks = linspace(1,length(colormap),length(hCb.TickLabels));
    hold off; grid on; view(0, 90); legend(hTx, 'TX');
    curTitleLabel = strrep(curFileName, '_', '-');
    title({'Path Losses on Map (Conti)', curTitleLabel});
    xlabel('Lon'); ylabel('Lat'); zlabel('Path Loss (dB)');

    % Save the plots.
    pathPathossesOnMapFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
        curFileName);
    saveas(hPathLossesOnMap, [pathPathossesOnMapFileToSave, '.png']);
    saveas(hPathLossesOnMap, [pathPathossesOnMapFileToSave, '.fig']);
end

disp('    Done!')

% EOF