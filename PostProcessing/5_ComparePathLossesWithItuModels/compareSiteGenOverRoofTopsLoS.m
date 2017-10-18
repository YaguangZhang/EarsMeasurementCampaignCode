%COMPARESITEGENOVERROOFTOPSLOS This script will plot the results from the
%ITU site-general model for LoS propagations over roof tops and compare it
%with what we got from the measurement campaign.
%
% Yaguang Zhang, Purdue, 10/17/2017

clear; clc; close all; 

%% Configurations

warning('on');

% Add libs to current path and set ABS_PATH_TO_EARS_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
% We will need later the category folder for .txt files.
ABS_PATH_TO_CATEGORY_TXTS = fullfile(pwd, ...
    '..', '4_1_PlotPathLossesByCategory', 'Categories');
cd('..'); setPath;

% Configure other paths accordingly.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'ComparePathLossesWithItuModels');

% Reuse results from loadMeasCampaignInfo.m, evalPathLosses.m.
ABS_PATH_TO_TX_INFO_LOGS_FILE= fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation', 'txInfoLogs.mat');
ABS_PATH_TO_PATH_LOSSES_FILE= fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation', ...
    'pathLossesWithGpsInfo.mat');

%% Before Processing the Data

disp(' ------------------------------- ')
disp('  compareSiteGenOverRoofTopsLoS ')
disp(' ------------------------------- ')

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

%% Get Info for Measurement Data Files and Calibration Polynomials

disp(' ')
disp('    Loading results from: ')
disp('      - loadMeasCampaignInfo.m')
disp('      - evalPathLosses.m')

assert(exist(ABS_PATH_TO_TX_INFO_LOGS_FILE, 'file')==2, ...
    'Couldn''t find txInfoLogs.mat! Please run PostProcessing/4_0_PathLossComputation/loadMeasCampaignInfo.m first.');
assert(exist(ABS_PATH_TO_PATH_LOSSES_FILE, 'file')==2, ...
    'Couldn''t find pathLossesWithGpsInfo.mat! Please run PostProcessing/4_0_PathLossComputation/evalPathLosses.m first.');

% The data have been processed before and the result files have been found.
disp('    Found all .mat files required.');
disp('        Loading the results...')

% Get records of the TxInfo.txt files (among other contant parameters for
% the measurement campaign, e.g. F_S, TX_LAT, TX_LON, and TX_POWER_DBM):
% 'TX_INFO_LOGS' and 'TX_INFO_LOGS_ABS_PAR_DIRS'.
load(ABS_PATH_TO_TX_INFO_LOGS_FILE);
% Get 'pathLossesWithGpsInfo', 'relPathsOutFilesUnderDataFolder', and
% 'maxMeasurablePathLossInfo'.
load(ABS_PATH_TO_PATH_LOSSES_FILE);

disp('    Done!')

%% Compute the ITU Path Losses

fInGHz = 28;

% Compute the Gaussian path loss random variables for a distance range.
dsInMStep = 20;
dsInM = 55:dsInMStep:1000;

pathLossesInDbGaussian = arrayfun(@(d) ...
    ituSiteGeneralOverRoofTopsLoS( fInGHz, d ), dsInM);

% Extract the means and variances for plotting.
pathLossesInDbMeans = arrayfun(@(p) ...
    p.pathLossInDbMean, pathLossesInDbGaussian)';
pathLossesInDbVars = arrayfun(@(p) ...
    p.pathLossInDbVar, pathLossesInDbGaussian)';

%% Extract Valid Path Losses

boolsValidPathlosses ...
    = checkValidityOfPathLossesWithGpsInfo(pathLossesWithGpsInfo, ...
    relPathsOutFilesUnderDataFolder);
validPathLossesWithValidGps ...
    = pathLossesWithGpsInfo(boolsValidPathlosses,:);
validRelPathsOutFilesUnderDataFolder ...
    = relPathsOutFilesUnderDataFolder(boolsValidPathlosses,:);

% Convert to cell for plotting.
validPLWithValidGPSCell = num2cell(validPathLossesWithValidGps, 2);
% Compute the TX and RX pair distances. Note that we will use the averaged
% GPS coordinates here.
distsFromTx = cellfun(@(s) ...
    norm([1000.*lldistkm([s(5) s(6)],[TX_LAT,TX_LON]), TX_HEIGHT_M-s(7)]), ...
    validPLWithValidGPSCell);

%% Load the Category Files

disp(' ')
disp('    Loading category .txt files ...')

catTxtsDirs = rdir(fullfile(ABS_PATH_TO_CATEGORY_TXTS, '*.txt'), ...
    '', false);
catTxts = arrayfun(@(d) loadCategoryTxt(d.name), catTxtsDirs);

disp('    Done!')

%% Find LoS Measurements

disp(' ')
disp('    Searching for LoS measurements ...')

% Here we will find the NLoS cases and exlude them from the whole
% measurement data set.
NLoSCategories = {...
    'Totally Blocked', ...
    'Totally Blocked by Buildings', ...
    'Totally Blocked by Foliage'};

catTxtsNLos = catTxts(ismember({catTxts.category}, NLoSCategories));
seriesNLoS = {catTxtsNLos.series}';
seriesNLoS = unique(vertcat(seriesNLoS{:}));

% This gives whether a row record in validPathLossesWithValidGps is LoS or
% not.
boolsLoSPathLossRecs = cellfun(@(p) ...
    ~contains(strrep(strrep(p, '/', '_'), '\', '_'), seriesNLoS), ...
    validRelPathsOutFilesUnderDataFolder);

disp('    Done!')

%% Plot

% The ITU result.
hResultItuSiteGenOverRoofsLoS = figure;
hold on;
hMean = plot(dsInM, pathLossesInDbMeans, 'k-');
h3Sigma = plot(dsInM, pathLossesInDbMeans + 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7);
plot(dsInM, pathLossesInDbMeans - 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7);
set(gca, 'XScale', 'log'); grid on;
curAxis = axis; axis([min(dsInM), max(dsInM), curAxis(3:4)]);
newXTicks = [10,100,200,500,1000];
set(gca, 'XTickLabels',newXTicks);
set(gca, 'XTick',newXTicks);
hLegend = legend([hMean, h3Sigma], 'Mean', '3 sigma range');
title('ITU Reference Model');
hold off;

% Plot path losses over distance from Tx (the averaged GPS coordinates are
% used), then add the ITU result onto it.
hCompBTLWithItuSiteGenOverRoofsLoS = figure; hold on; colormap jet;
% Our basic transission losses.
plot3k([distsFromTx, zeros(length(distsFromTx),1), ...
    validPathLossesWithValidGps(:,1)], 'Marker', {'.', 6});
curAxis = axis;
axis([min(distsFromTx)-10, max(distsFromTx)+100, curAxis(3:6)]);
% ITU model results.
yPlaneZeroPadding = zeros(length(dsInM),1);
hMean = plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans, 'k-');
h3Sigma = plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans + 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7);
plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans - 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7);
view(0, 0); 
set(gca, 'XScale', 'log'); grid on;
newXTicks = [10,100,200,500,1000];
set(gca, 'XTickLabels',newXTicks);
set(gca, 'XTick',newXTicks);
hLegend = legend([hMean, h3Sigma], 'Mean', '3 sigma range', ...
    'Location','southeast');
title('Path Losses over Distance (Large Scale & SIMO)');
title(hLegend, 'ITU Reference Model');
xlabel('Distance to Tx (m)'); ylabel(''); zlabel('Path Loss (dB)');
hold off;

% The same for only the LoS cases.
hCompBTLLoSOnlyWithItuSiteGenOverRoofsLoS = figure; hold on; colormap jet;
% Our basic transission losses.
plot3k([distsFromTx(boolsLoSPathLossRecs), zeros(sum(boolsLoSPathLossRecs),1), ...
    validPathLossesWithValidGps(boolsLoSPathLossRecs,1)], 'Marker', {'.', 6});
curAxis = axis;
axis([min(distsFromTx)-10, max(distsFromTx)+100, curAxis(3:6)]);
% ITU model results.
yPlaneZeroPadding = zeros(length(dsInM),1);
hMean = plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans, 'k-');
h3Sigma = plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans + 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7);
plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans - 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7);
view(0, 0); 
set(gca, 'XScale', 'log'); grid on;
newXTicks = [10,100,200,500,1000];
set(gca, 'XTickLabels',newXTicks);
set(gca, 'XTick',newXTicks);
hLegend = legend([hMean, h3Sigma], 'Mean', '3 sigma range', ...
    'Location','southeast');
title('Path Losses over Distance (Large Scale & SIMO, LoS Only)');
title(hLegend, 'ITU Reference Model');
xlabel('Distance to Tx (m)'); ylabel(''); zlabel('Path Loss (dB)');
hold off;

% Save the plots.
absPathToSavePlots = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'ituSiteGenOverRoofsLoS');
saveas(hResultItuSiteGenOverRoofsLoS, [absPathToSavePlots, '.fig']);
saveas(hResultItuSiteGenOverRoofsLoS, [absPathToSavePlots, '.png']);

absPathToSavePlots = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'compBTLWithItuSiteGenOverRoofsLoS');
saveas(hCompBTLWithItuSiteGenOverRoofsLoS, [absPathToSavePlots, '.fig']);
saveas(hCompBTLWithItuSiteGenOverRoofsLoS, [absPathToSavePlots, '.png']);

absPathToSavePlots = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'compBTLLoSOnlyWithItuSiteGenOverRoofsLoS');
saveas(hCompBTLLoSOnlyWithItuSiteGenOverRoofsLoS, [absPathToSavePlots, '.fig']);
saveas(hCompBTLLoSOnlyWithItuSiteGenOverRoofsLoS, [absPathToSavePlots, '.png']);

% EOF