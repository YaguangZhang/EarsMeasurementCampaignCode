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
% The defined range for the model.
dsInMRecRange = [55; 1000];
dsInM = dsInMRecRange(1):dsInMStep:dsInMRecRange(2);

ituModFctLoS = @(ds) arrayfun(@(d) ...
    ituSiteGeneralOverRoofTopsLoS( fInGHz, d ), ds);
ituPathLossesInDbGaussianLoS = ituModFctLoS(dsInM);

% Extract the means and variances for plotting.
pathLossesInDbMeans = arrayfun(@(p) ...
    p.pathLossInDbMean, ituPathLossesInDbGaussianLoS)';
pathLossesInDbVars = arrayfun(@(p) ...
    p.pathLossInDbVar, ituPathLossesInDbGaussianLoS)';

% Extend the model to 1 m (out of the defined range).
dsInMStepExt = 1;
dsInMExt = 1:dsInMStepExt:dsInMRecRange(1);

pathLossesInDbGaussianExt = arrayfun(@(d) ...
    ituSiteGeneralOverRoofTopsLoS( fInGHz, d ), dsInMExt);

% Extract the means and variances for plotting.
pathLossesInDbMeansExt = arrayfun(@(p) ...
    p.pathLossInDbMean, pathLossesInDbGaussianExt)';
pathLossesInDbVarsExt = arrayfun(@(p) ...
    p.pathLossInDbVar, pathLossesInDbGaussianExt)';

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
siteDistsFromTx = cellfun(@(s) ...
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
siteDistsFromTxLoS = siteDistsFromTx(boolsLoSPathLossRecs);
validPathLossesWithValidGpsLoS = ...
    validPathLossesWithValidGps(boolsLoSPathLossRecs, 1);
disp('    Done!')

%% Regression for Two Reference Models
% The close-in model:
%   PL(d) = PL(d0) + 10*n*log10(d/d0)
% And the Alpha/Beta/Gamma model:
%   PL(d) = 10*alpha*log10(d/d0) + beta + 10*gamma*log10(frequency in GHz)

[nLoS, closeInModFctLoS] = fitCloseIn(siteDistsFromTxLoS, ...
    validPathLossesWithValidGpsLoS, fInGHz*10^9);
% Use the ITU recommended value for site-general LoS propagation over
% roof-tops.
gamma0 = 1.96;
[alphaLoS, bettaLoS, ABGWithGivenGModFctLoS] = fitAlphaBetaWithGivenGamma(...
    siteDistsFromTxLoS, ...
    validPathLossesWithValidGpsLoS, fInGHz, gamma0);

% For plotting.
dsInMComb = [dsInM, dsInMExt];
pathLossesInDbMeansClo = closeInModFctLoS(dsInMComb);
pathLossesInDbMeansABG = ABGWithGivenGModFctLoS(dsInMComb);
% Compute the root mean squared error.

[rmseItuLoS, distsItuLoS] = computeRmseOverDist( ...
    siteDistsFromTxLoS, ...
    arrayfun(@(p) p.pathLossInDbMean, ituModFctLoS(siteDistsFromTxLoS)), ...
    validPathLossesWithValidGpsLoS);

[rmseCloLoS, distsCloLoS] = computeRmseOverDist( ...
    siteDistsFromTxLoS, ...
    closeInModFctLoS(siteDistsFromTxLoS), ...
    validPathLossesWithValidGpsLoS);

[rmseABGLoS, distsABGLoS] = computeRmseOverDist( ...
    siteDistsFromTxLoS, ...
    ABGWithGivenGModFctLoS(siteDistsFromTxLoS), ...
    validPathLossesWithValidGpsLoS);

%% Save the Results
absPathToSaveResults = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'resultsLoS');
save(absPathToSaveResults, ...
    'siteDistsFromTxLoS', 'validPathLossesWithValidGpsLoS', ...
    'ituModFctLoS', ...
    'nLoS', 'closeInModFctLoS', ...
    'alphaLoS', 'bettaLoS', 'gamma0', 'ABGWithGivenGModFctLoS', ...
    'rmseItuLoS', 'distsItuLoS', ...
    'rmseCloLoS', 'distsCloLoS', ...
    'rmseABGLoS', 'distsABGLoS');

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
newXTicks = [1,10,100,200,500,1000];
set(gca, 'XTickLabels',newXTicks);
set(gca, 'XTick',newXTicks);
hLegend = legend([hMean, h3Sigma], 'Mean', '3 sigma range');
title('ITU Reference Model');
hold off;

% Plot path losses over distance from Tx (the averaged GPS coordinates are
% used), then add the ITU result onto it.
hCompBTLWithItuSiteGenOverRoofsLoS = figure; hold on; colormap jet;
% Our basic transission losses.
plot3k([siteDistsFromTx, zeros(length(siteDistsFromTx),1), ...
    validPathLossesWithValidGps(:,1)], 'Marker', {'.', 6});
curAxis = axis;
axis([min([siteDistsFromTx; 1]), max(siteDistsFromTx)+100, curAxis(3:6)]);
% ITU model results.
yPlaneZeroPadding = zeros(length(dsInM),1);
hMean = plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans, 'k-');
h3Sigma = plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans + 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7);
plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans - 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7);
view(0, 0);
set(gca, 'XScale', 'log'); grid on;
newXTicks = [1,10,100,200,500,1000];
set(gca, 'XTickLabels',newXTicks);
set(gca, 'XTick',newXTicks);
hLegend = legend([hMean, h3Sigma], 'Mean', '3 sigma range', ...
    'Location','northwest');
title('Path Losses over Distance (Large Scale & SIMO)');
title(hLegend, 'ITU Reference Model');
xlabel('Distance to Tx (m)'); ylabel(''); zlabel('Path Loss (dB)');
hold off;

% The same for only the LoS cases.
hCompBTLLoSOnlyWithItuSiteGenOverRoofsLoS = figure; hold on; colormap jet;
% Our basic transission losses.
plot3k([siteDistsFromTxLoS, zeros(sum(boolsLoSPathLossRecs),1), ...
    validPathLossesWithValidGpsLoS], 'Marker', {'.', 6});
curAxis = axis;
axis([min([siteDistsFromTx; 1]), max(siteDistsFromTx)+100, curAxis(3:6)]);
% ITU model results.
yPlaneZeroPadding = zeros(length(dsInM),1);
hMean = plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans, 'k-');
h3Sigma = plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans + 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7);
plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans - 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7);
view(0, 0);
set(gca, 'XScale', 'log'); grid on;
newXTicks = [1,10,100,200,500,1000];
set(gca, 'XTickLabels',newXTicks);
set(gca, 'XTick',newXTicks);
hLegend = legend([hMean, h3Sigma], 'Mean', '3 sigma range', ...
    'Location','northwest');
title('Path Losses over Distance (Large Scale & SIMO, LoS Only)');
title(hLegend, 'ITU Reference Model');
xlabel('Distance to Tx (m)'); ylabel(''); zlabel('Path Loss (dB)');
hold off;

% Update: compare with two other models. The same for only the LoS cases.
hCompModelsLoS = figure; hold on; colormap jet;
% Our basic transission losses.
plot3k([siteDistsFromTxLoS, zeros(sum(boolsLoSPathLossRecs),1), ...
    validPathLossesWithValidGpsLoS], 'Marker', {'.', 6});
curAxis = axis;
axis([min([siteDistsFromTx; 1]), max(siteDistsFromTx)+100, curAxis(3:6)]);
% ITU model results.
ituLineWidth = 1.5;
yPlaneZeroPadding = zeros(length(dsInM),1);
hMean = plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans, 'k-', ...
    'LineWidth', ituLineWidth);
h3Sigma = plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans + 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7, ...
    'LineWidth', ituLineWidth);
plot3(dsInM, yPlaneZeroPadding, pathLossesInDbMeans - 3.*pathLossesInDbVars, '-.', ...
    'Color', ones(1,3).*0.7, ...
    'LineWidth', ituLineWidth);
% Extended ITU model results.
yPlaneZeroPaddingExt = zeros(length(dsInMExt),1);
hMeanExt = plot3(dsInMExt, yPlaneZeroPaddingExt, pathLossesInDbMeansExt, 'k-.');
h3SigmaExt = plot3(dsInMExt, yPlaneZeroPaddingExt, ...
    pathLossesInDbMeansExt + 3.*pathLossesInDbVarsExt, '-.', ...
    'Color', ones(1,3).*0.7);
plot3(dsInMExt, yPlaneZeroPaddingExt, pathLossesInDbMeansExt ...
    - 3.*pathLossesInDbVarsExt, '-.', ...
    'Color', ones(1,3).*0.7);
% Fitted close-in model.
yPlaneZeroPaddingComb = zeros(length(dsInMComb),1);
hClo = plot3(dsInMComb, yPlaneZeroPaddingComb, ...
    pathLossesInDbMeansClo, 'r--');
% Fitted Alpha/Beta/Gamma model (with a given gamma).
hABG = plot3(dsInMComb, yPlaneZeroPaddingComb, ...
    pathLossesInDbMeansABG, 'b:');
view(0, 0);
set(gca, 'XScale', 'log'); grid on;
newXTicks = [1,10,100,200,500,1000];
set(gca, 'XTickLabels',newXTicks);
set(gca, 'XTick',newXTicks);
hLegend = legend([hMean, h3Sigma, hMeanExt, hClo, hABG], ...
    'ITU Mean', 'ITU 3-sigma Range', 'ITU Extrapolated', ...
    'Close-in', 'Alpha/Beta/Gamma', ...
    'Location', 'northwest');
title('Path Losses over Distance (Large Scale & SIMO, LoS Only)');
title(hLegend, 'Reference Models');
xlabel('Distance to Tx (m)'); ylabel(''); zlabel('Path Loss (dB)');
hold off;

% Root mean errors for these models.
hCompModRmseLoS = figure; hold on;
hRmsqItu = plot(distsItuLoS, rmseItuLoS, 'k*');
hRmsqClo = plot(distsCloLoS, rmseCloLoS, 'ro');
hRmsqABG = plot(distsABGLoS, rmseABGLoS, 'b^');
hold off; grid on;
legend([hRmsqItu, hRmsqClo, hRmsqABG], 'ITU Mean', 'Close-in', ...
    'Alpha/Beta/Gamma');
% set(gca, 'XScale', 'log');
xlabel('Distance to Tx (m)');
ylabel('Mean Square Error for Path Loss (dB)');

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

absPathToSavePlots = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'compModelsLoS');
saveas(hCompModelsLoS, [absPathToSavePlots, '.fig']);
saveas(hCompModelsLoS, [absPathToSavePlots, '.png']);

absPathToSavePlots = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'compModRmseLoS');
saveas(hCompModRmseLoS, [absPathToSavePlots, '.fig']);
saveas(hCompModRmseLoS, [absPathToSavePlots, '.png']);

% EOF