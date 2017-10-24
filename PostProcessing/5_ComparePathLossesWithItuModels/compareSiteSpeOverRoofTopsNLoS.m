%COMPARESITEGENOVERROOFTOPSNLOS This script will plot the results from the
%ITU site-general model for NLoS propagations over roof tops and compare it
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

% Set this to be true to skip the input range check for the ITU NLoS model.
FLAG_IGNORE_OUT_OF_RANGE = false;

% Set this to be true to show each NLoS location on a map, interactively,
% for manual inspection.
FLAG_MANUAL_INSPECTION = false;

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

%% Load the Category Files

disp(' ')
disp('    Loading category .txt files ...')

catTxtsDirs = rdir(fullfile(ABS_PATH_TO_CATEGORY_TXTS, '*.txt'), ...
    '', false);
catTxts = arrayfun(@(d) loadCategoryTxt(d.name), catTxtsDirs);

disp('    Done!')

%% Extract Valid Path Losses

boolsValidPathlosses ...
    = checkValidityOfPathLossesWithGpsInfo(pathLossesWithGpsInfo, ...
    relPathsOutFilesUnderDataFolder);
validPathLossesWithValidGpsInfo ...
    = pathLossesWithGpsInfo(boolsValidPathlosses,:);
validRelPathsOutFilesUnderDataFolder ...
    = relPathsOutFilesUnderDataFolder(boolsValidPathlosses,:);

% Convert to cell for plotting.
validPLWithValidGPSCell = num2cell(validPathLossesWithValidGpsInfo, 2);
% Compute the TX and RX pair distances. Note that we will use the averaged
% GPS coordinates here.
siteDistsFromTx = cellfun(@(s) ...
    norm([1000.*lldistkm([s(5) s(6)],[TX_LAT,TX_LON]), TX_HEIGHT_M-s(7)]), ...
    validPLWithValidGPSCell);

%% Find NLoS Measurements

disp(' ')
disp('    Searching for NLoS measurements ...')

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
boolsNLoSPathLossRecs = cellfun(@(p) ...
    contains(strrep(strrep(p, '/', '_'), '\', '_'), seriesNLoS), ...
    validRelPathsOutFilesUnderDataFolder);
siteDistsFromTxNLoS = siteDistsFromTx(boolsNLoSPathLossRecs);
validPathLossesWithValidGpsNLoS = ...
    validPathLossesWithValidGpsInfo(boolsNLoSPathLossRecs, 1);
validRelPathsOutFilesNLoS ...
    = validRelPathsOutFilesUnderDataFolder(boolsNLoSPathLossRecs);

% Find all the unique series / site foulders.
validRelPathsSerNLoS = cellfun(@(p) ...
    regexp(p, ...
    '(\d+_[a-zA-Z]+[\\\/]Series_\d+)[\\\/]measureSignal_\d+.out$', ...
    'tokens'), ...
    validRelPathsOutFilesNLoS);
validRelPathsSerNLoS = [validRelPathsSerNLoS{:}]';
validRelPathsSerNLoSUnique = unique(validRelPathsSerNLoS);
validPathLossesWithValidGpsInfoNLoS ...
    = validPathLossesWithValidGpsInfo(boolsNLoSPathLossRecs,:);

disp(' ')
disp('All NLoS sites found:')
disp(' ')
disp(validRelPathsSerNLoSUnique);
disp(' ')
disp('    Done!')

%% Manually Inspection of the NLoS Sites

if FLAG_MANUAL_INSPECTION
    figure; hold on;
    % All the NLoS sites.
    disp(' ')
    disp('    NLoS sites geo info:')
    for idx = 1:length(validRelPathsSerNLoSUnique)
        curRelPath = validRelPathsSerNLoSUnique{idx};
        idxMeasRec = find(strcmp(validRelPathsSerNLoS, curRelPath), 1);
        curMeasRec = validPathLossesWithValidGpsInfoNLoS(idxMeasRec, :);
        curMeasSiteDist = siteDistsFromTxNLoS(idxMeasRec);
        disp([curRelPath, ', NLoSIndex = ', num2str(idx), ...
            ', latM = ', num2str(curMeasRec(5)), ...
            ', lonM = ', num2str(curMeasRec(6)), ...
            ', distTx = ', num2str(curMeasSiteDist)])
        plot(curMeasRec(6), curMeasRec(5), '*', 'Color', ones(1,3)*0.9);
    end
    disp(' ')
    
    plot_google_map('MapType', 'satellite');
    hSiteNLoS = nan;
    numSitesNLoS = length(validRelPathsSerNLoSUnique);
    for idx = 1:numSitesNLoS
        if isgraphics(hSiteNLoS) && isvalid(hSiteNLoS)
            delete(hSiteNLoS);
        end
        curRelPath = validRelPathsSerNLoSUnique{idx};
        idxMeasRec = find(strcmp(validRelPathsSerNLoS, curRelPath), 1);
        curMeasRec = validPathLossesWithValidGpsInfoNLoS(idxMeasRec, :);
        hSiteNLoS = plot(curMeasRec(6), curMeasRec(5), 'r*');
        title(curRelPath, 'Interpreter', 'none');
        disp('Press any key to continue...')
        pause;
    end
    
    disp('Done!')
end

%% Compute the ITU Path Losses

% For easy access, we will just hard code the parameters needed for each
% NLoS sites here. Ref:
%     function [ LNLoS1 ] ...
%        = ituSiteSpecificOverRoofTopsNLoS( fInGHz, dInM, ...
%           hRInM, wInM, ... % bInM,
%          phiInDegree, h1InM, h2InM ... %, lInM
%           )
fInGHz = 28;

% Just copied from the LoS version:
% -------------------------------
% Compute the Gaussian path loss random variables for a distance range.
dsInMStep = 20;
% The defined range for the model.
dsInMRecRange = [55; 1000];
dsInM = dsInMRecRange(1):dsInMStep:dsInMRecRange(2);
% Extend the model to 1 m (out of the defined range).
dsInMStepExt = 1;
dsInMExt = 1:dsInMStepExt:dsInMRecRange(1);
% -------------------------------

PARAMETERS = { ...
    ... % relPathSeries(For Date, type, series #), hRInM, wInM, phiInDegree
    {'20170617_LargeScale\Series_6',                   13,   40,         75}; ...
    {'20170617_LargeScale\Series_8',                   13,   40,         55}; ...
    {'20170617_LargeScale\Series_9',                   13,   40,         60}; ...
    {'20170619_LargeScale\Series_6',                   13,   40,         65}; ...
    {'20170619_LargeScale\Series_8',                   13,   40,         76}; ...
    {'20170620_LargeScale\Series_10',                  13,  100,         56}; ...
    {'20170620_LargeScale\Series_12',                  13,   50,         50}; ...
    {'20170620_LargeScale\Series_13',                  13,   50,         50}; ...
    {'20170620_LargeScale\Series_2',                   13,   85,         45}; ...
    {'20170620_LargeScale\Series_3',                   13,   40,         50}; ...
    {'20170620_LargeScale\Series_4',                   13,   40,       52.5}; ...
    {'20170620_LargeScale\Series_5',                   13,   40,         54}; ...
    {'20170620_LargeScale\Series_8',                   13,  100,         68}; ...
    {'20170621_SIMO\Series_1',                         13,  100,         70}; ...
    {'20170621_SIMO\Series_2',                         13,  100,         55}; ...
    {'20170621_SIMO\Series_3',                         13,  100,         65}; ...
    {'20170621_SIMO\Series_4',                         13,   50,         50}; ...
    {'20170622_LargeScale\Series_1',                   13,   40,         63}; ...
    {'20170622_SIMO\Series_2',                         13,   40,         50}; ...
    {'20170622_SIMO\Series_3',                         13,   85,         45}; ...
    {'20170622_SIMO\Series_4',                         13,   40,         65}; ...
    {'20170622_SIMO\Series_5',                         13,   40,         62}; ...
    {'20170623_SIMO\Series_1',                         13,   40,         78}; ...
    };

numSitesNLoS = length(validRelPathsSerNLoSUnique);
LNLoS1s = nan(numSitesNLoS, 1);
LNLoS1sDInM = nan(numSitesNLoS, 1);
for idx = 1:numSitesNLoS
    curRelPath = validRelPathsSerNLoSUnique{idx};
    assert(strcmp(PARAMETERS{idx}{1}, curRelPath), 'PARAMETERS does not match the found NLoS sites!')
    
    idxMeasRec = find(strcmp(validRelPathsSerNLoS, curRelPath), 1);
    curMeasRec = validPathLossesWithValidGpsInfoNLoS(idxMeasRec, :);
    
    dInM = siteDistsFromTxNLoS(idxMeasRec);
    hRInM = PARAMETERS{idx}{2};
    wInM = PARAMETERS{idx}{3};
    phiInDegree = PARAMETERS{idx}{4};
    h1InM = TX_HEIGHT_M;
    h2InM = curMeasRec(7);
    
    LNLoS1s(idx) = ituSiteSpecificOverRoofTopsNLoS( fInGHz, dInM, ...
        hRInM, wInM, ... % bInM,
        phiInDegree, h1InM, h2InM, ... %, lInM
        FLAG_IGNORE_OUT_OF_RANGE );
    LNLoS1sDInM(idx) = dInM;
end

% Sort the results according to dists.
results = sortrows([LNLoS1s, LNLoS1sDInM], 2);
LNLoS1s = results(:,1);
LNLoS1sDInM = results(:,2);

%% Compute the Root Mean Square

[indicesStarts, indicesEnds, dists] ...
    = findConsecutiveSubSeqs(siteDistsFromTxNLoS);

numResults = length(dists);
rmses = nan(numResults, 1);
errors = nan(length(validPathLossesWithValidGpsNLoS), 1);
for idx = 1:numResults
    curIndices = indicesStarts(idx):indicesEnds(idx);
    errors(curIndices) =  LNLoS1s(idx) ...
        - validPathLossesWithValidGpsNLoS(curIndices);
    rmses(idx) = sqrt(mean(...
        ( errors(curIndices) ).^2 ...
        ));
    
end

% Sort the results according to dists.
results = sortrows([rmses, dists], 2);
rmseItuNLoS = results(:,1);
distsItuNLoS = results(:,2);

rmseAggItuNLoS = sqrt(mean(...
        ( errors ).^2 ...
        ));

%% Regression for Two Reference Models
% The close-in model:
%   PL(d) = PL(d0) + 10*n*log10(d/d0)
% And the Alpha/Beta/Gamma model:
%   PL(d) = 10*alpha*log10(d/d0) + beta + 10*gamma*log10(frequency in GHz)

[nNLoS, closeInModFctNLoS] = fitCloseIn(siteDistsFromTxNLoS, ...
    validPathLossesWithValidGpsNLoS, fInGHz*10^9);
% Use the ITU recommended value for site-general NLoS propagation over
% roof-tops.
gamma0NLoS = 2.3;
[alphaNLoS, betaNLoS, ABGWithGivenGModFctNLoS] = fitAlphaBetaWithGivenGamma(...
    siteDistsFromTxNLoS, ...
    validPathLossesWithValidGpsNLoS, fInGHz, gamma0NLoS);

% For plotting.
dsInMComb = [dsInM, dsInMExt];
pathLossesInDbMeansClo = closeInModFctNLoS(dsInMComb);
pathLossesInDbMeansABG = ABGWithGivenGModFctNLoS(dsInMComb);
% Compute the root mean squared error.
[rmseCloNLoS, distsCloNLoS, rmseAggCloNLoS] = computeRmseOverDist( ...
    siteDistsFromTxNLoS, ...
    closeInModFctNLoS(siteDistsFromTxNLoS), ...
    validPathLossesWithValidGpsNLoS);

[rmseABGNLoS, distsABGNLoS, rmseAggABGNLoS] = computeRmseOverDist( ...
    siteDistsFromTxNLoS, ...
    ABGWithGivenGModFctNLoS(siteDistsFromTxNLoS), ...
    validPathLossesWithValidGpsNLoS);

%% Save the Results
absPathToLoadLoSResults = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'resultsNLoS.mat');
save(absPathToLoadLoSResults, ...
    'siteDistsFromTxNLoS', 'validPathLossesWithValidGpsNLoS', ...
    'validRelPathsSerNLoS', ...
    'nNLoS', 'closeInModFctNLoS', ...
    'alphaNLoS', 'betaNLoS', 'gamma0NLoS', 'ABGWithGivenGModFctNLoS', ...
    'LNLoS1sDInM', 'LNLoS1s', ...
    'rmseItuNLoS', 'distsItuNLoS', ...
    'rmseCloNLoS', 'distsCloNLoS', ...
    'rmseABGNLoS', 'distsABGNLoS');

%% Plot

% The ITU result.
hCompBTLNLoSOnlyWithItuSiteGenOverRoofsLoS = figure; hold on; colormap jet;
% Our basic transission losses.
plot3k([siteDistsFromTxNLoS, zeros(length(siteDistsFromTxNLoS),1), ...
    validPathLossesWithValidGpsNLoS], 'Marker', {'.', 6});
% ITU model results.
yPlaneOnePadding = -ones(length(LNLoS1sDInM),1);
hItuNLoS = plot3(LNLoS1sDInM, yPlaneOnePadding, LNLoS1s, 'k*');
% Fitted close-in model.
yPlaneZeroPaddingComb = zeros(length(dsInMComb),1);
hClo = plot3(dsInMComb, yPlaneZeroPaddingComb, ...
    pathLossesInDbMeansClo, 'r--');
% Fitted Alpha/Beta/Gamma model (with a given gamma).
hABG = plot3(dsInMComb, yPlaneZeroPaddingComb, ...
    pathLossesInDbMeansABG, 'b:');
% Put the parameters for the fitted models on the plots.
text(1000, 0, 70,...
    ['$PL_{Close-in}(d) = PL(d_0) + 10\times', num2str(nNLoS, '%.2f'), ...
    '\times log_{10}(\frac{d}{d_0})$'], ...
    'Rotation', 0, ...
    'Color', 'r', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom', 'Interpreter', 'latex');
text(1000, 0, 68,...
    ['$PL_{\alpha\beta\gamma}(d) =  10\times', num2str(alphaNLoS, '%.2f'), ...
    'log_{10}(\frac{d}{d_0}) +', ...
    num2str(betaNLoS, '%.2f'), ' + 10\times \gamma log_{10}(f)$'], ...
    'Rotation', 0, ...
    'Color', 'b', 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', 'Interpreter', 'latex');
% Put the aggregated root meas square errors on the plot.
text(1.3, 0, 155,...
    {'Root Mean Square Error:', ...
    ['  ITU Mean ', num2str(rmseAggItuNLoS, '%.2f'), ' dB'], ...
    ['  Close-in ', num2str(rmseAggCloNLoS, '%.2f'), ' dB'], ...
    ['  Alpha/Beta/Gamma ', num2str(rmseAggABGNLoS, '%.2f'), ' dB']}, ...
    'Rotation', 0, ...
    'Color', 'k', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', 'Interpreter', 'none');
view(0, 0);
set(gca, 'XScale', 'log'); grid on;
newXTicks = [1,10,100,200,500,1000];
set(gca, 'XTickLabels',newXTicks);
set(gca, 'XTick',newXTicks);
hLegend = legend([hItuNLoS, hClo, hABG], 'ITU NLoS', ...
    'Close-in', 'Alpha/Beta/Gamma', ...
    'Location','northwest');
title('Path Losses over Distance (Large Scale & SIMO, NLoS Only)');
xlabel('Distance to Tx (m)'); ylabel(''); zlabel('Path Loss (dB)');
hold off;

% Root mean errors for these models.
hCompModRmseLoS = figure; hold on;
hRmsqItu = plot(distsItuNLoS, rmseItuNLoS, 'k*');
hRmsqClo = plot(distsCloNLoS, rmseCloNLoS, 'ro');
hRmsqABG = plot(distsABGNLoS, rmseABGNLoS, 'b^');
hold off; grid on;
legend([hRmsqItu, hRmsqClo, hRmsqABG], 'ITU Mean', 'Close-in', ...
    'Alpha/Beta/Gamma');
% set(gca, 'XScale', 'log');
xlabel('Distance to Tx (m)');
ylabel('Mean Square Error for Path Loss (dB)');

%% Also Generate One Aggregate Plot for the RMSE
absPathToLoadLoSResults = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'resultsLoS.mat');
load(absPathToLoadLoSResults);

hCompModRmseAgg = figure; hold on;
hRmsqItuLoS = plot(distsItuLoS, rmseItuLoS, 'b*');
hRmsqCloLoS = plot(distsCloLoS, rmseCloLoS, 'bo');
hRmsqABGLoS = plot(distsABGLoS, rmseABGLoS, 'b^');
hRmsqItuNLoS = plot(distsItuNLoS, rmseItuNLoS, 'r*');
hRmsqCloNLoS = plot(distsCloNLoS, rmseCloNLoS, 'ro');
hRmsqABGNLoS = plot(distsABGNLoS, rmseABGNLoS, 'r^');
hold off; grid on; curAxis = axis;
axis([0,1000,0, min(1.5*curAxis(4), 55)]);
legend([hRmsqItuLoS, hRmsqCloLoS, hRmsqABGLoS, ...
    hRmsqItu, hRmsqClo, hRmsqABG], ...
    'ITU Mean LoS', 'Close-in LoS', 'Alpha/Beta/Gamma LoS', ...
    'ITU NLoS', 'Close-in NLoS', 'Alpha/Beta/Gamma NLoS', ...
    'Location','northwest');
transparentizeCurLegends;
% set(gca, 'XScale', 'log');
xlabel('Distance to Tx (m)');
ylabel('Mean Square Error for Path Loss (dB)');

% Save the plots.
absPathToSavePlots = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'compModelsNLoS');
saveas(hCompBTLNLoSOnlyWithItuSiteGenOverRoofsLoS, [absPathToSavePlots, '.fig']);
saveas(hCompBTLNLoSOnlyWithItuSiteGenOverRoofsLoS, [absPathToSavePlots, '.png']);

absPathToSavePlots = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'compModRmseNLoS');
saveas(hCompModRmseLoS, [absPathToSavePlots, '.fig']);
saveas(hCompModRmseLoS, [absPathToSavePlots, '.png']);

absPathToSavePlots = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'compModRmseAggregated');
saveas(hCompModRmseAgg, [absPathToSavePlots, '.fig']);
saveas(hCompModRmseAgg, [absPathToSavePlots, '.png']);
% EOF