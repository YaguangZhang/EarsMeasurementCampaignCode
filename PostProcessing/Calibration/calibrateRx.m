% CALIBRATERX Calibartion for the received signal power.
%
% Essentially, we need to find the relationship between the "power" of the
% received valtage from Gnu Radio and the true recieved signal power.
%
% Yaguang Zhang, Purdue, 08/14/2017

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

% Configure other paths accordingly.
ABS_PATH_TO_CALI_DATA = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'Data', '20170615_Calibration');
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'Calibration');
FLAG_USE_FILTERED_OUTPUT_FILES = true;

% We have two sets of reference data, with the Gnu Radio gain being set to
% 1dB and 76dB, respectively.
rxGains = [1; 76]; % In dB.
% Each set of data is stored in a folder named as "Gain_xxx".
calDataDirNamePrefix = 'Gain_';

% Reference received power measured by the spectrum analyzer.
measPowers = {[-19;-24;-29;-34;-39;-44;-49;-54;-59], ...
    [-39.6;-44.6;-49.6;-54.6;-59.6;-64.6;-69.6]};

%% Before Calibration

disp(' ------------- ')
disp('  calibrateRx')
disp(' ------------- ')

% Disable the tex interpreter in figures.
set(0,'DefaultTextInterpreter','none');

% Add libs to current path.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

% Make sure the calibration data are there.
disp(' ')
disp('    Checking calibration datasets...')

numDatasets = length(rxGains);
absPathsToCalDatasets = cell(numDatasets,1);
for idxRxGain = 1:numDatasets
    curRxGain = rxGains(idxRxGain);
    curExpFoldername = [calDataDirNamePrefix, num2str(curRxGain)];
    absPathsToCalDatasets{idxRxGain} = fullfile(ABS_PATH_TO_CALI_DATA, ...
        curExpFoldername);
    
    assert(exist( ...
        absPathsToCalDatasets{idxRxGain}, 'dir'...
        )==7, ['Expected folder ',curExpFoldername,' does not exist!']);
end

disp('    Done!')

%% Load the Calibration Datasets

disp(' ')
disp('    Loading calibration data...')

calData = cell(numDatasets,1);
for idxDataset = 1:numDatasets
    curDatasetAbsPath = absPathsToCalDatasets{idxDataset};
    
    if FLAG_USE_FILTERED_OUTPUT_FILES
        % Use the filtered (by GunRadio) Rx output files.
        curCalDataLogs = rdir(curDatasetAbsPath, ...
            'regexp(name, ''(_filtered\.out$)'')');
    else
        % Use the original output files.
        curCalDataLogs = rdir(curDatasetAbsPath, ...
            'regexp(name, ''(_\d+\.out$)'')');
    end
    % Make sure the log files are sorted by name before loading.
    [~, sortedIs] = sort({curCalDataLogs.name});
    curCalDataLogs = curCalDataLogs(sortedIs);
    
    % Load all the calibration data.
    curNumMeas = length(curCalDataLogs);
    curCalData = cell(curNumMeas,1);
    for idxCurMeas = 1:curNumMeas
        curCalData{idxCurMeas} = ...
            read_complex_binary(curCalDataLogs(idxCurMeas).name);
    end
    calData{idxDataset} = curCalData;
end

disp('    Done!')

%% Calibration

disp(' ')
disp('    Calibrating...')

% Sample rate used for GnuRadio.
Fs = 1.04 * 10^6;
% Low pass filter for the PSD.
maxFreqPassed = 46000; % In Hz.
% Number of samples to discard at the beginning.
numFirstSampsToDiscard = 1500;

% % Minimum valid calculated power.
% minValidCalPower = -140; % In dB.

% Minimum valid estimated SNR.
minValidEstSnr = 1.5;

% Compute the Rx Power (in dB) for each filtered output log file. We are
% doing the for loops again here instead of in the section above to better
% separate the code modules.
[calDataThresholded, estimatedSnrs, ... % For debugging.
    calculatedPowers] = deal(cell(numDatasets,1));
for idxDataset = 1:numDatasets
    curNumMeas = length(calData{idxDataset});
    curCalDataThr = cell(curNumMeas,1);
    curCalculatedP = nan(curNumMeas,1);
    curEstimatedSnrs = nan(curNumMeas,1);
    for idxCurMeas = 1:curNumMeas
        curSeries = calData{idxDataset}{idxCurMeas};
        % Discard the first numFirstSampsToDiscard samples.
        curSeries = curSeries((numFirstSampsToDiscard+1):end);
        % First of all, threshold the I&Q waveforms of each signal vector
        % to eliminate corss-correlation and system noise.
        [signalReal, ~, hNoiseSigmaReal] = ...
            thresholdWaveform(real(curSeries), true);
        [signalImag, ~, hNoiseSigmaImag] = ...
            thresholdWaveform(imag(curSeries), true);
        % Compute the complex FFT of the resulted signal; Store it in the
        % cell curCalDataThr for debugging.
        curCalDataThr{idxCurMeas} = signalReal+1i.*signalImag;
        
        % Signal to process.
        X = curCalDataThr{idxCurMeas};
        L = length(X);
        % FFT results.
        Y = fft(X);
        % Two-sided spectrum.
        P2 = abs(Y/L);
        % Single-sided spectrum.
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        % Frequency domain.
        f = Fs*(0:(L/2))/L;
        % Plot the result for debugging.
        hPSD = figure; hold on;
        hP1 = plot(f,P1);
        curAxis = axis; curAxis(2) = f(end);
        hLPF = plot([maxFreqPassed, maxFreqPassed], ...
            [curAxis(3),curAxis(4)], '-.r');
        x = [maxFreqPassed maxFreqPassed f(end) f(end)];
        y = [curAxis(3) curAxis(4) curAxis(4) curAxis(3)];
        patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
        
        % Compute the power.
        powerSpectralDen = P2.^2;
        maxIdxPassed = find(f>maxFreqPassed, 1)-1;
        maxIdxFiltered = find(f>2*maxFreqPassed, 1)-1;
        % Discard DC component.
        curCalculatedP(idxCurMeas) = ...
            trapz(powerSpectralDen(2:maxIdxPassed));
        % As the reference noise power, compute the power from
        % maxFreqPassed to 2*maxFreqPassed.
        powerFiltered = trapz(powerSpectralDen(...
            (maxIdxPassed+1):maxIdxFiltered) ...
            );
        curEstimatedSnrs(idxCurMeas) = ...
            curCalculatedP(idxCurMeas)/powerFiltered;
        
        text(maxFreqPassed, mean([curAxis(3) curAxis(4)]), ...
            ['Estimated SNR = ', ...
            num2str(curEstimatedSnrs(idxCurMeas), '%.2f')]);
        hold off;
        title(['Single-Sided Amplitude Spectrum P1 - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)]);
        xlabel('f (Hz)'); ylabel('|P1(f)|'); axis(curAxis);
        legend([hP1, hLPF], 'P1', 'LPF'); grid on;
        
        % Save the plots.
        pathNewPsdFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), '-psd']);
        saveas(hPSD, [pathNewPsdFileToSave, '.png']);
        pathNewCompNoiseSigmaToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), '-noise-sigma-']);
        saveas(hNoiseSigmaReal, [pathNewCompNoiseSigmaToSave, 'real.png']);
        saveas(hNoiseSigmaImag, [pathNewCompNoiseSigmaToSave, 'imag.png']);
        %         % Also a .fig copy. saveas(hPSD, [pathNewPsdFileToSave,
        %         '.fig']); saveas(hNoiseSigmaReal,
        %         [pathNewCompNoiseSigmaToSave, 'real.fig']);
        %         saveas(hNoiseSigmaImag, [pathNewCompNoiseSigmaToSave,
        %         'imag.fig']);
        
    end
    calDataThresholded{idxDataset} = curCalDataThr;
    % Change to dB and remove the gain from the Gnu Radio.
    calculatedPowers{idxDataset} = 10.*log(curCalculatedP)./log(10) ...
        - rxGains(idxDataset);
    estimatedSnrs{idxDataset} = curEstimatedSnrs;
end

disp('    Done!')

%% Plot

disp(' ')
disp('    Plotting...')

seriesColors = colormap(parula);
[numSeriesColors, ~] = size(seriesColors);
rng(2);
indicesColorToUse = randi([1 numSeriesColors],1,numDatasets);

hFigCalibration = figure; hold on;
calPs = vertcat(calculatedPowers{:});
calPs = calPs(~isinf(calPs));
meaPs = vertcat(measPowers{:});
axisToSet = [min(meaPs) max(meaPs) ...
    min(calPs) max(calPs)];
% Fitted least-squares line.
[hLsLines, lsLinesPolys] = deal(cell(numDatasets,1));
for idxDataset = 1:numDatasets
    xs = measPowers{idxDataset};
    ys = calculatedPowers{idxDataset};
    
    % For plotting, only avoid points with inf as calculated power.
    xsToShow = xs(~isinf(ys));
    ysToShow = ys(~isinf(ys));
    
    % Plot the results.
    colorToUse = seriesColors(indicesColorToUse(idxDataset),:);
    % Non-inf points.
    scatter(xsToShow, ysToShow, '*', 'MarkerEdgeColor', colorToUse, ...
        'LineWidth',1.5);
end
% Set the visible area of the plot now according to the data points shown.
axis(axisToSet); axis equal; finalAxis = axis; axis manual;
% Add the lslines.
for idxDataset = 1:numDatasets
    xs = measPowers{idxDataset};
    ys = calculatedPowers{idxDataset};
    
    %     % For fitting, remove points with too low calculated power.
    %     xsToFit = xs(ys>=minValidCalPower); ysToFit =
    %     ys(ys>=minValidCalPower);
    
    % For fitting, remove points with too low estimated SNR.
    boolsPtsToFit = estimatedSnrs{idxDataset}>=minValidEstSnr;
    xsToFit = xs(boolsPtsToFit);
    ysToFit = ys(boolsPtsToFit);
    
    % Cover the unused points in the plot.
    hIgnoredPts = plot(xs(~boolsPtsToFit),ys(~boolsPtsToFit), ...
        'r*', 'LineWidth',1.5);
    
    % Linear fitting.
    lsLinePoly = polyfit(xsToFit, ysToFit, 1);
    
    % Plot the fitted line.
    xRangeToShow = linspace(finalAxis(1),finalAxis(2));
    valuesLsLine = polyval(lsLinePoly,xRangeToShow);
    hLsLines{idxDataset} = ...
        plot(xRangeToShow,valuesLsLine, ...
        'Color',colorToUse,'LineStyle', '--');
    
    % Show the polynomial on the plot.
    if lsLinePoly(2)>0
        strPoly=['y = ',num2str(lsLinePoly(1)),'x+',num2str(lsLinePoly(2))];
    elseif lsLinePoly(2)<0
        strPoly=['y = ',num2str(lsLinePoly(1)),'x',num2str(lsLinePoly(2))];
    else % lsLinePoly(2)==0
        strPoly=['y = ',num2str(lsLinePoly(1)),'x'];
    end
    idxMiddlePtToFit = floor(length(xsToFit)/2);
    hPolyText = text(xsToFit(idxMiddlePtToFit), ...
        ysToFit(idxMiddlePtToFit), strPoly);
    set(hPolyText,'Rotation',rad2deg(atan(lsLinePoly(1))));
    
    lsLinesPolys{idxDataset} = lsLinePoly;
end
% plot([finalAxis(1) finalAxis(2)], [minValidCalPower minValidCalPower],
% ...
%     'r--');
% x = [finalAxis(1) finalAxis(1) finalAxis(2) finalAxis(2)]; minY = -200; y
% = [minY minValidCalPower minValidCalPower minY];
% patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
title('Calibration results');
xlabel('Measured Power (dB)');
ylabel('Calculated Power (dB)');
grid on; hold off;
pathCalFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'Calibration');
saveas(hFigCalibration, [pathCalFileToSave, '.fig']);
saveas(hFigCalibration, [pathCalFileToSave, '.png']);
disp('    Done!')
% EOF