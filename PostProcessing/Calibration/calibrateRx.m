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

% We have two sets of reference data, with the Gnu Radio gain being set to
% 1dB and 76dB, respectively.
rxGains = [1; 76]; % In dB.
% Each set of data is stored in a folder named as "Gain_xxx".
calDataDirNamePrefix = 'Gain_';

% Gnu Radio sample rate.


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
    
    % We will use the filtered Rx output files.
    curCalDataLogs = rdir(curDatasetAbsPath, ...
        'regexp(name, ''(_filtered\.out$)'')');
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

% Sample rate used.
Fs = 1.04 * 10^6;
% Low pass filter for the PSD.
maxFreqPassed = 46000; % In Hz.
% Minimum valid calculated power.
minValidCalPower = -140; % In dB.

pathCalFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'Calibration');

% Compute the Rx Power (in dB) for each filtered output log file. We are
% doing the for loops again here instead of in the section above to better
% separate the code modules.
[calDataThresholded, ... % For debugging.
    calculatedPowers] = deal(cell(numDatasets,1));
for idxDataset = 1:numDatasets
    curNumMeas = length(calData{idxDataset});
    curCalDataThr = cell(curNumMeas,1);
    curCalculatedP = nan(curNumMeas,1);
    for idxCurMeas = 1:curNumMeas
        curSeries = calData{idxDataset}{idxCurMeas};
        % First of all, threshold the I&Q waveforms of each signal vector
        % to eliminate corss-correlation and system noise.
        signalReal = thresholdWaveform(real(curSeries));
        signalImag = thresholdWaveform(imag(curSeries));
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
        hold off;
        title(['Single-Sided Amplitude Spectrum P1 - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)]);
        xlabel('f (Hz)'); ylabel('|P1(f)|'); axis(curAxis);
        legend([hP1, hLPF], 'P1', 'LPF'); grid on;
        % Save the plot.
        pathNewPsdFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['psd-set-',num2str(idxCurMeas),'-pt-',num2str(idxCurMeas)]);
        saveas(hPSD, [pathNewPsdFileToSave, '.fig']);
        saveas(hPSD, [pathNewPsdFileToSave, '.png']);
        
        % Compute the power.
        powerSpectralDen = P2.^2;
        maxIdxPassed = find(f>maxFreqPassed, 1)-1;
        curCalculatedP(idxCurMeas) = ...
            trapz(powerSpectralDen(1:maxIdxPassed));
    end
    calDataThresholded{idxDataset} = curCalDataThr;
    % Change to dB and remove the gain from the Gnu Radio.
    calculatedPowers{idxDataset} = 10.*log(curCalculatedP)./log(10) ...
        - rxGains(idxDataset);
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
    
    % For fitting, remove points with too low calculated power.
    xsToFit = xs(ys>=minValidCalPower);
    ysToFit = ys(ys>=minValidCalPower);
    
    lsLinePoly = polyfit(xsToFit, ysToFit, 1);
    
    xRangeToShow = linspace(finalAxis(1),finalAxis(2));
    valuesLsLine = polyval(lsLinePoly,xRangeToShow);
    hLsLines{idxDataset} = plot(xRangeToShow,valuesLsLine, ...
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
plot([finalAxis(1) finalAxis(2)], [minValidCalPower minValidCalPower], ...
    'r--');
x = [finalAxis(1) finalAxis(1) finalAxis(2) finalAxis(2)];
minY = -200;
y = [minY minValidCalPower minValidCalPower minY];
patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
title('Calibration results');
xlabel('Measured Power (dB)');
ylabel('Calculated Power (dB)');
grid on; hold off;
pathCalFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'Calibration');
saveas(hFigCalibration, [pathCalFileToSave, '.fig']);
saveas(hFigCalibration, [pathCalFileToSave, '.png']);
disp('    Done!')
% EOF