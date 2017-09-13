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
% Set this to be true to use the filtered (by GunRadio) outputs; Otherwise,
% the original samples without post-processing will be used.
FLAG_USE_FILTERED_OUTPUT_FILES = true;

% We have two sets of reference data, with the Gnu Radio gain being set to
% 1dB and 76dB, respectively.
rxGains = [1; 76]; % In dB.
% Each set of data is stored in a folder named as "Gain_xxx".
calDataDirNamePrefix = 'Gain_';

% Reference received power measured by the spectrum analyzer.
measPowers = {[-19;-24;-29;-34;-39;-44;-49;-54;-59], ...
    [-39.6;-44.6;-49.6;-54.6;-59.6;-64.6;-69.6]};

% Set this to be true to eliminate noise (i.e. compute noise sigma) using
% signal amplitude; Otherwise, we will eliminate the noise in the real and
% imaginary parts of the signal separately.
FLAG_NOISE_ELI_VIA_AMP = true;

% Manually ignore some of the measurements.
BOOLS_MEAS_TO_FIT = {[1 1 1 1 1 1 0 0 0], ...
    [1 1 1 1 1 1 1]};

% Sample rate used for GnuRadio.
Fs = 1.04 * 10^6;
% Low pass filter for the PSD. Tried before: 46000; 39500.
maxFreqPassed = 30000; % In Hz.

% Minimum valid calculated power.
minValidCalPower = -inf; % In dB. Before: -140.

% Minimum valid estimated SNR.
minValidEstSnr = 0; % Before: 1.5.

% Number of samples to discard at the beginning.
numStartSampsToDiscard = 100000; % ~0.1s
% After discarding these samples, furthermore only keep the middle part of
% the signal for calibration.
timeLengthAtCenterToUse = 3; % In second.

% The 76dB gain dataset (set #1 needs more noise elimination).
NUMS_SIGMA_FOR_THRESHOLD = [3.5, 3.5];

% Set this to true if you want figures to be generated silently.
FLAG_GEN_PLOTS_SILENTLY = true;

% For any figure generated, a .pgn screen shot is always saved; Set this to
% be true to also save a .fig version (which may dramatically slow down the
% program).
FLAG_SAVE_FIG_COPY = false;

%% Before Calibration

disp(' ------------- ')
disp('  calibrateRx')
disp(' ------------- ')

% Disable the tex interpreter in figures.
set(0,'DefaultTextInterpreter','none');
% Hide figures if it's required.
if FLAG_GEN_PLOTS_SILENTLY
    set(0,'DefaultFigureVisible','off');
else
    set(0,'DefaultFigureVisible','on');
end

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
disp('    Verifying calibration datasets...')

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
    
    NUM_SIGMA_FOR_THRESHOLD = NUMS_SIGMA_FOR_THRESHOLD(idxDataset);
    
    for idxCurMeas = 1:curNumMeas
        % Display the progress.
        disp(['        Set: ', num2str(idxDataset),'/',...
            num2str(numDatasets), ...
            '    Measurement: ', num2str(idxCurMeas),'/',...
            num2str(curNumMeas)]);
        
        % Prepare the samples for calibration.
        curSeries = calData{idxDataset}{idxCurMeas};
        % Discard the first numStartSampsToDiscard of samples.
        curSeries = curSeries((numStartSampsToDiscard+1):end);
        % Further more, only keep the middle part for calibration.
        numSampsToKeep = ceil(timeLengthAtCenterToUse*Fs);
        numSampsCurSeries = length(curSeries);
        if numSampsToKeep > numSampsCurSeries
            warning('There are not enough samples to keep. We will use all remaining ones.');
        else
            idxRangeToKeep = floor(0.5.*numSampsCurSeries ...
                + [-1,1].*numSampsToKeep./2);
            curSeries = curSeries(idxRangeToKeep(1):idxRangeToKeep(2));
        end
        % Make sure we end up with even number of samples.
        if mod(length(curSeries),2)==1
            curSeries = curSeries(1:(end-1));
        end
        
        % Before calibration, plot the PSD of the input samples (after
        % range limitation but before noise elimination).
        X0 = curSeries;
        L0 = length(X0);
        % Apply a L-point minimum 4-term Blackman-Harris window to the
        % signal before fft.
        X0Win = blackmanharris(L0).*X0;
        % FFT results.
        Y0 = fftshift(fft(X0));
        Y0Win = fftshift(fft(X0Win));
        % Frequency domain.
        f0 = (-L0/2:L0/2-1)*(Fs/L0);
        idxDC0 = L0/2+1;
        % PSD.
        powerSpectralDen0 = abs(Y0).^2/L0;
        powerSpectralDen0Win = abs(Y0Win).^2/L0;
        
        % Plot the result in dB for debugging.
        hPSDInputInDb = plotPsdInDbWithLPF(powerSpectralDen0, f0, ...
            maxFreqPassed);
        title({['Estimated Power Spectrum Density (dB) - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)], ...
            'Before Noise Eliminiation'});
        hPSDInputInDbWin = plotPsdInDbWithLPF(powerSpectralDen0Win, f0, ...
            maxFreqPassed);
        title({['Estimated Power Spectrum Density (dB) - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)], ...
            'Before Noise Eliminiation, Windowed'});
        
        % For calibration, first of all, threshold the I&Q waveforms of
        % each signal vector to eliminate corss-correlation and system
        % noise.
        [signalReal, ~, hNoiseSigmaReal] = ...
            thresholdWaveform(real(curSeries), true);
        [signalImag, ~, hNoiseSigmaImag] = ...
            thresholdWaveform(imag(curSeries), true);
        % Update: it makes more sense to eliminate pts according to the
        % amplitude of the signal, instead of doing it separately for the
        % real and image parts.
        [~, boolsEliminatedPts, hNoiseSigmaAmp] = ...
            thresholdWaveform(abs(curSeries), true);
        curSeriesEliminated = curSeries;
        curSeriesEliminated(boolsEliminatedPts) = 0;
        
        % Depending on the flag, choose to use which version of
        % noise-eliminated signal (via real and imaginary parts separately,
        % or via the amplitude as a whole).
        if FLAG_NOISE_ELI_VIA_AMP
            curCalDataThr{idxCurMeas} = curSeriesEliminated;
        else
            % Compute the complex FFT of the resulted signal; Store it in
            % the cell curCalDataThr for debugging.
            curCalDataThr{idxCurMeas} = signalReal+1i.*signalImag;
        end
        
        % Signal (noise eliminiated) to process.
        X = curCalDataThr{idxCurMeas};
        L = length(X);
        
        % Apply a L-point minimum 4-term Blackman-Harris window to the
        % signal before fft.
        %  X = blackmanharris(L).*X;
        % FFT results.
        Y = fftshift(fft(X));
        % Frequency domain.
        f = (-L/2:L/2-1)*(Fs/L);
        idxDC = L/2+1;
        % PSD.
        powerSpectralDen = abs(Y).^2/L;
        
        % Plot the result for debugging.
        hPSD = figure; hold on;
        hPowerSpectralDen = plot(f,powerSpectralDen);
        curAxis = axis; curAxis(1) = f(1); curAxis(2) = f(end);
        if ~isinf(maxFreqPassed)
            hLPF = plot([maxFreqPassed, -maxFreqPassed; ...
                maxFreqPassed, -maxFreqPassed], ...
                [curAxis(3),curAxis(3); ...
                curAxis(4),curAxis(4)], '-.r');
            x = [maxFreqPassed maxFreqPassed f(end) f(end)];
            y = [curAxis(3) curAxis(4) curAxis(4) curAxis(3)];
            patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
            x = [-maxFreqPassed -maxFreqPassed f(1) f(1)];
            y = [curAxis(3) curAxis(4) curAxis(4) curAxis(3)];
            patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
        end
        
        % Compute the power.
        indicesFPassed = find(abs(f)<=maxFreqPassed);
        % Compute the power by integral. Note that we will discard the DC
        % component here (will be passed by the LPF).
        psdPassed = powerSpectralDen;
        psdPassed(idxDC) = 0;
        psdPassed = psdPassed(indicesFPassed);
        curCalculatedP(idxCurMeas) = trapz(f(indicesFPassed), psdPassed);
        % For the noise, we compute the power right outside of the LPF but
        % limit the integral range to be as wide as the LPF. That is, as
        % the reference noise power, compute the power from maxFreqPassed
        % to 2*maxFreqPassed in both the positive and negative parts.
        indicesFFiltered = setdiff(find(abs(f)<=2*maxFreqPassed), ...
            indicesFPassed);
        psdFiltered = powerSpectralDen(indicesFFiltered);
        powerFiltered = trapz(f(indicesFFiltered), psdFiltered);
        curEstimatedSnrs(idxCurMeas) = ...
            curCalculatedP(idxCurMeas)/powerFiltered;
        
        text(min(maxFreqPassed, 50000), mean([curAxis(3) curAxis(4)]), ...
            ['Estimated SNR = ', ...
            num2str(curEstimatedSnrs(idxCurMeas), '%.2f')]);
        hold off;
        title(['Estimated Power Spectrum Density - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)]);
        xlabel('f (Hz)'); ylabel('Estimated PSD (V^2/Hz)'); axis(curAxis);
        if ~isinf(maxFreqPassed)
            legend([hPowerSpectralDen, hLPF(1)], 'PSD', 'LPF');
        else
            legend(hPowerSpectralDen, 'PSD');
        end
        grid minor;
        
        % The same PSD plot in dB.
        hPSDdB = plotPsdInDbWithLPF(powerSpectralDen, f, maxFreqPassed, ...
            curEstimatedSnrs(idxCurMeas));
        title(['Estimated Power Spectrum Density (dB) - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)]);
        
        % Paths to save the plots.
        pathInputPSDdBFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), ...
            '-psd-input-dB']);
        pathWindowedInputPSDdBFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), ...
            '-psd-input-windowed-dB']);
        pathNewPSDFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), '-psd']);
        pathNewPSDdBFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), '-psd-dB']);
        pathNewCompNoiseSigmaToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), '-noise-sigma-']);
        % A .png figure for easy access.
        saveas(hPSDInputInDb, [pathInputPSDdBFileToSave, '.png']);
        saveas(hPSDInputInDbWin, [pathWindowedInputPSDdBFileToSave, ...
            '.png']);
        saveas(hPSD, [pathNewPSDFileToSave, '.png']);
        saveas(hPSDdB, [pathNewPSDdBFileToSave, '.png']);
        % Also a .fig copy.
        if FLAG_SAVE_FIG_COPY
            saveas(hPSDInputInDb, [pathInputPSDdBFileToSave, '.fig']);
            saveas(hPSDInputInDbWin, [pathWindowedInputPSDdBFileToSave, ...
            '.fig']);
            saveas(hPSD, [pathNewPSDFileToSave, '.fig']);
            saveas(hPSDdB, [pathNewPSDdBFileToSave, '.fig']);
        end
        
        % Plot the noise elimination for the real & imaginary parts only if
        % we use them for calculating the signal power.
        if ~FLAG_NOISE_ELI_VIA_AMP
            saveas(hNoiseSigmaReal, [pathNewCompNoiseSigmaToSave, 'real.png']);
            saveas(hNoiseSigmaImag, [pathNewCompNoiseSigmaToSave, 'imag.png']);
            if FLAG_SAVE_FIG_COPY
                saveas(hNoiseSigmaReal, [pathNewCompNoiseSigmaToSave, 'real.fig']);
                saveas(hNoiseSigmaImag, [pathNewCompNoiseSigmaToSave, 'imag.fig']);
            end
        end
        
        % Always plot the amplitude noise elimination for comparsion.
        saveas(hNoiseSigmaAmp, [pathNewCompNoiseSigmaToSave, 'amp.png']);
        if FLAG_SAVE_FIG_COPY
            saveas(hNoiseSigmaAmp, [pathNewCompNoiseSigmaToSave, 'amp.fig']);
        end
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
    
    % For fitting, remove points with too low calculated power.
    xsToFit = xs(ys>=minValidCalPower);
    ysToFit = ys(ys>=minValidCalPower);
    
    % For fitting, remove points with too low estimated SNR.
    boolsPtsToFit = estimatedSnrs{idxDataset}>=minValidEstSnr;
    % Also get rid of measurements to use in line fitting.
    boolsPtsToFit = boolsPtsToFit & BOOLS_MEAS_TO_FIT{idxDataset}';
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
if ~isinf(minValidCalPower)
    plot([finalAxis(1) finalAxis(2)], [minValidCalPower minValidCalPower], ...
        'r--');
    x = [finalAxis(1) finalAxis(1) finalAxis(2) finalAxis(2)];
    minY = -200;
    y = [minY minValidCalPower minValidCalPower minY];
    patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
end
title('Calibration results');
xlabel('Measured Power (dB)');
ylabel('Calculated Power (dB)');
grid minor; hold off;
pathCalFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'Calibration');
saveas(hFigCalibration, [pathCalFileToSave, '.png']);
if FLAG_SAVE_FIG_COPY
    saveas(hFigCalibration, [pathCalFileToSave, '.fig']);
end
disp('    Done!')

if FLAG_GEN_PLOTS_SILENTLY
    set(0,'DefaultFigureVisible','on');
end
% EOF