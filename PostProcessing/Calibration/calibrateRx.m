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
    '/Users/zhan1472/Google Drive/Annapolis Measurement Campaign';

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
        % Compute the complex FFT of the resulted signal.
        curCalDataThr{idxCurMeas} = signalReal+1i.*signalImag;
        powerSpectralDen = fft(curCalDataThr{idxCurMeas});
        curCalculatedP(idxCurMeas) = trapz(abs(powerSpectralDen));
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
hFigCalibration = figure; hold on;
for idxDataset = 1:numDatasets
    % For plotting, remove points with inf as calculated power.
    xs = measPowers{idxDataset};
    ys = calculatedPowers{idxDataset};
    xs = xs(~isinf(ys));
    ys = ys(~isinf(ys));
    scatter(xs, ys, '*');
end
lsline;
hold off;
pathFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'Calibration');
saveas(hFigCalibration, [pathFileToSave, '.fig']);
saveas(hFigCalibration, [pathFileToSave, '.png']);
disp('    Done!')
% EOF