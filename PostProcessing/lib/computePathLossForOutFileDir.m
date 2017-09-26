function [ pathLossInDb ] ...
    = computePathLossForOutFileDir(curOutFileDir, rxGain, ...
    noiseEliminationFct, powerShiftsForCali)
%COMPUTEPATHLOSSFOROUTFILEDIR Load the Gnu Radio samples stored in the .out
%file specified by the input dir struct outFileDir, and compute the path
%loss for it.
%
% We will consider both the TX calibration and the antenna normalization.
%
% Inputs:
%   - curOutFileDir
%     A dir struct specifying which .out file will be processed. Note: curOutFileDir should at least has the field name, which contains
%     the full absolute path for the .out file to be processed, or only the
%     file name; For the second case, another field folder is required for the
%     full abosolute path of the parent folder.
%   - rxGain
%     A scalar. The Gnu Radio gain (in dB).
%   - noiseEliminationFct
%     A function to run the loaded signal through for a noise-eliminated
%     copy of the signal. The syntax of it should be like:
%         [signalNoiseEliminated, boolsEliminatedPts] = ...
%             noiseEliminationFct(signalComplexArray)
%   - powerShiftsForCali
%     Essentially the parameter b in y=ax+b for the calibration line
%     corresponding to the current rxGain.
%
% Yaguang Zhang, Purdue, 09/26/2017

%% Parameters

% TX power into upconverter in dBm.
txPower = -8;

% Low pass filter for the PSD. Tried before: 46000; 39500.
maxFreqPassed = 20000; % In Hz.
% High pass filter to remove the DC component.
minFreqPassed = 1; % In Hz.

% Number of samples to discard at the beginning.
numStartSampsToDiscard = 100000; % ~0.1s
% After discarding these samples, furthermore only keep the middle part of
% the signal for calibration.
timeLengthAtCenterToUse = 1; % In second.

% Sample rate used for GnuRadio.
try
    Fs = evalin('base', 'Fs');
catch
    warning('GnuRadio sample frequency not found in the base workspace.')
    warning('Will use the default value 1.04 * 10^6.')
    Fs = 1.04 * 10^6;
end

%% Load Data

if isfield(curOutFileDir, 'folder')
    [~, outFileName] = fileparts(curOutFileDir.name);
    absPathOutFile = fullfile(curOutFileDir.folder, [outFileName, '.out']);
else
    absPathOutFile = curOutFileDir.name;
end
curSignal = read_complex_binary(absPathOutFile);

%% Calculate Power

% If any of the input

% Discard the first numStartSampsToDiscard of samples.
curSignal = curSignal((numStartSampsToDiscard+1):end);
% Further more, only keep the middle part for calibration.
numSampsToKeep = ceil(timeLengthAtCenterToUse*Fs);
numSampsCurSeries = length(curSignal);
if numSampsToKeep > numSampsCurSeries
    warning('There are not enough samples to keep. We will use all remaining ones.');
else
    idxRangeToKeep = floor(0.5.*numSampsCurSeries ...
        + [-1,1].*numSampsToKeep./2);
    curSignal = curSignal(idxRangeToKeep(1):idxRangeToKeep(2));
end
% Make sure we end up with even number of samples.
if mod(length(curSignal),2)==1
    curSignal = curSignal(1:(end-1));
end

% Noise elimination.
[~, boolsEliminatedPts] = ...
    noiseEliminationFct(curSignal);
curSignalEliminated = curSignal;
curSignalEliminated(boolsEliminatedPts) = 0;

% For the signal to process, i.e. the noise eliminiated signal, compute the
% PSD.
X = curSignalEliminated;
L = length(X);
% FFT results.
Y = fftshift(fft(X));
% Frequency domain.
f = (-L/2:L/2-1)*(Fs/L);
idxDC = L/2+1;
% PSD.
powerSpectralDen = abs(Y).^2/L;

% Compute the power.
boolsFPassed = abs(f)<=maxFreqPassed ...
    & abs(f)>=minFreqPassed;
% Compute the power by integral. Note that we will always discard the DC
% component here (although it may be passed by the filters).
psdPassed = powerSpectralDen;
psdPassed(~boolsFPassed) = 0;
psdPassed(idxDC) = 0;
calcP = trapz(f, psdPassed);


%% Rx Calibration

% Change to dB and remove the gain from the Gnu Radio.
calcPInDbShifted = 10.*log10(calcP) - rxGain;
measPInDb = calcPInDbShifted + powerShiftsForCali;

%% Antenna Nomalization
antGain = 0;

%% Final Result
pathLossInDb = txPower + antGain - measPInDb;

end
% EOF