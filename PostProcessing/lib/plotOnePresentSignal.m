function [ hFig ] = plotOnePresentSignal( signal, ...
    numPreSamples, numPostSamples, figureName, Fs)
%PLOTONEPRESENTSIGNAL Plot the tallest bump of the signal.
%
% Inputs:
%   - signal
%     A vector of complex numbers.
%   - numPreSamples, numPostSamples
%     Integers. The range of signal to be plotted can be adjusted using
%     numPreSamples and numPostSamples. The samples before and after the
%     first tallest signal sample found will be plotted accordingly if
%     there are enough samples available.
%   - figureName
%     Optional. A string to specify the figure's name. 
%   - Fs
%     Optional. A number to specify the sample rate for the input signal.
%     Used to properly set the x axis / time line labels. Default
%     1.04MSamples/s.
%
% Yaguang Zhang, Purdue, 06/18/2017

if nargin<5
    Fs = 1.04*10^6;
end

RATIO_VS_TALLEST = 0.9;

signalReal = real(signal);
signalImag = imag(signal);

% Construct data for the x axis (time line).
numSamples = length(signal);
xTimes = [1:numSamples];


[valueMax, ~] = max(real(signal));
% Find the first peak that is tall enough.
idxRefSample = find(signalReal>RATIO_VS_TALLEST.*valueMax, 1);

minIdxToPlot = idxRefSample - numPreSamples;
if(minIdxToPlot<1)
    minIdxToPlot = 1;
end
maxIdxToPlot = idxRefSample + numPostSamples;
if(maxIdxToPlot>length(signal))
    minIdxToPlot = length(signal);
end

if nargin>3
    hFig = figure('Name',figureName);
else
    hFig = figure;
end

subplot(4,1,1); hold on;
hR = plot(signalReal(minIdxToPlot:maxIdxToPlot), 'b-');
hI = plot(signalImag(minIdxToPlot:maxIdxToPlot), 'r-.');
title('First Detected Signal (Volt)')
hold off; legend([hR, hI], 'Real', 'Image'); axis tight;
subplot(4,1,2); hold on;
hR = plot(10.*log10(signalReal(minIdxToPlot:maxIdxToPlot).^2), 'b-');
hI = plot(10.*log10(signalImag(minIdxToPlot:maxIdxToPlot).^2), 'r-.');
title('First Detected Signal (in dB Volt)')
hold off; legend([hR, hI], 'Real(dB)', 'Image(dB)'); axis tight;
subplot(4,1,3); hold on;
signalAmp = abs(signalImag(minIdxToPlot:maxIdxToPlot));
hAmp = plot(signalAmp, 'b-');
title('First Detected Signal (Amplitude in Volt)')
hold off; legend([hAmp], 'Amplitude'); axis tight;
subplot(4,1,4); hold on;
hAmp = plot(10.*log10(signalAmp), 'b-');
title('First Detected Signal (Amplitude in dB Volt)')
hold off; legend([hAmp], 'dB Amplitude'); axis tight;
if nargin>3
    suptitle(figureName);
end
% EOF