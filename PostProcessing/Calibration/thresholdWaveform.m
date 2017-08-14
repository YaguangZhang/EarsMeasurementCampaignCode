function [ waveformThresholded ] = thresholdWaveform( waveform )
%THRESHOLDWAVEFORM Thresholded the waveform (a real vector) to eliminate
%corss-correlation and system noise.
%
% We have implemented the algorithm provided by Professor Chris Anderson:
%
%   Find the standard deviation in the first 5-10% of samples before the
%   first peak, set the threshold at 3*sigma above 0, and remove everything
%   below the threshold.
%
% For complex signals, you may want to apply this function twice to its
% real and imaginary parts, separately.
%
% Yaguang Zhang, Purdue, 08/14/2017

% Find the first peak that is tall enough.
RATIO_VS_TALLEST = 0.9;
[valueMax, ~] = max(waveform);
idxRefSample = find(waveform>RATIO_VS_TALLEST.*valueMax, 1);

% Find the standard deviation of the start part of the sample segment
% before the peak, and set the threshold accordingly.
RATIO_SEGMENT_FOR_STD = 0.1;
threshold = 3*std(waveform( ...
    1:floor(RATIO_SEGMENT_FOR_STD*(idxRefSample-1)) ...
    ));

% Set everything below the threshold to 0.
waveformThresholded = waveform;
waveformThresholded(waveform<threshold) = 0;
end
% EOF