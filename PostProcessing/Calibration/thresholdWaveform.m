function [ waveformThresholded, boolsEleminatedPts, hDebugFig ] ...
    = thresholdWaveform( waveform, flagDebug )
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
% Optionally, if flagDebug is set and its value is true, we will generate
% plots for debugging.
%
% Yaguang Zhang, Purdue, 08/14/2017

hDebugFig = nan;
if nargin < 2
    flagDebug = false;
end

% Find the peak that is "tall enough":
%   1. Use amplitude of the samples
%    2. For all amplitudes, minus 2*median
%   3. Eleminate negative points
%    4. Find the first sample higher than RATIO_VS_TALLEST of the tallest
%   one.
RATIO_VS_TALLEST = 0.75;
waveformAmps = abs(waveform);

twiceMedian = 2*median(waveformAmps);
shiftedWaveformAmps = waveformAmps - twiceMedian;
shiftedWaveformAmps(shiftedWaveformAmps<0) = 0;
[ampMax, ~] = max(shiftedWaveformAmps);

idxRefSample = find(shiftedWaveformAmps>RATIO_VS_TALLEST.*ampMax, 1);

% Find the standard deviation of the middle part of the sample segment
% before the peak, and set the threshold accordingly.
RATIO_SEGMENT_FOR_STD = 0.1;
numSampsToConsider = floor(RATIO_SEGMENT_FOR_STD*(idxRefSample-1));
idxRangeSampsToCons = floor(idxRefSample/2-numSampsToConsider/2): ...
    floor(idxRefSample/2+numSampsToConsider/2);
threshold = 3*std(waveform( idxRangeSampsToCons ));

% Set everything below the threshold to 0.
waveformThresholded = waveform;
boolsEleminatedPts = abs(waveform)<threshold;
waveformThresholded(boolsEleminatedPts) = 0;

if flagDebug
    % Generate plots for debugging.
    hDebugFig = figure;
    
    subplot(3,1,1); hold on;
    hShiftedWaveformAmp = plot(1:idxRefSample, ...
        shiftedWaveformAmps(1:idxRefSample), '.', ...
        'Color', [1,1,1]*0.4);
    maxShifedWaveformAmp = plot([0 0;1 1]*idxRefSample, ...
        [1, -1;1 -1]*ampMax, 'r-.');
    hold off; grid on; axis tight;
    legend([hShiftedWaveformAmp, maxShifedWaveformAmp(1)], ...
        'Used to compute delta', 'Other points before peak #1'); 
    title('Samples for computing noise sigma');
    
    subplot(3,1,2); hold on;
    hPtsNotCons = plot(1:idxRefSample, ...
        waveform(1:idxRefSample), '.', ...
        'Color', [1,1,1]*0.4);
    hPtsConsidered = plot(idxRangeSampsToCons, ...
        waveform(idxRangeSampsToCons), '*b');
    hold off; grid on; axis tight;
    legend([hPtsConsidered, hPtsNotCons], ...
        'Used to compute delta', 'Other points before peak #1'); 
    title('Samples for computing noise sigma');
    
    subplot(3,1,3); hold on;
    idxMaxToShow = min(floor(idxRefSample*2),length(waveform));
    waveformToShow = waveform(1:idxMaxToShow);
    plot(1:idxMaxToShow, ...
        waveformToShow, '.', ...
        'Color', 'b');
    indicesDiscared = find(abs(waveformToShow)<threshold);
    hEleminatedPts = plot(indicesDiscared, ...
        waveformToShow(indicesDiscared), '.', ...
        'Color', [1,1,1]*0.5);
    hTwiceMedian = ...
        plot([0 0;1 1]*(2*idxRefSample+1), ...
        [1, -1;1 -1]*twiceMedian, 'r-.');
    grid on; axis tight; finalAxis = axis;
    hDiscardRegion = ...
        patch([0,0,2*idxRefSample+1,2*idxRefSample+1], ...
        [-threshold,threshold,threshold,-threshold], ...
        [1,1,1]*0.9);
    set(hDiscardRegion, 'LineStyle', 'None');
    uistack(hDiscardRegion,'bottom')
    hold off;
    axis(finalAxis);
    legend([hDiscardRegion, hEleminatedPts, hTwiceMedian(1)], ...
        'Region to eleminate', 'eleminated Samples', '+/- 2*Median'); 
    title('Eleminated samples (Set to 0)');
end

end
% EOF