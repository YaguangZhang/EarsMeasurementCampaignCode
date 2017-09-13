function [ hPsdInDb ] ...
    = plotPsdInDbWithLPF( powerSpectralDen, f, maxFreqPassed, ...
    estimatedSnr )
%PLOTPSDINDBWITHLPF Plot the PSD (in dB) for the input power spectral
%density with frequency range f.
%
% AnD maxFreqPassed is the cutoff frequency of the input LPF.
% Correspondingly, the signal frequency range that can pass will also be
% plotted.
%
% Optionally, if estimatedSnr is specified, it will be shown on the plot,
% too.
%
% Yaguang Zhang, Purdue, 09/13/2017

hPsdInDb = figure; hold on;
powerSpectralDenIndB = 10*log10(powerSpectralDen);
hPowerSpectralDenIndB = plot(f, powerSpectralDenIndB);
curAxis = axis; curAxis(1) = f(1); curAxis(2) = f(end);
if ~isinf(maxFreqPassed)
    hLPFIndB = plot([maxFreqPassed, -maxFreqPassed; ...
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
if nargin>3
    text(min(maxFreqPassed, 50000), mean([curAxis(3) curAxis(4)]), ...
        ['Estimated SNR = ', ...
        num2str(estimatedSnr, '%.2f')]);
end
hold off;
xlabel('f (Hz)'); ylabel('Estimated PSD (V^2/Hz in dB)');
axis(curAxis);
if ~isinf(maxFreqPassed)
    legend([hPowerSpectralDenIndB, hLPFIndB(1)], 'PSD (dB)', 'LPF');
else
    legend(hPowerSpectralDenIndB, 'PSD (dB)');
end
grid minor;

end
% EOF