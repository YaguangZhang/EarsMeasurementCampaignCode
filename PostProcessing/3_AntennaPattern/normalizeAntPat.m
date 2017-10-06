function [ normalizedAntPatAz, normalizedAntPatEl ] ...
    = normalizeAntPat( antPatAz, antPatEl, ...
    maxAntGainInDb )
%NORMALIZEANTPAT Normalize the antenna pattern specified.
%
% Inputs:
%   - antPatAz, antPatEl
%     The antenna patterns, for the Azimuth and the Elevation sweep,
%     respectively; Each of them is a struct containing fields:
%       - azs
%         The azimuth angles in degree.
%       - els
%         The elevation angles in degree.
%       - amps
%         The linear amplitudes of the samples.
%       - phases
%         The phases of the samples.
%     All of these fields contains a column vector with each row
%     corresponding to a sweep sample.
%   - maxAntGainInDb
%     The maximum antenna gain in dB.
% 
% Yaguang Zhang, Purdue, 10/02/2017

normalizedAntPatAz = antPatAz;
normalizedAntPatEl = antPatEl;

ampsInDb = 10.*log10([antPatAz.amps; antPatEl.amps]);
ampsInDb = ampsInDb-max(ampsInDb)+maxAntGainInDb;

amps = 10.^(ampsInDb./10);

idxD = length(antPatAz.amps);
normalizedAntPatAz.amps = amps(1:idxD);
normalizedAntPatEl.amps = amps((idxD+1):end);

end
% EOF