function [ normalizedAntPat ] = normalizeAntAmp( antPat, maxAntGainInDb )
%PLOTANTPATTERN Plot the the antenna plane pattern in a 2D polar coordinate
%system.
%
% Inputs:
%   - antPat 
%     The antenna patterns, for the Azimuth or the Elevation sweep; It is a
%     struct containing fields:
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

normalizedAntPat = antPat;
ampsInDb = 10.*log10(antPat.amps);
ampsInDb = ampsInDb-max(ampsInDb)+maxAntGainInDb;

normalizedAntPat.amps = 10.^(ampsInDb./10);

end
% EOF