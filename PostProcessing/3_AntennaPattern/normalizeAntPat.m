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

% Supported method:
%    - 'asAWhole'
%      This method uses the maxium amplitude of the combined data set of
%      the data from both the aximuth and elevation sweep, and normalizes
%      that data set as a whole;
%    - 'separately'
%      This method normalizes the two sweep patterns separately, to make
%      sure the maximum gain for both sweep planes are the same.
NORM_METHOD = 'separately';

normalizedAntPatAz = antPatAz;
normalizedAntPatEl = antPatEl;

switch NORM_METHOD
    case 'asAWhole'
        ampsInDb = 10.*log10([antPatAz.amps; antPatEl.amps]);
        ampsInDb = ampsInDb-max(ampsInDb)+maxAntGainInDb;
        
        amps = 10.^(ampsInDb./10);
        
        idxD = length(antPatAz.amps);
        normalizedAntPatAz.amps = amps(1:idxD);
        normalizedAntPatEl.amps = amps((idxD+1):end);
    case 'separately'
        ampsInDbAz = 10.*log10(antPatAz.amps);
        normalizedAmpsInDbAz = ampsInDbAz-max(ampsInDbAz)+maxAntGainInDb;
        normalizedAntPatAz.amps = 10.^(normalizedAmpsInDbAz./10);
        
        ampsInDbEl = 10.*log10(antPatEl.amps);
        normalizedAmpsInDbEl = ampsInDbEl-max(ampsInDbEl)+maxAntGainInDb;
        normalizedAntPatEl.amps = 10.^(normalizedAmpsInDbEl./10);
    otherwise
end

end
% EOF