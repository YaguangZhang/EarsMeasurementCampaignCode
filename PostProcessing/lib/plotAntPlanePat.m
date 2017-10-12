function [ polarAx, hpbw ] ...
    = plotAntPlanePat( polarAx, anglesInDegree, ampsInDb )
%PLOTANTPATTERN Plot the the antenna plane pattern in a 2D polar coordinate
%system.
%
% Inputs:
%   - anglesInDegree
%     A real vector specifying the angles in degree.
%   - ampsInDb
%     The amplitude in dB.
%
% Yaguang Zhang, Purdue, 10/02/2017

try
    % If we can find the value for normalizing the antenna patter, we will
    % use it in the HPBW computation.
    maxPowerInDb = evalin('base', 'maxAntGainInDb');
    [ hpbw, anglesHpbw, powerInDbHpbw ] ...
        = computeHPBW( anglesInDegree, ampsInDb, maxPowerInDb );
catch
    [ hpbw, anglesHpbw, powerInDbHpbw ] ...
        = computeHPBW( anglesInDegree, ampsInDb );
end

% Plot.
hold on;
polarplot(deg2rad(anglesInDegree), ampsInDb-min(ampsInDb));
polarAx.RAxis.Label.String = 'dB';
polarplot(deg2rad([0 0; anglesHpbw']), ...
    [0, 0; powerInDbHpbw'-min(ampsInDb)], '--k');
hold off;

end
% EOF