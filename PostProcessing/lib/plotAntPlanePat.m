function [ hPolar ] = plotAntPlanePat( polarAx, anglesInDegree, ampsInDb )
%PLOTANTPATTERN Plot the the antenna plane pattern in a 2D polar coordinate
%system.
%
% Inputs:
%   - polarAx 
%     The polar axes where we should make the plot, e.g. polarAx =
%     subplot(2,1,1,polaraxes).
%   - anglesInDegree
%     A real vector specifying the angles in degree.
%   - ampsInDb
%     The amplitude in dB.
% 
% Yaguang Zhang, Purdue, 10/02/2017

polarplot(deg2rad(anglesInDegree), ampsInDb-min(ampsInDb));
polarAx.RAxis.Label.String = 'dB';

end
% EOF