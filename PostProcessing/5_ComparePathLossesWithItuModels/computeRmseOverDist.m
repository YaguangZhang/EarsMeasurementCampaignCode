function [ rmses, dists ] = computeRmseOverDist(dsFromTx, modelLs, measLs)
%COMPUTERMSEOVERDCIST Computes the Root Mean Square Errors over distance. 
%
% The samples from the same distance will be used for each RMSE computation
% (for one specific distance). 
%
% Inputs:
%   - dsFromTx
%     The distances from the TX.
%   - modelLs
%     The paths losses from the model.
%   - measLs
%     The measured reference paths losse.
% Output:
%   - rmses, dists
%     The the resulted root mean square errors and the corresponding
%     distances.
%
% Yaguang Zhang, Purdue, 10/19/2017

% First we need to get groups of measLs with the same dsFromTx.
[indicesStarts, indicesEnds, dists] ...
    = findConsecutiveSubSeqs(dsFromTx);

numResults = length(dists);
rmses = nan(numResults, 1);
for idx = 1:numResults
    curIndices = indicesStarts(idx):indicesEnds(idx);
    rmses(idx) = sqrt(mean(...
        ( modelLs(curIndices) - measLs(curIndices) ).^2 ...
        ));
end

% Sort the results according to dists.
results = sortrows([rmses, dists], 2);
rmses = results(:,1);
dists = results(:,2);
end
% EOF