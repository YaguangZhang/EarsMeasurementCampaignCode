function [ outFilesDirs, gpsFilesDirs ] ...
    = loadFilesFromSeriesDir( seriesDir )
%LOADFILESFROMSERIESDIR Load the Gnu Radio out files, as well as the
% GPS log files, stored in one _Series folder specified by the input dir
% struct seriesDir.
%
% The outputs, outFilesDirs and gpsFilesDirs, are the resulted dir struct
% arrays, for Gnu Radio out files and the GPS log files, respectively.
%
% Yaguang Zhang, Purdue, 09/26/2017

% Check the folder's name.
[~, seriesFolderName] = fileparts(seriesDir.name);
idxSeriesTokens = regexp(seriesFolderName, ...
    '^Series_(\d+)$', 'tokens');
assert(length(idxSeriesTokens)==1, ...
    'The name for the input folder does not start with `Series_`!');

outFilesDirs = rdir(fullfile(seriesDir.folder, seriesFolderName, ...
    '*.out'), 'regexp(name, ''measureSignal_\d+.out$'')', true);
gpsFilesDirs = rdir(fullfile(seriesDir.folder, seriesFolderName, ...
    '*_GPS.log'));

end
% EOF