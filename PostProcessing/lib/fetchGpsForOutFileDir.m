function [ lat, lon, gpsLog ] ...
    = fetchGpsForOutFileDir( outFileDir )
%FETCHGPSFOROUTFILEDIR Find the cooresponding GPS log file and fetch the
%coordinates stored in it for the .out file specified by the input dir
%struct outFileDir.
%
% Outputs (lat, lon) is the resulted coordinate and gpsLog is the GPS
% information struct containing all the fields stored in the GPS log file.
%
% Yaguang Zhang, Purdue, 09/26/2017

% Get the filename for the .out file.
[~, outFileName] = fileparts(outFileDir.name);

% Find the GPS log file with the same time stamp (specified in the
% filename).
gpsFileDir = rdir(fullfile(outFileDir.folder, [outFileName, '_GPS.log']), '', false);

% Parse the GPS log.
gpsLog = parseGpsLog(gpsFileDir.name);
gpsLogSample = nmealineread(gpsLog.gpsLocation);

lat = gpsLogSample.latitude;
lon = gpsLogSample.longitude;

% Add a minus sign if it is W or S.
if(isW(gpsLog.gpsLocation))
    lon = -lon;
end
if(isS(gpsLog.gpsLocation))
    lat = -lat;
end
end
% EOF