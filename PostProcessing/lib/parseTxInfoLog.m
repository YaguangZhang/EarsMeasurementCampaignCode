function [ txInfoLog, absParDir ] = parseTxInfoLog( logAbsDir )
%PARSETXINFOLOG Parse a TxInfo.txt log file specified by the input absolute
%path.
%
% Yaguang Zhang, Purdue, 10/04/2017

% Start line pattern.
startIndicator = '----';
% Sample line pattern:
%   Series #,   Location,   TXAZ,   TXEL,   RXAZ,   RXEL,   note
sampleLinePat = '%d %s %f %f %f %f %s';
% Sample line delimiter.
sampleLineDel = ',';

logFileID = fopen(logAbsDir, 'r');
% Find the start line.
curLine = '';
while(~strcmp(curLine(1:min(length(curLine),4)), startIndicator))
    curLine = fgetl(logFileID);
end
% Get rid of the title line for the table.
fgetl(logFileID);

% Fetch all the log info.
logs = textscan(logFileID, sampleLinePat, ...
    'Delimiter', sampleLineDel, 'CollectOutput', 0);
fclose(logFileID);

% Convert the log info to a struct array.
numSeries = length(logs{1});
txInfoLog = repmat(struct('seriesNum', nan, ...
    'location', '', ...
    'txAz', nan, ...
    'txEl', nan, ...
    'rxAz', nan, ...
    'rxEl', nan, ...
    'note', nan), numSeries, 1);
for idxSeries = 1:numSeries
    % Always available: seriesNums and locations.
    txInfoLog(idxSeries).seriesNum = logs{1}(idxSeries);
    txInfoLog(idxSeries).location = logs{2}{idxSeries};
    % May not be available: txAzs, txEls, rxRzs, txEls, and notes.
    if length(logs{3})>=idxSeries
        txInfoLog(idxSeries).txAz = logs{3}(idxSeries);
    end
    if length(logs{4})>=idxSeries
        txInfoLog(idxSeries).txEl = logs{4}(idxSeries);
    end
    if length(logs{5})>=idxSeries
        txInfoLog(idxSeries).rxAz = logs{5}(idxSeries);
    end
    if length(logs{6})>=idxSeries
        txInfoLog(idxSeries).rxEl = logs{6}(idxSeries);
    end
    if length(logs{7})>=idxSeries
        txInfoLog(idxSeries).note = logs{7}{idxSeries};
    end
end

% Store the absolute path to the folder containing the corresponding series
% data.
absParDir = fileparts(logAbsDir);

end
% EOF