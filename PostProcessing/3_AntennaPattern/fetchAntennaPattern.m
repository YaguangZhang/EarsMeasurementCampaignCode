% FETCHANTENNAPATTERN Fetch the azimuth and elevation S21 sweeps for the
% frequency needed.
%
% We will also save the results into a .mat file for future use, as well as
% generating plots for debugging. Currently, we only load the amplitude
% (i.e. phases are ignored).
%
% Yaguang Zhang, Purdue, 09/29/2017

clear; clc; close all;

%% Configurations

% Add libs to current path and set ABS_PATH_TO_EARS_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% Configure other paths accordingly.
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'AntennaPattern');

% Frequency desired. We were doing the 28GHz measurements, with the
% built-in upconverter and downconverter working for 2.5 GHz intermediate
% frequency.
inFreq  = 2.5*10^9;

% For antenna pattern normalization.
maxAntGainInDb = 20;

%% Before Processing the Data

disp(' --------------------- ')
disp('  fetchAntennaPattern')
disp(' --------------------- ')

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

% Set the path to the antenna pattern log files.
absPathAntPatFolder = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'Data', 'Antenna Pattern');
absPathAntPatAz = fullfile(absPathAntPatFolder, 'NSF_28GHZ_AZ_Cut.txt');
absPathAntPatEl = fullfile(absPathAntPatFolder, 'NSF_28GHZ_EL_Cut.txt');

%% Load the Antenna Data

disp(' ')
disp('    Loading antenna pattern data from the log files ...')

assert(exist(absPathAntPatAz, 'file')==2, ...
    ['Couldn''t find the log file: ', absPathAntPatAz]);
assert(exist(absPathAntPatEl, 'file')==2, ...
    ['Couldn''t find the log file: ', absPathAntPatEl]);

% 1601 points were swept from 2-3 GHz. These inputs coorespond to 27.5-28.5
% GHz for the TX signal.
scanndedFreqs = 2000000000:625000:3000000000;
intFreq = 2500000000; % 2.5 GHz => 25.5+2.5 = 28 GHz.

pat28GAz = parseAntPatLog( absPathAntPatAz, ...
    scanndedFreqs, intFreq );
pat28GEl = parseAntPatLog( absPathAntPatEl, ...
    scanndedFreqs, intFreq );
% The measurement system was set up for an Azimuth cut in both cases, even
% though the second is an elevation pattern, so here we need to switch the
% az and el for pat28GEl.
els = pat28GEl.azs;
pat28GEl.azs = pat28GEl.els;
pat28GEl.els = els;

% At last, we need to normalize the results.
pat28GAzNorm = normalizeAntAmp(pat28GAz, maxAntGainInDb);
pat28GElNorm = normalizeAntAmp(pat28GEl, maxAntGainInDb);

disp('    Done!')

%% Save the Results

disp(' ')
disp('    Saving antenna pattern data ...')

pathToSaveLoadedPattern = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'antennaPattern.mat');
save(pathToSaveLoadedPattern, 'pat28GAzNorm', 'pat28GElNorm');

disp('    Done!')

%% Plot

disp(' ')
disp('    Plotting...')

[hPat2DRef, hPat3DRef, ...
    hInterPat3DOnLineLinear, hInterPat3DOnLineDb, ...
    hInterPat3DWeightedSumLinear, hInterPat3DWeightedSumDb] ...
    = plotAntPattern(pat28GAzNorm, pat28GElNorm);

% Save the plots.
absPathCurFile = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'pat2DRef');
saveas(hPat2DRef, [absPathCurFile, '.png']);
saveas(hPat2DRef, [absPathCurFile, '.fig']);

absPathCurFile = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'pat3DRef');
saveas(hPat3DRef, [absPathCurFile, '.png']);
saveas(hPat3DRef, [absPathCurFile, '.fig']);

absPathCurFile = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'interPat3DOnLineLinear');
saveas(hInterPat3DOnLineLinear, [absPathCurFile, '.png']);
saveas(hInterPat3DOnLineLinear, [absPathCurFile, '.fig']);

absPathCurFile = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'interPat3DOnLineDb');
saveas(hInterPat3DOnLineDb, [absPathCurFile, '.png']);
saveas(hInterPat3DOnLineDb, [absPathCurFile, '.fig']);

absPathCurFile = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    'interPat3DWeightedSumLinear');
saveas(hInterPat3DWeightedSumLinear, [absPathCurFile, '.png']);
saveas(hInterPat3DWeightedSumLinear, [absPathCurFile, '.fig']);

absPathCurFile = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'interPat3DWeightedSumDb');
saveas(hInterPat3DWeightedSumDb, [absPathCurFile, '.png']);
saveas(hInterPat3DWeightedSumDb, [absPathCurFile, '.fig']);

disp('    Done!')

% EOF