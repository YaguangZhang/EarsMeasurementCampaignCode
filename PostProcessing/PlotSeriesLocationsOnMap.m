% PLOTSERIESLOCATIONSONMAP Plot the GPS points for the series data in a
% specified folder on a Google map.
%
% All data for each series should be organized in its own folder with name
% "Series_#" (e.g. "Series_1"). The GPS data should be contained in files
% named like "measureSignal_1497709963_GPS.log".
%
% Yaguang Zhang, Purdue, 06/07/2017

% Load data and set the current Matlab directory.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
addpath(fullfile(fullfile(pwd, 'lib')));

PATH_FOLDER_TO_PROCESS = fullfile(pwd, '..', '..', 'Data', '20170617_LargeScale');

% For each folder, read in all the GPS log files.

v = read_complex_binary(filePath);

vToPlot = v(1:100000);
figure;
hold on;
plot(real(vToPlot), 'b-');
plot(imag(vToPlot), 'r-.');
legend('Real','Image');
% EOF