% SETPATH Add lib folders into Matlab path.
%
% Yaguang Zhang, Purdue, 07/11/2017

cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
addpath(fullfile(pwd, '..', 'gnuradio-tools', 'matlab'));
addpath(genpath(fullfile(pwd, 'lib')));

% EOF