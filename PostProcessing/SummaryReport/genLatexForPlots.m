% GENLATEXFORPLOTS Generate a Latex snippet for embedding in the summary
% report.
%
% Note the script has only been tested on a Windows machine.
%
% Yaguang Zhang, Purdue, 07/12/2017

clear; clc; close all;

%% Configurations

% The absolute path to the shared Google Drive folder "Annapolis
% Measurement Campaign". Please make sure it is correct for the machine
% which will run this script.
%  - On Mac Lemma:
%    '/Users/zhan1472/Google Drive/Annapolis Measurement Campaign';
%  - On Windows Dell:
%    '\\LEMMA\Google Drive\Annapolis Measurement Campaign';
%  - Local copy on Windows Dell:
%    'C:\Users\Zyglabs\Documents\MEGAsync\EARS';
ABS_PATH_TO_EARS_SHARED_FOLDER = ...
    'C:\Users\Zyglabs\Documents\MEGAsync\EARS';

% Configure other paths accordingly.
ABS_PATH_TO_DATA = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, 'Data');
ABS_PATH_TO_PLOTS = fullfile(ABS_PATH_TO_EARS_SHARED_FOLDER, ...
    'PostProcessingResults', 'SummaryReport', 'plots');

%% Before Generating the Snippet

disp(' ------------------ ')
disp('  genLatexForPlots')
disp(' ------------------ ')

% Add libs to current path.
cd(fileparts(mfilename('fullpath')));
ABS_PATH_TO_SAVE_SNIPPET = fullfile(pwd, 'photosAndPlots.tex');
addpath(fullfile(pwd));
cd('..'); setPath;

% The plots should already be there.
assert(exist(ABS_PATH_TO_PLOTS, 'dir')==7, ...
    'No plots were found! Please run genPlots first.');


% Load the plot information generated by genPlots.
disp(' ')
disp('    Loading information for the plots...')
% allSeriesParentDirs and allSeriesDirs will be available after this.
load(fullfile(ABS_PATH_TO_PLOTS, 'plotInfo.mat'));
disp('    Done!')

%% Generate the Snippet Site by Site

disp(' ')
disp('  => Generating Latex snippet for plots...')

% Write handler.
fID = fopen(ABS_PATH_TO_SAVE_SNIPPET,'wt');
% Set Latex graphic path.
fprintf(fID, ...
    '\\graphicspath{{../../../PostProcessingResults/SummaryReport/plots/}}\n');

% For each day and a specific type of measurements (e.g. LargeScale, SIMO
% and Conti)...
for idxPar=1:length(allSeriesParentDirs)
    disp([num2str(idxPar), '/', num2str(length(allSeriesParentDirs))])
    
    % Note that depending on the machine where genPlots was excuted, the
    % exact paths to the plots may change, so we will only use the folder
    % names.
    [~, curParFolderName] = fileparts(allSeriesParentDirs(idxPar).name);
    
    % Find the GPS on map figure.
    imageGpsOnMap = rdir(fullfile(ABS_PATH_TO_PLOTS, ...
        [curParFolderName, '_gpsSamplesOnMap.png']),'',1);
    
    % Find the .txt log for measurement sites.
    txInfo = rdir(fullfile(ABS_PATH_TO_DATA, curParFolderName, ...
        'TxInfo.txt'));
    
    % Find all photos for the measurement sites.
    photosAntenJpg = rdir(fullfile(ABS_PATH_TO_PLOTS, ...
        [curParFolderName, '_antenna*.jpg']));
    photosSetupJpg = rdir(fullfile(ABS_PATH_TO_PLOTS, ...
        [curParFolderName, '_setup*.jpg']));
    % Convert to .png for Latex if necessary.
    jpgFiles = [photosAntenJpg; photosSetupJpg];
    for idxJpgFile = 1:length(jpgFiles)
        pngFile = regexprep(jpgFiles(idxJpgFile).name,'.jpg$','.png');
        if (exist(pngFile, 'file')~=2)
            jpgIm = imread(jpgFiles(idxJpgFile).name);
            imwrite(jpgIm, pngFile);
        end
    end
    photosAnten = rdir(fullfile(ABS_PATH_TO_PLOTS, ...
        [curParFolderName, '_antenna*.png']),'',1);
    photosSetup = rdir(fullfile(ABS_PATH_TO_PLOTS, ...
        [curParFolderName, '_setup*.png']),'',1);
    
    % Output to the Latex snippet.
    fprintf(fID, '\\section{%s}\n', replace(curParFolderName, '_', ' - '));
    % GPS on map and txInfo.
    fprintf(fID, '\\subsection{Overview}\n');
    fprintf(fID, '\\begin{figure}[ht] \\caption{GPS samples on Google map}\n');
    fprintf(fID, '\\includegraphics[width=\\textwidth]{%s}\\centering\\end{figure}\n', ...
        imageGpsOnMap.name);
    fprintf(fID, '\\begin{minipage}{\\textwidth} Measurement logs:\n\n');
    fprintf(fID, '\\lstinputlisting[basicstyle=\\ttfamily\\scriptsize, breaklines]{%s} \\end{minipage}\n', ...
        replace(txInfo.name, '\', '/'));
    fprintf(fID, '\\clearpage\n');
    
    % Illustration for sample measurements.
    fprintf(fID, '\\subsection{Plots for Sample Measurements}\n');
    curSeriesDirs = allSeriesDirs{idxPar};
    for idxSeries = 1:length(curSeriesDirs)
        partFolder = regexp(curSeriesDirs(idxSeries).name,'Data\/(201\d+_\w+)\/','tokens');
        assert(length(partFolder)==1, 'There should be exactly 1 date & type pair as the name for the parent data folder!');
        partFolder = partFolder{1}{1};
        fprintf(fID, '\\subsubsection{Site label: %s}\n', replace([partFolder,'_Series_',num2str(idxSeries)], '_', '\_'));
        % Photos.
        fprintf(fID, '\\begin{figure}[ht] \\caption{Photo from the antenna}\n');
        fprintf(fID, '\\includegraphics[width=0.9\\textwidth]{%s}\\centering\\end{figure}\n', ...
            ['"',replace(photosAnten(idxSeries).name, '.png', ''),'"']);
        fprintf(fID, '\\begin{figure}[ht] \\caption{Photo for the setup}\n');
        fprintf(fID, '\\includegraphics[width=0.9\\textwidth]{%s}\\centering\\end{figure}\n', ...
            ['"',replace(photosSetup(idxSeries).name, '.png', ''),'"']);
        % Signals. Find the signal over time plots for this series.
        imagesSig = rdir(fullfile(ABS_PATH_TO_PLOTS, ...
            [curParFolderName, '_oneSigPerMeas_series_', num2str(idxSeries), '_*.png']),'',1);
        for idxMea = 1:length(imagesSig)
            imgSigTitle = replace(replace(...
                imagesSig(idxMea).name, '.png', ''), ...
                '_', ' ');
            fprintf(fID, '\\begin{figure}[ht] \\caption{%s}\n', imgSigTitle);
            fprintf(fID, '\\includegraphics[width=0.9\\textwidth]{%s}\\centering\\end{figure}\n', ...
                imagesSig(idxMea).name);
        end
        fprintf(fID, '\\clearpage\n');
    end
end

fclose(fID);
disp('     Done!')

% EOF