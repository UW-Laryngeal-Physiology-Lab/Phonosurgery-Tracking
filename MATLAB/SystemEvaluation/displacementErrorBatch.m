function displacementErrorBatch(actDisp,endOffset,numFrames)
% DISPLACEMENTERRORBATCH Computes displacement error for batch of trials
%
% displacementErrorBatch(ACTDISP,ENDOFFSET,NUMFRAMES) Prompts user for
% directory with folders corresponding to individual displacement trials.
% Instrument displacement is calculated for each trial.  Error statistics
% are calculated between the actual instrument displacement (in mm) given
% by scalar ACTDISP and the calculated displacements.  The user is prompted
% to save an .xls containing the trial statistics.  ENDOFFSET and NUMFRAMES
% determine the tracked frame subset used to calculate the statistics.
% 
% Expected Folder Structure
% - MainDirectory
%       - TrialFolder
%       - TrialFolder
%       - TrialFolder
% The user is prompted for the MainDirectory.  All folder within the
% MainDirectory are assumed trialfolders.  The program looks for files
% named TrialFolder_left.mat and TrialFolder_right.mat in the TrialFolder.
% These files should be UGT generated tracking data for the left and right
% view respectively. If these files are not found, the folder is not
% included as a trial.
%
% Statistics are calculated for a subset of the tracked frames.  Assume M
% frames are available.  Then, statistics are calculated over frames [M -
% ENDOFFSET,NUMFRAMES + (M - ENDOFFSET)].  Displacement is calculated as
% the difference between the instrument track point 3D position in these
% frames and the first tracked frame.  A frame subset should be chosen such
% that the instrument has been completely displaced and is no longer
% vibrating.

% Prompt User For Directory Location
dName = uigetdir(cd,'Select Directory with Trial Folders');

% Prompt User for Calibration File
[calName,calPath] = uigetfile('*.mat','Load Camera Calibration File');
calFile = fullfile(calPath,calName);

% Get Folders
dirStruct = dir(dName);
dirStruct = dirStruct([dirStruct.isdir]);

%% Setup Data Cell
dataCell = {'Source','Folder','Max','Min','Range','Mean','RMS','xs_l',...
            'ys_l','xs_r','ys_r'};

%% Loop Through Analysis
for k = 1:numel(dirStruct)
    foldName = dirStruct(k).name;
    
    % Don't do anything if directory is root or parent
    if(any(strcmp(foldName,{'.','..'})))
        continue;
    end
    
    disp(foldName);
    
    % Grab Analysis Files
    leftFile = fullfile(dName,foldName,[foldName '_left.mat']);
    rightFile = fullfile(dName,foldName,[foldName '_right.mat']);
    
    % Get Displacment Error Analysis Data
    [e3d,p2d] = displacementErrorAnalysis(calFile,leftFile,rightFile,...
                                endOffset,numFrames,actDisp);
    
    % Append to Data Cell
    dataCell = [dataCell;
                {dName,foldName,e3d.maxError,e3d.minError,...
                e3d.rangeError,e3d.meanError,e3d.rmsError,...
                p2d.xs_l,p2d.ys_l,p2d.xs_r,p2d.ys_r}];
end 

%% Save Data Cell
[fName,pName] = uiputfile('*.xls','Save Displacement Error Stats');
if(fName ~= 0)
    xlswrite(fullfile(pName,fName),dataCell);
end