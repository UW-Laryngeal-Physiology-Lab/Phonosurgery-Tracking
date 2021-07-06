function genResultsViewData(camMatFile,leftAnalysisFile,rightAnalysisFile)
% GENRESULTSVIEWDATA Generates .mat for Results Viewer
%
% genResultsViewData(CAMMATFILE,LEFTANALYSISFILE,RIGHTANALYSISFILE)
% Generates .mat that contains stereo tracking information regarding
% left-right view data LEFTANALYSISFILE & RIGHTANALYSISFILE generated with
% UGA.  CAMMATFILE corresponds to the stereo calibration data for the two
% views.  A prompt will appear for the save location of a .mat file.  This
% file contains the data needed to display tracking data using
% ResultsViewer_console.

%% Input Arguments
if(nargin == 0)
    [camName,camPath] = uigetfile('*.mat','Load Camera Cal File');
    camMatFile = fullfile(camPath,camName);
end

% Left Analysis File
if(nargin < 2)
   [leftName,leftPath] = uigetfile('*.mat','Load Left View Analysis File');
   leftAnalysisFile = fullfile(leftPath,leftName);
end

% Right Analysis File
if(nargin < 3)
    [rightName,rightPath] = uigetfile('*.mat','Load Right View Analysis File');
    rightAnalysisFile = fullfile(rightPath,rightName);
end

%% Load Stereo Calibration Parameters
% Camera Calibration
cameraParameters = load(camMatFile,'fc_left','fc_right',...
                            'cc_left','cc_right','om','T','kc_left',...
                            'kc_right');
[~,~,F] = calibMats(cameraParameters); 

%% Load Analysis Files
s_l = load(leftAnalysisFile); s_l = s_l.saveStruct;
s_r = load(rightAnalysisFile); s_r = s_r.saveStruct;

%% Grab Needed Structures
numFrames = s_l.numFrames;
i1_l = s_l.inst; i1_r = s_r.inst;

% Build Tracked Vector
tracked = (s_l.label == 0) & (s_r.label == 0);

% Setup Structures
%i1_l.epiLines = nan*ones(numFrames,3);
i1_l.corPt = nan*ones(numFrames,2);

%% Generate Point Correspondences
disp('Generating Point Correspondences'); 
%i1_l.epiLines(tracked,:) = (F * [i1_l.trackPt(tracked,:)';...
%                      ones(1,size(i1_l.trackPt(tracked,:),1))])';

% Compute Corresponding Point (Line Intersection)
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_l = getMid(i1_l);
midline1_r = getMid(i1_r);

% Compute Intersection of Midline & Epipolar Line
%C1 = cross(i1_l.epiLines(tracked,:),...
%           [cos(midline1_r(tracked,2)),sin(midline1_r(tracked,2)),...
%            -midline1_r(tracked,1)]);
%i1_l.corPt(tracked,:) = C1(:,1:2) ./ repmat(C1(:,3),1,2);
[dummy1,corPt,dummy] = orientPtCorr(i1_l.trackPt(tracked,:),...
    midline1_l(tracked,:),midline1_r(tracked,:),F,10);
i1_l.corPt(tracked,:) = corPt;

%% Save File
tracked = double(tracked);
[sName,sPath] = uiputfile('*.mat','Save Results View Data');
save(fullfile(sPath,sName),'i1_l','i1_r','tracked'); 



