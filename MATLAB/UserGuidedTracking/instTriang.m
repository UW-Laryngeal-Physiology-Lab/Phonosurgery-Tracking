function [startPt,endPt,startOr,endOr,sp2D_l,sp2D_r,ep2D_l,ep2D_r] = ...
                                                       instTriang(varargin)
% INSTTRIANG Triangulates Instrument position from UGA data
%
% [STARTPT,ENDPT,STARTOR,ENDOR.SP2D_L,SP2D_R,EP2D_L,EP2D_R] = instTriang()
% Triangulates the position of a phonosurgery instrument that was tracked
% using UGA.  Prompts the user for the location of .mat files  containing
% the stereo camera calibration data, the left view tracking data and the
% right view tracking data.  Four (N x 3) matrices are returning consisting
% of row vectors of the form [X,Y,Z].  Movement of the instrument is given
% by STARTPT and ENDPT. Where ENDPT - STARTPT gives the movement in 3D
% between frameS.  For example the displacement between frame N and N+1 is
% given by ENDPT(N) - STARTPT(N). STARTOR and ENDOR is the normalized
% instrument orientation in 3D. NaN is used as a placeholder for data at
% frames at which there is no tracking data.  SP2D_L and SP2D_R are image
% start points in the left(L) and right(R) view.  EP2D_L and EP2D_R are the
% corresponding endpoints.  For example the 2D displacement of the
% instrument between frame N and N+1 in the left view is given by 
% EP2D_L(N) - SP2D_L(N).
%
% [...]= instTriang(CAMCALFILE,LEFTANALYSISFILE,
%                                            RIGHTANALYSISFILE)
% Uses locations of stereo camera calibration file CAMCALFILE, and left and
% right view UGA analysis files LEFTANALYSISFILE and RIGHTANALYSISFILE.

%% Input Parser
if(nargin == 3)
    camCalFile = varargin{1};
    leftAnalysisFile = varargin{2};
    rightAnalysisFile = varargin{3};
else
    % Load Camera Calibration File
    [calName,calPath] = uigetfile('*.mat',...
        'Load Camera Calibration Parameters');
    camCalFile = fullfile(calPath,calName);
    
    % Load Analysis Files
    [fName_l,pName_l] = uigetfile('*.mat','Load Left View Analysis File',calPath);
    [fName_r,pName_r] = uigetfile('*.mat','Load Right View Analysis File',pName_l);
    
    leftAnalysisFile = fullfile(pName_l,fName_l);
    rightAnalysisFile = fullfile(pName_r,fName_r);
end

% Read Camera Parameters
cameraParameters = load(camCalFile);
[PL,PR,F] = calibMats(cameraParameters);

% Read Analysis Parameters
s_l = load(leftAnalysisFile); s_l = s_l.saveStruct;
s_r = load(rightAnalysisFile); s_r = s_r.saveStruct;

% Grab Inst Structs
numFrames = s_l.numFrames;
m1_l = s_l.mark; m1_r = s_r.mark;
i1_l = s_l.inst; i1_r = s_r.inst;

% Utilize joint endIndex if the number of frames in the left and right
% video are different.  This was added due to an issue with a single frame
% discrepency between the left and right video.  By using the joint end the
% the smaller subset of frames is used.
jointEnd = min([numel(s_l.label),numel(s_r.label)]);

gFrames = (s_l.label(1:jointEnd-1) == 0) & ...
          (s_r.label(1:jointEnd-1) == 0) & ...
          (s_l.label(2:jointEnd) == 0) & ...
          (s_r.label(2:jointEnd) == 0);
startFrames = find(gFrames);
frames = startFrames;
endFrames = startFrames + 1;

% Setup Structures
startPt = nan*ones(numFrames,3);
startOr = nan*ones(numFrames,3);
endPt = nan*ones(numFrames,3);
endOr = nan*ones(numFrames,3);

startPt_l = nan*ones(numFrames,2);
endPt_l = nan*ones(numFrames,2);
startOrient_l = nan*ones(numFrames,2);
endOrient_l = nan*ones(numFrames,2);

startPt_r = nan*ones(numFrames,2);
endPt_r = nan*ones(numFrames,2);
startOrient_r = nan*ones(numFrames,2);
endOrient_r = nan*ones(numFrames,2);

%% Generate Point Correspondences
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_l = getMid(i1_l);
midline1_r = getMid(i1_r);

% Left View Start & End
startPt_l(frames,:) = i1_l.trackPt(startFrames,:);
endPt_l(frames,:) = startPt_l(startFrames,:) + ...
    (m1_l.subT(endFrames,:) - m1_l.subT(startFrames,:)); 

% Compute Correspondences
[start_op_l,start_tp_r,start_op_r] = orientPtCorr(startPt_l(...
    startFrames,:),midline1_l(startFrames,:),midline1_r(startFrames,:),...
    F,10);
[end_op_l,end_tp_r,end_op_r] = orientPtCorr(endPt_l(...
    frames,:),midline1_l(endFrames,:),midline1_r(endFrames,:),F,10);

% Left View Orient & End Point
startOrient_l(frames,:) = start_op_l;
endOrient_l(frames,:) = end_op_l;

% Right View Start, End & Orient
startPt_r(frames,:) = start_tp_r;
endPt_r(frames,:) = end_tp_r;
startOrient_r(frames,:) = start_op_r;
endOrient_r(frames,:) = end_op_r;

%% Run Triangulation
triangFun = @(left,right) stereo_triangulation(left',right',...
    cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

[startPt_sub,dummy] = triangFun(startPt_l(frames,:),...
                                startPt_r(frames,:));
[startOr_sub,dummy] = triangFun(startOrient_l(frames,:),...
                                startOrient_r(frames,:));
[endPt_sub,dummy] = triangFun(endPt_l(frames,:),...
                              endPt_r(frames,:));
[endOr_sub,dummy] = triangFun(endOrient_l(frames,:),...
                              endOrient_r(frames,:));

% Get 3D Points
startPt(frames,:) = startPt_sub';
endPt(frames,:) = endPt_sub';

% Get 3D Orientation Vectors
getOrient = @(or_sub,pt_sub) (or_sub' - pt_sub') ./ ...
    repmat(sqrt(sum((or_sub' - pt_sub').^2,2)),1,3);
startOr(frames,:) = getOrient(startOr_sub,startPt_sub);
endOr(frames,:) = getOrient(endOr_sub,endPt_sub);

% 2D Outputs
sp2D_l = startPt_l; sp2D_r = startPt_r;
ep2D_l = endPt_l; ep2D_r = endPt_r;
