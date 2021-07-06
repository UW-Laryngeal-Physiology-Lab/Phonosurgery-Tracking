function [pos3D,orient3D,pos2D_L,pos2D_R,orient2D_L,orient2D_R] = ...
                                                        instTraj(varargin)
% INSTTRARJ Estimates instrument trajecotry from UGT data
%
% [POS3D,ORIENT3D,POS2D_L,POS2D_R,ORIENT2D_L,ORIENT2D_R] = instTraj()
% Estimates trajectory of a phonomicrosurgery instrument whose features
% were tracked using UGT. Prompts the user for the location of .mat files
% containing the stereo camera calibration data, the left view tracking
% data and the right view tracking data.  3D instrument position and
% orientation vector per frame are returned in POS3D and ORIENT3D.  POS2D_L
% and POS2D_R are the left video track point and its correspondence point
% in the right video.  ORIENT2D_L and ORIENT2D_R are the points in the left
% and right view used to form the orientation vector with the left track
% point and correspondence point. See below for detailed descriptions of
% all returned arrays.
%
% [...]= instTraj(CAMCALFILE,LEFTANALYSISFILE,RIGHTANALYSISFILE)
% Uses locations of stereo camera calibration file CAMCALFILE, and left and
% right view UGA analysis files LEFTANALYSISFILE and RIGHTANALYSISFILE.
%
% 3D Position :
% The returned array POS3D represents the 3D trajectory of the instrument.
% Specifically, it is the 3D position of the track point in the left camera
% video.  POS2D_L gives the 2D position of the track point in the left
% camera video.  POS2D_R gives the correspondence point in the right video
% used to triangulate 3D position.  This correspondence point is found by
% intersecting the epipolar line associated with the left video trackpoint
% with the instrument's imaged midline in the right video.
%
% 3D Orientation Vector :
% The returned array ORIENT3D represents the instrument's orientation
% vector.  This is a 3D unit vector aligned with the axis of the
% instrument. The vector's origin is the 3D position point of the
% instrument returned in POS3D.  The image of this vector is given by the
% two points POS2D_L and ORIENT2D_L in the left video and POS2D_R and
% ORIENT2D_R in the right video.  Note the images of this vector are not
% of unit length.  
%
% RETURNS
% POS3D : (N x 3) array where N is number of video frames tracking was
% performed over.  Each row has the form [X Y Z], where [X Y Z] is the 3D
% position of the instrument in the left camera frame at an individual
% frame.  Position is only reported for frames corresponding to the
% instrument being correctly tracked (in both views).  If the instrument is
% not correctly tracked in a frame, the position coordinates are valued
% NaN.
%
% ORIENT3D : (N x 3) array.  Each row has the form [dX dY dZ], where 
% [dX dY dZ] is the 3D orientation vector relative to POS3D in the left
% camera frame.  A value is only reported for frames corresponding to the
% instrument being correctly tracked (in both views).  If the instrument is
% not correctly tracked in a frame, the coordinates are valued NaN.
%
% POS2D_L : (N x 2) array.  Each row has the form [x y], where [x y] is the
% 2D position of the instrument's track point in the left video.  Position
% is only reported for frames corresponding to the instrument being
% correctly tracked.  If the instrument is not correctly tracked in a
% frame, the position coordinates are valued NaN.
%
% POS2D_R : (N x 2) array.  Each row has the form [x y], where [x y] is the
% 2D position of the point in the right video that corresponds to the track
% point in the left video.Position is only reported for frames
% corresponding to the instrument being correctly tracked.  If the
% instrument is not correctly tracked in a frame, the position coordinates
% are valued NaN.
%
% ORIENT2D_L : (N x 2) array.  Each row has the form [dX dY], where [dX dY]
% is the image of the orientation vector relative to POS2D_L in the left
% video.  Note, this is not a unit vector.  A value is only reported for
% frames corresponding to the instrument being correctly tracked.  If the
% instrument is not correctly tracked in a frame, the coordinates
% are valued NaN.
%
% ORIENT2D_R : (N x 2) array.  Each row has the form [dX dY], where [dX dY]
% is the image of the orientation vector relative to POS2D_R in the right
% video.  Note, this is not a unit vector.  A value is only reported for
% frames corresponding to the instrument being correctly tracked.  If the
% instrument is not correctly tracked in a frame, the coordinates
% are valued NaN.

%% Check Inputs
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

%% Load Data
% Read Camera Parameters
cameraParameters = load(camCalFile);
[PL,PR,F] = calibMats(cameraParameters);

% Read Tracking Parameters
s_l = load(leftAnalysisFile); s_l = s_l.saveStruct;
s_r = load(rightAnalysisFile); s_r = s_r.saveStruct;

% Grab Inst Structs
numFrames = s_l.numFrames;
m1_l = s_l.mark; m1_r = s_r.mark;
i1_l = s_l.inst; i1_r = s_r.inst;

%% Initialize Data Structures
% Utilize joint endIndex if the number of frames in the left and right
% video are different.  This was added due to an issue with a single frame
% discrepency between the left and right video.  By using the joint end the
% the smaller subset of frames is used.
jointEnd = min([numel(s_l.label),numel(s_r.label)]);

% Determine frames with correct tracking parameters
gFrames = ((s_l.label(1:jointEnd) == 0) & (s_r.label(1:jointEnd) == 0)); 
frames = find(gFrames);

% Initialize Structures %
% 3D Position & Orientation 
pos3D = nan*ones(numFrames,3);
orient3D = nan*ones(numFrames,3);

% 2D Position & Orientation
pos2D_L = nan*ones(numFrames,2);
orient2D_L = nan*ones(numFrames,2);
pos2D_R = nan*ones(numFrames,2);
orient2D_R = nan*ones(numFrames,2);

%% Generate Point Correspondences
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_l = getMid(i1_l);
midline1_r = getMid(i1_r);

% Left View Start & End
pos2D_L(frames,:) = i1_l.trackPt(frames,:);

% Compute Orient and Correspondence Points WRT Left Track Point
[op_l,tp_r,op_r] = orientPtCorr(pos2D_L(frames,:),midline1_l(frames,:),...
                                                midline1_r(frames,:),F,10);
% Left View Orient
orient2D_L(frames,:) = op_l;

% Right View Track Point Correspondence & Orient
pos2D_R(frames,:) = tp_r;
orient2D_R(frames,:) = op_r;

%% Run Triangulation
triangFun = @(left,right) stereo_triangulation(left',right',...
    cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

[pos3D_sub,dummy] = triangFun(pos2D_L(frames,:),...
                                pos2D_R(frames,:));
[orient3D_sub,dummy] = triangFun(orient2D_L(frames,:),...
                                orient2D_R(frames,:));

% Get 3D Points
pos3D(frames,:) = pos3D_sub';

% Get 3D Orientation Vectors
getOrient = @(or_sub,pt_sub) (or_sub' - pt_sub') ./ ...
    repmat(sqrt(sum((or_sub' - pt_sub').^2,2)),1,3);
orient3D(frames,:) = getOrient(orient3D_sub,pos3D_sub);
