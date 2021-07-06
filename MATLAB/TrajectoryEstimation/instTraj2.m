function [pos3D,orient3D,pos2D_L,pos2D_R,orient2D_L,orient2D_R] = ...
                    instTraj2(inst_L,label_L,inst_R,label_R,camCalFile)
% Estimates instrument trajectory using UGT structs and cal file
% [POS3D,ORIENT3D,POS2D_L,POS2D_R,ORIENT2D_L,ORIENT2D_R] = ...
%                       instTraj2(INST_L,LABEL_L,INST_R,LABEL_R,CAMCALFILE)
% Estimates trajectory of a phonomicrosurgery instrument whose features
% were tracked using UGT. INST_L and INST_R are the inst structs (contains
% trackPt, rho, and theta fields) from the left and right view UGT tracking
% data files.  LABEL_L and LABEL_R are the label arrays from the left and
% right view UGT tracking data files.  These arrays must be the same size
% otherwise the function will return empty output arguments.  CAMCALFILE is
% a the location of a .mat file containing the stereo camera calibration
% parameters.  This file should be generated using the
% cameraCalibrationToolbox function stereo_gui.  3D instrument position and
% orientation vector per frame are returned in POS3D and ORIENT3D.  POS2D_L
% and POS2D_R are the left video track point and its correspondence point
% in the right video.  ORIENT2D_L and ORIENT2D_R are the points in the left
% and right view used to form the orientation vector with the left track
% point and correspondence point. See below for detailed descriptions of
% all returned arrays.
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

%% Load Calibration Data
cameraParameters = load(camCalFile);
[PL,PR,F] = calibMats(cameraParameters);

%% Init Data Structures
% Verify label arrays are the same length
if(numel(label_L) ~= numel(label_R))
    disp(['Error.  Label Structures are not same size.  ',...
          'Returning empty output arguments']);
      pos3D = []; orient3D = []; 
      pos2D_L = []; pos2D_R = [];
      orient2D_L = []; orient2D_R = [];
      return;
end

numFrames = numel(label_L);

% Determine frames with correct tracking parameters
gFrames = (label_L == 0) & (label_R == 0);
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

%% Generate Point Correspondence
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_L = getMid(inst_L);
midline1_R = getMid(inst_R);

% Left View Start & End
pos2D_L(frames,:) = inst_L.trackPt(frames,:);

% Compute Orient and Correspondence Points WRT Left Track Point
[op_L,tp_R,op_R] = orientPtCorr(pos2D_L(frames,:),midline1_L(frames,:),...
                                                midline1_R(frames,:),F,10);
% Left View Orient
orient2D_L(frames,:) = op_L;

% Right View Track Point Correspondence & Orient
pos2D_R(frames,:) = tp_R;
orient2D_R(frames,:) = op_R;
%% Run Midpoint Triangulation
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