function [error3d,startPos2d] = displacementErrorAnalysis(calFile,...
                            leftFile,rightFile,endOffset,numFrames,actDisp)
% DISPLACEMENTERRORANALYSIS Computes dispalcement error statistics
%
% [ERROR3D,STARTPOS2D] = displacementErrorAnalysis(CALFILE,LEFTFILE,
%                                   RIGHTFILE,ENDOFFSET,NUMFRAMES,ACTDISP)
% Computes displacement error statistics.  Uses stereo calibration data
% in CALFILE and UGT generated tracking data in LEFTFILE and RIGHTFILE to
% estimate instrument 3D trajectory.  3D instrument displacement is found
% between the first tracked frame and subset of frames near the end.  Let M
% be the number of tracked frames.  The subset over which statistics are
% calculated is defined as [M - ENDOFFSET,NUMFRAMES + (M - ENDOFFSET)].
% Statistics are calculated with respect to the actual 3D displacement (in
% mm) of the instrument given by ACTDISP.  Error statistics are returned
% in structure ERROR3D.  Track point starting position in the left and
% right trackfile are given in the structure STARTPOS2D.  See below for
% detailed description of both structures.  
%
% ERROR3D Structure
% Error is defined as the difference between the ACTDISP and calculated 3D
% displacement from the first tracked frame to frames in the subset defined
% by ENDOFFSET and NUMFRAMES.  For example, assume 1000 tracked frames,
% ENDOFFSET = 100, NUMFRAMES 50.  Then the instrument displacement will be
% calculated for frames [900,950].  The displacement is calculated as the
% position at these frames minus the position at the first tracked frame.
% 
% Structure Elements : maxError, minError, rangeError, meanError, rmsError
%
% STARTPOS2D Structure
% 3D trajectory is calculated from the tracked position of a single point
% on the instrument's imaged midlines in the left and right view video.
% STARTPOS2D gives information on the position of this point in the first
% tracked frame.  This is important for understanding the effect of in
% image position on accuracy
%
% Structure Elements : xs_l, ys_l, xs_r, ys_r

%% Get Transformation to Rectified Centered View
s = load(calFile,'R','T'); R = s.R; T = s.T;
[R_view,T_view] = rectifyCenterView(R,T);

%% Get 3D Trajectory Subset to be analyzed
% Triangulate 3D Trajectory
[startPt,~,~,~,sp2d_l,sp2d_r] = ...
    instTriang(calFile,leftFile,rightFile);

% Calculate Frame Subset over which trajectory is available %
startFrame = find(~isnan(startPt(:,1)),1,'first');

% The end frame is either first nan after start or last frame
firstNanOffset = find(isnan(startPt((startFrame+1):end,1)),1,'first');

if(isempty(firstNanOffset))
    % End frame is last frame
    endFrame = size(startPt,1);
else
    % End frame is (firstNan - 1) after start
    endFrame = startFrame + firstNanOffset - 1;
end

% Grab Signal Subsets
f = startFrame:endFrame;
startPt = startPt(f,:);
sp2d_l = sp2d_l(f,:);
sp2d_r = sp2d_r(f,:);

%% Displacement Error Analysis
% Assume the instrument has been completely displaced asnd is static at
% frames : [N - 2*numFrames,N - numFrames].  Displacement accuracy
% statistics will be calculated over this subset.  The displacement is
% be calculated between the first frame and every frame in this subset.  
N = size(startPt,1);
frameSub = (0:numFrames) + (N - endOffset);
[error3d,startPos2d] = analyzeSubset(frameSub,startPt,sp2d_l,sp2d_r,...
                                     R_view,T_view,actDisp);

function [error3d,startPos2d] = ...
    analyzeSubset(frameSub,startPt,sp2d_l,sp2d_r,R,T,actualDisp)
% ANALYZESUBSET Performs Statistical Analysis on 3D Signal
%
% [ERROR3D,STARTPOS2D] = analyzeSubset(FRAMESUB,STARTPT, SP2D_L,EP2D_L,
%                                      SP2D_R,EP2D_R,R,T)
% Computes 3D displacement statistics.  Uses (Nx3) matrix STARTPT
% representing the instrument position.  All statistics are computed over
% FRAMESUB.  Which is a vector of continuous frame numbers.  For example
% FRAMESUB = 300:400.  A displacement signal as found as the difference in
% instrument position between the first frame and those in FRAMESUB.  An
% error signal is computed as the displacement signal minus the actual
% signal.  This signal is the deviation from the actual instrument
% displacement.  See below for a description of statistics computed.
% STARTPOS2D is a structure that has the 2D starting coordinates of the
% instrument in the left and right view.

% Grab Subsets
initPt = startPt(1,:); % Initial Instrument Position
startPt = startPt(frameSub,:);

% Convert 3D Data into rectified centered view
initPt = ((R * initPt.') + T).';
startPt = ((R * startPt.') + repmat(T,1,size(startPt,1))).';

% 3D Stats
dispSig = sqrt(sum(...
                    (startPt - repmat(initPt,size(startPt,1),1)).^2,2));
errorSig = dispSig - actualDisp;
maxError = max(errorSig);
minError = min(errorSig);
rangeError = maxError - minError;
meanError = mean(errorSig,1);
rmsError = sqrt(mean(errorSig.^2));

% 3D Error Statistics
error3d = struct('maxError',maxError,'minError',minError,'rangeError',...
    rangeError,'meanError',meanError,'rmsError',rmsError);

% Starting Position Information
startPos2d = struct('xs_l',sp2d_l(1,1),'ys_l',sp2d_l(1,2),...
    'xs_r',sp2d_r(1,1),'ys_r',sp2d_r(1,2));