function [mStruct,hStruct,rho,theta] = imProEval(mark,inst,instVid,backVid,frames,frameDex,params)
% IMPROEVAL Returns ImPro Algorithm Masks for a frame after track/detect
%
% Usage:
% [MSTRUCT,HSTRUCT,RHO,THETA] = imProEval(MARK,INST,INSTVID,BACKVID,
%                                                               FRAMEDEX)
% [MSTRUCT,HSTRUCT,RHO,THETA] = imProEval(MARK,INST,INSTVID,BACKVID,
%                                                       FRAMES,FDEX,PARAMS)
% Computes data related to single frame of instrument detection/tracking.
% Uses outputs of detectTrackFun corresponding to a single instrument in
% single view: MARK (marker struct) & INST (instrument struct).  INSTVID is
% the videoReader object corresponding to the video the tracking was
% performed on from FRAMES(1) to FRAMES(2).  FDEX is the relative index @
% which the data is generated for.  The first frame of BACKVID is used for
% background subtraction.  Algorithm parameters can be set using structure
% PARAMS.  See below for a description of the structure.  MSTRUCT is is
% struct with masks corresponding to different image processing steps of
% instrument tracking/detection.  HSTRUCT contains hough transform data and
% {RHO,THETA} are the boundary line parameters.  See below for a
% description of MSTRUCT.
%
% PARAMS (by category)
% Boundary Line Fitting : {thetaResolution,rhoResolution}
% Image Processing : {backThresh,wBack,wDark,instThresh,diskR,fattenR,
%                     instThresh}
%
% MSTRUCT fields
% blobMask : Instrument blobs detected in frame
% pMask : Prediction mask for specific instrument
% candMask : Boundary Line Candidate Mask
% lineCandMask : Line Candidates prior to Hough based boundary line fitting
% lineMask: Lines corresponding to fitted boundaries and midline

%% Input Arguments
if nargin == 6
    % No Parameter Structure 
    params = struct();
end

pStruct = checkParameters(params);
trackParams = pStruct;
trackParams.suppress = getSuppress(21,9.9,pStruct);

%% Run ImPro Routines
k = frameDex;

% Grab Frames
backIm = rgb2gray(backVid.read(1)); %backIm(:,:,2:3) = [];
frameIm = rgb2gray(instVid.read(frames(1) + (k-1))); %frameIm(:,:,2:3) = [];
imSize = size(backIm);

[candIm,backMask] = genCandIm(frameIm,backIm,pStruct.backThresh,...
                              pStruct.wBack,pStruct.wDark);
blobIm = instBlobDetect(candIm,backMask,...
                            pStruct.instThresh,pStruct.diskR);

% Prediction Mask
markDisp = mark.mCorn(k,:) - mark.mCorn(k-1,:);
pMask = predictMask(mark.mCorn(k,:),markDisp,...
                    inst.trackPt(k-1,:),inst.rho(k-1,:),...
                    inst.theta(k-1,:),imSize);

% Boundary Line Estimation
param1 = trackParams;
param1.thetaVals = (-5:trackParams.thetaResolution:5) + ...
                    round((180/pi) * mean(inst.theta(k-1,:)));
[~,rho,theta,~,tStruct] = blobToBound(and(blobIm,pMask),candIm,param1);

% Build Structures
hStruct = tStruct.houghStruct;
mStruct = struct('blobMask',blobIm,'pMask',pMask,'candMask',...
                 tStruct.candMask,'lineCandMask',tStruct.lineCandMask,...
                 'lineMask',drawLineMask(imSize,rho,theta));
end

%% Utility Functions

function pStruct = checkParameters(paramStruct)
% CHECKPARAMETERS Checks parameter structure
%
% PSTRUCT = checkParameters(PARAMSTRUCT) Checks PARAMSTRUCT for fields
% utilized by algorithms in this function.  If the field is missiing in
% PARAMSTRUCT it is set to the default values.

pStruct = paramStruct;

% Check for Needed Marker Fields
pStruct = checkField(pStruct,'Update','None');
pStruct = checkField(pStruct,'Prediction','Constant Velocity');

% Check for Needed Boundary Line Fitting Fields
pStruct = checkField(pStruct,'thetaResolution',1);
pStruct = checkField(pStruct,'rhoResolution',1);

% Check for Needed Image Processing Fields
pStruct = checkField(pStruct,'backThresh',30);
pStruct = checkField(pStruct,'wBack',1.5);
pStruct = checkField(pStruct,'wDark',1);
pStruct = checkField(pStruct,'instThresh',280);
pStruct = checkField(pStruct,'diskR',5);
end

function paramStruct = checkField(paramStruct,fieldName,defaultVal)
% CHECKFIELD Checks for field in structure
%
% PARAMSTRUCT = checkField(PARAMSTRUCT,FIELDNAME,DEFAULTVAL) Checks if
% PARAMSTRUCT has field FIELDNAME.  If it does not it is added to
% PARAMSTRUCT with value DEFAULTVAL.

if(~isfield(paramStruct,fieldName))
    paramStruct.(fieldName) = defaultVal;
end
end

function sup = getSuppress(rhoSize,thetaSize,paramStruct)
% GETSUPPRESS Computes Hough Peaks Suppression Mask Size
%
% SUP = getSuppress(RHOSIZE,THETASIZE,PARAMSTRUCT) Computes suppresion mask
% for Hough Peak detection.  RHOSIZE and THETASIZE are the mask dimension
% size in pixels and degrees respectively.  The appropriate resolution
% fields of PARAMSTRUCT is utilized to calcuate the actual mask size [rho
% theta] that is returned as SUP.

oddRho = getOdd(round(rhoSize/paramStruct.rhoResolution));
oddTheta = getOdd(round(thetaSize/paramStruct.thetaResolution));
sup = [oddRho,oddTheta];

end

function oddVal = getOdd(val)
% GETODDVAL Finds nearest odd integer
%
% ODDVAL = getOdd(VAL)

if(val/2 == round(val/2))
    % It's even, convert it
    oddVal = val-1;
else
    % It's Odd
    oddVal = val;
end
end
