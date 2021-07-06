function [theta] = orientFromMark(frameIm,wCorn,wStats,weighting,pStruct)
% ORIENTFROMMARK Estimates Instrument Orientation from the Marker Region
%
% THETA = orientFromMark(FRAMEIM,WCORN,WSTATS,PC,WEIGHTING,PSTRUCT)  Uses
% horizontal gradient information from the marker region to estimate the
% orientation of the instrument where theta is the orientation in radians.
% In order to find the angle of the hough line associated with the
% orientation subtract pi/2 from theta.  PSTRUCT is an optional structure
% of parameters see below for description.
%
% PSTRUCT fields:
% gmiThresh : Gradient Magnitude Threshold used to locate interest edge
% points on the marker (default = 10)
% sElem : Structuring element applied duriong interest pixel detection.
% (default = [1,1,1;1,0,1;1,1,1])

%% Input Arguments
if(nargin == 4)
    pStruct = struct();
end

% Default Parameters
if(~isfield(pStruct,'gmiThresh'))
    pStruct.gmiThresh = 5;
end

if(~isfield(pStruct,'sElem'))
    pStruct.sElem = [1,1,1;1,0,1;1,1,1];
end

%% Run Algorithm
regMask = false(size(frameIm));
[gmi,gx,gy] = computeGradientMagIm(frameIm);

% Compute Region Mask
regMask((0:wStats.wSize-1) + wCorn(1,2),...
        (0:wStats.wSize-1) + wCorn(1,1)) = weighting;

% Get Interest Pixels
iMask = (abs(gx) > abs(gy)) & regMask & (abs(gmi) > pStruct.gmiThresh);
iMask2 = imopen(iMask,pStruct.sElem);

% Compute Orientation & Estimate Theta
orient = atan(gx(iMask2) ./ gy(iMask2));
o2 = orient;
o2(orient < 0) = pi + orient(orient < 0);
theta = median(o2);
