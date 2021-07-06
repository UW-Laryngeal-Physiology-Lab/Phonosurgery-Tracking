function [leftEdge,rightEdge] = boundEdgeDetect(frameIm,backIm,pStruct,subIm)
% BOUNDEDGEDETECT Detects points that look like instrument boundary edges
%
% [LEFTEDGE,RIGHTEDGE] = boundEdgeDetect(FRAMEIM,BACKIM,PSTRUCT,SUBIM)
% Detects pixels in FRAMEIM with properties like left and right
% instrument edges.  LEFTEDGE and RIGHTEDGE are binary images.  Background
% subtraction is performed using static background image BACKIM.  Algorithm
% behavior is controlled by parameters in PSTRUCT, see below for
% description.  Pixels can be detected within a rectangular subimage of
% FRAMEIM by provided SUBIM.  It should be an array with the form [x_left,
% x_right,y_top,y_bottom].  
%
% First, gradient is computed on frame image.  Then, background pixels are
% thresholded out.  Next, horizontal edges are thresholded. Following this
% weak edges are thresholded out using pStruct.edgeThresh.  Non-maxima
% suppression is applied to thin remaining vertical edges. Final set of
% edges are grouped into left and right edges depending on gradient
% orientation.  Algorithm expects instrument pixels to be dark and
% background pixels to be light.
%
% PSTRUCT Fields :
% backThresh : Threshold used for background subtraction.
% edgeThresh : Magnitude edge threshold

%% Generate Background Mask
% Either use entire image or subimage
if(nargin == 3)
    backMask = (double(backIm) - double(frameIm)) < pStruct.backThresh;
    %[~,backMask] = genCandIm(frameIm,backIm,pStruct.backThresh,...
    %                        pStruct.w_back,pStruct.w_dark);
elseif(nargin == 4)
    % Operate on a SubImage
    backMask = (double(backIm(subIm(3):subIm(4),subIm(1):subIm(2))) - ...
        double(frameIm(subIm(3):subIm(4),subIm(1):subIm(2)))) < pStruct.backThresh;
    
    %[~,backMask] = genCandIm(...
    %    frameIm(subIm(3):subIm(4),subIm(1):subIm(2)),...
    %    backIm(subIm(3):subIm(4),subIm(1):subIm(2)),...
    %    pStruct.backThresh,pStruct.w_back,pStruct.w_dark);
end

%% Use Thresholding techniques to locate possible instrument edge pixels
% Compute the Gradient
if(nargin == 3)
    % Use Full Image
    [gmi,gx,gy] = computeGradientMagIm(frameIm);
else
    % Use Subimage
    [gmi,gx,gy] = computeGradientMagIm(frameIm(subIm(3):subIm(4),subIm(1):subIm(2)));
end

% Remove background pixels
gmi(backMask) = 0;

% Threshold out horizontal edges
gmi(abs(gx) > abs(gy)) = 0;

% Threshold out weak edges 
edgeMask = abs(gmi) > pStruct.edgeThresh;
imSize = size(edgeMask);

%% Non-Maxima Suppression
% Clear Vertical Boundaries
edgeMask(1:imSize(1),1) = 0;
edgeMask(1:imSize(1),end) = 0;
lMaskIdx = find(edgeMask);

% Compute Neighbors
neighborsMat = repmat([-imSize(1),imSize(1)],numel(lMaskIdx),1) + ...
               repmat(lMaskIdx,1,2);
trueNeighbor = (neighborsMat > 0 & neighborsMat <= numel(edgeMask) & edgeMask(neighborsMat));
idxVec = [1;2];

% Perform NMS
notMax = true(numel(lMaskIdx),1);
for k = 1:numel(lMaskIdx)
    gradVal = gmi(lMaskIdx(k));
    if(gradVal >= max(gmi(neighborsMat(k,idxVec(trueNeighbor(k,:))))))
        notMax(k) = 0;
    end
end

% DEBUG NMS
%{
removeMask(lMaskIdx(notMax)) = 1;
%}

% Remove Non-Maxima
edgeMask(lMaskIdx(notMax)) = 0;

%% Group Pixels into left and right edges
if(nargin == 3)
    leftEdge = edgeMask & (gy > 0);
    rightEdge = edgeMask & (gy < 0);
elseif(nargin == 4)
    % Fill the SubIm
    leftEdge = false(size(frameIm));
    rightEdge = leftEdge;
    
    leftEdge(subIm(3):subIm(4),subIm(1):subIm(2)) = edgeMask & (gy > 0);
    rightEdge(subIm(3):subIm(4),subIm(1):subIm(2)) = edgeMask & (gy < 0);
end
