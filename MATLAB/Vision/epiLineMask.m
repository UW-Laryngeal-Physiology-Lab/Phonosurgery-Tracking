function epiMask = epiLineMask(sourcePt,F,imSize)
% EPILINEMASK Creates Mask of Bresenham Drawn Epipolar Line
%
% EPIMASK = epiLineMask(SOURCEPT,F,IMSIZE) Creates mask image EPIMASK with
% Epiplor line corresponding to SOURCEPT and fundamental matrix F.  F takes
% a point in the source image into an line in the destination image. The
% returned mask has dimensions IMSIZE.  The line is drawn using Brenham's
% Line Drawing Algo.

% Compute Epipolar Line
eline = F * [sourcePt(:);1];

% Determine Line End Points
if (abs(eline(1) / eline(2)) <= 1)
    % Use X Based Parameterization
    X = [1,imSize(2)];
    Y = -(eline(1)*X + eline(3))/eline(2);
else
    Y = [1,imSize(1)];
    X = -(eline(2)*Y + eline(3))/eline(1);
end

% Get Valid Line Indexes into Mask Image
ils = intLineSeg([X(1),Y(1)],[X(2),Y(2)]);
[xIdx,yIdx] = ils.getIntLocs();
goodIdx = and(yIdx >= 1,yIdx <= imSize(1));

% Draw Mask
epiMask = false(imSize);
epiMask(sub2ind(imSize,yIdx(goodIdx),xIdx(goodIdx))) = 1;

    
    