function cornerDex = getCornerDex(matSize)
% GETCORNERDEX Gets the Indices of the Corners of a 2D Matrix
%
% CORNERDEX = getCornerDex(MATSIZE) Returns the linear indices
% corresponding to the corners of a 2D matrix of size MATSIZE
% [numRows,numCols].  CORNERDEX is (4x1) array corresponding to corners
% [upperLeft,lowerLeft,upperRight,lowerRight].

% 2011-02-28 KS
cornerDex = [1;matSize(1);1+matSize(1)*(matSize(2)-1);...
    matSize(1)*matSize(2)];
