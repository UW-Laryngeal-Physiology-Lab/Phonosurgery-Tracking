function [x1_op,x2,x2_op] = orientPtCorr(x1,mid1,mid2,F,alpha)
% ORIENTPTCORR Generates two point correspondences used for 3D orient est.
%
% [X1_OP,X2,X2_OP] = orientPtCorr(X1,MID1,MID2,F,ALPHA) Finds two view pt
% correspondence pair (X1,X2) and (X1_OP,X2_OP).  The pt variables are (N x
% 2) matrices of the form [x,y].  MID1 and MID2 are (N x 2) matrices of the
% form [rho,theta] that quantify the 2D orientation associated with each
% point.  F is the fundamental matrix related the two view calibration.
% Alpha is a parameter that defines the distance between the source pts
% (X1,X2) and their orientation points.  The orientation points are
% calculated such that they exist on their associate midline at a vertical
% position lower than that of the source point.

% View 1 Orientation Point
x1_op = x1 + ...
        repmat(alpha*sign(cos(mid1(:,2))),1,2) .* ...
        [-sin(mid1(:,2)),cos(mid1(:,2))];

% View 2 Source Point
x2 = midEpiIntersect(x1,F,mid2);

% View 2 Orientation Point
x2_op = midEpiIntersect(x1_op,F,mid2);