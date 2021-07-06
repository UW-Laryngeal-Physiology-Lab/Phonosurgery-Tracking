function [R_new,T_new] = rectifyCenterView(R,T)
% RECTIFYCENTERVIEW Computes transform from left cam to centerRect view
%
% [RNEW,TNEW] = rectifyCenterView(R,T) Rigid body transform parameters R
% and T are taken to be the transformation from the left camera position to
% the right camera position of a stereo pair. An "intermediate" view exists
% at the midpoint of the baseline between cameras.  The orientation of this
% view is that of the left and right cameras after applying a rectifying
% transformation (assuming equal intrinsics).  This function computes the
% rigid body transformation from the left camera frame to the intermidate
% camera frame in RNEW and TNEW.  That is : 
% X_intermedite = RNEW *X_leftCamera + TNEW.

% See Fusiello et al.  A compact algorithm for rectification of stereo
% pairs for the algorithm used to compute the rectifying transformation.

% Find right camera center in left camera frame (world)
cRight = -R.' * T;

% Find New Camera Axis in left camera frame
r1 = cRight ./ norm(cRight);
r2 = cross([0;0;1],r1);
r3 = cross(r1,r2);

R_new = [r1.';r2.';r3.'];

% Baseline midpoint in left camera frame is translation
T_new = -R_new * (cRight/2);