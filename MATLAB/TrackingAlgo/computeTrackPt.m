function trackPt = computeTrackPt(rho,theta,mCorn,halfSize)
% COMPUTETRACKPT Computes tracking point in marker window
%
% TRACKPT = computeTrackPt(RHO,THETA,MCORN,HALFSIZE) Computes tracking
% point in marker window as the intersection of the windows vertical point
% and instrument midline.  The instrument midline is computed as the mean
% of the instrument's boundary line parameters RHO & THETA.  MCORN is the
% corner of the marker window and HALFSIZE the halfSize of the marker
% window.  

% Note the line parameterization (rho,theta) treats image indices in (0,0)
% origin coordinate system.  Therefore (-1) must be added to matrix 
% indices and (+1) must be added to coordinates to get indices.

midRho = mean(rho); midTheta = mean(theta);

yPt = (mCorn(1,2) + halfSize)-1;
trackPt = [1 + ((midRho - yPt*sin(midTheta))/cos(midTheta)),yPt+1]; 