function pMask = predictMask(mCorn,markDisp,prevTrackPt,prevRho,...
                                                        prevTheta,imSize)
% PREDICTMASK Predicts instrument region mask position
%
% PMASK = predictMask(MCORN,MARKDISP,PREVTRACKPT,PREVRHO,PREVTHETA,IMSIZE)
% Uses the upper left hand corner of the marker in the current frame MCORN,
% the marker window displacement between the current and previous frame MARKDISP, 
% the trackPt position in the previous frame PREVTRACKPT, the previous
% frame boundary line esimates PREVRHO and PREVTHETA to predict the postion
% of a mask containing the instrument's cylindrical region in the current
% frame.  A binary matrix of IMSIZE is returned, where 1's indicate mask
% pixels.
%
% An inverse transformation is used to determine mask pixels.  


% Binary Mat that will contain the mask
pMask = false(imSize);

% Propagate TrackPt using Marker Window Displacement
newTrackPt = prevTrackPt + markDisp;

% Determine coordinate transformation such that midline is oriented
% vertically and newTrackPt position does not change.
theta = -mean(prevTheta);
R = [cos(theta),-sin(theta);sin(theta),cos(theta)];
translation = -((R*newTrackPt') - newTrackPt');

% Define the subregion in Current Frame for Mask %
% The top of the marker window bounds the mask position vertically.  The
% horizontal frame edges bound the mask position horizontally.
pSub = false([mCorn(2),imSize(2)]);
[X,Y] = meshgrid(1:imSize(2),1:mCorn(2));

% Transform the subregion coordinate system such that the midline is
% oriented vertically and the newTrackPt position does NOT change.
transformVals = R*[X(:)';Y(:)'] + repmat(translation,1,numel(X));

% Points whose transformation falls within xBounds define the mask
xBounds = newTrackPt(1) + [-0.75,0.75]*abs(diff(prevRho));
pSub(and(transformVals(1,:) > xBounds(1),transformVals(1,:) < xBounds(2))) = 1;
pMask(1:mCorn(2),1:imSize(2)) = pSub;
