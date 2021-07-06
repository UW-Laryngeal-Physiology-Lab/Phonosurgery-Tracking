function [rhoFit,thetaFit,refitMSE,numInliers] = ...
                refitBoundLine(rho,theta,x,y,inlierThresh)
% REFITBOUNDLINE Refits hough based boundary lines using TLS
%
% [RHOFIT,THETAFIT] = refitBoundLine(RHO,THETA,X,Y,INLIERTHRESH) Uses
% coarse hough based estimation of boundary line {RHO,THETA} (THETA in
% radians) to fit a boundary line using Linear Total Least Squares.
% Boundary line candidate points in image (X,Y) are tested against
% INLIERTHRESH (squared distance from coarse line).  All found inlier
% points are used to fit the boundary line given by (RHOFIT,THETAFIT).

% Note the line parameterization (rho,theta) treats image indices in (0,0)
% origin coordinate system.  Therefore (-1) must be added to matrix 
% indices and (+1) must be added to coordinates to get indices.

% Hough Line Parameterization & Inlier Finding
numPts = numel(x);
l = [cos(theta),sin(theta),-rho];
e = (l * [(x(:)-1).';(y(:)-1).';ones(1,numPts)]).^2;
inliers = e < inlierThresh;
numInliers = sum(inliers);

% TLS
[rhoFit,thetaFit,refitMSE] = fitLineTLS(x(inliers)-1,y(inliers)-1);
%{
A = [x(inliers)-1,y(inliers)-1,ones(numInliers,1)];
[~,~,V] = svd(A);
l_fit = V(:,3) / sqrt(sum(V(1:2,3).^2));

% Return Result with Proper Sign
if(l_fit(3) > 0)
    l_fit = -l_fit;
end
rhoFit = -l_fit(3);

thetaFit = acos(abs(l_fit(1)));
if(l_fit(2) < 0)
    thetaFit = -thetaFit;
end

if(nargout >= 3)
    % Refitting MSE
    refitMSE = mean((A*l_fit).^2);
end
%}
    
    


