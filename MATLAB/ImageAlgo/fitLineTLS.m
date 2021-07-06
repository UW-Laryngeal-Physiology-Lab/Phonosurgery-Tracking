function [rhoFit,thetaFit,fitMSE] = fitLineTLS(x,y)
% FITLINETLS Fits line using total least squares
%
% [RHOFIT,THETAFIT,REFITMSE] = fitLineTLS(X,Y) Fits a line using foot of
% normal notation (rho,theta) to a set of x-y pairs.  X and Y are (Nx1)
% vectors whose indices correspond to x-y pairs.  The total least squares
% error criterion (perpendicular distance) is used to fit the line.  RHOFIT
% and THETAFIT are set S.T. -pi/2 < THETAFIT <= pi/2.

% Solve For Line Parameterization using SVD
A = [x(:),y(:),ones(numel(x),1)];
[~,~,V] = svd(A);
l_fit = V(:,3) / sqrt(sum(V(1:2,3).^2));

% Return Result with Proper Sign
if(l_fit(3) > 0)
    l_fit = -l_fit;
end
rhoFit = -l_fit(3);

% Get Value of Angle
thetaFit = acos(abs(l_fit(1)));

% Put Angle in Proper Quadrant
if(l_fit(1) > 0)
    if(l_fit(2) < 0)
        thetaFit = -thetaFit;
    end
else
    if(l_fit(2) > 0)
        thetaFit = -thetaFit;
        rhoFit = -rhoFit;
    else
        rhoFit = -rhoFit;
    end
end

if(nargout >= 3)
    % Refitting MSE
    fitMSE = mean((A*l_fit).^2);
end