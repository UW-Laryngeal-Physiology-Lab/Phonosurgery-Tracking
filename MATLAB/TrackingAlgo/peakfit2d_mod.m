function [P,fitMSE] = peakfit2d_mod(K)
% PEAKFIT2D_MOD Uses quadratic fitting to find subpixel 2D peak
%
% [P,FITMSE] = peakfit2d_mod(K) Fits subpixel peak location to data in (3 x
% 3) matrix K.  The max value of the matrix should be at entry K(2,2).  The
% data in K is fit to a 2D quadratic to find the peak.  P is 2 element
% array [x_offset,y_offset] which is an offset from location (2,2) where
% the fitted peak is located.  FITMSE is the MSE of the quadratic fitting
% to K.

% Code based on : http://www.mathworks.com/matlabcentral/fileexchange/26504-sub-sample-peak-fitting-2d


% approximate polynomial parameter
% F(x,y) => z = a*x^2+b*x*y+c*x+d+e*y^2+f*y
% a = 6*a_true
% b = 4*b_true
% c = 6*c_true
% d = 9*d_true
% e = 6*e_true
% f = 6*f_true
a = (K(2,1)+K(1,1)-2*K(1,2)+K(1,3)-2*K(3,2)-2*K(2,2)+K(2,3)+K(3,1)+K(3,3));
b = (K(3,3)+K(1,1)-K(1,3)-K(3,1));
c = (-K(1,1)+K(1,3)-K(2,1)+K(2,3)-K(3,1)+K(3,3));
d = (2*K(2,1)-K(1,1)+2*K(1,2)-K(1,3)+2*K(3,2)+5*K(2,2)+2*K(2,3)-K(3,1)-K(3,3));
e = (-2*K(2,1)+K(1,1)+K(1,2)+K(1,3)+K(3,2)-2*K(2,2)-2*K(2,3)+K(3,1)+K(3,3));
f = (-K(1,1)-K(1,2)-K(1,3)+K(3,1)+K(3,2)+K(3,3));

% (ys,xs) is subpixel shift of peak location relative to point (2,2)
ys = (6*b*c-8*a*f)/(16*e*a-9*b^2);
xs = (6*b*f-8*e*c)/(16*e*a-9*b^2);
P = [xs,ys];

if(nargout == 2)
    % Compute Fit Error
    x_i = repmat([-1;0;1],3,1);
    y_i = [-1;-1;-1;0;0;0;1;1;1];
    B = [x_i.^2,x_i.*y_i,x_i,ones(numel(x_i),1),y_i.^2,y_i];
    fitMSE = mean((B*([a/6,b/4,c/6,d/9,e/6,f/6].') - ...
        [K(1,:).';K(2,:).';K(3,:).']).^2);
end
    