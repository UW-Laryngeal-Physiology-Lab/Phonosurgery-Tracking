function gIm = gaussian2dIm(X0,Y0,sigmaX,sigmaY,nRows,nCols)
% GAUSSIAN2DIM Generates Image of 2D Gaussian
% 
% GIM = gaussian2dIm(X0,Y0,SIGMAX,SIGMAY,NROWS,NCOLS) Generates uint8 image
% of 2D gaussian centered at (X0,Y0) with standard deviations SIGMAX and
% SIGMAY.  The image is a matrix of dimension (NROWS x NCOLS).

% Amplitude
A = 150;

% Generate Image
[X,Y] = meshgrid(1:nCols,1:nRows);
 
gIm = uint8(A * exp(-(...
        ((X - X0).^2)/(2*sigmaX.^2) + ...
        ((Y - Y0).^2)/(2*sigmaY.^2) ...
        )));