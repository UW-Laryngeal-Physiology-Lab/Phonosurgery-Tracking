function g = shifted2dGaussSeq(delX,delY,nRows,nCols,params)
% SHIFTED2DGAUSSSEQ Generates Greyscale images sequence of shifted gaussian
%
% G = shifted2dGaussSeq(DELX,DEL,NROWS,nCols,PARAMS)
% Inputs:
% DELX - Delta X Sequence in pixels (Nx1)
% DELY - Delta Y Sequence in pixels (Nx1)
% NROWS - Number of image rows
% nCols - Number of image columns
% PARAMS - Optional Parameters Structure.  See below for description of
% parameters
% OUTPUTS:
% g - (NROWS X nCols X N) uint8 sequence, containing a 2D gaussing that
% moves around according to DELX and DELY.  
% DESCRIPTION
% Generates an image sequence with a 2D gaussian.  At each frame the
% gaussian is shifted in the x and y direction based on DELX and DELY.

%% Input Parameters
if nargin == 4
    params = struct();
end

if(isfield(params,'X0'))
    X0 = params.X0;
else
    X0 = round(nCols/2);
end

if(isfield(params,'Y0'))
    Y0 = params.Y0;
else
    Y0 = round(nRows/2);
end

A = 150;
sigmax = 15/2;
sigmay = 15/2;
%% Generate Shifted Gaussian
g = zeros([nRows,nCols,numel(delX)],'uint8');
[X,Y] = meshgrid(1:nCols,1:nRows);

shiftX = [0 ; cumsum(delX(:))];
shiftY = [0 ; cumsum(delY(:))];

for k = 1:numel(shiftX)
    g(:,:,k) = uint8(A * exp(-(((X - shiftX(k) - X0).^2)/(2*sigmax.^2) + ...
                         ((Y - shiftY(k) - Y0).^2)/(2*sigmay.^2))));
end
    