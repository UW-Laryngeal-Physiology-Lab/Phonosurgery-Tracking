function [X_RMSE,Y_RMSE] = testHarness1(delX,delY,nRows,nCols,sigma)
% TESTHARNESS1 Test harness for subPixTempTrack testing
%
% [X_RMSE,Y_RMSE] = testHarness1(DELX,DELY,NROWS,NCOLS,SIGMA) Generates
% video of shifted 2D Gaussian with standard deviation SIGMA.  DELX and
% DELY determines the interframe shift.  The video has dimensions NROWS by
% NCOLS.  subPixTempTrack is ran on the video.  The RMSE error between the
% actual x and y displacement and that estimated by subPixTempTrack is
% returned in X_RMSE and Y_RMSE.

if nargin == 4
    sigma = 0;
end

nFrames = numel(delX);

g = shifted2dGaussSeq(delX,delY,nRows,nCols);
for k = 1:size(g,3)
    g(:,:,k) = g(:,:,k) + uint8(sigma*randn(size(g,1),size(g,2)));
end

vWriter = VideoWriter('shiftGauss.avi','Uncompressed AVI');
vWriter.open();
vWriter.writeVideo(reshape(g,[nRows,nCols,1,nFrames+1]));
vWriter.close();

vidObj = VideoReader('shiftGauss.avi');
[tCorn,nCorn,subT,wStats] = subPixTempTrack(vidObj,[1,nFrames+1],20,80,'Update','None');
clear vidObj;

del_est = diff(subT,1,1);
X_RMSE = sqrt(mean((delX - del_est(:,1)).^2));
Y_RMSE = sqrt(mean((delY - del_est(:,2)).^2));

fprintf('XRMSE : %3.6f\nYRMSE : %3.6f\n',X_RMSE,Y_RMSE);