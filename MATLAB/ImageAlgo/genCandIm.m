function [candIm,backMask] = genCandIm(frameIm,backIm,backThresh,...
    w_back,w_dark)
% GENCANDIM Generates Candidate Image and Background Mask
%
% [CANDIM,BACKMASK] = genCandIm(FRAMEIM,BACKIM,BACKTHRESH,W_BACK,W_DARK)
% Generates "Candidate Image" CANDIM using simulated microsurgery frame
% image FRAMEIM and image w/o instruments BACKIM.

% Background Threshold Mask
backMask = (double(backIm) - double(frameIm)) < backThresh;

% Candidate Image
candIm = w_dark * (255 - double(frameIm)) + ...
         w_back * (double(backIm) - double(frameIm));