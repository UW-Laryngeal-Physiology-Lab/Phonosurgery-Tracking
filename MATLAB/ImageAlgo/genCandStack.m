function [candStack,backMaskStack] = genCandStack(imStack,backIm,...
    backThresh,w_back,w_dark)
% GENCANDIM Generates Candidate Image and Background Mask Stack
%
% [CANDSTACK,BACKMASKSTACK] = genCandIm(FRAMEIM,BACKIM,BACKTHRESH,W_BACK,
% W_DARK) Generates "Candidate Image" CANDIM using simulated microsurgery
% frame image FRAMEIM and image w/o instruments BACKIM.

% Background Threshold Mask
backMaskStack = (double(repmat(backIm,[1,1,1,size(imStack,4)])) - ...
    double(imStack)) < backThresh;

% Candidate Image
candStack = w_dark * (255 - double(imStack)) + ...
         w_back * (double(repmat(backIm,[1,1,1,size(imStack,4)])) - ...
         double(imStack));