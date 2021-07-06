function [lineMask,fitBlobs,midRho,midTheta] = ...
                                        segFun(frameIm,backFrame,tCurr,thetaVals,rhoResolution,suppress)
% SEGFUN Runs Instrument Boundary Line Segmentation
%
% [LINEMASK,FITBLOBS,MIDRHO,MIDTHETA] = segFun(FRAMEIM,BACKFRAME,TCURR)
                                    
% Generate Candidate Image
[candIm,backMask] = genCandIm(frameIm,backFrame,20,1.5,1);
                                    
% Generate Instrument Blobs
instBlobs = instBlobDetect(candIm,backMask,280,5);

% Fit To Largest Blob above Template
fitBlobs = instBlobs;
fitBlobs(tCurr(2):end,:) = 0;
cc = bwconncomp(fitBlobs);
blobSizes = cellfun(@numel,cc.PixelIdxList);
[~,maxIdx] = max(blobSizes);
fitBlobs = zeros(size(fitBlobs)); fitBlobs(cc.PixelIdxList{maxIdx}) = 1;

% Boundary Candidates
candMask = findBoundCand(cc.PixelIdxList{maxIdx},size(fitBlobs),3);

% Find Midline
%[rho,theta] = boundLines(...
%    and(edge(candIm,'canny'),lineCands(candIm,candMask,20)),0.05);
[rho,theta] = boundLines(...
    lineCands(candIm,candMask,20),thetaVals,rhoResolution,suppress);


midRho = mean(rho); midTheta = mean(theta);
lineMask = drawLineMask(size(frameIm),[rho,midRho],[theta,midTheta]);


