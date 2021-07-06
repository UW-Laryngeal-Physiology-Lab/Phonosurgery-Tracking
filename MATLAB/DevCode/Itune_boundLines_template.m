function lineMask = Itune_boundLines_template(frameIm,pStruct)

imSize = size(frameIm);
[candIm,backMask] = genCandIm(frameIm,pStruct.backIm,pStruct.backThresh,...
                        pStruct.w_back,pStruct.w_dark);
blobMask = instBlobDetect(candIm,backMask,pStruct.instThresh,...
                            pStruct.diskR);
                        
% Fit To Largest Blob above Template
fitBlobs = blobMask;
fitBlobs(pStruct.tCurr(2):end,:) = 0;
cc = bwconncomp(fitBlobs);
blobSizes = cellfun(@numel,cc.PixelIdxList);
[~,maxIdx] = max(blobSizes);

% Boundary Candidates
candMask = findBoundCand(cc.PixelIdxList{maxIdx},size(fitBlobs),3);

% Find Midline
[rho,theta] = boundLines(...
    lineCands(candIm,candMask,20),0.05);

%midRho = mean(rho); midTheta = mean(theta);
lineMask = drawLineMask(imSize,[rho,midRho],[theta,midTheta]);