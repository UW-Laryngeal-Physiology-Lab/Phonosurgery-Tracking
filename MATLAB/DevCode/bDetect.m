function [rho,theta] = bDetect(vidObj,frames,pL,pR,m1,m2)

%% Internal Parameters
numFrames = 1 + ([-1,1]*frames(:));

delRho = [30,40];
delTheta = [-3,3];
voteThresh = 10;
suppress = [21,101];

rho = zeros(numFrames,4);
theta = zeros(numFrames,4);

tic();
for k = 1:numFrames
    frameNum = frames(1) + (k-1);
    frameIm = vidObj.read(frameNum); frameIm = frameIm(:,:,1);
    
    maxPt = max([m1(k,2),m2(k,2)]);
    
    % Detect Edges
    leftEdge = Itune_edgeDetect(frameIm,pL); leftEdge(maxPt:end,:) = 0;
    rightEdge = Itune_edgeDetect(frameIm,pR); rightEdge(maxPt:end,:) = 0;
    
    [H,match,T,R] = parAccum(leftEdge,rightEdge,delRho,delTheta,voteThresh);
    peaks = houghpeaks(H,2,'NHoodSize',suppress,'Threshold',0);
    
    peakMatchRhoIdx = sub2ind(size(match),peaks(:,1),peaks(:,2),2*ones(size(peaks,1),1));
    peakMatchThetaIdx = sub2ind(size(match),peaks(:,1),peaks(:,2),ones(size(peaks,1),1));
    
    rho(k,:) = [R(peaks(:,1)),match(peakMatchRhoIdx).'];
    theta(k,:) = (pi/180) * [T(peaks(:,2)),match(peakMatchThetaIdx).'];
        
end
tiempo = toc(); 
fprintf('%4.3f ms per Frame\n',1000*(tiempo/(numFrames)));
    
