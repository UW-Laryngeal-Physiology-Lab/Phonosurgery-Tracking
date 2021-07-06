%instTrackingEval_script.m
%% Video Load
% Video Objects
errorCheck = 1;
instVid = VideoReader('exercise_1_cap_right.avi');
backVid = VideoReader('testMove_back_right.avi');
backIm = rgb2gray(backVid.read(1));

%% Setup
% Frames to Analyze
imSize = [instVid.Height,instVid.Width];
frames = [5300,5320]; numFrames = 1 + (frames(2)-frames(1));

%% Detection & Tracking
params = struct('Update','None','Prediction','Constant Velocity',...
                'rhoResolution',1,'thetaResolution',1,'instThresh',450);
            
[m1,m2,i1,i2] = detectTrackFun2(instVid,backVid,frames);

%% Generate Result Video
disp('Generating Result Video'); tic();
resultVid = zeros([imSize,3,numFrames],'uint8');
getMidLineMask = @(iStruct,k) drawLineMask(imSize,mean(iStruct.rho(k,:)),mean(iStruct.theta(k,:)));
getLineMask = @(iStruct,k) drawLineMask(imSize,...
    [iStruct.rho(k,:),mean(iStruct.rho(k,:))],...
    [iStruct.theta(k,:),mean(iStruct.theta(k,:))]);

for k = 1:numFrames
    temp = rgb2gray(instVid.read((k-1) + frames(1)));
    
    lineMask1 = getLineMask(i1,k); 
    lineMask2 = getLineMask(i2,k);
    
    lineMask1(round(i1.trackPt(k,2)) + (-5:5),round(i1.trackPt(k,1)) + (-5:5)) = 1;
    lineMask2(round(i2.trackPt(k,2)) + (-5:5),round(i2.trackPt(k,1)) + (-5:5)) = 1;
    
    resultVid(:,:,:,k) = genOverlayIm(temp,or(lineMask1,lineMask2));
end
toc();
implay(resultVid);

%% Instrument "width"
inst = i2;

instWidth = abs(diff(inst.rho,1,2));
figure(); subplot(3,1,1);
plot(1:numel(instWidth),inst.rho(:,1),'gx',...
     1:numel(instWidth),inst.rho(:,2),'gx');
xlabel('Frame #');ylabel('\rho');title('Bound Line \rho values');
 
subplot(3,1,2);
plot(1:numel(instWidth),instWidth); xlabel('Frame');
ylabel('|\rho_1 - \rho_2|'); title('Instrument Width');

errors = find(abs(instWidth(1) - instWidth) > 10);
disp('Error Frames');
disp(errors);

%{
%% Track Point Interframe Displacement
trackInWindow = abs((mark.mCorn(2:end,1) - inst.trackPt(2:end,1)) - (mark.mCorn(1:end-1,1) - inst.trackPt(1:end-1,1)));
figure();plot(trackInWindow);
xlabel('Frame #');
ylabel('x_{trackInWindow}(n) - x_{trackInWindow}(n-1)');
title('Interframe X Displacement of Track Point');
%}
% Theta Values
%figure();
subplot(3,1,3);
plot(1:numFrames,(180/pi)*inst.theta(:,1).','bx',...
     1:numFrames,(180/pi)*inst.theta(:,2).','bx');%,...
     %1:numFrames,(180/pi)*mean(inst.theta,2));
xlabel('Frame #'); ylabel('\theta_{1},\theta_{2}');
title('Boundary Line Theta Values');

%% Reproduction
% Edge Parameter Structure
edgeStruct = struct('backThresh',30,'w_back',1.5,'w_dark',1,...,
    'edgeThresh',20,'backIm',backIm);
edgeStructL = edgeStruct; edgeStructL.edgeName = 'left';
edgeStructR = edgeStruct; edgeStructR.edgeName = 'right';

% Parameter Struct
pStruct = struct();
pStruct.delRho = [35,45];
pStruct.delTheta = [-1.5,1.5];
pStruct.Update = 'None';
pStruct.Prediction = 'Constant Velocity';
pStruct.thetaRes = 0.1;

pStruct_1 = pStruct;
pStruct_2 = pStruct;


k = 5;
frameIm = rgb2gray(instVid.read((k-1) + frames(1)));

% Prediction Mask
mark1Disp = m1.mCorn(k,:) - m1.mCorn(k-1,:);
pMask1 = predictMask(m1.mCorn(k,:),mark1Disp,i1.trackPt(k-1,:),...
                [0,45],i1.theta(k-1,:),imSize);

mark2Disp = m2.mCorn(k,:) - m2.mCorn(k-1,:);
pMask2 = predictMask(m2.mCorn(k,:),mark2Disp,i2.trackPt(k-1,:),...
                [0,45],i2.theta(k-1,:),imSize);

pStruct_1.thetaBounds = (180/pi * mean(i1.theta(k-1,:))) + ...
                        [-5,5];
pStruct_2.thetaBounds = (180/pi * mean(i2.theta(k-1,:))) + ...
                        [-5,5];
                    
% Detect Edges
leftEdge = Itune_edgeDetect(frameIm,edgeStructL); 
leftEdge_1 = leftEdge .* pMask1;
leftEdge_2 = leftEdge .* pMask2;

rightEdge = Itune_edgeDetect(frameIm,edgeStructR); 
rightEdge_1 = rightEdge .* pMask1;
rightEdge_2 = rightEdge .* pMask2;

% Compute Initial Estimates
[r_1,t_1,ais_1] = parLineDetect(leftEdge_1,rightEdge_1,pStruct_1);
[r_2,t_2,ais_2] = parLineDetect(leftEdge_2,rightEdge_2,pStruct_2);

% Sort [left,right]
[r_1,sortIdx] = sort(r_1,'ascend'); t_1 = t_1(sortIdx);
[r_2,sortIdx] = sort(r_2,'ascend'); t_2 = t_2(sortIdx);

% Refit Boundary Lines
[yL_2,xL_2] = find(leftEdge_2); [yR_2,xR_2] = find(rightEdge_2);
[rL_2,tL_2] = refitBoundLine(r_2(1),t_2(1),xL_2,yL_2,3);
[rR_2,tR_2] = refitBoundLine(r_2(2),t_2(2),xR_2,yR_2,3);
rho2 = [rL_2,rR_2]; theta2 = [tL_2,tR_2];

