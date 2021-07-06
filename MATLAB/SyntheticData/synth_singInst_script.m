% synth_singInst_script.m
%% Simulation Parameters
x = [160:200,200*ones(1,50)];
y = [repmat(200,[1,20]),200:230,230*ones(1,40)];
%ang = ((0.2/numel(x)) * (0:numel(x)));
ang = 0.2*ones(size(x));


imSize = [600,400];
frames = [1,numel(x)];
numFrames = 1 + (frames(end) - frames(1));

vidData = zeros([imSize numel(x)],'uint8');

backInfo = struct('mean',220','std',5);
instInfo = struct('mean',40,'std',5,'size',[400,40]);
markerInfo = struct('offset',5,'numStripes',4,'thickness',10,...
    'darkMean',30,'darkStd',5,'lightMean',220,'lightStd',5);
posInfo = struct('x',[],'y',[],'ang',[]);

%% Generate Videos
for k = 1:numel(x)
    posInfo.x = x(k); posInfo.y = y(k); posInfo.ang = ang(k);
    
    [instIm,instMask] = instGen(imSize,instInfo,markerInfo,posInfo);
    backIm = simpBackGen(imSize,backInfo);
    vidData(:,:,k) = synthIm(instIm,instMask,backIm,10,1);
end

clear iVid bVid;

vw = VideoWriter('instVid.avi','Uncompressed AVI');
vw.open();
vw.writeVideo(reshape(vidData,[imSize,1,numel(x)]));
vw.close();

backSample = simpBackGen(imSize,backInfo);
vw = VideoWriter('backVid.avi','Uncompressed AVI');
vw.open();
vw.writeVideo(reshape(backSample,[imSize,1,1]));
vw.close();

%% Detect & Tracking
iVid = VideoReader('instVid.avi');
bVid = VideoReader('backVid.avi');

errorCheck = 1;
params = struct('Update','None','Prediction','Constant Velocity',...
                'rhoResolution',1,'thetaResolution',5);

[m1_l,i1_l,~,dtfl_ais] = ...
    detectTrackFun_1(iVid,bVid,frames,params,'Error Check',errorCheck);

%% View Results
resultsVid = zeros([imSize,3,frames(2)],'uint8');
getLineMask = @(iStruct,k) drawLineMask(imSize,...
    [iStruct.rho(k,:),mean(iStruct.rho(k,:))],...
    [iStruct.theta(k,:),mean(iStruct.theta(k,:))]);

for k = 1:numFrames
   temp = rgb2gray(iVid.read((k-1) + frames(1))); 
   
   lineMask1 = getLineMask(i1_l,k); 
   
   lineMask1(round(i1_l.trackPt(k,2)) + (-5:5),round(i1_l.trackPt(k,1)) +...
       (-5:5)) = 1;
   
   resultsVid(:,:,:,k) = genOverlayIm(temp,lineMask1);
end

implay(resultsVid);

