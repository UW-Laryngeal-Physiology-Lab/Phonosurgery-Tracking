% staticNoise_dataGen_Script.m
% Detect Track Script for single instrument
%% Video Load
% Video Objects
suffix = '157';
errorCheck = 1;
lVid = mmreader(['leftInst_' suffix '.avi']); 
lBack = mmreader(['leftBgnd_' suffix '.avi']);
rVid = mmreader(['rightInst_' suffix '.avi']); 
rBack = mmreader(['rightBgnd_' suffix '.avi']);

%% Load Stereo Calibration
% Camera Calibration
cameraParameters = load('Calib_Results_stereo.mat','fc_left','fc_right',...
                            'cc_left','cc_right','om','T','kc_left',...
                            'kc_right');
[PL,PR,F] = calibMats(cameraParameters);  

%% Setup
% Frames to Analyze
imSize = [lVid.Height,lVid.Width];
frames = [1,90]; numFrames = 1 + (frames(2)-frames(1));

%% Detection & Tracking
params = struct('Update','None','Prediction','Constant Velocity',...
                'rhoResolution',1,'thetaResolution',0.1);

disp('Left View Detection & Tracking');
[m1_l,i1_l,~,dtfl_ais] = ...
    detectTrackFun_1(lVid,lBack,frames,params,'Error Check',errorCheck);

disp('Right View Detection & Tracking');
[m1_r,i1_r,~,dtfr_ais] = ...
    detectTrackFun_1(rVid,rBack,frames,params,'Error Check',errorCheck);

%% Generate Point Correspondences
disp('Generating Point Correspondences'); 
tic();
i1_l.epiLines = (F * [i1_l.trackPt';ones(1,size(i1_l.trackPt,1))])';

% Compute Corresponding Point (Line Intersection)
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_r = getMid(i1_r);

C1 = cross(i1_l.epiLines,...
           [cos(midline1_r(:,2)),sin(midline1_r(:,2)),-midline1_r(:,1)]);
i1_l.corPt = C1(:,1:2) ./ repmat(C1(:,3),1,2);

toc();

%% Run the Triangulation
disp('Running Triangulation');tic();
% Instrument Midlines
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_r = getMid(i1_r);

i1_l.startPt = i1_l.trackPt(1:end-1,:);
i1_r.startPt = midEpiIntersect(i1_l.startPt,F,midline1_r(1:end-1,:));

%i1_l.endPt = i1_l.startPt + diff(m1_l.mCorn,1,1);
%i1_l.endPt = i1_l.startPt + m1_l.delT;
i1_l.endPt = i1_l.startPt + diff(m1_l.subT,1,1);
i1_r.endPt = midEpiIntersect(i1_l.endPt,F,midline1_r(2:end,:));

[XL1_start,XR1_start] = stereo_triangulation(i1_l.startPt',i1_r.startPt',cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

[XL1_end,XR1_end] = stereo_triangulation(i1_l.endPt',i1_r.endPt',cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

displacement_1 = sqrt(sum((XL1_start - XL1_end).^2,1));

toc();

figure();
plot(2:numFrames,displacement_1); xlabel('Frame'); ylabel('Displacement (mm)');
title('Instrument_1 3D Displacement');

%% Generate Result Video
disp('Generating Result Stereo Video'); tic();
stereoVid = zeros([[1,2].*imSize,3,numFrames],'uint8');
getMidLineMask = @(iStruct,k) drawLineMask(imSize,mean(iStruct.rho(k,:)),mean(iStruct.theta(k,:)));
getLineMask = @(iStruct,k) drawLineMask(imSize,...
    [iStruct.rho(k,:),mean(iStruct.rho(k,:))],...
    [iStruct.theta(k,:),mean(iStruct.theta(k,:))]);

for k = 1:numFrames
    temp_l = rgb2gray(lVid.read((k-1) + frames(1)));
    temp_r = rgb2gray(rVid.read((k-1) + frames(1)));
    
    lineMask1_l = getLineMask(i1_l,k); 

    lineMask1_l(round(i1_l.trackPt(k,2)) + (-5:5),round(i1_l.trackPt(k,1)) + (-5:5)) = 1;
    
    lineMask1_r = getLineMask(i1_r,k); 
    
    lineMask1_r(round(i1_l.corPt(k,2)) + (-5:5),round(i1_l.corPt(k,1)) + (-5:5)) = 1;
    
    stereoVid(:,:,:,k) = cat(2,genOverlayIm(temp_l,lineMask1_l),...
                               genOverlayIm(temp_r,lineMask1_r));
end
toc();
implay(stereoVid);

%% Save a .mat with all info
save('scriptDataTemp','m1_l','i1_l','dtfl_ais',...
     'm1_r','i1_r','dtfr_ais','frames','numFrames',...
     'lVid','rVid','imSize','PL','PR','F');
