% detectTrackScript.  Monocular Tracking Multiple Instruments
%% Video Load
% Video Objects
errorCheck = 1;
lVid = VideoReader('exercise_cap_left.avi'); lBack = VideoReader('exercise_back_left.avi');
rVid = VideoReader('exercise_cap_right.avi'); rBack = VideoReader('exercise_back_right.avi');
%% Load Stereo Calibration Parameters
% Camera Calibration
cameraParameters = load('Calib_Results_stereo.mat','fc_left','fc_right',...
                            'cc_left','cc_right','om','T','kc_left',...
                            'kc_right');
[PL,PR,F] = calibMats(cameraParameters);    

%% Setup
% Frames to Analyze
imSize = [lVid.Height,lVid.Width];
frames = [2950,3000]; numFrames = 1 + (frames(2)-frames(1));

%% Detection & Tracking
params = struct('Update','None','Prediction','Constant Velocity',...
                'rhoResolution',1,'thetaResolution',1,'instThresh',450);

disp('Left View Detection & Tracking');
%[m1_l,m2_l,i1_l,...
% i2_l,~,dtfl_ais] = detectTrackFun(lVid,lBack,frames,params,'Error Check',errorCheck);
[m1_l,m2_l,i1_l,i2_l] = detectTrackFun2(lVid,lBack,frames);

disp('Right View Detection & Tracking');
%[m1_r,m2_r,i1_r,...
% i2_r,~,dtfr_ais] = detectTrackFun(rVid,rBack,frames,params,'Error Check',errorCheck);
[m1_r,m2_r,i1_r,i2_r] = detectTrackFun2(rVid,rBack,frames);

%% Generate Point Correspondences
disp('Generating Point Correspondences'); 
tic();
i1_l.epiLines = (F * [i1_l.trackPt';ones(1,size(i1_l.trackPt,1))])';
i2_l.epiLines = (F * [i2_l.trackPt';ones(1,size(i2_l.trackPt,1))])';

% Compute Corresponding Point (Line Intersection)
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_r = getMid(i1_r);
midline2_r = getMid(i2_r);

C1 = cross(i1_l.epiLines,...
           [cos(midline1_r(:,2)),sin(midline1_r(:,2)),-midline1_r(:,1)]);
i1_l.corPt = C1(:,1:2) ./ repmat(C1(:,3),1,2);

C2 = cross(i2_l.epiLines,...
           [cos(midline2_r(:,2)),sin(midline2_r(:,2)),-midline2_r(:,1)]);
i2_l.corPt = C2(:,1:2) ./ repmat(C2(:,3),1,2);
toc();

%{
[XL_1,XR_1] = stereo_triangulation(i1_l.trackPt',i1_l.corPt',cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

[XL_2,XR_2] = stereo_triangulation(i2_l.trackPt',i2_l.corPt',cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

displacementL_1 = sqrt(sum(diff(XL_1,1,2).^2,1));
displacementL_2 = sqrt(sum(diff(XL_2,1,2).^2,1)); toc();

figure(); subplot(1,2,1);
plot(displacementL_1); xlabel('Frame'); ylabel('Displacement (mm)');
title('Instrument_1 3D Position');
subplot(1,2,2);
plot(displacementL_2); xlabel('Frame'); ylabel('Displacement (mm)');
title('Instrument_2 3D Position');
%plot3(XL(1,:),XL(2,:),XL(3,:));
%plot3(XR(1,:),XR(2,:),XR(3,:));
%comet3(XL(1,:),XL(2,:),XL(3,:));
%}
%% Run the Triangulation
disp('Running Triangulation');tic();
% Instrument Midlines
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_r = getMid(i1_r);
midline2_r = getMid(i2_r);

i1_l.startPt = i1_l.trackPt(1:end-1,:);
i1_r.startPt = midEpiIntersect(i1_l.startPt,F,midline1_r(1:end-1,:));

%i1_l.endPt = i1_l.startPt + diff(m1_l.mCorn,1,1);
%i1_l.endPt = i1_l.startPt + m1_l.delT;
i1_l.endPt = i1_l.startPt + diff(m1_l.subT,1,1);
i1_r.endPt = midEpiIntersect(i1_l.endPt,F,midline1_r(2:end,:));

i2_l.startPt = i2_l.trackPt(1:end-1,:);
i2_r.startPt = midEpiIntersect(i2_l.startPt,F,midline2_r(1:end-1,:));

%i2_l.endPt = i2_l.startPt + diff(m2_l.mCorn,1,1);
%i2_l.endPt = i2_l.startPt + m2_l.delT;
i2_l.endPt = i2_l.startPt + diff(m2_l.subT,1,1);
i2_r.endPt = midEpiIntersect(i2_l.endPt,F,midline2_r(2:end,:));

[XL1_start,XR1_start] = stereo_triangulation(i1_l.startPt',i1_r.startPt',cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

[XL1_end,XR1_end] = stereo_triangulation(i1_l.endPt',i1_r.endPt',cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

[XL2_start,XR2_start] = stereo_triangulation(i2_l.startPt',i2_r.startPt',cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

[XL2_end,XR2_end] = stereo_triangulation(i2_l.endPt',i2_r.endPt',cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

displacement_1 = sqrt(sum((XL1_start - XL1_end).^2,1));
displacement_2 = sqrt(sum((XL2_start - XL2_end).^2,1));

toc();

figure(); subplot(1,2,1);
plot(2:numFrames,displacement_1); xlabel('Frame'); ylabel('Displacement (mm)');
title('Instrument_1 3D Displacement');
subplot(1,2,2);
plot(2:numFrames,displacement_2); xlabel('Frame'); ylabel('Displacement (mm)');
title('Instrument_2 3D Displacement');

pos_1 = cumsum(XL1_start - XL1_end,2);
pos_2 = cumsum(XL2_start - XL2_end,2);
numPosPts = size(pos_1,2);

figure('Name','Instrument 1 Position');
subplot(3,1,1);plot(1:numPosPts,pos_1(1,:));
xlabel('Frames');ylabel('Position (mm)');title('X Position');
subplot(3,1,2);plot(1:numPosPts,pos_1(2,:));
xlabel('Frames');ylabel('Position (mm)');title('Y Position');
subplot(3,1,3);plot(1:numPosPts,pos_1(3,:));
xlabel('Frames');ylabel('Position (mm)');title('Z Position');

figure('Name','Instrument 2 Position');
subplot(3,1,1);plot(1:numPosPts,pos_2(1,:));
xlabel('Frames');ylabel('Position (mm)');title('X Position');
subplot(3,1,2);plot(1:numPosPts,pos_2(2,:));
xlabel('Frames');ylabel('Position (mm)');title('Y Position');
subplot(3,1,3);plot(1:numPosPts,pos_2(3,:));
xlabel('Frames');ylabel('Position (mm)');title('Z Position');

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
    lineMask2_l = getLineMask(i2_l,k);
    
    lineMask1_l(round(i1_l.trackPt(k,2)) + (-5:5),round(i1_l.trackPt(k,1)) + (-5:5)) = 1;
    lineMask2_l(round(i2_l.trackPt(k,2)) + (-5:5),round(i2_l.trackPt(k,1)) + (-5:5)) = 1;
    
    lineMask1_r = getLineMask(i1_r,k); 
    lineMask2_r = getLineMask(i2_r,k);
    
    lineMask1_r(round(i1_l.corPt(k,2)) + (-5:5),round(i1_l.corPt(k,1)) + (-5:5)) = 1;
    lineMask2_r(round(i2_l.corPt(k,2)) + (-5:5),round(i2_l.corPt(k,1)) + (-5:5)) = 1;
    
    stereoVid(:,:,:,k) = cat(2,genOverlayIm(temp_l,or(lineMask1_l,lineMask2_l)),...
                               genOverlayIm(temp_r,or(lineMask1_r,lineMask2_r)));
end
toc();
implay(stereoVid);

%% Save a .mat with all info
save('scriptDataTemp','m1_l','m2_l','i1_l','i2_l',...
     'm1_r','m2_r','i1_r','i2_r','frames','numFrames',...
     'lVid','rVid','imSize','PL','PR','F');
