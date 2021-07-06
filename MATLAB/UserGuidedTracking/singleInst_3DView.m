%% Load Stereo Calibration Parameters
% Camera Calibration
cameraParameters = load('Calib_Results_stereo.mat','fc_left','fc_right',...
                            'cc_left','cc_right','om','T','kc_left',...
                            'kc_right');
[PL,PR,F] = calibMats(cameraParameters);    

%% Analysis File Names
[fName_l,pName_l] = uigetfile('*.mat','Load Left View Analysis File');
[fName_r,pName_r] = uigetfile('*.mat','Load Right View Analysis File');

%% Load Analysis Files
s_l = load(fullfile(pName_l,fName_l)); s_l = s_l.saveStruct;
s_r = load(fullfile(pName_r,fName_r)); s_r = s_r.saveStruct;

% Grab Inst Structs
numFrames = s_l.numFrames;
m1_l = s_l.mark; m1_r = s_r.mark;
i1_l = s_l.inst; i1_r = s_r.inst;
gFrames = (s_l.label(1:end-1) == 0) & ...
          (s_r.label(1:end-1) == 0) & ...
          (s_l.label(2:end) == 0) & ...
          (s_r.label(2:end) == 0);
startFrames = find(gFrames);
endFrames = startFrames + 1;

% Setup Structures
i1_l.epiLines = nan*ones(numFrames,3);
i1_l.corPt = nan*ones(numFrames,2);
i1_l.startPt = nan*ones(numFrames-1,2);
i1_r.startPt = nan*ones(numFrames-1,2);
i1_l.endPt = nan*ones(numFrames-1,2);
i1_r.endPt = nan*ones(numFrames-1,2);

%% Generate Point Correspondences
disp('Generating Point Correspondences'); 
tic();
i1_l.epiLines(gFrames,:) = (F * [i1_l.trackPt(gFrames,:)';...
                      ones(1,size(i1_l.trackPt(gFrames,:),1))])';

% Compute Corresponding Point (Line Intersection)
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_r = getMid(i1_r);

% Compute Intersection of Midline & Epipolar Line
C1 = cross(i1_l.epiLines(gFrames,:),...
           [cos(midline1_r(gFrames,2)),sin(midline1_r(gFrames,2)),...
            -midline1_r(gFrames,1)]);
i1_l.corPt(gFrames,:) = C1(:,1:2) ./ repmat(C1(:,3),1,2);

%% Run the Triangulation
disp('Running Triangulation');tic();
% Instrument Midlines
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_r = getMid(i1_r);

% Start Points
i1_l.startPt = i1_l.trackPt;
i1_r.startPt(startFrames,:) = midEpiIntersect(...
    i1_l.startPt(startFrames,:),F,midline1_r(startFrames,:));

% End Points
i1_l.endPt(endFrames,:) = i1_l.startPt(startFrames,:) + (m1_l.subT(endFrames,:) - m1_l.subT(startFrames,:)); 
i1_r.endPt(endFrames,:) = midEpiIntersect(i1_l.endPt(endFrames,:),F,midline1_r(endFrames,:));

% Perform Triangulation
[XL1_start,XR1_start] = stereo_triangulation(...
    i1_l.startPt(startFrames,:)',i1_r.startPt(startFrames,:)',...
    cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

% Perform Triangulation
[XL1_end,XR1_end] = stereo_triangulation(...
    i1_l.endPt(endFrames,:)',i1_r.endPt(endFrames,:)',...
    cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

displacement_1 = nan*ones(numFrames,1);
displacement_1(startFrames,:) = sqrt(sum((XL1_start - XL1_end).^2,1));

toc();

figure();
plot(1:numFrames,displacement_1); xlabel('Frame'); ylabel('Displacement (mm)');
title('Instrument_1 3D Displacement');

pos_1 = nan*ones(3,numFrames);
pos_1(:,startFrames) = cumsum(XL1_end - XL1_start,2);

figure('Name','Instrument 1 Position (displacement)');
subplot(3,1,1);plot(1:numFrames,pos_1(1,:));
xlabel('Frames');ylabel('Position (mm)');title('X Position');
subplot(3,1,2);plot(1:numFrames,pos_1(2,:));
xlabel('Frames');ylabel('Position (mm)');title('Y Position');
subplot(3,1,3);plot(1:numFrames,pos_1(3,:));
xlabel('Frames');ylabel('Position (mm)');title('Z Position');

figure('Name','Instrument 1 Position (start point)');
subplot(3,1,1);plot(XL1_start(1,:));
xlabel('Frames');ylabel('Position (mm)');title('X Position');
subplot(3,1,2);plot(XL1_start(2,:));
xlabel('Frames');ylabel('Position (mm)');title('Y Position');
subplot(3,1,3);plot(XL1_start(3,:));
xlabel('Frames');ylabel('Position (mm)');title('Z Position');

figure('Name','Instrument 1 Position (end point)');
subplot(3,1,1);plot(XL1_end(1,:));
xlabel('Frames');ylabel('Position (mm)');title('X Position');
subplot(3,1,2);plot(XL1_end(2,:));
xlabel('Frames');ylabel('Position (mm)');title('Y Position');
subplot(3,1,3);plot(XL1_end(3,:));
xlabel('Frames');ylabel('Position (mm)');title('Z Position');

%% Track Point Based
[XL1_TP,XR1_TP] = stereo_triangulation(...
    i1_l.trackPt(gFrames,:)',i1_r.trackPt(gFrames,:)',...
    cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

figure('Name','Instrument 1 Position (end point)');
subplot(3,1,1);plot(XL1_TP(1,:));
xlabel('Frames');ylabel('Position (mm)');title('X Position');
subplot(3,1,2);plot(XL1_TP(2,:));
xlabel('Frames');ylabel('Position (mm)');title('Y Position');
subplot(3,1,3);plot(XL1_TP(3,:));
xlabel('Frames');ylabel('Position (mm)');title('Z Position');