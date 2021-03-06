% orientPtCorr_ver.m
% Verification fo orientPtCorr.m  Verification is based on two metrics
% point correspondence and positioing on midlines.  For both the track
% point and orientation point correspondence metrics (@ each frame
% resulting in an error signal) are computed that evalute that the two
% points correspond (in a stereo configuration geometric sense using
% Fundamental matrix) and that within each view they fall on their
% respective midline.  A value of zero indicates no error.  A dataset of
% scissor data for a subject performing the six square simulated
% phonomicrosurgery task is used along with its stereo calibration data.
% KS 2011_09_16
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
mid1 = [mean(i1_l.rho,2),mean(i1_l.theta,2)];
mid2 = [mean(i1_r.rho,2),mean(i1_r.theta,2)];

f = (s_l.label == 0) & (s_r.label == 0);

%% Run orientPtCorr
[op_l,tp_r,op_r] = orientPtCorr(i1_l.trackPt(f,:),mid1(f,:),mid2(f,:),F,10);

%% View 1
[c_v1,ls1_v1,ls2_v1] = orientPtCorr_ver_computeMetrics(i1_l.trackPt(f,:),...
    mid1(f,:),tp_r,mid2(f,:),F);
[c_v2,ls1_v2,ls2_v2] = orientPtCorr_ver_computeMetrics(op_l,...
    mid1(f,:),op_r,mid2(f,:),F);

figure('Name','Track Point Error Signals');
subplot(3,1,1); plot(c_v1); title('Correspondence');
subplot(3,1,2); plot(ls1_v1); title('View 1 Line');
subplot(3,1,3); plot(ls2_v1); title('View 2 Line');

figure('Name','Orientation Point Error Signals');
subplot(3,1,1); plot(c_v2); title('Correspondence');
subplot(3,1,2); plot(ls1_v2); title('View 1 Line');
subplot(3,1,3); plot(ls2_v2); title('View 2 Line');

