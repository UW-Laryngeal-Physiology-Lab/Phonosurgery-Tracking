% tremorAnalysis_3D.m
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
i1_l.orientPt = nan*ones(numFrames,2);
i1_l.startPt = nan*ones(numFrames-1,2);
i1_r.orientPt = nan*ones(numFrames,2);
i1_r.startPt = nan*ones(numFrames-1,2);
i1_l.endPt = nan*ones(numFrames-1,2);
i1_r.endPt = nan*ones(numFrames-1,2);

%% Generate Point Correspondences
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_l = getMid(i1_l);
midline1_r = getMid(i1_r);

% Get TrackPt and Orientation Pt Correspondence
[op_l,tp_r,op_r] = orientPtCorr(i1_l.trackPt(gFrames,:),...
    midline1_l(gFrames,:),midline1_r(gFrames,:),F,10);

i1_l.corPt(gFrames,:) = tp_r;
i1_l.orientPt(gFrames,:) = op_l;

%% Run the Triangulation
disp('Running Triangulation')

% Start Points
i1_l.startPt = i1_l.trackPt;
i1_r.startPt = i1_l.corPt;

% End Points
i1_l.endPt(endFrames,:) = i1_l.startPt(startFrames,:) + (m1_l.subT(endFrames,:) - m1_l.subT(startFrames,:)); 
i1_r.endPt(endFrames,:) = midEpiIntersect(i1_l.endPt(endFrames,:),F,midline1_r(endFrames,:));

% Orient Points
i1_r.orientPt(gFrames,:) = op_r;

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

%{
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
%}

%% Linear Interpolation for Missing Points
allStartFrames = startFrames(1):startFrames(end);
allEndFrames = endFrames(1):endFrames(end);

X_start = zeros([3,numel(allStartFrames)]);
X_start(1,:) = interp1(startFrames,XL1_start(1,:),allStartFrames);
X_start(2,:) = interp1(startFrames,XL1_start(2,:),allStartFrames);
X_start(3,:) = interp1(startFrames,XL1_start(3,:),allStartFrames);

X_end = zeros([3,numel(allEndFrames)]);
X_end(1,:) = interp1(endFrames,XL1_end(1,:),allEndFrames);
X_end(2,:) = interp1(endFrames,XL1_end(2,:),allEndFrames);
X_end(3,:) = interp1(endFrames,XL1_end(3,:),allEndFrames);

% View Interpolated Start & End Points
figure('Name','Instrument 1 Position (start point)');
subplot(3,1,1);plot(X_start(1,:));
xlabel('Frames');ylabel('Position (mm)');title('X Position');
subplot(3,1,2);plot(X_start(2,:));
xlabel('Frames');ylabel('Position (mm)');title('Y Position');
subplot(3,1,3);plot(X_start(3,:));
xlabel('Frames');ylabel('Position (mm)');title('Z Position');

%{
figure('Name','Instrument 1 Position (end point)');
subplot(3,1,1);plot(X_end(1,:));
xlabel('Frames');ylabel('Position (mm)');title('X Position');
subplot(3,1,2);plot(X_end(2,:));
xlabel('Frames');ylabel('Position (mm)');title('Y Position');
subplot(3,1,3);plot(X_end(3,:));
xlabel('Frames');ylabel('Position (mm)');title('Z Position');
%}

%% FFT Analysis
f = (90/numel(allStartFrames)) * (0:(numel(allStartFrames)-1));
f_sIdx = find(f >= 0,1,'first');
f_eIdx = find(f >= 15,1,'first');
f_idx = f_sIdx:f_eIdx;

w1 = hann(numel(X_start(1,:)));

fft_start_x = abs(fft(w1' .* (X_start(1,:) - mean(X_start(1,:)))));
fft_start_y = abs(fft(w1' .* (X_start(2,:) - mean(X_start(2,:)))));
fft_start_z = abs(fft(w1' .* (X_start(3,:) - mean(X_start(3,:)))));

figure('Name','Instrument Position FFT (dB)');
subplot(3,1,1);plot(f(f_idx),fft_start_x(f_idx));
xlabel('Frame');title('X Position FFT Mag');
subplot(3,1,2);plot(f(f_idx),fft_start_y(f_idx));
xlabel('Frame');title('Y Position FFT Mag');
subplot(3,1,3);plot(f(f_idx),fft_start_z(f_idx));
xlabel('Frame');title('Z Position FFT Mag');

%% STFT Analysis
wLen = 3;
NWindow = 90 * wLen; wdw = hann(NWindow);
Nfft = round(90/0.1);
f = (90/Nfft)*(0:(Nfft-1));
idx1 = find(f >= 1/wLen,1,'first'); 
idx2 = find(f >= 15,1,'first');
iSub = idx1:idx2;

[S_x,F_x,T_x,P_x] = spectrogram(X_start(1,:) - mean(X_start(1,:)),wdw,round(NWindow/2),Nfft,90);
[S_y,F_y,T_y,P_y] = spectrogram(X_start(2,:) - mean(X_start(2,:)),wdw,round(NWindow/2),Nfft,90);
[S_z,F_z,T_z,P_z] = spectrogram(X_start(3,:) - mean(X_start(3,:)),wdw,round(NWindow/2),Nfft,90);

figure();
surf(1 + 90*T_x,F_x(iSub),sqrt(P_x(iSub,:)),'EdgeAlpha',0); view(0,90);
title('X STFT');
figure();
surf(1 + 90*T_y,F_y(iSub),sqrt(P_y(iSub,:)),'EdgeAlpha',0);view(0,90);
title('YSTFT');
figure();
surf(1 + 90*T_z,F_z(iSub),sqrt(P_z(iSub,:)),'EdgeAlpha',0);view(0,90);
title('ZSTFT');

% PSD Estimate
PSD_x = mean(P_x,2);
PSD_y = mean(P_y,2);
PSD_z = mean(P_z,2);

figure();
subplot(3,1,1);
plot(F_x(iSub),PSD_x(iSub));
xlabel('Frequency (Hz)');title('X Direction Periodogram');
subplot(3,1,2);
plot(F_y(iSub),PSD_y(iSub));
xlabel('Frequency (Hz)');title('Y Direction Periodogram');
subplot(3,1,3);
plot(F_z(iSub),PSD_z(iSub));
xlabel('Frequency (Hz)');title('Z Direction Periodogram');