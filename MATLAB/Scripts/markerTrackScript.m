% markerTrackScript.m Monocular Marker Tracking
%% Video Load
% Video Objects
vid = VideoReader('exercise_1_cap_right.avi');

%% Setup
% Frames to Analyze
imSize = [vid.Height,vid.Width];
frames = [2250,2500]; numFrames = 1 + (frames(2)-frames(1));

wSize = 59;
nSize = 3;
update = 'None';
prediction = 'Constant Velocity';

%% Get Window 
f1 = vid.read(frames(1)); f1 = f1(:,:,1);
[centerPt,tSize] = squareDraw(f1,wSize);
windowStruct = struct('centerPt',centerPt,'tSize',tSize);
%% Template Tracking
[tCorn,nCorn,subT,wStats,ais] = subPixTempTrack(vid,frames,wSize,...
    nSize,'Update','None','Prediction','Constant Velocity',...
    'Template Window',windowStruct);

% View Results
tVid = trackVid(vid,frames,tCorn,nCorn,wStats);
implay(tVid);

figure();
subplot(2,1,1);
plot(ais.maxNCCScore);
xlabel('Frame'); ylabel('Max NCC Score');
subplot(2,1,2);
plot(ais.fitMSE);
xlabel('Frame'); ylabel('Sub Pix Fit MSE');