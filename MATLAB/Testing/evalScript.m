%% Evaluate a Detection
fNum = 6672;

% Load Old Tracking Data Needed
s = load(trackFile_1);
vidObj = mmreader(s.saveStruct.vidName);
frameIm = rgb2gray(vidObj.read(fNum));
backIm = s.saveStruct.backIm;

% Run the Detector
pStruct = struct('edgeThresh',10);
tCurr = newTrackData_1.mark.tCorn(fNum,:);
[r,t,dbl_ais] = detection_bLines(frameIm,backIm,...
    round(tCurr),newTrackData_1.mark.wStats,newTrackData_1.pc,pStruct);

% Display the Results of the new tracking
lineMask = drawLineMask(size(frameIm),r,t);
oIm = genOverlayIm(frameIm,lineMask);
figure(); imshow(oIm);

%% Evaluate a Tracking
fNum = 10444;

% Load Old Tracking Data Needed
s = load(trackFile_1);
vidObj = mmreader(s.saveStruct.vidName);
frameIm = rgb2gray(vidObj.read(fNum));
backIm = s.saveStruct.backIm;

% Run the Tracker
pStruct = struct('edgeThresh',10);
tPrev = newTrackData_1.mark.tCorn(fNum-1,:);
if(fNum > 2 && s.saveStruct.label(fNum - 2) == 0)
    tPrev2 = newTrackData_1.mark.tCorn(fNum-2,:);
else
    tPrev2 = tPrev;
end

trackPtPrev = newTrackData_1.inst.trackPt(fNum-1,:);
thetaPrev = newTrackData_1.inst.theta(fNum-1,:);
nPrev = newTrackData_1.mark.nCorn(fNum-1,:);
nSize = newTrackData_1.mark.wStats.nx-1;

tCurr = newTrackData_1.mark.tCorn(fNum,:);

                [r,t,dbl_ais] = track_bLines(frameIm,backIm,...
                                        round(tCurr),tPrev,trackPtPrev,thetaPrev,pStruct);
                                    
% Display the Results of the new tracking
lineMask = drawLineMask(size(frameIm),r,t);
oIm = genOverlayIm(frameIm,lineMask);
figure(); imshow(oIm);                                   
                                    