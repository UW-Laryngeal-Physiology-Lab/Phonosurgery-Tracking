function newTrackData = trackWithOldData(tFile)
% TESTTRACKINGAGAINSTDATA Tests latest tracking algo against old data
%
% NEWTRACKDATA = testTrackingAgainstData() Prompts user for tracking .mat
% file.  The latest tracking data is run on the video data corresponding to
% the tracking file.  Tracking/Detection is ran based on frame labels in
% tracking file.  The resulting tracking data is returned in NEWTRACKDATA.
%
% NEWTRACKDATA = testTrackingAgainstData(TRACKFILE) Uses the tracking data
% in TRACKFILE.  

%% Setup the Tracking using the old data
% Display Tracking File Name
fprintf('Tracking File : %s\n',tFile);

% Load the old data
s = load(tFile); tData = s.saveStruct;
oldLabel = tData.label;

% Load Video & Store Video Parameters
vidName = tData.vidName;
backName = tData.backName;
vidObj = mmreader(vidName);
backObj = mmreader(backName);
imSize = [vidObj.Height,vidObj.Width];
numFrames = vidObj.NumberOfFrames;

% Store Backframe
backIm = rgb2gray(backObj.read(1));

% Init Structs
%label = nan*ones(numFrames,1);
mark = struct('tCorn',nan*ones(numFrames,2),'wStats',[],...
    'nCorn',nan*ones(numFrames,2),'subT',nan*ones(numFrames,2));
inst = struct('rho',nan*ones(numFrames,2),'theta',...
    nan*ones(numFrames,2),'trackPt',nan*ones(numFrames,2));
algoInfo = struct('nccScore',nan*ones(numFrames,1),'fitMSE',...,
    nan*ones(numFrames,1),'votesLeft',nan*ones(numFrames,1),'votesRight',...
    nan*ones(numFrames,1),'leftNumInliers',nan*ones(numFrames,1),...
    'leftRefitMSE',nan*ones(numFrames,1),'rightNumInliers',...
    nan*ones(numFrames,1),'rightRefitMSE',nan*ones(numFrames,1));

% Template Related
neighborhoodSize = 2;
pc = tData.pc;
templateFrame = tData.templateFrame;

% Get Starting and Endframe
startFrame = tData.templateFrame;
endFrame = find(~isnan(oldLabel),1,'last');
fprintf('Tracking from Frame %u to %u\n',startFrame,endFrame);

%% Run Tracking Algorithm
pStruct = struct('edgeThresh',20);

frameIm = rgb2gray(vidObj.read(startFrame));

% Marker Window
mark.tCorn(startFrame,:) = tData.mark.tCorn(startFrame,:);
mark.subT(startFrame,:) = tData.mark.subT(startFrame,:);
mark.wStats = tData.mark.wStats;
mark.nCorn(startFrame,:) = tData.mark.nCorn(startFrame,:);
tempWindow = mark.tCorn(startFrame,:);
tSize = mark.wStats.wSize;

% Marker Window Algo Info
algoInfo.nccScore(startFrame,1) = 1;
algoInfo.fitMSE(startFrame,1) = 0;

% Boundary Lines
[r,t,dbl_ais] = detection_bLines(frameIm,backIm,tempWindow,mark.wStats,pc,pStruct);
inst.rho(startFrame,:) = r;
inst.theta(startFrame,:) = t;
inst.trackPt(startFrame,:) = computeTrackPt(r,t,tempWindow,...
                                                    round((tSize-1)/2));

% Boundary Line Algo Info                                               
algoInfo.votesLeft(startFrame,:) = dbl_ais.votesLeft;
algoInfo.votesRight(startFrame,:) = dbl_ais.votesRight;
algoInfo.leftNumInliers(startFrame,:) = dbl_ais.leftNumInliers;
algoInfo.leftRefitMSE(startFrame,:) = dbl_ais.leftRefitMSE;
algoInfo.rightNumInliers(startFrame,:) = dbl_ais.rightNumInliers;
algoInfo.rightRefitMSE(startFrame,:) = dbl_ais.rightRefitMSE;
    
                                                
for frameNum = 1+startFrame:endFrame
    % Status Indicator
    if(mod(frameNum,100) == 0)
        fprintf('Status : %u / %u\n',frameNum,endFrame);
    end
    
    % Grab the Current Frame
    frameIm = rgb2gray(vidObj.read(frameNum));
    
    fail_bLine = 0; % Flag used to identify algo error
    switch(oldLabel(frameNum-1))
        case 0
            % Perform Tracking
            % Template
            tPrev = mark.tCorn(frameNum-1,:);
            if(frameNum > 2 && oldLabel(frameNum - 2) == 0)
                tPrev2 = mark.tCorn(frameNum-2,:);
            else
                tPrev2 = tPrev;
            end
            trackPtPrev = inst.trackPt(frameNum-1,:);
            thetaPrev = inst.theta(frameNum-1,:);
            nPrev = mark.nCorn(frameNum-1,:);
            nSize = mark.wStats.nx-1;

            [tCurr,nCurr,maxScore,fitMSE] = track_temp(frameIm,pc,tPrev,...
                                                tPrev2,nPrev,nSize);

            % Boundary Line Estimates
            try
                [r,t,dbl_ais] = track_bLines(frameIm,backIm,...
                                        round(tCurr),tPrev,trackPtPrev,thetaPrev,pStruct);
            catch ME
                fail_bLine = 1;
            end   
        otherwise
            % Perform Detection
            % Template
            [tCurr,nCurr,maxScore,fitMSE] = detection_temp(frameIm,...
                pc,mark.wStats);

            % Boundary Lines
            % Use Marker Orientation
            %lastCorr = find(oldLabel(1:frameNum-1) == 0,1,'last');
            %tCorn_lastCorr = mark.tCorn(lastCorr,:);
            %trackPt_lastCorr = inst.trackPt(lastCorr,:);
            %thetaPrev = mean(inst.theta(lastCorr,:));
            %trackPtDelta = trackPt_lastCorr - tCorn_lastCorr;
            try
                %[r,t,dbl_ais] = track_bLines(frameIm,backIm,...
                %    round(tCurr),tCurr,tCurr + trackPtDelta,thetaPrev);
                [r,t,dbl_ais] = detection_bLines(frameIm,backIm,...
                    round(tCurr),mark.wStats,pc,pStruct);
            catch ME
                fail_bLine = 1;
            end   
    end
    
    % Store Results
    mark.tCorn(frameNum,:) = round(tCurr);
    mark.subT(frameNum,:) = tCurr;
    mark.nCorn(frameNum,:) = round(nCurr);

    % Template Related Algo Info
    algoInfo.nccScore(frameNum,1) = maxScore;
    algoInfo.fitMSE(frameNum,1) = fitMSE;

    if(fail_bLine == 0)
        % Boundary Line Estimates
        inst.rho(frameNum,:) = r;
        inst.theta(frameNum,:) = t;
        halfSize = round((mark.wStats.wSize-1)/2);
        inst.trackPt(frameNum,:) = computeTrackPt(r,t,tCurr,halfSize);

        % Boundary Line Algo Info
        algoInfo.votesLeft(frameNum,:) = dbl_ais.votesLeft;
        algoInfo.votesRight(frameNum,:) = dbl_ais.votesRight;
        algoInfo.leftNumInliers(frameNum,:) = dbl_ais.leftNumInliers;
        algoInfo.leftRefitMSE(frameNum,:) = dbl_ais.leftRefitMSE;
        algoInfo.rightNumInliers(frameNum,:) = dbl_ais.rightNumInliers;
        algoInfo.rightRefitMSE(frameNum,:) = dbl_ais.rightRefitMSE;
    end 
end

%% Return the Results
newTrackData = struct('mark',mark,'pc',pc,...
    'neighborhoodSize',neighborhoodSize,'templateFrame',templateFrame,...
    'inst',inst,'algoInfo',algoInfo,'vidName',vidName,'backName',...
    backName,'imSize',imSize,'numFrames',numFrames,'backIm',backIm,...
    'oldLabel',oldLabel);

