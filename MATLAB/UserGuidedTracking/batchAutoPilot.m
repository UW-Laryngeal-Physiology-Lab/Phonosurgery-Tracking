function batchAutoPilot(batchInfoFile,resultFileName,threshVec)
% BATCHAUTOPILOT Runs UGA Auto-Pilot for batch of analysis files
%
% batchAutoPilot(BATCHINFOFILE,RESULTFILENAME) Runs auto-pilot on a batch
% of UGA analysis files.  BATCHINFOFILE is .txt file of the form:
% 
% fileName_1 stopFrame
% fileName_2 stopFrame
% fileName_3 stopFrame
%
% Where fileName_x corresponds to a UGA_1 generate analysis file in which
% the template has been initialized.  This program runs the auto-pilot
% routine of UGA_1 until the associated stopFrame.  If the algorithm hits a
% frame at which it is underconfident the auto-pilot stops prematurely at
% that frame.  RESULTFILENAME is .txt file that store the result of the
% batch run.  It is of the form:
%
% fileName_1 startFrame desiredStopFrame actualStopFrame
 

%% Input Arguments
% Check if user passes confidence threshold vector
if(nargin < 3)
    threshVec = [16.83,0.09,3.61,0.018,-44.5,-43.5,1.78,20,20,0.4];
end

%% Read batchInfoFile
fileList = {}; stopFrameList = [];
fid_info = fopen(batchInfoFile,'r');

eof = 0;
while(eof == 0)
   tLine = fgetl(fid_info); 
   if(tLine == -1)
       eof = 1;
   else
       % Parse Line into fileList & Number
       C = textscan(tLine,'%s %d');
        
       fileList = [fileList;C{1}];
       stopFrameList = [stopFrameList;C{2}];
   end
end
fclose(fid_info);

%% Init
disp('Start Batch Auto-Pilot');
numFiles = numel(fileList);

% Check threshVec
typeVec = [0,0,0,0,1,1,0,1,1,1];
if(numel(threshVec) ~= numel(typeVec))
    disp('Error: Incorrect Number of Elements in threshVec');
    disp('End Batch Auto-Pilot');
    return;
end

saveFields = {'mark','pc','neighborhoodSize','templateFrame','inst',...
'algoInfo','vidName','backName','imSize','numFrames',...
'backIm','state','currentLabel','label','currFrame','analysisFileName'};

% Create .txt file
fid = fopen(resultFileName,'w+');
fprintf(fid,'Batch Auto Pilot\n');

%% Batch Loop
for m = 1:numFiles
    % Load Analysis File
    disp(['Running: ' fileList{m}]);
    
    load(fileList{m},'saveStruct');
    stopFrame = stopFrameList(m);
    
    s = struct(); % Struct used for holding variables
    s.threshVec = threshVec;
    s.typeVec = typeVec;
    
    % Load all Save Fields
    for n = 1:numel(saveFields)
        s.(saveFields{n}) = saveStruct.(saveFields{n});
    end
    
    % Setup Video Objects
    s.vidObj = mmreader(s.vidName);
    s.backObj = mmreader(s.backName);
    
    % Write File Name to .txt
    fprintf(fid,'%s\t',fileList{m});
    
    % Check for stopFrame special condition
    if(stopFrame == 0)
        stopFrame = s.vidObj.NumberOfFrames - 1;
    end
    
    % Auto-Pilot
    startFrame = s.currentLabel;
    s.currFrame = startFrame;
    
    if(stopFrame > s.currentLabel && stopFrame <= s.numFrames)
        % Write Start Frame and Desired Stop Frame to .txt
        fprintf(fid,'%d\t%d\t',s.currentLabel,stopFrame);
        
        % Label the First Frame Correct
        s = labelCurr(s,0);
        
        for k = s.currFrame:stopFrame
            % Form Feature Vector
            f = [k-1,k];
            fVec = [sqrt(sum(diff(s.mark.subT(f,:),1,1).^2,2)),...
            abs(diff(s.algoInfo.nccScore(f),1,1)),...
            abs(diff(diff(s.inst.rho(f,:),1,2),1,1)),...
            abs(diff(mean(s.inst.theta(f,:),2),1,1)),...
            diff(s.algoInfo.leftNumInliers(f),1,1),...
            diff(s.algoInfo.rightNumInliers(f),1,1),...
            abs(diff(s.inst.trackPt(f,1) - s.mark.tCorn(f,1),1,1)),...
            s.algoInfo.leftNumInliers(k),...
            s.algoInfo.rightNumInliers(k),...
            s.algoInfo.nccScore(k)];
        
            % Check Algorithm Confidence
            result = threshClassifyFun(fVec,s.threshVec,s.typeVec);

            if(any(result))
                % Algorithm has low confidence stop
                disp('Low Confidence Stop Early');
                break;
            end

            s = labelCurr(s,0);
        end
    else
        disp('Error: Incorrect Stop Frame');
    end
    
    % Update Analysis File
    saveAnalysis(fileList{m},s);
    
    % Write Stop Frame to .txt
    fprintf(fid,'%d\n',s.currentLabel-1);
end
disp('End Batch Auto Pilot');
fclose(fid);

function saveAnalysis(fileName,dStruct)
% SAVEANALYSIS Saves .mat analysis file corresponding to UGA
%
% saveAnalysis(FILENAME,DSTRUCT) Saves the UGA data in DATASTRUCT to .mat
% file FILENAME.

saveFields = {'mark','pc','neighborhoodSize','templateFrame','inst',...
'algoInfo','vidName','backName','imSize','numFrames',...
'backIm','state','currentLabel','label','currFrame','analysisFileName'};

saveStruct = dStruct;
allFields = fieldnames(saveStruct);

% Remove all non save fields
for k = 1:numel(allFields)
    cName = allFields{k};
    if(~any(strcmp(cName,saveFields)))
        saveStruct = rmfield(saveStruct,cName);
    end
end

disp(['Saving: ' fileName]);
save(fileName, 'saveStruct');

function dStruct = labelCurr(dStruct,labelNum)
% LABELCURR Labels the current frame and runs tracker on next frame
%
% DSTRUCT = labelCurr(DSTRUCT,LABELNUM) Labels the current frame of DSTRUCT
% as LABELNUM.  Then runs tracker on proceeding frame and updates current
% frame and label in DSTRUCT.

% Label The Unlabeled Frames
dStruct.label(dStruct.currFrame) = labelNum;

% Run Tracker
dStruct = performAnalysis(dStruct,dStruct.currFrame + 1);

% Update Current Frame and Label
dStruct.currFrame = dStruct.currFrame + 1;
dStruct.currentLabel = dStruct.currFrame;


function dStruct = performAnalysis(dStruct,frameNum)
% PERFORMANALYSIS Runs Instrument Tracking/Detection on Single Frame
%
% DSTRUCT = performAnalysis(DSTRUCT,FRAMENUM) Instrument Tracking/Detection
% algorithm.  The previous frame label determines whether the detection or
% tracking is performed.

% Grab Frame Data
frameIm = rgb2gray(dStruct.vidObj.read(frameNum));

fail_bLine = 0; % Flag used to identify algo error
switch(dStruct.label(frameNum-1))
    case 0
        % Perform Tracking
        % Template
        tPrev = dStruct.mark.tCorn(frameNum-1,:);
        if(frameNum > 2 && dStruct.label(frameNum - 2) == 0)
            tPrev2 = dStruct.mark.tCorn(frameNum-2,:);
        else
            tPrev2 = tPrev;
        end
        trackPtPrev = dStruct.inst.trackPt(frameNum-1,:);
        thetaPrev = dStruct.inst.theta(frameNum-1,:);
        nPrev = dStruct.mark.nCorn(frameNum-1,:);
        nSize = dStruct.mark.wStats.nx-1;
        
        [tCurr,nCurr,maxScore,fitMSE] = track_temp(frameIm,dStruct.pc,tPrev,...
                                            tPrev2,nPrev,nSize);
                                        
        % Boundary Line Estimates
        try
            [r,t,dbl_ais] = track_bLines(frameIm,dStruct.backIm,...
                                    round(tCurr),tPrev,trackPtPrev,thetaPrev);
        catch ME
            fail_bLine = 1;
        end
    otherwise
        % Perform Detection
        % Template
        [tCurr,nCurr,maxScore,fitMSE] = detection_temp(frameIm,...
            dStruct.pc,dStruct.mark.wStats);
        
        % Boundary Lines
        % Determine Last Correct Labeling
        % Use Last Correct to Seed Detection
        lastCorr = find(dStruct.label(1:frameNum-1) == 0,1,'last');
        tCorn_lastCorr = dStruct.mark.tCorn(lastCorr,:);
        trackPt_lastCorr = dStruct.inst.trackPt(lastCorr,:);
        thetaPrev = mean(dStruct.inst.theta(lastCorr,:));
        trackPtDelta = trackPt_lastCorr - tCorn_lastCorr;
        try
            [r,t,dbl_ais] = track_bLines(frameIm,dStruct.backIm,...
                round(tCurr),tCurr,tCurr + trackPtDelta,thetaPrev);
            %[r,t,dbl_ais] = detection_bLines(frameIm,dStruct.backIm,...
            %    round(tCurr),dStruct.mark.wStats,dStruct.pc);
        catch ME
            fail_bLine = 1;
        end
end

% Store Results
dStruct.mark.tCorn(frameNum,:) = round(tCurr);
dStruct.mark.subT(frameNum,:) = tCurr;
dStruct.mark.nCorn(frameNum,:) = round(nCurr);

% Template Related Algo Info
dStruct.algoInfo.nccScore(frameNum,1) = maxScore;
dStruct.algoInfo.fitMSE(frameNum,1) = fitMSE;

if(fail_bLine == 0)
% Boundary Line Estimates
dStruct.inst.rho(frameNum,:) = r;
dStruct.inst.theta(frameNum,:) = t;
halfSize = round((dStruct.mark.wStats.wSize-1)/2);
dStruct.inst.trackPt(frameNum,:) = computeTrackPt(r,t,tCurr,halfSize);

% Boundary Line Algo Info
dStruct.algoInfo.votesLeft(frameNum,:) = dbl_ais.votesLeft;
dStruct.algoInfo.votesRight(frameNum,:) = dbl_ais.votesRight;
dStruct.algoInfo.leftNumInliers(frameNum,:) = dbl_ais.leftNumInliers;
dStruct.algoInfo.leftRefitMSE(frameNum,:) = dbl_ais.leftRefitMSE;
dStruct.algoInfo.rightNumInliers(frameNum,:) = dbl_ais.rightNumInliers;
dStruct.algoInfo.rightRefitMSE(frameNum,:) = dbl_ais.rightRefitMSE;
end
