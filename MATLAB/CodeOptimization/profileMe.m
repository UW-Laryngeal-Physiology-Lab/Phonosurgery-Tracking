for k = frames(1):frames(end)
    frameIm = rgb2gray(vidObj.read(k));
    fail_bLine = 0;
    
    switch(label(k-1))
        case 0
            % Perform Tracking
            % Template
            tPrev = mark.tCorn(k-1,:);
            if(label(k - 2) == 0)
                tPrev2 = mark.tCorn(k-2,:);
            else
                tPrev2 = tPrev;
            end
            trackPtPrev = inst.trackPt(k-1,:);
            thetaPrev = inst.theta(k-1,:);
            nPrev = mark.nCorn(k-1,:);
            nSize = mark.wStats.nx-1;

            [tCurr,nCurr,maxScore,fitMSE] = track_temp(frameIm,pc,tPrev,...
                                                tPrev2,nPrev,nSize);

            % Boundary Line Estimates
            try
                [r,t,dbl_ais] = track_bLines(frameIm,backIm,...
                                        round(tCurr),tPrev,trackPtPrev,thetaPrev);
            catch ME
                fail_bLine = 1;
            end
    otherwise
        % Perform Detection
        % Template
        [tCurr,nCurr,maxScore,fitMSE] = detection_temp(frameIm,...
            pc,mark.wStats);
        
        % Boundary Lines
        % Determine Last Correct Labeling
        % Use Last Correct to Seed Detection
        lastCorr = find(label(1:k-1) == 0,1,'last');
        tCorn_lastCorr = mark.tCorn(lastCorr,:);
        trackPt_lastCorr = inst.trackPt(lastCorr,:);
        thetaPrev = mean(inst.theta(lastCorr,:));
        trackPtDelta = trackPt_lastCorr - tCorn_lastCorr;
        try
            [r,t,dbl_ais] = track_bLines(frameIm,backIm,...
                round(tCurr),tCurr,tCurr + trackPtDelta,thetaPrev);
            %[r,t,dbl_ais] = detection_bLines(frameIm,backIm,...
            %    round(tCurr),mark.wStats,pc);
        catch ME
            fail_bLine = 1;
        end
    end

% Store Results
mark.tCorn(k,:) = round(tCurr);
mark.subT(k,:) = tCurr;
mark.nCorn(k,:) = round(nCurr);

% Template Related Algo Info
algoInfo.nccScore(k,1) = maxScore;
algoInfo.fitMSE(k,1) = fitMSE;

if(fail_bLine == 0)
% Boundary Line Estimates
inst.rho(k,:) = r;
inst.theta(k,:) = t;
halfSize = round((mark.wStats.wSize-1)/2);
inst.trackPt(k,:) = computeTrackPt(r,t,tCurr,halfSize);

% Boundary Line Algo Info
algoInfo.votesLeft(k,:) = dbl_ais.votesLeft;
algoInfo.votesRight(k,:) = dbl_ais.votesRight;
algoInfo.leftNumInliers(k,:) = dbl_ais.leftNumInliers;
algoInfo.leftRefitMSE(k,:) = dbl_ais.leftRefitMSE;
algoInfo.rightNumInliers(k,:) = dbl_ais.rightNumInliers;
algoInfo.rightRefitMSE(k,:) = dbl_ais.rightRefitMSE;
end
end