function viewInstTrackData(saveStruct)
% VIEWINSTTRACKDATA View Data generated from monocular inst tracking

%% Determine Temporal Bounds
startIdx = find(~isnan(saveStruct.label),1,'first');
endIdx = find(~isnan(saveStruct.label),1,'last');
t = startIdx:endIdx;
incorrect = (saveStruct.label(t,:) ~= 0);
%incorrect = false(numel(t),1);

%% Template
x = saveStruct.mark.subT(t,1); x(incorrect) = nan;
y = saveStruct.mark.subT(t,2); y(incorrect) = nan;
nccScore = saveStruct.algoInfo.nccScore(t); nccScore(incorrect) = nan;
fitMSE = saveStruct.algoInfo.fitMSE(t); fitMSE(incorrect) = nan;

figure('Name','Template Position');
subplot(2,1,1); plot(t,x); title('X Position');
subplot(2,1,2); plot(t,y); title('Y Position');

figure('Name','Template Tracking Algo Info');
subplot(2,1,1);plot(t,nccScore); title('NCC Score');
subplot(2,1,2);plot(t,fitMSE); title('Fit MSE');

%% Boundary Line
tp_x = saveStruct.inst.trackPt(t,1); tp_x(incorrect) = nan;
tp_y = saveStruct.inst.trackPt(t,2); tp_y(incorrect) = nan;
r_L = saveStruct.inst.rho(t,1); r_L(incorrect) = nan;
r_R = saveStruct.inst.rho(t,2); r_R(incorrect) = nan;
t_L = saveStruct.inst.theta(t,1); t_L(incorrect) = nan;
t_R = saveStruct.inst.theta(t,2); t_R(incorrect) = nan;

width = r_L - r_R; angDiff = t_L - t_R;

figure('Name','Instrument Track Point');
subplot(2,1,1); plot(t,tp_x); title('Track Pt. X');
subplot(2,1,2); plot(t,tp_y); title('Track Pt. Y');

figure('Name','Boundary Line Parameterization');
subplot(3,1,1); plot(t,r_L); title('\rho_L');
subplot(3,1,2); plot(t,r_R); title('\rho_R');
subplot(3,1,3); hold all; plot(t,(180/pi)*t_L); plot(t,(180/pi)*t_R); hold off; title('\theta');

figure('Name','Boundary Line Characteristics');
subplot(2,1,1); plot(t,width); title('Width');
subplot(2,1,2); plot(t,(180/pi)*angDiff); title('Angular Difference');

votesLeft = saveStruct.algoInfo.votesLeft(t); votesLeft(incorrect) = nan;
votesRight = saveStruct.algoInfo.votesRight(t); votesRight(incorrect) = nan;
leftRefit = saveStruct.algoInfo.leftRefitMSE(t); leftRefit(incorrect) = nan;
rightRefit = saveStruct.algoInfo.rightRefitMSE(t); rightRefit(incorrect) = nan;
leftInliers = saveStruct.algoInfo.leftNumInliers(t); leftInliers(incorrect) = nan;
rightInliers = saveStruct.algoInfo.rightNumInliers(t); rightInliers(incorrect) = nan;

figure('Name','Boundary Line Algo Info');
subplot(3,1,1);
plot(t,votesLeft,t,votesRight); title('Line Votes');
subplot(3,1,2);
plot(t,leftRefit,t,rightRefit); title('Refitting MSE');
subplot(3,1,3);
plot(t,leftInliers,t,rightInliers); title('Number of Inliers');