% errorEval_script.m
%% Load Analysis File
[fName,pName] = uigetfile('*.mat','Load Analysis File');

s = load(fullfile(pName,fName));
s = s.saveStruct;
label = s.label;

errorLabel = label;
errorLabel(label == 1) = 1;
errorLabel(label == 2) = 1;
errorLabel(label == 3) = 1;

%% Locate Correct Tracking (Non-Detection)
corrTrackResp = filter([1,1],1,errorLabel);
corrTrack = find(corrTrackResp == 0);
numCorrect = numel(corrTrack);
disp([num2str(numCorrect) ' correct tracking points']);

%% Locate Mistracking Points
errorTrackResp = filter([1,-1],1,errorLabel);
errorTrack = find(errorTrackResp == 1);
%errorTrack = find(label == 1);
numErrorTrack = numel(errorTrack);
disp([num2str(numErrorTrack) ' error track points']);

%% Grab Estimates and Measurements
disp_sig = sqrt(sum(diff(s.mark.subT,[],1).^2,2));
nccDisp_sig = abs(diff(s.algoInfo.nccScore,[],1));
width_sig = diff(s.inst.rho,[],2);
orient_sig = mean(s.inst.theta,2);
deltaWidth_sig = diff(width_sig,1);
deltaOrient_sig = diff(orient_sig,1);
deltaRho_sig = abs(diff(s.inst.rho,[],1));
deltaTheta_sig = abs(diff(s.inst.theta,[],1));
deltaVotesL_sig = diff(s.algoInfo.votesLeft,[],1);
deltaVotesR_sig = diff(s.algoInfo.votesRight,[],1);
deltaInliersL_sig = diff(s.algoInfo.leftNumInliers,[],1);
deltaInliersR_sig = diff(s.algoInfo.rightNumInliers,[],1);
trackPtOffset = s.inst.trackPt(:,1) - s.mark.tCorn(:,1);
deltaTrackPtOffset_sig = abs(diff(trackPtOffset,[],1));

% Correct Template Tracking
disp_corr = disp_sig(corrTrack - 1);
ncc_corr = s.algoInfo.nccScore(corrTrack);
nccDisp_corr = nccDisp_sig(corrTrack - 1);
subFitMSE_corr = s.algoInfo.fitMSE(corrTrack);

% Correct Instrument Pose
delta_rhoL_corr = deltaRho_sig(corrTrack - 1,1);
delta_rhoR_corr = deltaRho_sig(corrTrack - 1,2);
delta_thetaL_corr = deltaTheta_sig(corrTrack - 1,1);
delta_thetaR_corr = deltaTheta_sig(corrTrack - 1,2);
width_corr = width_sig(corrTrack);
orient_corr = orient_sig(corrTrack);
delta_width_corr = deltaWidth_sig(corrTrack - 1);
delta_orient_corr = deltaOrient_sig(corrTrack - 1);
delta_trackPtOffset_corr = deltaTrackPtOffset_sig(corrTrack - 1);

votesL_corr = s.algoInfo.votesLeft(corrTrack);
votesR_corr = s.algoInfo.votesRight(corrTrack);
delta_votesL_corr = deltaVotesL_sig(corrTrack - 1);
delta_votesR_corr = deltaVotesR_sig(corrTrack - 1);

inliersL_corr = s.algoInfo.leftNumInliers(corrTrack);
inliersR_corr = s.algoInfo.rightNumInliers(corrTrack);
delta_inliersL_corr = deltaInliersL_sig(corrTrack - 1);
delta_inliersR_corr = deltaInliersR_sig(corrTrack - 1);

refitL_corr = s.algoInfo.leftRefitMSE(corrTrack);
refitR_corr = s.algoInfo.rightRefitMSE(corrTrack);

% Error Template Tracking
disp_error = disp_sig(errorTrack - 1);
ncc_error = s.algoInfo.nccScore(errorTrack);
nccDisp_error = nccDisp_sig(errorTrack - 1);
subFitMSE_error = s.algoInfo.fitMSE(errorTrack);

% Error Instrument Pose
delta_rhoL_error = deltaRho_sig(errorTrack - 1,1);
delta_rhoR_error = deltaRho_sig(errorTrack - 1,2);
delta_thetaL_error = deltaTheta_sig(errorTrack - 1,1);
delta_thetaR_error = deltaTheta_sig(errorTrack - 1,2);
width_error = width_sig(errorTrack);
orient_error = orient_sig(errorTrack);
delta_width_error = deltaWidth_sig(errorTrack - 1);
delta_orient_error = deltaOrient_sig(errorTrack - 1);
delta_trackPtOffset_error = deltaTrackPtOffset_sig(errorTrack - 1);

votesL_error = s.algoInfo.votesLeft(errorTrack);
votesR_error = s.algoInfo.votesRight(errorTrack);
delta_votesL_error = deltaVotesL_sig(errorTrack - 1);
delta_votesR_error = deltaVotesR_sig(errorTrack - 1);

inliersL_error = s.algoInfo.leftNumInliers(errorTrack);
inliersR_error = s.algoInfo.rightNumInliers(errorTrack);
delta_inliersL_error = deltaInliersL_sig(errorTrack - 1);
delta_inliersR_error = deltaInliersR_sig(errorTrack - 1);

refitL_error = s.algoInfo.leftRefitMSE(errorTrack);
refitR_error = s.algoInfo.rightRefitMSE(errorTrack);
%% Visualize Data Template Tracking
y_corr = ones(numel(corrTrack),1);
y_error = ones(numel(errorTrack),1);

figure();
plot(disp_corr,y_corr,'bx',disp_error,y_error,'ro');
title('Template Displacement');

figure();
subplot(2,1,1);
plot(ncc_corr,y_corr,'bx',ncc_error,y_error,'ro');
title('NCC Score');
subplot(2,1,2);
plot(nccDisp_corr,y_corr,'bx',nccDisp_error,y_error,'ro');
title('NCC Delta Score');

%{
figure();
plot(subFitMSE_corr,y_corr,'bx',subFitMSE_error,y_error,'ro');
xlabel('Sub Fit MSE');
title('Sub Fit MSE');
%}

%% Visualize Data Instrument Pose
y_corr = ones(numel(corrTrack),1);
y_error = ones(numel(errorTrack),1);

% Boundary Parameter Deltas
figure();
subplot(4,1,1);
plot(delta_rhoL_corr,y_corr,'bx',delta_rhoL_error,y_error,'ro');
title('\rho_L Displacement');
subplot(4,1,2);
plot(delta_rhoR_corr,y_corr,'bx',delta_rhoR_error,y_error,'ro');
title('\rho_R Displacement');
subplot(4,1,3);
plot(delta_thetaL_corr,y_corr,'bx',delta_thetaL_error,y_error,'ro');
title('\theta_L Displacement');
subplot(4,1,4);
plot(delta_thetaR_corr,y_corr,'bx',delta_thetaR_error,y_error,'ro');
title('\theta_R Displacement');

% Boundary Line Measurments
figure();
subplot(4,1,1);
plot(votesL_corr,y_corr,'bx',votesL_error,y_error,'ro');
title('Left Line Votes');
subplot(4,1,2);
plot(votesR_corr,y_corr,'bx',votesR_error,y_error,'ro');
title('Right Line Votes');
subplot(4,1,3);
plot(delta_votesL_corr,y_corr,'bx',delta_votesL_error,y_error,'ro');
title('Left Line Votes Displacement');
subplot(4,1,4);
plot(delta_votesR_corr,y_corr,'bx',delta_votesR_error,y_error,'ro');
title('Right Line Votes Displacement');

figure();
subplot(4,1,1);
plot(inliersL_corr,y_corr,'bx',inliersL_error,y_error,'ro');
title('Left Line inliers');
subplot(4,1,2);
plot(inliersR_corr,y_corr,'bx',inliersR_error,y_error,'ro');
title('Right Line inliers');
subplot(4,1,3);
plot(delta_inliersL_corr,y_corr,'bx',delta_inliersL_error,y_error,'ro');
title('Left Line inliers Displacement');
subplot(4,1,4);
plot(delta_inliersR_corr,y_corr,'bx',delta_inliersR_error,y_error,'ro');
title('Right Line inliers Displacement');

% Instrument Width & Orientation
figure();
subplot(5,1,1);
plot(width_corr,y_corr,'bx',width_error,y_error,'ro');
title('Width');
subplot(5,1,2);
plot(orient_corr,y_corr,'bx',orient_error,y_error,'ro');
title('Orient');
subplot(5,1,3);
plot(delta_width_corr,y_corr,'bx',delta_width_error,y_error,'ro');
title('Width Displacement');
subplot(5,1,4);
plot(delta_orient_corr,y_corr,'bx',delta_orient_error,y_error,'ro');
title('Orient Displacement');
subplot(5,1,5);
plot(delta_trackPtOffset_corr,y_corr,'bx',delta_trackPtOffset_error,y_error,'ro');
title('Track Pt. Displacement Delta');
