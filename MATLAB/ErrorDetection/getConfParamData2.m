function [corr,err] = ...
                        getConfParamData2(trackFile,startFrame,endFrame)
% GETCONFPARAMDATA2 Gets parameter vals used by algo. confidence classifer
%
% [CORR,ERR] = getConfParamData(TRACKFILE) Extracts confidence parameter
% values for correctly tracked and tracking error frames in TRACKFILE.
% TRACKFILE is the file saved by the user-guided tracking application UGT.
% Two structures, CORR and ERR are returned.  These structures contain the
% parameter values for correctly tracked frames in CORR and tracking error
% frames in ERR.
%
% [CORR,ERR] = getConfParamData(TRACKFILE,STARTFRAME) Starts extraction of
% parameter values from STARTFRAME.
%
% [CORR,ERR] = getConfParamData(TRACKFILE,STARTFRAME,ENDFRAME) Extracts
% confidence parameter values from STARTFRAME to ENDFRAME.  
%
% At each frame, the phonomicrosurgery instrument tracking algorithm
% extracts a set of features related to the instrument's visible boundaries
% and a marker placed on the instrument's cylindrical rod.  After finding
% these features, a set of confidence parameters are calculated.  These
% parameters determine whether or not the algorithm is confident about the
% features it has extracted.  If it is not confident, it asks for
% user-intervention.  The confidence parameters are:
% 1) delta_NCC : The difference between the marker window NCC coefficient
% at the current frame and the previous frame {smoothness parameter}
% 2) delta_width : The difference between the instrument width at the
% current frame and previous frame {smoothness paramter}.  Width is
% calculated using the rho parameter of the left and right boundary line
% estimate.  Width = (rho_{left} - rho_{right}).
% 3) delta_orientation : The difference between the midline orientation at
% the current frame and previous frame {smoothness parameter}.  Orientation
% is calculated using the theta parameter of the left and right boundary
% line estimate.  Orientation = (theta_{left} + theta_{right})/2.
% 4) inliersL and inliersR : The number of inlier points used to fit the
% left and right boundary line estimate respectively.
% 5) delta_inliersL and delta_inliersR : The diffrence betwen the number of
% inliers in the current frame and the previous frame {smoothness
% parameter}.

%% Load the Tracking Data
s = load(trackFile);

%% Get Label and Tracked Parmeter subsets
% Get Label Data
s = s.saveStruct;
label = s.label;

% Determine Start and End Frame
if(nargin == 1)
    % Use First Labelled Frame if not provided by user
    startFrame = find(~isnan(label),1); % First Labelled Frame
end
    
if(nargin < 3)
    % End Frame is first unlabelled frame after start frame
    % if not defined by the user
    nanFrames = find(isnan(label));
    endFrame = nanFrames(find(nanFrames > startFrame,1));
end

frames = startFrame:endFrame;

% Grab Signals from startFrame to endFrame. 
label = label(frames);
nccScore = s.algoInfo.nccScore(frames,:);
rho = s.inst.rho(frames,:);
theta = s.inst.theta(frames,:);
leftNumInliers = s.algoInfo.leftNumInliers(frames,:);
rightNumInliers = s.algoInfo.leftNumInliers(frames,:);
markY = s.mark.tCorn(frames,2);

% Use all Error Types
errorLabel = label;
errorLabel(label > 1) = 1;

%% Locate Correct Tracking Frames (Non-Detection)
% All frames that have been correctly tracked need to be extracted from the
% label data.  This is done by noting that the frame prior to correctly
% tracked frame is labelled zero.

% Compute the sum between adjacent labels.  A correctly tracked frame is 
% indicated by values that sum to zero.
corrTrackResp = filter([1,1],1,errorLabel); 
corrTrackResp(1) = 1; % Exclude First Frame

% Find all adjacent labels that summed to zero.  The second point in that
% adjacent set the correctly tracked frame.  
corrTrack = find(corrTrackResp == 0);
numCorrect = numel(corrTrack);
disp([num2str(numCorrect) ' correct tracking points']);

%% Locate Tracking Error Frames
% All frames that are tracking errors need to be extracted from the label
% data.  In general, error frames are labelled one. The frame previous to a
% tracking error frame is labelled correct (zero).  

% Compute difference between adjacent labels.  Tracking error is indicated
% by adjacent labels whose difference is 1.  
errorTrackResp = filter([1,-1],1,errorLabel); errorTrackResp(1) = 0;

% Find all adjacent labels that have a difference of one.  The second label
% in that adjacent set is a tracking error.
errorTrack = find(errorTrackResp == 1);
numErrorTrack = numel(errorTrack);
disp([num2str(numErrorTrack) ' blur/occlusion points']);

%% Find Correct Tracking Frame Parameter Values
% Store parameter values for correctly tracked frames in corr struct.  The
% parameters are : 
% delta_NCC, delta_width, delta_orient, inliersL, inliersR, 
% delta_inliersL, delta_inliersR
%
% delta_[parameterName] are smoothness parameters.  They are calculated
% by taking the difference between the current and previous frame value.

% NCC Coefficient
diffNCC_sig = abs(diff(nccScore,[],1));
delta_NCC_corr = diffNCC_sig(corrTrack - 1);

% Instrument Width
width_sig = diff(rho,[],2);
diffWidth_sig = abs(diff(width_sig,1));
delta_width_corr = diffWidth_sig(corrTrack - 1);

% Instrument Orientation
orient_sig = mean(theta,2);
diffOrient_sig = abs(diff(orient_sig,1));
delta_orient_corr = diffOrient_sig(corrTrack - 1);

% Boundary Line Fitting Inliers
inliersL_corr = leftNumInliers(corrTrack);
inliersR_corr = rightNumInliers(corrTrack);


%diffInliersL_sig = diff(leftNumInliers,[],1);
%diffInliersR_sig = diff(rightNumInliers,[],1);
%delta_inliersL_corr = diffInliersL_sig(corrTrack - 1);
%delta_inliersR_corr = diffInliersR_sig(corrTrack - 1);
diffInliersDistL_sig = abs(diff(markY - leftNumInliers,[],1));
diffInliersDistR_sig = abs(diff(markY - rightNumInliers,[],1));
delta_inliersDistL_corr = diffInliersDistL_sig(corrTrack - 1);
delta_inliersDistR_corr = diffInliersDistR_sig(corrTrack - 1);

% Correct Structure Stores all the parameters for correct tracking frames
corr = struct('delta_NCC_corr',delta_NCC_corr,...
                 'delta_width_corr',delta_width_corr,...
                 'delta_orient_corr',delta_orient_corr,...
                 'inliersL_corr',inliersL_corr,...
                 'inliersR_corr',inliersR_corr,...
                 'delta_inliersDistL_corr',delta_inliersDistL_corr,...
                 'delta_inliersDistR_corr',delta_inliersDistR_corr,...
                 'frames_corr',frames(corrTrack));

%% Find Tracking Error Frame Parameter Values
% Store parameter values for tracking error frames in err struct.  See
% above section for description of parameters.

% NCC Coefficient
diffNCC_sig = abs(diff(nccScore,[],1));
delta_NCC_error = diffNCC_sig(errorTrack - 1);

% Instrument Width
width_sig = diff(rho,[],2);
diffWidth_sig = abs(diff(width_sig,1));
delta_width_error = diffWidth_sig(errorTrack - 1);

% Instrument Orientation
orient_sig = mean(theta,2);
diffOrient_sig = abs(diff(orient_sig,1));
delta_orient_error = diffOrient_sig(errorTrack - 1);

% Boundary Line Fitting Inliers
inliersL_error = leftNumInliers(errorTrack);
inliersR_error = rightNumInliers(errorTrack);
%diffInliersL_sig = diff(leftNumInliers,[],1);
%diffInliersR_sig = diff(rightNumInliers,[],1);
%delta_inliersL_error = diffInliersL_sig(errorTrack - 1);
%delta_inliersR_error = diffInliersR_sig(errorTrack - 1);
diffInliersDistL_sig = abs(diff(markY - leftNumInliers,[],1));
diffInliersDistR_sig = abs(diff(markY - rightNumInliers,[],1));
delta_inliersDistL_error = diffInliersDistL_sig(errorTrack - 1);
delta_inliersDistR_error = diffInliersDistR_sig(errorTrack - 1);

% err Structure Stores all the parameters for tracking error frames
err = struct('delta_NCC_error',delta_NCC_error,...
                 'delta_width_error',delta_width_error,...
                 'delta_orient_error',delta_orient_error,...
                 'inliersL_error',inliersL_error,...
                 'inliersR_error',inliersR_error,...
                 'delta_inliersDistL_error',delta_inliersDistL_error,...
                 'delta_inliersDistR_error',delta_inliersDistR_error,...
                 'frames_error',frames(errorTrack));

