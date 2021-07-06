% frameDropTraj_script.m
%% File Info
calFile = 'Calib_Results_stereo.mat';
trackingDirectory = 'capture_1';
leftFile = 'F_G1_L.mat';
rightFile = 'F_G1_R.mat';

%% Frame Drop Info
disp('----frameDropTraj_script.m------');
missingFrames_L = [];
missingFrames_R = [];

%% Get Data Files
disp('STATUS : Getting Data Files');

% Calibration File
[calName,calPath] = uigetfile('*.mat','Load Camera Calibration File');
if(calName == 0)
    disp('ERROR : User Cancelled Calibration File Selection');
    disp('STATUS : Exit');
    return;
end

camCalFile = fullfile(calPath,calName);

% Tracking Files
[tFileName_L,tPath_L] = uigetfile('*.mat','Load Left Tracking File',...
                                  calPath);
if(tFileName_L == 0)
    disp('ERROR : User Cancelled Left Tracking File Selection');
    disp('STATUS : Exit');
    return;
end
                                                            
[tFileName_R,tPath_R] = uigetfile('*.mat','Load Right Tracking File',...
                                  tPath_L);
if(tFileName_R == 0)
    disp('ERROR : User Cancelled Right Tracking File Selection');
    disp('STATUS : Exit');
    return;
end

leftTrackFile = fullfile(tPath_L,tFileName_L);
rightTrackFile = fullfile(tPath_R,tFileName_R);

%% Load Data
disp('STATUS : Loading Data');

% Read Camera Parameters
cameraParameters = load(camCalFile);
[PL,PR,F] = calibMats(cameraParameters);

% Read Tracking Parameters
s_L = load(leftTrackFile); s_L = s_L.saveStruct;
s_R = load(rightTrackFile); s_R = s_R.saveStruct;

% Grab label vectors
label_L = s_L.label; label_R = s_R.label;

% Grab Inst Structs
numFrames_L = s_L.numFrames; numFrames_R = s_R.numFrames;
i1_L = s_L.inst; i1_R = s_R.inst;

%% Expand structures to account for dropped frames
disp('STATUS : GENERATING 3D Data');
[i1_L,label_L] = frameDropExpand(i1_L,label_L,missingFrames_L);
[i1_R,label_R] = frameDropExpand(i1_R,label_R,missingFrames_R);

% Verify Both Views Have Same Number of Frame
if(numel(label_L) ~= numel(label_R))
    disp('ERROR : Number of frames do not match.  Check that you entered correct missing frames.');
    return;
end

numFrames = numel(label_L);

%% Init Data Structures
% Determine frames with correct tracking parameters
gFrames = (label_L == 0) & (label_R == 0);
frames = find(gFrames);

% Initialize Structures %
% 3D Position & Orientation 
pos3D = nan*ones(numFrames,3);
orient3D = nan*ones(numFrames,3);

% 2D Position & Orientation
pos2D_L = nan*ones(numFrames,2);
orient2D_L = nan*ones(numFrames,2);
pos2D_R = nan*ones(numFrames,2);
orient2D_R = nan*ones(numFrames,2);

%% Generate Point Correspondence
getMid = @(iStruct) [mean(iStruct.rho,2),mean(iStruct.theta,2)];
midline1_L = getMid(i1_L);
midline1_R = getMid(i1_R);

% Left View Start & End
pos2D_L(frames,:) = i1_L.trackPt(frames,:);

% Compute Orient and Correspondence Points WRT Left Track Point
[op_L,tp_R,op_R] = orientPtCorr(pos2D_L(frames,:),midline1_L(frames,:),...
                                                midline1_R(frames,:),F,10);
% Left View Orient
orient2D_L(frames,:) = op_L;

% Right View Track Point Correspondence & Orient
pos2D_R(frames,:) = tp_R;
orient2D_R(frames,:) = op_R;
%% Run Midpoint Triangulation
triangFun = @(left,right) stereo_triangulation(left',right',...
    cameraParameters.om,...
    cameraParameters.T,cameraParameters.fc_left,cameraParameters.cc_left,...
    cameraParameters.kc_left,0,cameraParameters.fc_right,cameraParameters.cc_right,...
    cameraParameters.kc_right,0);

[pos3D_sub,dummy] = triangFun(pos2D_L(frames,:),...
                                pos2D_R(frames,:));
[orient3D_sub,dummy] = triangFun(orient2D_L(frames,:),...
                                orient2D_R(frames,:));

% Get 3D Points
pos3D(frames,:) = pos3D_sub';

% Get 3D Orientation Vectors
getOrient = @(or_sub,pt_sub) (or_sub' - pt_sub') ./ ...
    repmat(sqrt(sum((or_sub' - pt_sub').^2,2)),1,3);
orient3D(frames,:) = getOrient(orient3D_sub,pos3D_sub);
%% Interpolation
disp('STATUS : Interpolating 3D Data');

% Build Joint label array
% Either unlabeled (nan) --> unlabeled (nan)
% Both Correct (0) --> Correct (0)
% Correct & Error --> Error
% Blur & non-blur Error --> non-blur error (because we interp blur)
label = nan*ones(numel(label_L),1);

for k = 1:numel(label)
    if(isnan(label_L(k)) || isnan(label_R(k)))
        % Either unlabeled
        label(k) = nan;
    elseif(label_L(k) == 0 && label_R(k) == 0)
        % Both labelled correct
        label(k) = 0;
    elseif(label_L(k) >= 0 || label_R(k) >= 0)
        % Both labeled.  Either labelled Error.  Assign the higher value so
        % blur and non-blur are assigned non-blur.
        label(k) = max([label_L(k),label_R(k)]);
    end
end

% Find Error Labels Proceeding Correct Labels (atartErrors)
errorLabel = ~isnan(label) & label > 0;
startError = 1 + find(diff(errorLabel) == 1);

% Identify interpFrames
% 1) Single error label sandwiched between correct labels
% 2) Multiple blur labels sandwiched between correct labels
interpFrames = findInterpFrames(label);

% Run Interpolation
frames_i = find(gFrames | interpFrames);
X = nan*ones(numel(gFrames),1); Y = X; Z = X; % Position
X_or = X; Y_or = X; Z_or = X; % Orientation

X(frames_i) = interp1(frames,pos3D(frames,1),frames_i);
Y(frames_i) = interp1(frames,pos3D(frames,2),frames_i);
Z(frames_i) = interp1(frames,pos3D(frames,3),frames_i);

X_or(frames_i) = interp1(frames,orient3D(frames,1),frames_i);
Y_or(frames_i) = interp1(frames,orient3D(frames,2),frames_i);
Z_or(frames_i) = interp1(frames,orient3D(frames,3),frames_i);

%% Plot Signals
disp('STATUS : Plotting 3D Data');

% Determine Frame Subset
f1 = find(~isnan(X),1,'first'); f2 = find(~isnan(X),1,'last');
fSub = f1:f2;

figure();
subplot(3,1,1); plot(fSub,X(fSub));
xlabel('Frame'); title('X');
subplot(3,1,2); plot(fSub,Y(fSub));
xlabel('Frame'); title('Y');
subplot(3,1,3); plot(fSub,Z(fSub));
xlabel('Frame'); title('Z');

%% Compute Surgical Metrics
disp('STATUS : Computing Metrics');

% Path Length
X_vel = diff(X); Y_vel = diff(Y); Z_vel = diff(Z);
pathLenIdx = ~isnan(X_vel) & ~isnan(Y_vel) & ~isnan(Z_vel);
pathLen = sum(sqrt(X_vel(pathLenIdx).^2 + Y_vel(pathLenIdx).^2 + ...
              Z_vel(pathLenIdx).^2));

% Depth Perception
dp_sig = sum([X_vel,Y_vel,Z_vel] .* [X_or(1:end-1),Y_or(1:end-1),...
                                     Z_or(1:end-1)],2);
depthPerception = sum(abs(dp_sig(~isnan(dp_sig))));
          
% Smoothness %
% X Smoothness
jerk_X = diff(X_vel,2,1).^2;
smoothness_X = 1/(numel(~isnan(jerk_X))) * ...
                            sqrt(sum(jerk_X(~isnan(jerk_X))));
% Y Smoothness
jerk_Y = diff(Y_vel,2,1).^2;
smoothness_Y = 1/(numel(~isnan(jerk_Y))) * ...
                            sqrt(sum(jerk_Y(~isnan(jerk_Y))));
                        
% Z Smoothness
jerk_Z = diff(Z_vel,2,1).^2;
smoothness_Z = 1/(numel(~isnan(jerk_Z))) * ...
                            sqrt(sum(jerk_Z(~isnan(jerk_Z))));
                        
% Net Smoothness
netJerk = sqrt(sum([jerk_X,jerk_Y,jerk_Z],2));
netSmoothness = 1/(numel(~isnan(netJerk))) * ...
                            sqrt(sum(netJerk(~isnan(netJerk))));
                        
% Display Metrics
disp(['Path Length (mm) : ' num2str(pathLen) ' mm'])
disp(['Depth Perception (mm) : ' num2str(depthPerception)]);
disp(['X Smoothness : ' num2str(smoothness_X)]);
disp(['Y Smoothness : ' num2str(smoothness_Y)]);
disp(['Z Smoothness : ' num2str(smoothness_Z)]);
disp(['Net Smoothness : ' num2str(netSmoothness)]);

%% Save XLS File