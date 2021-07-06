% TrajectoryEst_script.m
% Script used to generate instrument 3D trajectory data and motion metrics
% from UGT tracking data files and camera calibration file.  The trajectory
% consists of the 3D instrument position over the frame sequence.  This
% position is triangulated for all correctly tracked frames (in both
% views). Additionally, it is interpolated for error frames sandwiched
% between two correct frames and consecutive blur frames sandwiched between
% correct frames.  The axes X,Y,Z refer to the coordinate frame of the left
% camera.  Calculated motion metrics are : Path Length, Depth Percpetion, X
% Smoothness, Y Smoothness, Z Smoothness, Net Motion Smoothness, and Net
% Orientation Smoothness.
%
% The following 6 variables need to be set in the 'Script Setup' portion
% below these comments
%
% 1) camCalFile :
% 
% 2) leftTrackFile :
% 
% 3) rightTrackFile : J:\Phonosurgery\Simulated_Surgery_Software
% 4) missingFrames_L : 
% 5) missingFrames_R : 
% 6) plotFlag : [1]
%
% The tracking data files should correspond to the same instrument over the
% same video frame subset.  
%% Script Setup
disp('------------------------------------------------------------------');
disp('TrajectoryEst_script.m');
% Calibration File
% The location of the cameraCalibrationToolbox generated stereo calibration
% file.
camCalFile = 'J:\Phonosurgery\Simulated_Surgery_Software\Calibration\2015_6_18_3\Calib_Results_stereo.mat';

% UGT Tracking data files
% Tracking data files generated from the left and right video for the same
% surgical instrument over the same video frame subset.  These files
% contain feature data (estimated midline and track point) that is used to
% estimate 3D instrument position at each frame.
leftTrackFile = 'J:\Phonosurgery\Learning Curve Automatic Tracking\C003\C003_S05T2L\C003_S05T2L_left_scissors';
rightTrackFile ='J:\Phonosurgery\Learning Curve Automatic Tracking\C003\C003_S05T2L\C003_S05T2L_right_scissors';

% Frame Drop/Failed Grab Arrays
% Frames that were dropped or unsuccessfully grabbed during video
% acquisition.  As of 12/15/2011 this data is saved in a .txt file along
% with the video files.  Enter the dropped frames numbers into the array.
% If no frames were dropped set the variable equal to an empty array []
missingFrames_L = [];
missingFrames_R = [2835];

% Plot Flag
% Controls if the 3D position is plotted in a new figure window when the
% script is run.
plotFlag = 1;

%% Load Data
disp('STATUS : Loading Data');

% Grab Structures from Tracking Files
s_L = load(leftTrackFile); s_L = s_L.saveStruct; 
disp(['Loaded : ' leftTrackFile]);

s_R = load(rightTrackFile); s_R = s_R.saveStruct;
disp(['Loaded : ' rightTrackFile]);

% Grab label vectors
label_L = s_L.label; label_R = s_R.label;

% Grab Numframes and Inst Structs
numFrames_L = s_L.numFrames; numFrames_R = s_R.numFrames;
inst_L = s_L.inst; inst_R = s_R.inst;

%% Expand structures to account for dropped frames
disp('STATUS : Accounting for dropped frames');
%label_L = label_L(find(label_L==0));
%label_R = label_R(find(label_R==0));
[inst_L,label_L] = frameDropExpand(inst_L,label_L,missingFrames_L);
[inst_R,label_R] = frameDropExpand(inst_R,label_R,missingFrames_R);

% Verify Both Views Have Same Number of Frame
if(numel(label_L) ~= numel(label_R))
    disp('ERROR : Number of frames do not match.  Check that you entered correct missing frames.');
    return;
end

numFrames = numel(label_L);

%% Run Trajectory Estimation
disp('STATUS : Triangulating 3D Trajectory');
gFrames = (label_L == 0) & (label_R == 0);
frames = find(gFrames);

[pos3D,orient3D,pos2D_L,pos2D_R,orient2D_L,orient2D_R] = ...
                       instTraj2(inst_L,label_L,inst_R,label_R,camCalFile);
                   
%% Interpolation
disp('STATUS : Interpolating 3D Data');
% Build Joint label array
% How label_L and label_R are used to generate label :
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

% If you want to transform the coordinate frame from the left camera to
% something else, I would put that code here.

%% Plot Signals
if(plotFlag)
    disp('STATUS : Plotting 3D Data');

    % Determine Frame Subset
    f1 = find(~isnan(X),1,'first'); f2 = find(~isnan(X),1,'last');
    fSub = f1:f2;

    figure();
    subplot(3,1,1); plot(fSub,X(fSub));
    xlabel('Frame'); title('X Axis Trajectory');
    subplot(3,1,2); plot(fSub,Y(fSub));
    xlabel('Frame'); title('Y Axis Trajectory');
    subplot(3,1,3); plot(fSub,Z(fSub));
    xlabel('Frame'); title('Z Axis Trajectory');
    
    figure;                            
    plot3(X(fSub),Y(fSub),Z(fSub));     % View 3D plot if desired
    xlabel('X');                        % Added by A. Vamos 
    ylabel('Y');
    zlabel('Z');
    title('3D Motion trajectory');
end

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
% First the jerk signal is calculated using the 3rd order foward difference
% (estimate of 3rd derivative).  Smoothness is calculated as the square
% root of (1/2) times the sum of this signal squared, normalized by
% the number of frames used to compute the sum.

% X Smoothness
jerk_X = diff(X,3,1);
smoothness_X = (1/sum(~isnan(jerk_X))) * ...
                            sqrt((1/2)*sum(jerk_X(~isnan(jerk_X)).^2));
% Y Smoothness
jerk_Y = diff(Y,3,1);
smoothness_Y = (1/sum(~isnan(jerk_Y))) * ...
                            sqrt((1/2)*sum(jerk_Y(~isnan(jerk_Y)).^2));
                        
% Z Smoothness
jerk_Z = diff(Z,3,1);
smoothness_Z = (1/sum(~isnan(jerk_Z))) * ...
                            sqrt((1/2)*sum(jerk_Z(~isnan(jerk_Z)).^2));
                        
% Net Smoothness
netJerk = sqrt(sum([jerk_X.^2,jerk_Y.^2,jerk_Z.^2],2));
netSmoothness = 1/(sum(~isnan(netJerk))) * ...
                            sqrt((1/2)*sum(netJerk(~isnan(netJerk)).^2));
                        
% Orientation Smoothness
orientJerk_X = diff(X_or,3,1);
orientJerk_Y = diff(Y_or,3,1);
orientJerk_Z = diff(Z_or,3,1);
netOrientJerk = sqrt(sum([orientJerk_X.^2,orientJerk_Y.^2,...
                                                    orientJerk_Z.^2],2));
netOrientSmoothness = (1/sum(~isnan(netOrientJerk))) * ...
           sqrt((1/2) * sum(netOrientJerk(~isnan(netOrientJerk)).^2));
                        
% Display Metrics
disp(['Path Length (mm) : ' num2str(pathLen) ' mm'])
disp(['Depth Perception (mm) : ' num2str(depthPerception)]);
disp(['X Smoothness : ' num2str(smoothness_X)]);
disp(['Y Smoothness : ' num2str(smoothness_Y)]);
disp(['Z Smoothness : ' num2str(smoothness_Z)]);
disp(['Net Motion Smoothness : ' num2str(netSmoothness)]);
disp(['Net Orientation Smoothness : ' num2str(netOrientSmoothness)]);

%% Generate Excel Spreadsheet 
disp('STATUS : Saving Spreadsheet');
% Generate Cell Arrays to be dumped into XLSX Spreadshset %
infoSheet = {'camCalFile',camCalFile;
             'Left Tracking File',s_L.analysisFileName;
             'Right Tracking File',s_R.analysisFileName};

data3DSheet = {'Frame','X (mm)','Y (mm)','Z (mm)'};
data3DSheet = [data3DSheet;num2cell([(1:numel(X)).',X,Y,Z])];

metricsSheet = {'Path Length (mm)',pathLen;
                 'Depth Perception (mm)',depthPerception;
                 'X Smoothness',smoothness_X;
                 'Y Smoothness',smoothness_Y;
                 'Z Smoothness',smoothness_Z;
                 'Net Motion Smoothness',netSmoothness;
                 'Net Orientation Smoothness',netOrientSmoothness};

% Prompt User to Save %           
% Get Starting Save Directory
[pathStr,name] = fileparts(leftTrackFile);
suggestName = fullfile(pathStr,[name '.xlsx']);

[saveName,savePath] = uiputfile('*.xlsx','Save Data to Spreadsheet',suggestName);

if(saveName ~= 0)
    fullSaveName = fullfile(savePath,saveName);
    disp(['STATUS : Saving ' fullSaveName]);
    xlswrite(fullSaveName,infoSheet,'Info');
    xlswrite(fullSaveName,data3DSheet,'Trajectory');
    xlswrite(fullSaveName,metricsSheet,'Metrics');
else
    disp('STATUS : File Save Cancelled');
end

disp('STATUS : Done');
disp('------------------------------------------------------------------');

