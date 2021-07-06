%surgicalMetrics.m
%% Load Raw Data
[pos3D,orient3D] = instTraj();
%% Find Frame Subset
frameStart = find(~isnan(pos3D(:,1)),1,'first');
frameEnd = find(~isnan(pos3D(:,1)),1,'last');
f = frameStart:frameEnd;

%% Get Interpolated Signal
pos3D_sub = pos3D(f,:);
orient3D_sub = orient3D(f,:);

gFrames = ~isnan(pos3D_sub(:,1));

interpFun = @(sig,idx) interp1(f(gFrames).',sig(gFrames,idx),f.','linear');

X_interp = interpFun(pos3D_sub,1);
Y_interp = interpFun(pos3D_sub,2);
Z_interp = interpFun(pos3D_sub,3);

dX_interp = interpFun(orient3D_sub,1);
dY_interp = interpFun(orient3D_sub,2);
dZ_interp = interpFun(orient3D_sub,3);

pos = [X_interp,Y_interp,Z_interp];
orient = [dX_interp,dY_interp,dZ_interp];

%% Basic Plots
% Plot Raw Position
figure('Name','Instrument Position');
subplot(3,1,1);
plot(f,pos(:,1)); title('X Position (mm)');
subplot(3,1,2);
plot(f,pos(:,2)); title('Y Position (mm)');
subplot(3,1,3);
plot(f,pos(:,3)); title('Z Position (mm)');

figure('Name','Instrument Orientation');
subplot(3,1,1);
plot(f,orient(:,1)); title('X Orientation');
subplot(3,1,2);
plot(f,orient(:,2)); title('Y Orientation');
subplot(3,1,3);
plot(f,orient(:,3)); title('Z Orientation');

% Plot Velocity
vel = sqrt(sum(diff(pos,1,1).^2,2));
figure('Name','3D Velocity');
plot(f(2:end),vel);
xlabel('Frame Number');
ylabel('Velocity (mm/s)');
title('Velocity');

%% Time
disp(['Frames : ' num2str(numel(f))]);

%% Path Length
pathLen = sum(vel);
disp(['Path Length(mm): ' num2str(pathLen)]);

%% Depth Perception
dispSig = diff(pos,1,1);
dp_sig = sum(dispSig .* orient(1:end-1,:),2);

figure();
plot(f(1:end-1),dp_sig);title('Depth Perception Signal (mm)');
disp(['Net Depth Perception (mm): ' num2str(sum(abs(dp_sig)))]);

%% Smoothness
jerk = diff(vel,2,1);
jerk_sig = sqrt(sum(jerk.^2,2));

figure();
plot(jerk_sig);title('Jerk Signal');
disp(['Smoothness : ' num2str(sum(jerk_sig))]);