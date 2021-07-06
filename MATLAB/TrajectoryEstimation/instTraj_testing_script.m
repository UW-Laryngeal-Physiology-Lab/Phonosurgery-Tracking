% instTraj_testing_script.m
% Tests that new trajectory estimation routine instTraj generates the same
% data as old routine instTriang.  

%% Generate Trajectory Data
% New Routine
[pos3D,orient3D,pos2D_L,pos2D_r,orient2D_L,orient2D_R] = instTraj();

% Old Routine
[startPt,~,startOr,~,sp2D_l,sp2D_r] = instTriang();

%% Verify
% Determine Common Frames
cFrames = ~isnan(pos3D(:,1)) & ~isnan(startPt(:,1));

error2 = @(sig1,sig2) sum(sum((sig1 - sig2).^2,2),1);

% Compute error
pos3D_err = error2(pos3D(cFrames,:),startPt(cFrames,:));
fprintf('3D Position Error : %4.4f\n',pos3D_err);

or3D_err = error2(orient3D(cFrames,:),startOr(cFrames,:));
fprintf('3D Orientation Vector Error : %4.4f\n',or3D_err);






