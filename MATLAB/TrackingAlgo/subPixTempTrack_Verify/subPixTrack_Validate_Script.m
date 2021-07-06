%subPixTrack_Validate_Script.m
%% Test 1
disp('Test 1: +X Unit Steps');
nFrames = 50;
nRows = 200; nCols = 200;
delX = ones(nFrames,1);
delY = zeros(nFrames,1);

[X_RMSE,Y_RMSE] = testHarness1(delX,delY,nRows,nCols);

%% Test 2
disp('Test 2: -X Unit Steps');
nFrames = 50;
nRows = 200; nCols = 200;
delX = ones(nFrames,1);
delY = zeros(nFrames,1);

[X_RMSE,Y_RMSE] = testHarness1(delX,delY,nRows,nCols);

%% Test 3
disp('Test 3: Y Unit Steps');
nFrames = 50;
nRows = 200; nCols = 200;
delX = zeros(nFrames,1);
delY = ones(nFrames,1);

[X_RMSE,Y_RMSE] = testHarness1(delX,delY,nRows,nCols);

%% Test 4
disp('Test 4: -Y Unit Steps');
nFrames = 50;
nRows = 200; nCols = 200;
delX = zeros(nFrames,1);
delY = -ones(nFrames,1);

[X_RMSE,Y_RMSE] = testHarness1(delX,delY,nRows,nCols);

%% Test 5
disp('Test 5: +X 0.5 Steps');
nFrames = 100;
nRows = 200; nCols = 200;
delX = 0.5*ones(nFrames,1);
delY = zeros(nFrames,1);

[X_RMSE,Y_RMSE] = testHarness1(delX,delY,nRows,nCols);

%% Test 6 Random
disp('Test 6: Random Steps');
nFrames = 100;
nRows = 400; nCols = 400;
delX = 3*randn(nFrames,1); delY = 3*randn(nFrames,1);

[X_RMSE,Y_RMSE] = testHarness1(delX,delY,nRows,nCols,5);
