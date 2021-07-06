% simpleSynth_1.m
%% Setup
% Image Size
numRows = 640;
numCols = 480;

% Background
backMean = 220;
backStd = 10;
backSynth = uint8(repmat(backMean,numRows,numCols) + ...
    backStd*randn([numRows,numCols]));

% Foreground Object
objMean = 50;
objStd = 10;
objPattern = uint8(repmat(objMean,numRows,numCols) + ...
    objStd*randn([numRows,numCols]));

rectWidth = 40;
rectHeight = 300;
rotDeg = -40; 

rectMask1 = false([numRows,numCols]);
rectMask1(100:(100+rectHeight),(numCols/2):(numCols/2)+rectWidth) = 1;
rectMask = imrotate(rectMask1,rotDeg,'crop');

% Create Synthestic Image
noiseStd = 10;
synthIm = backSynth;
synthIm(rectMask) = objPattern(rectMask);
h_gauss = fspecial('gaussian',4,2);
synthIm = imfilter(synthIm,h_gauss,'symmetric') + ...
    uint8(noiseStd*randn([numRows,numCols]));

figure();
imshow(synthIm); title('Synthetic Image');
%% Simple Intensity Threshold
instBlob = synthIm < 100;
instBlob = imclose(instBlob,strel('disk',5));

figure();
imshow(instBlob); title('Instrument Blob');

%% Edge Detection
[gradMag,gx,gy] = computeGradientMagIm(synthIm);
gradThresh = 10;
edgePoints = gradMag > gradThresh;

edgeMask = and(instBlob,and((abs(gy) > abs(gx)),edgePoints));
figure();
imshow(edgeMask); title('Boundary Candidates');
%% Line Fitting
% Hough Transform
[H,T,R] = hough(edgeMask);
peaks = houghpeaks(H,2); % Finds Top 2 Line Peaks
lines = houghlines(edgeMask,T,R,peaks);

figure();
imshow(synthIm);
for k = 1:numel(lines)
    hold on;
    line([lines(k).point1(1) lines(k).point2(1)],...
         [lines(k).point1(2) lines(k).point2(2)],'Color','r','LineWidth',1.5);
end
