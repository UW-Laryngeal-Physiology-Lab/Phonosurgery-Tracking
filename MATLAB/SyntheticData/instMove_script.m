%% Image Parameters
% Image Size
numRows = 640;
numCols = 480;
imSize = [numRows,numCols];

% Background Texture
backMean = 210;
backStd = 10;
backStats = struct('mean',backMean,'std',backStd);

% Instrument Texture
instMean = 50;
instStd = 10;
instStats = struct('mean',instMean,'std',instStd);

% Marker Texture
mStats = struct('lowMean',70,'lowStd',10,'highMean',200,'highStd',10);

% Instrument & Marker Size
rectWidth = 40;
rectHeight = 300;
markerHeight = 50;

% Texture Patterns & Marker Image
[instPattern,backPattern] = genPatterns(instStats,backStats,imSize);
markerIm = genMarkerIm(rectWidth,markerHeight,mStats);

%% Instrument Displacement
% Instrument Starting Location
X0 = [100,10];

% Instrument Velocity
nPts = 22; % Choose So Instrument Stays in Image
dx = 10; dy = 8;
delT = [0,0;dx*ones(nPts,1),dy*ones(nPts,1)];

% Compute Upper Left Corner of Instrument Versus Time
rectUL = repmat(X0,size(delT,1),1) + cumsum(delT);
%% Generate Rectangle Movement Movie
instVid = zeros([imSize,1,size(rectUL,1)],'uint8');

for k = 1:size(rectUL,1)
    instVid(:,:,1,k) = ...
        drawVertMarkerIm(rectUL(k,:),rectHeight,rectWidth,markerIm,...
            markerHeight,instPattern,backPattern);
end
clear vidObj;
vWriter = VideoWriter('tempVid.avi','Uncompressed AVI');
vWriter.open();
vWriter.writeVideo(instVid);
vWriter.close();
vidObj = VideoReader('tempVid.avi');

%% Track It!
[tCorner,nCorner,wStats] = simpTempTrack(vidObj,[1,nPts],71,2);
tVid = trackVid(vidObj,[1,nPts],tCorner,nCorner,wStats);
implay(tVid);
