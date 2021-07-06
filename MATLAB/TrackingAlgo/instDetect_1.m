function [inst,ais] = instDetect_1(blobIm,mark,param)
%% INSTDETECT_1 Detects single instrument from blobImage using marker
%
% [INST,AIS] = instDetect_1(BLOBIM,MARK,PARAM) Locates blob region in
% BLOBIM representative of surgical instrument.  Instrument marker struct
% MARK guides the detection algorithm.  PARAM is an optional structure
% (with optional fields) with parameters that control the underlying
% algorithm.  See below for description of MARK and PARAM.  INST is a
% binary image containing the single instrument blob.
%
% MARK fields 
% mCorn : corner [x,y] of marker window associated with
% a single instrument 
% mSize : Side length in pixels of marker window
%
% PARAM fields (all optional) 
% axisRatioThresh : (default=1/3) Threshold used to determine whether a
% connected component in BLOBIM has a shape that is instrument like.  It is
% the ratio (minorAxisLength/majorAxisLength) of an ellipse fit to the
% connected component.

%% Input Arguments
% Parameter Defaults
axisRatioThresh = 1/2;

% Assign User Defined Parameters
if nargin == 4
    fNames = fieldnames(param);
    for k = 1:numel(fNames)
        switch(param.(fNames{1}))
            case 'axisRatioThresh'
                axisRatioThresh = param.axisRatioThresh;
        end
    end
end
%% Find Blobs with Instrument Shape Features
% Blob Splitting Based on Vertical Cutoff (Lowest Marker Point)
% Clear Blob Pixels below top of Markers
blobIm(mark.mSize + mark.mCorn(1,2):end,:) = 0;

cc = bwconncomp(blobIm);

% Generate Feature Vector for Blobs
stats = regionprops(cc,'MajorAxisLength','MinorAxisLength','Perimeter');
fVector = [[stats.MinorAxisLength]./[stats.MajorAxisLength];...
           [stats.Perimeter]];

% Sort by Blob Perimeter
% Only Look @ Blobs w/ "axis ratio" below threshold
axisRatio = fVector(1,:) < axisRatioThresh;
[~,sortDex] = sort(fVector(2,:),2,'descend');
sortedCCs = sortDex(axisRatio(sortDex));

% Update the ais
if(nargout == 2)
    % Perimeter of Instrument-Like Blobs
    ais.instLikePerims = fVector(2,sortedCCs);
    
    fprintf('%u Instrument-Like Blobs Detected\n',numel(sortedCCs));
end

%% Assign Blob to Instrument based on Connectivity
imSize = size(blobIm);
inst = [];
for k = 1:numel(sortedCCs)
    connPix = checkConn(mark.mCorn,mark.mSize,...
                                    cc.PixelIdxList{sortedCCs(k)},imSize);
    if(connPix > 0)
        inst = false(imSize);
        inst(cc.PixelIdxList{sortedCCs(k)}) = 1;
        inst(mark.mCorn(1,2):end,:) = 0;
        break;
    end
end

function connPix = checkConn(mCorn,mSize,blobIdx,imSize)
% CHECKCONN Checks if blob pixels fall within marker window
%
% CONNPIX = checkConn(MCORN,MSIZE,BLOB) Checks if binary image BLOB has
% pixels that fall within a window with ul-corner MCORN [x,y] and
% sidelengths MSIZE.  The number of pixels that fall within this window is
% returned by CONNPIX.

[row,col] = ind2sub(imSize,blobIdx);
connPix = sum((row >= mCorn(1,2)) & (row <= (mSize-1) + mCorn(1,2)) & ...
              (col >= mCorn(1,1)) & (col <= (mSize-1) + mCorn(1,1)));