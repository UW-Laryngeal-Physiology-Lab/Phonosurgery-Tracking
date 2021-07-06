function [inst1,inst2,ais] = instDetect_2(blobIm,mark1,mark2,param)
% INSTDETECT Detects 2 instruments in blob image based on markers
%
% [INST1,INST2,AIS] = instDetect(BLOBIM,MARK1,MARK2,PARAM) Locates blob
% regions in BLOBIM representative of two instruments.  Instrument marker
% structs MARK1 & MARK2 guide the detection algorithm.  PARAM is an
% optional structure (with optional fields) with parameters that control
% the underlying algorithm.  See below for descriptions of MARK1/MARK2 and
% PARAM.  INST1 & INST2 are binary images containing single instrument
% blobs.
%
% MARK1/MARK2 fields 
% mCorn : corner [x,y] of marker window associated with a single instrument
% mSize : Side length in pixels of marker window
%
% PARAM fields (all optional) 
% axisRatioThresh : (default=1/3) Threshold used to determine whether a
% connected component in BLOBIM has a shape that is instrument like.  It is
% the ratio (minorAxisLength/majorAxisLength) of an ellipse fit to the
% connected component.

%% Input Arguments
% Parameter Defaults
axisRatioThresh = 1/3;

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

% Blob Splitting Based on Vertical Connectivity ?
%% Find Blobs with Instrument Shape Features
% Blob Splitting Based on Vertical Cutoff (Lowest Marker Point)
% Clear Blob Pixels below top of Markers
blobIm(max(mark1.mSize,mark2.mSize) + ...
       max(mark1.mCorn(1,2),mark2.mCorn(1,2)):end,:) = 0;

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
if(nargout == 3)
    % Perimeter of Instrument-Like Blobs
    ais.instLikePerims = fVector(2,sortedCCs);
    
    fprintf('%u Instrument-Like Blobs Detected\n',numel(sortedCCs));
end

%% Assign 2 Largest Blobs to instruments based on connectivity
blob_1 = false(size(blobIm)); blob_1(cc.PixelIdxList{sortedCCs(1)}) = 1;
blob_2 = false(size(blobIm)); blob_2(cc.PixelIdxList{sortedCCs(2)}) = 1;

connMat = [checkConn(mark1.mCorn,mark1.mSize,blob_1),...
           checkConn(mark1.mCorn,mark1.mSize,blob_2);...
           checkConn(mark2.mCorn,mark2.mSize,blob_1),...
           checkConn(mark2.mCorn,mark2.mSize,blob_2)];

% Update ais
if(nargout == 3)
    ais.connMat = connMat;
end
           
% Check for Valid Connectivity
check1 = sum(connMat > 0,1); check2 = sum(connMat > 0,2);
checkVal = sum([check1,check2']);
if(checkVal ~= 4)
    disp('Detection Error.  Invalid Connectivity to Markers');
    inst1 = [];inst2 = []; return;
end

%% Return Data Structures
switch (connMat(1,1) > 0)
    case 1
        inst1 = blob_1;inst1(mark1.mCorn(1,2):end,:) = 0;
        inst2 = blob_2;inst2(mark2.mCorn(1,2):end,:) = 0;
    case 0
        inst1 = blob_2;inst1(mark2.mCorn(1,2):end,:) = 0;
        inst2 = blob_1;inst2(mark1.mCorn(1,2):end,:) = 0;
end

function connPix = checkConn(mCorn,mSize,blob)
% CHECKCONN Checks if blob pixels fall within marker window
%
% CONNPIX = checkConn(MCORN,MSIZE,BLOB) Checks if binary image BLOB has
% pixels that fall within a window with ul-corner MCORN [x,y] and
% sidelengths MSIZE.  The number of pixels that fall within this window is
% returned by CONNPIX.

connPix = sum(sum(blob((0:mSize-1)+mCorn(1,2),(0:mSize-1)+mCorn(1,1))));
