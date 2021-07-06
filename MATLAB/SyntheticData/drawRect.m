function rectIm = drawRect(rectCorner,instPattern,backPattern,imSize)
% DRAWRECT Draws Instrument Rectangle 
%
% RECTIM = drawRect(RECTCORNER,INSTSTATS,BACKSTATS,IMSIZE)

rectIm = backPattern;
%% Generate Mask of Instrument
instMask = false(imSize);

% Boundary Lines

bCell = cell(4,1);
bCell{1,1} = intLineSeg(rectCorner(1,:),rectCorner(2,:));
bCell{2,1} = intLineSeg(rectCorner(2,:),rectCorner(3,:));
bCell{3,1} = intLineSeg(rectCorner(3,:),rectCorner(4,:));
bCell{4,1} = intLineSeg(rectCorner(1,:),rectCorner(4,:));

bLocs = [];
for k = 1:4
    [x,y] = bCell{k,1}.getIntLocs();
    bLocs = [bLocs;sub2ind(imSize,y,x)];
end

instMask(bLocs) = 1;
instMask = imfill(instMask,'holes');

%% Generate RectIm
rectIm(instMask) = instPattern(instMask);

% Create Synthestic Image
noiseStd = 5;
h_gauss = fspecial('gaussian',4,2);
rectIm = imfilter(rectIm,h_gauss,'symmetric') + ...
    uint8(noiseStd*randn(imSize));


