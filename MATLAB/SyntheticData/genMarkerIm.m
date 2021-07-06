function markerIm = genMarkerIm(mWidth,mHeight,mStats)
% GENMARKERIM Generates Image of Vertical Marker
%
% MARKERIM = genMarkerIm(MWIDTH,MHEIGHT)  Generates image of marker with
% rectangular dimensions (pixels) MWIDTH x MHEIGHT.

markerIm = randi([0,1],[mHeight,mWidth]);
highImage = markerIm == 1;
numHigh = sum(sum(highImage));
numLow = numel(markerIm) - numHigh;


markerIm(highImage) = mStats.highMean + mStats.highStd*randn(numHigh,1);
markerIm(~highImage) = mStats.lowMean + mStats.lowStd*randn(numLow,1);
markerIm = uint8(markerIm);











