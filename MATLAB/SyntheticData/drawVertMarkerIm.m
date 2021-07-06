function rectIm = drawVertMarkerIm(rectUl,rectHeight,rectWidth,...
    markerIm,markerHeight,instPattern,backPattern)

rectIm = backPattern;
rows = (0:(rectHeight-1)) + rectUl(2);
cols = (0:(rectWidth-1)) + rectUl(1);

rectIm(rows,cols) = instPattern(rows,cols);

markerRows = (1:markerHeight) + (rows(end)-round(1.5*markerHeight));
rectIm(markerRows,cols) = markerIm;

% Simulate Noise & PSF
noiseStd = 10;
h_gauss = fspecial('gaussian',4,2);
rectIm = imfilter(rectIm,h_gauss,'symmetric') + ...
    uint8(noiseStd*randn(size(rectIm)));
