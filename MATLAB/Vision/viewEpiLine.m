function epiIm = viewEpiLine(destIm,sourcePt,F)
% VIEWEPILINE Draws Image with Epipolar line overlayed
%
% EPIIM = viewEpiLine(DESTIM,SOURCEPT,F)  Draws Epipolar line in grayscale
% image DESTIM corresponding to SOURCEPT.  F is the fundamental matrix
% describing the transformation of a point in the source image to a line in
% the destination image.  If no given no output arguments a figure with the
% overlay image is automatically generated.

imSize = size(destIm);
lineMask = epiLineMask(sourcePt,F,imSize);
epiIm = genOverlayIm(destIm,lineMask);

if nargout == 0
    figure();
    imshow(epiIm);
end