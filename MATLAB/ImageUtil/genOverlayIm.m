function overlayIm = genOverlayIm(grayIm,maskIm)
% GENOVERLAY Generates Color Overlay Image of Mask on Grayscale Im
%
% OVERLAYIM = genOverlayIm(GRAYIM,MASKIM)  Generates an image that is a red
% overlay of binary mask MASKIM on grayscale image GRAYIM.

overlayIm = repmat(grayIm,[1,1,3]);

% Mask In Red
overlayIm(maskIm) = 255;
overlayIm([false(numel(grayIm),1);repmat(maskIm(:),2,1)])= 0;
