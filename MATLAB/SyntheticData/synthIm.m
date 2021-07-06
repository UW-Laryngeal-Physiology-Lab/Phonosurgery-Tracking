function sIm = synthIm(instIm,instMask,backIm,hSize,sigma)
% SYNTHIM Synthesizes Simulated Instrument Image 
%
% SIM = synthIm(INSTIM,INSTMASK,BACKIM,HSIZE,SIGMA) Sythesizes simulated
% instrument image using instrument image and mask (location of instrument
% in image) INSTIM and INSTMASK.  The instrument is placed on the
% background image BACKIM.  Applies gaussian smoothing with kernel of HSIZE
% and width SIGMA to model PSF.

% Place Instruemtn on background
backIm(instMask) = instIm(instMask);

% Model PSF
h_gauss = fspecial('gaussian',hSize,sigma);
sIm = imfilter(backIm,h_gauss,'symmetric') + ...
    uint8(sigma*randn(size(backIm)));
