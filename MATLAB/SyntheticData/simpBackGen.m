function backIm = simpBackGen(imSize,backInfo)
% SIMPBACKGEN Generates Simple Synthetic Background Image
%
% BACKIM = simpBackGen(IMSIZE,BACKINFO) Generates background image of size 
% IMSIZE ([rows,cols]) that follows parameters in BACKINFO (gaussian
% model).
%
% BACKINFO Fields
% Mean : Mean of background intensity (gaussian model)
% Std : Intensity Std (gaussian model)

backIm = uint8(backInfo.mean + backInfo.std*randn(imSize));

