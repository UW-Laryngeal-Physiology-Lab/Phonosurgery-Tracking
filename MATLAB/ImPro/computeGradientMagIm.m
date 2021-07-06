function [gradMagIm,gx,gy] = computeGradientMagIm(imData)
% COMPUTEGRADIENTMAGIM Computes Image of Gradient Magnitude
%
% [GRADMAG,GX,GY] = computeGradientMagIm(IMDATA) Computes Gradient
% Magnitude GRADMAGIM, and directional components GX,GY of grayscale image
% IMDATA.

H = (1/8)*fspecial('sobel');
gx = imfilter(double(imData),H,'symmetric');
gy = imfilter(double(imData),H','symmetric');
gradMagIm = sqrt(gx.^2 + gy.^2);


