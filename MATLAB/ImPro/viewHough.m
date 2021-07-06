function h_fig = viewHough(H,T,R)
% VIEWHOUGH Visualizes Hough Line Parameter Space
%   Detailed explanation goes here
h_fig = figure();

% Display the Hough matrix.
maxVal = max(max(H));
imshow(imadjust(mat2gray(H,[0.5*maxVal,maxVal])),'XData',T,'YData',R,...
      'InitialMagnification','fit');
title('Hough Space');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(hot)

end

