function errStats = errorStats3D(errSig)
% ERRORSTATS3D Computes stats for 3D error signal
%
% ERRSTATS = errorStats3D(ERRSIG) Computes statistics for 3D error signal
% (Nx3) ERRSIG.  The stats are returned in structure ERRSTATS for each
% dimension of the signal (x,y,z).  See below for description fo ERRSTATS
% structure.
%
% ERRSTATS Fields
% x : Structure with fields Max, Mean, Std, & RMS corresponding to x dim
% y : Structure with fields Max, Mean, Std, & RMS corresponding to y dim
% z : Structure with fields Max, Mean, Std, & RMS corresponding to z dim

xErr = errSig(:,1); yErr = errSig(:,2); zErr = errSig(:,3);
%figure();
%subplot(3,1,1); plot(xErr); ylabel('x_{err}');
%subplot(3,1,2); plot(yErr); ylabel('y_{err}');
%subplot(3,1,3); plot(zErr); ylabel('z_{err}');

xErrMean = mean(xErr);
[xAbsMax,maxDex] = max(abs(xErr)); xErrMax = xErr(maxDex);
xErrStd = std(xErr);
xErrRMS = sqrt(mean(xErr.^2));

yErrMean = mean(yErr);
[yAbsMax,maxDex] = max(abs(yErr)); yErrMax = yErr(maxDex);
yErrStd = std(yErr);
yErrRMS = sqrt(mean(yErr.^2));

zErrMean = mean(zErr);
[zAbsMax,maxDex] = max(abs(zErr)); zErrMax = zErr(maxDex);
zErrStd = std(zErr);
zErrRMS = sqrt(mean(zErr.^2));

x = struct('Max',xErrMax,'Mean',xErrMean,'Std',xErrStd,'RMS',xErrRMS);
y = struct('Max',yErrMax,'Mean',yErrMean,'Std',yErrStd,'RMS',yErrRMS);
z = struct('Max',zErrMax,'Mean',zErrMean,'Std',zErrStd,'RMS',zErrRMS);

errStats = struct('x',x,'y',y,'z',z);

% disp('X'); disp([xErrMax;xErrMean;xErrStd]);
% disp('Y'); disp([yErrMax;yErrMean;yErrStd]);
% disp('Z'); disp([zErrMax;zErrMean;zErrStd]);
% 
% figure();hist(xErr);title('3D X Error Histogram');
% figure();hist(yErr);title('3D Y Error Histogram');
% figure();hist(zErr);title('3D Z Error Histogram');
% 
% results = [xErrMax,yErrMax,zErrMax;...
%            xErrMean,yErrMean,zErrMean;...
%            xErrStd,yErrStd,zErrStd];