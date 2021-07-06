function errStats = errorStats2D(errSig)
% ERRORSTATS2D Computes Stats for 2D error signal
% 
% ERRSTATS = errorStats2D(ERRSIG) Computes statistics for 2D error signal
% (Nx2) ERRSIG.  The stats are returned in structure ERRSTATS for each
% dimension of the signal (x,y).  See below for description of ERRSTATS
% structure.
%
% ERRSTATS Fields
% x : Structure with fields Max, Mean, Std, RMS corresponding to x dim
% y : Structure with fields Max, Mean, Std, RMS corresponding to y dim

xErr = errSig(:,1); yErr = errSig(:,2);

% figure();
% subplot(2,1,1); plot(xErr); ylabel('x_{err}');
% subplot(2,1,2); plot(yErr); ylabel('y_{err}');

xErrMean = mean(xErr);
[xAbsMax,maxDex] = max(abs(xErr)); xErrMax = xErr(maxDex);
xErrStd = std(xErr);
xErrRMS = sqrt(mean(xErr.^2));

yErrMean = mean(yErr);
[yAbsMax,maxDex] = max(abs(yErr)); yErrMax = yErr(maxDex);
yErrStd = std(yErr);
yErrRMS = sqrt(mean(yErr.^2));

x = struct('Max',xErrMax,'Mean',xErrMean,'Std',xErrStd,'RMS',xErrRMS);
y = struct('Max',yErrMax,'Mean',yErrMean,'Std',yErrStd,'RMS',yErrRMS);
errStats = struct('x',x,'y',y);

%disp('2D x'); disp([xErrMax;xErrMean;xErrStd]);
%disp('2D y');disp([yErrMax;yErrMean;yErrStd]);

%figure();hist(xErr);title('2D X Error Histogram');
%figure();hist(yErr);title('2D Y Error Histogram');