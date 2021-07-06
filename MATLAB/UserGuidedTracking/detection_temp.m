function [tCurr,nCurr,maxScore,fitMSE] = detection_temp(frameIm,pc,wStats)
% DETECTION_TEMP Detects weighted template in entire frame image

[tCurr,~,maxScore,fitMSE] = findTempMatch_w(frameIm,pc,[1,1],...
                                                    max(size(frameIm)));
nCurr = tCurr - (round(((wStats.nx-1) - (wStats.wSize))/2)*[1,1]);
