function [N,X] = trackPeakNeighbors(mark,vidObj,frames)
% TRACKPEAKNEIGHBORS Histogram of tracking peak spatial neighbors
%
% [N,X] = trackPeakNeighbors(MARK,VIDOBJ,FRAMES) Uses tracking results
% structure MARK to find the NCC values of all spatial neighbors to
% detected peaks.  Returned histogram data [N,X] represents the fraction of
% peak value obtained by the neighbors.  X is the domain of the histogram
% in fraction of peak value.  N is the fraction of all neighbors that
% obtained the peak value @ X.  VIDOBJ is a videoReader object
% corresponding to the monocular sequence tracked and FRAMES is a 2 element
% vector corresponding to the subset over VIDOBJ that was tracked.

numFrames = 1 + (frames(2) - frames(1));

%% Iterate Through Peaks
rowAdd = [-ones(3,1);zeros(2,1);ones(3,1)];
colAdd = [(-1:1)';[-1;1];(-1:1)'];

peakNeighbs = zeros(numFrames-1,8);

for k = 1:numFrames-1
    [~,nccScore] = trackEval(mark,vidObj,frames,k+1);
    [maxVal,maxIdx] = max(nccScore(:));
    neighbIdx = (maxIdx + colAdd*size(nccScore,1)) + rowAdd;
    peakNeighbs(k,:) = nccScore(neighbIdx)/maxVal;
end

%% Create Histogram
[N,X] = hist(peakNeighbs(:),(0.0:0.01:1.1));
N = N/numel(peakNeighbs);
figure();stem(X,N,'Marker','None');
xlabel('Fraction of Peak Value');
ylabel('Fraction of Neighbors Over Frame Sequence');
title('NCC Peak Neighbors Histogram');