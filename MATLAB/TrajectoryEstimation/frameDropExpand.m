function [instNew,labelNew] = frameDropExpand(inst,label,droppedFrames)
% FRAMEDROPEXPAND Expands instStruct and labelVec to deal w/ dropped frames
%
% [INSTNEW,LABELNEW] = frameDropExpand(INST,LABEL,DROPPEDFRAMES) Expands
% tracking structure INST and array LABEL to account for dropped frames
% during video acquisition.  DROPPEDFRAMES is an array indicating the frame
% numbers with respect to the first frame that were dropped.  See below for
% a more detailed  description of dropped frame numbers.  NaN values are
% pushed into the arrays within INST at indices corrresponding to dropped
% frames.  This structure with expanded elements is returned in INSTNEW.
% Other labeled frames (label = 4) are pushed into LABEL at indices
% corresponding to dropped frames.  The expanded array is returned in
% LABELNEW.
%
% Dropped Frames
% 

% Determine Number of Frame Types %
% Grabbed Frame : A video frame that was grabbed and saved
% Dropped Frame : A video frame that was unsuccessfully grabbed and not
% saved
numDroppedFrames = numel(droppedFrames);
numGrabbedFrames = numel(label);
totalFrames = numDroppedFrames + numGrabbedFrames;

% Build index arrays %
idx = 1:totalFrames;
grabbedIdx = idx; grabbedIdx(droppedFrames) = [];
droppedIdx = idx(droppedFrames);

% Build Output Structures
instNew = struct('rho',NaN*ones(totalFrames,2),...
                 'theta',NaN*ones(totalFrames,2),...
                 'trackPt',NaN*ones(totalFrames,2));
labelNew = NaN*ones(totalFrames,1);

% Fill grabbed data
instNew.rho(grabbedIdx,:) = inst.rho;
instNew.theta(grabbedIdx,:) = inst.theta;
instNew.trackPt(grabbedIdx,:) = inst.trackPt;

labelNew(grabbedIdx) = label;

% Label Dropped Frames as other error
labelNew(droppedIdx) = 4;




