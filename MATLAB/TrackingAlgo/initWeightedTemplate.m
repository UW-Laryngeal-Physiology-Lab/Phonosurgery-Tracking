function [tPos,tSize,pc] = initWeightedTemplate(frameIm,tempStartSize)
% INITWEIGHTEDTEMPLATE Initializes Weighted Template
%
% [TPOS,TSIZE,PC] = initWeightedTemplate(FRAMEIM,TEMPSTARTSIZE) Initializes
% template for binary weighted matching interactively.  FRAMEIM is
% displayed and the user must interactively select the template window
% position and size.  The upper left hand corner of the window in [x,y] form is returned
% in TPOS.  The window has size TSIZE.  Following this a window zoomed in
% on the template is displayed.  The user must interactively define a
% bounding polygon which defines the binary weighting mask.  All points
% within the polygon will be considered part of the template (assigned a
% weighting of 1) while all other points are considered background
% (assigned weighting of 0).  The precomputed struct used by template
% matching/tracking routines is returned as PC.  For more information on
% this struct see the function wMatchBuildPcStruct.  

% Get Square Window Center
[centerPt,tSize] = squareDraw(frameIm,tempStartSize);
tPos = centerPt - ((tSize-1)/2);

% Get Weighting Window
temp = frameIm((0:tSize-1)+tPos(2),(0:tSize-1)+tPos(1),:);
h_weightWindow = figure('Name','Define Interest Region');imshow(temp);
h_poly = impoly(gca,'Closed',1);
wait(h_poly);
weights = double(h_poly.createMask());
close(h_weightWindow);

% Build Weighted Matching Precompute Struct
pc = wMatchBuildPcStruct(temp,weights);