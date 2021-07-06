function [tempNew,neighbNew,maxScore,fitMSE,simSurf] = ...
             findTempMatch_sub(frameIm,template,tCurr,nCurr,nSize,varargin)
% FINDTEMPMATCH_SUB Finds best match to template in frame neighborhood
%
% [TEMPNEW,NEIGHBNEW,MAXSCORE,SIMSURF] = ...
%       findTempMatch_sub(FRAMEIM,TEMPLATE,TCURR,NCURR,NSIZE,MEASURE='ncc')
% Given an image template TEMPLATE and search neighborhood window (upper
% left hand corner [x,y]->NCURR, NSIZE) finds best match in image FRAMEIM
% using similarity measure.  Utilizes location of template window in
% previous image [x,y]->TCURR (window upper left hand corner) to find new
% template window upper left hand corner [x,y]->TEMPNEW and search
% neighborhood upper left hand corner [x,y]->NEIGHBNEW.  Optional variable
% MEASURE can be 'ncc' or 'ssd' and determines the similarity measure.
%
% Input Variables
% FRAMEIM Image to search for template match in
% TEMPLATE Template image to be matched
% TCURR 2 element array [x,y].  Defines location of upper left hand corner 
% of template image (window)
% NCURR 2 element array [x,y].  Define location of upper left hand corner
% of search window.
% NSIZE The width and height of search window.
% MEASURE Optional variable that determines similarity metric.  When not
% set NCC is used, can be defined also as 'ssd'.
%
% Output Variables
% TEMPNEW Upper left hand corner of best template match in FRAMEIM of same
% size as TEMPLATE
% NEIGHBNEW Upper left hand corner of next search window.  A window with
% width and height NSIZE with upper left corner NEIGHBNEW is centered on the
% window defined by TEMPNEW.
% MAXSCORE The similairty score associated with the match.
% SIMSURF Similarity Surface used to find best match location.  Either NCC
% or SSD depending on user input.

%% Input Arguments
% Determine Similarity Measure
switch nargin
    case 6
        % User Inputted Similarity Type
        switch varargin{1}
            case 'ssd'
                sim = 'ssd';
            otherwise
                sim = 'ncc';
        end
    otherwise
        sim = 'ncc';
end

tSize = size(template,1); % Assumes Square Template

% Define Pixels in search neighborhood.  Prevent Searching outside the
% image
nx = max([1,nCurr(1)]):min(nCurr(1)+(nSize-1),size(frameIm,2));
ny = max([1,nCurr(2)]):min(nCurr(2)+(nSize-1),size(frameIm,1));

[ssdScore,nccScore] = template_matching(template,...
                        frameIm(ny,nx,:));

% Neighborhood Mask
% Compute nccMask to exclude windows that don't fit in neighborhood
nMask = false(numel(ny),numel(nx));
nMask(1+((tSize-1)/2):(end-(tSize-1)/2),...
      (1+(tSize-1)/2):(end-(tSize-1)/2)) = 1; 

switch sim
    case 'ncc'
        score = nMask .* nccScore;
        [maxScore,maxDex] = max(score(:));
        [maxRow,maxCol] = ind2sub(size(score),maxDex);
        if(nargout == 5); simSurf = score; end;
    case 'ssd'
        score = nMask .* ssdScore;
        [maxScore,maxDex] = max(score(:));
        [maxRow,maxCol] = ind2sub(size(score),maxDex);
        if(nargout == 5); simSurf = score; end;
end

% Subpixel Estimation
[subAdjust,fitMSE] = peakfit2d_mod(score( (-1:1)+maxRow, (-1:1)+maxCol));
maxCol = subAdjust(1) + maxCol;
maxRow = subAdjust(2) + maxRow;


% Upper Left Corner Best Match Window in Absolute Coordinates
%maxCoords = [nx(1),ny(1)] + [maxCol-1,maxRow-1] +  ...
%    -(tSize - 1)/2 * ones(1,2);
maxCoords = [nx(1),ny(1)] + ...
            (([maxCol,maxRow] - ((tSize - 1)/2 * ones(1,2))) - ones(1,2));
%dispVec = maxCoords - tCurr;
%disp(dispVec);

% Update Template & Neighborhood
%neighbNew = nCurr + dispVec;
%tempNew = tCurr + dispVec;
tempNew = maxCoords;
neighbNew = tempNew - (round((nSize - tSize)/2)*[1,1]);

