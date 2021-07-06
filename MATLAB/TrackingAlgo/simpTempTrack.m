function [tCorner,nCorner,wStats,ais] = simpTempTrack(vidObj,frames,...
                            templateSize,neighborhoodSize,varargin)
% SIMPTEMPTRACK Tracks Template Region in Video
%
% [TCORNER,NCORNER,WSTATS,AIS] = simpTempTrack(...
%    VIDOBJ,FRAMES,TEMPLATESIZE,NEIGHBORHOODSIZE,...
%                                           PropertyName,PropertyValue,...)
%
% Performs simple template tracking from FRAMES(1) to FRAMES(2) of
% VideoReader object VIDOBJ.  The first frame is used for initialziation.
% The user must place a sqaure window of side length TEMPLATESIZE on the
% region to be tracked.  A neighborhood centered at the template of side
% length NEIGHBORHOODSIZE is used to search for the template match in the
% next frame.  The matching process is repeated until FRAMES(2) is reached.
% TCORNER and NCORNER are (numFrames x 2) matrices.  Each row corresponds
% to a frame and represents the upper left hand corner [x,y] of the best
% template match or searched neighborhood for that frame.  WSTATS is a
% structure with fields wSize, nx, and ny representing the template windows
% size and the neighborhood window side lengths (nx = ny).  The 'Update'
% property determines how the template is updated after a best match has
% been found, the 'Color' property determines whether to convert video
% frames to grayscale, and the 'Prediction' property determines how the
% algorithm predicts where the next template will be.
%
% Properties
% 'Update' : {'Replace','None'}
% 'Replace' (default value) - After finding the best match in a frame to
% the template the template is replaced with the best match.
% 'None' - The template initialized from the first frame is used in every
% frame to find the best match.
% 'Color' : {0,1}
% 0 (default) : The Video is grayscale.  Frames will be converted to 
% grayscale after loading.
% 1 : The Video is color nothing is done after loading
% 'Prediction' : {'Constant Position','Constant Velocity'}  The algorithm
% uses a simple motion model and gating to reduce the search space over
% which the template is searched for.  The size of the square gating region
% depends on the parameter NEIGHBORHOODSIZE, where the neighborhood is
% located is determined by this parameter.
% 'Constant Position' : The gating neighborhood is centered on the location
% of the previous template match.
% 'Constant Velocity' : Velocity is estimated from the previous two
% template match locations and the gating neighborhood is displaced from
% the location centered on the previous template match an amount based on
% the velocity estimate.

%% Input Arguments
% Default Arguments
update = 'Replace';
color = 0;
prediction = 'Constant Position';

if nargin > 4
    for k = 1:2:(numel(varargin)-1)
        pName = varargin{k};
        switch pName
            case 'Update'
                update = varargin{k+1};
            case 'Color'
                color = varargin{k+1};
            case 'Prediction'
                prediction = varargin{k+1};
        end
    end
end

% Template Size should be odd
if(templateSize/2 ~= round(templateSize/2))
    % It is odd
    tSize = templateSize;
else
    % It is even
    tSize = templateSize - 1;
end

nSize = round(tSize * neighborhoodSize)-1;
wStats = struct('wSize',[],'nx',nSize+1,'ny',nSize+1);

% Video Related
numFrames = 1 + (frames(2) - frames(1));

% Saves Corner of Search Neighborhood & Best Template Match
tCorner = zeros(numFrames,2);
nCorner = zeros(numFrames,2);

%% Manually Initialize Template
% Take Template as Neighborhood Center
switch color
    case 0
        prevFrame = rgb2gray(vidObj.read(frames(1)));
    case 1
        prevFrame = vidObj.read(frames(1));
end

% User Defined Neighborhood Region
[centerPt,tSize] = squareDraw(prevFrame,tSize);
tempWindow = centerPt - ((tSize-1)/2);
wStats.wSize = tSize;

%{
h_fig1 = figure('Name','Select Matching Neighborhood');
imshow(prevFrame);title('Select Matching Neighborhood');
h_template = imrect(gca,[1,1,tSize-1,tSize-1]);
tempWindow = round(wait(h_template));
delete(h_template); close(h_fig1);
%}

%% Init Windows
tCurr = tempWindow(1:2);
nCurr = tCurr - repmat(round((nSize - tSize)/2),1,2);
temp = prevFrame((0:tSize-1)+tempWindow(2),(0:tSize-1)+tempWindow(1),:);

nCorner(1,:) = nCurr;
tCorner(1,:) = tCurr;

%{
% Place Template Window & Center of Neighborhood
nCrop = [tempWindow(1) - round((nSize-tSize)/2) ...
         tempWindow(2) - round((nSize-tSize)/2) ...
         nSize nSize];

% Compute nccMask to exclude windows that don't fit in neighborhood
nMask = false(nCrop(4)+1,nCrop(3)+1);
nMask((tSize-1)/2:1+(end-(tSize-1)/2),...
      (tSize-1)/2:1+(end-(tSize-1)/2)) = 1; 


% Store Template Neighborhood Corner for 1st Frame
nCorner(1,:) = nCrop(1:2);
tCorner(1,:) = tempWindow(1:2);
%}

%% Run Template Matching
ais = struct('maxNCCScore',ones(numFrames,1));

tic();
for k = 2:numFrames
    % Store Search Neighborhood
    nCorner(k,:) = nCurr;
    
    switch color
        case 0
            nextFrame = rgb2gray(vidObj.read(frames(1)+(k-1)));
        case 1
            nextFrame = vidObj.read(frames(1)+(k-1));
    end
    
    [tNext,nNext,score] = ...
        findTempMatch(nextFrame,temp,tCurr,nCurr,nSize,'ncc');
    
    % Update Template
    switch update
        case 'None'
            % Keep Template Constant
        case 'Replace'
            % Replace Template with Best Match
            temp = nextFrame((0:tSize-1)+tNext(2),(0:tSize-1)+tNext(1),:);
    end
    
    % Update Template
    tCurr = tNext;
    tCorner(k,:) = tCurr;
    
    % Update Search Neighborhood
    if (k > 2)
        switch prediction
            case 'Constant Position'
                nCurr = nNext;
            case 'Constant Velocity'
                nCurr = nNext + (tCorner(k,:) - tCorner(k-1,:));
        end
    else
        nCurr = nNext;
    end
    
    if (nargout == 4)
        % Update the ais
        ais.maxNCCScore(k,1) = score;
    end
    %{
    % Grab Next Frame, Neighborhood, & Template
    nextFrame = vidObj.read(frames(1) + (k-1));
    nextNeighb = imcrop(nextFrame,nCrop);
    temp = imcrop(prevFrame,tempWindow);
    
    [~,nccScore] = template_matching(temp,...
                                            nextNeighb);
    nccScore = nMask .* nccScore;
    [~,maxDex] = max(nccScore(:));
    [maxRow,maxCol] = ind2sub(size(nccScore),maxDex);
    
    % Upper Left Corner Best Match Window in Absolute Coordinates
    maxCoords = nCrop(1,1:2) + [maxCol-1,maxRow-1] +  ...
        -(tSize - 1)/2 * ones(1,2);
    dispVec = maxCoords - tempWindow(1,1:2);
    %disp(dispVec);
    
    % Update Template & Neighborhood
    nCrop(1,1:2) = nCrop(1,1:2) + dispVec;
    tempWindow(1,1:2) = tempWindow(1,1:2) + dispVec;
    
    
    % Store New Corner Locations
    tCorner(k,:) = tempWindow(1:2);
    prevFrame = nextFrame; 
    %}
end
tiempo = toc(); 
fprintf('%4.3f ms per Frame\n',1000*(tiempo/(numFrames-1)));
