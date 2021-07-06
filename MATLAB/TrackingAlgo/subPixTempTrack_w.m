function [tCorner,nCorner,subTCorner,wStats,ais] = ...
    subPixTempTrack_w(vidObj,frames,templateSize,neighborhoodSize,varargin)
% SUBPIXTEMPTRACK Tracks weighted template in Video @ subpixel resolution
%
% [TCORNER,NCORNER,SUBTCORNER,WSTATS,AIS] = subPixTempTrack(
%                           VIDOBJ,FRAMES,TEMPLATESIZE,NEIGHBORHOODSIZE,
%                           PropertyName,PropertyValue)
% Performs template tracking from FRAMES(1) to FRAMES(2) of videoReader
% object VIDOBJ.  The first frame is used for initialization.  The user
% must set the center of a square window of side length TEMPLATESIZE on the
% region to be tracked.  After selecting the window a polygon over the
% window is drawn by the used to define the interest region of the
% template.  This is used to remove the influence of background regions on
% the template.  A neighborhood centered at the template of side length
% NEIGHBORHOODSIZE is used to search for the template match in the next
% frame.  The matching process is repeated until FRAMES(2) is reached.
% TCORNER, SUBT, NCORNER are (numFrames x 2) matrices.  Each row
% corresponds to a frame and represents the upper left hand corner (x,y) of
% the best template match or searched neighborhood for that frame.  TCORNER
% and NCORNER are at pixel resolution.  SUBT is the subpixel resolution
% match of the template.  WSTATS is a structure with fields wSize, nx, and
% ny representing the template window's size and the neighborhood window
% side lengths (nx = ny).  The 'Update' property determines how the
% template is updated after a best match has been found, the 'Color'
% property determines whether to convert video frames to grayscale, and the
% 'Prediction' property determines how the algorithm predicts where the
% next template will be.  See below for description of properties.
%
% Properties
% 'Update' : {'Replace','None'}
% 'Replace' (default value) After finding the best match in a frame to
% the template the template is replaced with the best match.
% 'None' - The template initialized from the first frame is used in every
% frame to find the best match.
%
% 'Color' : {0,1}
% 0 (default) : The Video is grayscale.  Frames will be converted to 
% grayscale after loading.
% 1 : The Video is color nothing is done after loading
%
% 'Prediction' : {'Constant Position','Constant Velocity'}  The algorithm
% uses a simple motion model and gating to reduce the search space over
% which the template is searched for.  The size of the square gating region
% depends on the parameter NEIGHBORHOODSIZE.  Where the neighborhood is
% located is determined by this parameter.
% 'Constant Position' : The gating neighborhood is centered on the location
% of the previous template match.
% 'Constant Velocity' : Velocity is estimated from the previous two
% template match locations and the gating neighborhood is displaced from
% the location centered on the previous template match an amount based on
% the velocity estimate.
%
% 'Template Window' : {WINDOWSTRUCT} Normally the algorithm uses an
% interactive window to have the user select the template within the first
% frame.  This option allows the location of the template window to instead
% be passed as an argument.  A structure windowStruct is used to pass the
% template window center point and size.  See below for a description of
% WINDOWSTRUCT.
% WINDOWSTRUCT Fields :
% centerPt : [x,y] center of template window location
% tSize : Side length in pixels of template window.   The value is rounded
% up to the nearest odd integer to force the center pt. to the center of
% the template window.


%% Input Arguments
% Default Arguments
update = 'Replace';
color = 0;
prediction = 'Constant Position';
useWindowStruct = 0;

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
            case 'Template Window'
                windowStruct = varargin{k+1};
                useWindowStruct = 1;
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
subTCorner = zeros(numFrames,2);

%% Manually Initialize Template
% Take Template as Neighborhood Center
switch color
    case 0
        prevFrame = rgb2gray(vidObj.read(frames(1)));
    case 1
        prevFrame = vidObj.read(frames(1));
end

if(useWindowStruct == 0)
    % User Defined Template Window
    [centerPt,tSize] = squareDraw(prevFrame,tSize);
else
    % User Passed Template Window
    centerPt = windowStruct.centerPt;
    tSize = round(windowStruct.tSize);
    if(tSize/2 == round(tSize/2))
        % Force Odd Size
        tSize = tSize + 1;
    end
end    
tempWindow = centerPt - ((tSize-1)/2);
wStats.wSize = tSize;

%% Init Windows
tCurr = tempWindow(1:2);
nCurr = tCurr - repmat(round((nSize - tSize)/2),1,2);
temp = prevFrame((0:tSize-1)+tempWindow(2),(0:tSize-1)+tempWindow(1),:);

% Determine Template Weighting Mask
h_weightWindow = figure('Name','Define Interest Region');imshow(temp);
h_poly = impoly(gca,'Closed',1);
wait(h_poly);
weights = double(h_poly.createMask());
close(h_weightWindow);

pc = wMatchBuildPcStruct(double(temp),weights);

nCorner(1,:) = nCurr;
tCorner(1,:) = tCurr;
subTCorner(1,:) = tCurr;

%% Run Template Matching
ais = struct('maxNCCScore',ones(numFrames,1),'fitMSE',zeros(numFrames,1));

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
    
    try
        [tNext,nNext,score,fitMSE] = ...
            findTempMatch_w(nextFrame,pc,nCurr,nSize);
    catch ME
        disp(['Error w/ subpixel tracking: ' num2str(k)]);
        throw(ME);
    end
    
    subTCorner(k,:) = tNext;
    tNext = round(tNext);
    nNext = round(nNext);

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
    
    if (nargout == 5)
        % Update the ais
        ais.maxNCCScore(k,1) = score;
        ais.fitMSE(k,1) = fitMSE;
    end
end
tiempo = toc(); 
fprintf('%4.3f ms per Frame\n',1000*(tiempo/(numFrames-1)));