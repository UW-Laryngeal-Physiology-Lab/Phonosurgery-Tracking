function [centerPt,wSize] = squareDraw(imData,wSize)
% SQUAREDRAW Interactive square draw on image based on midpoint
%
% [CENTERPT,WSIZE] = SQUAREDRAW(IMDATA,WSIZE)  Opens a new figure and
% displays IMDATA.  The user selects a point in the image and a square
% centered on the point with sidelength WSIZE is drawn.  The user is
% prompted to change the center point or window size until they are
% satisfied.  The center point in image coordinates and square window side
% length (forced to be odd) are returned in CENTERPT [x,y] AND WSIZE.  If
% the user closes the window before setting point and window size empty
% arguments are returned.

% Force Window Size to Odd Value
wSize = round(wSize);
if(wSize/2 == round(wSize/2))
    % Change to Odd Size
    wSize = wSize - 1;
end

% Draw Image
h_imFig = figure(); h_axes = axes();
imshow(imData,'Parent',h_axes);

% Init Center Point
centerPt = getCp(h_axes);

% Check if user did not select a point
if(isempty(centerPt))
    % Return empty output variables
    centerPt = []; wSize = [];
    
    % Close figure if still open
    if(ishghandle(h_imFig));
        close(h_imFig);
    end
    return; % Leave Function
end

% Draw Box
corners = boxCorners(centerPt,wSize);
h_box = drawBox(corners,h_axes);

% Let the user modify the center point and window size until they press the
% done button.
doneFlag = 0;
while(doneFlag == 0)   
    % Verify Figure is still open.  If not return empty output arguments
    if(~ishghandle(h_imFig))
        centerPt = []; wSize = [];
        return;
    end
    
    % Request Users Next Action
    ButtonName = questdlg('Template Window Setup.  Close window to cancel.', ...
                    'Template Window Setup', ...
                    'Change Window Size','Change Center Point','Done','Done');
    
    switch ButtonName
        case 'Change Window Size'
            delete(h_box);
            wSize = changeWindowSize(wSize);

            % Draw Box
            corners = boxCorners(centerPt,wSize);
            h_box = drawBox(corners,h_axes);

        case 'Change Center Point'
            delete(h_box);
            % Get Center Point
            cpTemp = getCp(h_axes);
            
            % Check if user did not select a point
            if(isempty(cpTemp))
                if(~ishghandle(h_imFig))
                    % If figure was closed exit function with empty output
                    % arguments
                    centerPt = []; wSize = [];
                    return; % Leave Function
                end
            else
                % Valid centerPt
                centerPt = cpTemp;
            end

            % Draw Box
            corners = boxCorners(centerPt,wSize);
            h_box = drawBox(corners,h_axes);

        case 'Done'
            doneFlag = 1;
            
       otherwise
            % User Closed Dialog.  Treat as Cancel.
            centerPt = []; wSize = [];
            doneFlag = 1;
    end
end
close(h_imFig);


function cp = getCp(h_axes)
% GETCP User Interactive Selection of Window Center Point
% User Select Point
h_point = impoint(h_axes);

% Get Position of Point.  Return empty array if user closed window prior to
% selecting point.
cp = [];
if(~isempty(h_point))
    cp = round(h_point.getPosition());
end

delete(h_point);

function corners = boxCorners(squareCenter,wSize)
% BOXCORNERS Corners of Box in [x,y]
corners = (wSize - 1)/2 * [-1,-1;1,-1;1,1;-1,1] + repmat(squareCenter(:).',4,1);

function h_box = drawBox(corners,h_axes)
% DRAWBOX Draws box on axes based on corners
db = [1,2,4,3,2,4,1,3];
h_box = line(corners(db,1),corners(db,2),'Color','b','Parent',h_axes);

function wSize = changeWindowSize(currSize)
% CHANGEWINDOWSIZE Interactive Dialog for Changing Window Size
inputAns = inputdlg('Window Size','Change Window Size',1,{num2str(currSize)});

if(isempty(inputAns))
    wSize = currSize;
else
    wSize = round(str2num(inputAns{1}));
    if(wSize/2 == round(wSize/2))
        % Change to Odd Size
        wSize = wSize - 1;
    end
end

