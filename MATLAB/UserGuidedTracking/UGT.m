function varargout = UGT(varargin)
% UGT MATLAB code for UGT.fig
%      UGT, by itself, creates a new UGT or raises the existing
%      singleton*.
%
%      H = UGT returns the handle to a new UGT or the handle to
%      the existing singleton*.
%
%      UGT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UGT.M with the given input arguments.
%
%      UGT('Property','Value',...) creates a new UGT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UGT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UGT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singlet on)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UGT

% Last Modified by GUIDE v2.5 21-Oct-2015 12:01:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UGT_OpeningFcn, ...
                   'gui_OutputFcn',  @UGT_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before UGT is made visible.
function UGT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UGT (see VARARGIN)

% Choose default command line output for UGT
handles.output = hObject;

% Setup Finite State Machine %
% States
% idle : Application has opened, waiting for new tracking or tracking to be
% loaded
% waitForTemplate : Video has been loaded for new tracking.  Waiting for
% user to initialize the template.
% waitForLabel : Template has been initialized and tracking has been
% perormed on current frame.  Waiting for user to label frame.
% automation : An automation routine (autoPilot, haveFaith, or multiError)
% is running
handles.state = 'idle';
set(hObject,'Units','pixels');
set(handles.statusEdit,'String','Waiting for File');

% Automation Stop Flag
% Used to stop automation prior to reaching the stop frame
handles.autoStopFlag = 0;

% Setup Auto-Pilot
% threshVec = [delta_NCC,delta_width,delta_orient,delta_inliersDistL,...
%              delta_inliersDistR,inliersL,inliersR]
% Threshold Parameters Set 2011-11-25
handles.threshVec = [0.0332,1.72,0.00987,10,10,10,10];
handles.typeVec = [0,0,0,0,0,1,1];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UGT wait for user response (see UIRESUME)
% uiwait(handles.UGT);


function confVec = getConfParams(mark,inst,algoInfo,frameNum)
% GETCONFPARAMS Computes Vector of Confidence Parameters
%
% CONFVEC = getConfParams(INST,ALGOINFO,FRAMENUM) Computes vector of
% confidence parameters for a given frame number FRAMENUM.  Tracking data
% structures MARK, INST, and ALGOINFO are used to compute individual
% confidence parameters.  Tracking data should be available at FRAMENUM and
% FRAMENUM-1 in these data structures.  A vector CONFVEC is returned, see
% below for description.
%
% CONFVEC = [delta_NCC, delta_width, delta_Orient, delta_inliersDistL, 
%            delta_inliersDistR, leftInliers,rightInliers]
% 2011-11-25 KS

% Previous & Current Frame
f = [frameNum-1,frameNum];

% Compute Confidence Parameters %
% NCC
delta_NCC = abs(diff(algoInfo.nccScore(f),1,1));

% Inst Width
widthSig = diff(inst.rho(f,:),1,2);
delta_widthSig = abs(diff(widthSig,1,1));

% Inst Orientation
orientSig = mean(inst.theta(f,:),2);
delta_orientSig = abs(diff(orientSig,1,1));

% inliersDist
markY = mark.tCorn(f,2);
inliersDistL = markY - algoInfo.leftNumInliers(f);
inliersDistR = markY - algoInfo.rightNumInliers(f);
delta_inliersDistL = abs(diff(inliersDistL,1,1));
delta_inliersDistR = abs(diff(inliersDistR,1,1));

% Build Confidence Parameter Vector
confVec = [
           delta_NCC,...
           delta_widthSig,...
           delta_orientSig,...
           delta_inliersDistL,...
           delta_inliersDistR,...
           algoInfo.leftNumInliers(frameNum),...
           algoInfo.rightNumInliers(frameNum)
           ];
           
function result = checkAlgoConf(handles,frameNum)
% CHECKALGOCONF Checks Algorithm Confidence Parameters
%
% RESULT = CHECKALGOCONF(HANDLES,FRAMENUM) Checks algorithm tracking
% confidence at FRAMENUM.  HANDLES is the gui HANDLES structure containing
% mark, inst, algoInfo data structures with tracking data and threshVec and
% typeVec vectors.  RESULT is a binary vector, value = 0 means that the
% algorithm is confident about the parameter.  A value = 1 means the
% algorithm underconfident about the parameter.

% Get Confidence Parameter Vector
fVec = getConfParams(handles.mark,handles.inst,handles.algoInfo,frameNum);

% Get Result of Confidence Classification
result = threshClassifyFun(fVec,handles.threshVec,handles.typeVec);

% --- Outputs from this function are returned to the command line.
function varargout = UGT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in previousFrameButton.
function previousFrameButton_Callback(hObject, eventdata, handles)
% hObject    handle to previousFrameButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(strcmp(handles.state,{'automation','idle'}))
    % Do Nothing if tracking file has not been started or
    % automation is running
    return;
end

if(handles.currFrame > 1)
    handles.currFrame = handles.currFrame - 1;
    displayFun(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in nextFrameButton.
function nextFrameButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextFrameButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(strcmp(handles.state,{'automation','idle'}))
    % Do Nothing if tracking file has not been started or
    % automation is running
    return;
end

if(handles.currFrame < handles.numFrames)
    handles.currFrame = handles.currFrame + 1;
    displayFun(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in initializeTemplateButton.
function initializeTemplateButton_Callback(hObject, eventdata, handles)
% hObject    handle to initializeTemplateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% If idle or running automation do nothing
if(any(strcmp(handles.state,{'idle','automation'})))
    return;
end

% If waiting for label, verify with user that re-initializing the template
% will result in a loss of tracking data.
if(strcmp(handles.state,'waitForLabel'))
    % Prompt User
    msgAns = questdlg('Re-initializing the template will close the current tracking file.  Any unsaved progress will be lost.  Continue?',...
    'Re-initialize Template','Yes','No','No');
    
    % Do nothing if user does not want to reinitialize template
    if(~strcmp('Yes',msgAns))
        return;
    end

    % Clear Analysis File Name
    handles.analysisFileName = [];
    set(handles.analysisFileEdit,'String','');
    
    % Clear Data Structures
    handles = initNonVideoStructures(handles);
    guidata(hObject,handles);
    displayFun(handles);
    
end

currFrame = handles.currFrame;

% Interactively Get Template Window
defStartSize = 99; % Default Window Starting Size
frameIm = rgb2gray(handles.vidObj.read(handles.currFrame));
[centerPt,tSize] = squareDraw(frameIm,defStartSize);

% Check for user cancel
if(isempty(centerPt))
    errordlg('Cancelled Template Initializiation Process','Cancelled Init','modal');
    return;
end

% Upper Left Corner of Template
tempWindow = centerPt - ((tSize-1)/2);

nSize = round(tSize * handles.neighborhoodSize)-1;
wStats = struct('wSize',tSize,'nx',nSize+1,'ny',nSize+1);
nCurr = tempWindow - repmat(round((nSize - tSize)/2),1,2);
temp = frameIm((0:tSize-1)+tempWindow(2),(0:tSize-1)+tempWindow(1),:);

% Interactively get Template Weighting Mask
h_weightWindow = figure('Name','Define Interest Region');imshow(temp);
h_poly = impoly(gca,'Closed',1);
wait(h_poly);

% If User closed window do not continue template initialization process
if(isempty(h_poly))
    if(ishghandle(h_weightWindow))
        close(h_weightWindow);
    end
    errordlg('Cancelled Template Initializiation Process','Cancelled Init','modal');
    return;
end

% Store Weighting Mask and Close Window
weights = double(h_poly.createMask());
close(h_weightWindow);

% Template Related Data Structures
handles.pc = wMatchBuildPcStruct(double(temp),weights);
handles.mark.tCorn(currFrame,:) = tempWindow;
handles.mark.subT(currFrame,:) = tempWindow;
handles.mark.wStats = wStats;
handles.mark.nCorn(currFrame,:) = nCurr;
handles.templateFrame = currFrame;

% Template Related Algo Info
handles.algoInfo.nccScore(currFrame,1) = 1;
handles.algoInfo.fitMSE(currFrame,1) = 0;

% Estimate Boundary Lines
[r,t,dbl_ais] = detection_bLines(frameIm,handles.backIm,tempWindow,wStats,handles.pc);
handles.inst.rho(currFrame,:) = r;
handles.inst.theta(currFrame,:) = t;
handles.inst.trackPt(currFrame,:) = computeTrackPt(r,t,tempWindow,...
                                                    round((tSize-1)/2));
% Boundary Line Algo Info
handles.algoInfo.votesLeft(currFrame,:) = dbl_ais.votesLeft;
handles.algoInfo.votesRight(currFrame,:) = dbl_ais.votesRight;
handles.algoInfo.leftNumInliers(currFrame,:) = dbl_ais.leftNumInliers;
handles.algoInfo.leftRefitMSE(currFrame,:) = dbl_ais.leftRefitMSE;
handles.algoInfo.rightNumInliers(currFrame,:) = dbl_ais.rightNumInliers;
handles.algoInfo.rightRefitMSE(currFrame,:) = dbl_ais.rightRefitMSE;

% Set State to waitForLabel
handles.currentLabel = currFrame;
handles.state = 'waitForLabel';
set(handles.statusEdit,'String','Waiting for Label');

displayFun(handles);
guidata(hObject,handles);

function currFrameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to currFrameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currFrameEdit as text
%        str2double(get(hObject,'String')) returns contents of currFrameEdit as a double

% Do Nothing if idle or automation
if(strcmp(handles.state,{'idle','automation'}))
    return;
end

frameNum = round(str2double(get(hObject,'String')));
if(frameNum >= 1 && frameNum <= handles.numFrames)
    handles.currFrame = frameNum;
    set(handles.currFrameEdit,'String',num2str(frameNum));
    displayFun(handles);
else
    set(handles.currFrameEdit,'String',num2str(handles.currFrame));
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function currFrameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currFrameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function totalFramesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to totalFramesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of totalFramesEdit as text
%        str2double(get(hObject,'String')) returns contents of totalFramesEdit as a double


% --- Executes during object creation, after setting all properties.
function totalFramesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to totalFramesEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in labelCorrectButton.
function labelCorrectButton_Callback(hObject, eventdata, handles)
% hObject    handle to labelCorrectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~strcmp(handles.state,'waitForLabel'))
    return;
end

handles = labelCurr(handles,0);
displayFun(handles);
guidata(hObject,handles);


% --- Executes on button press in labelBlurButton.
function labelBlurButton_Callback(hObject, eventdata, handles)
% hObject    handle to labelBlurButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~strcmp(handles.state,'waitForLabel'))
    return;
end

handles = labelCurr(handles,1);
displayFun(handles);
guidata(hObject,handles);

% --- Executes on button press in labelOcclusionButton.
function labelOcclusionButton_Callback(hObject, eventdata, handles)
% hObject    handle to labelOcclusionButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~strcmp(handles.state,'waitForLabel'))
    return;
end

handles = labelCurr(handles,2);
displayFun(handles);
guidata(hObject,handles);

% --- Executes on button press in labelOutOfFrameButton.
function labelOutOfFrameButton_Callback(hObject, eventdata, handles)
% hObject    handle to labelOutOfFrameButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~strcmp(handles.state,'waitForLabel'))
    return;
end

handles = labelCurr(handles,3);
displayFun(handles);
guidata(hObject,handles);

% --- Executes on button press in labelOtherButton.
function labelOtherButton_Callback(hObject, eventdata, handles)
% hObject    handle to labelOtherButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~strcmp(handles.state,'waitForLabel'))
    return;
end

handles = labelCurr(handles,4);
displayFun(handles);
guidata(hObject,handles);

function handles = labelCurr(handles,labelNum)
% LABELCURR Labels the Current Frame & Runs tracking on next frame
%
% HANDLES = labelCurr(HANDLES,LABELNUM)  Labels the current frame as
% integer LABELNUM -> {0:Correct,1:Blur,2:Occlusion,3:OtherError}.  After
% labeling the tracking is run on the next frame and the current frame and
% label are updated.  If the current frame is not the current label the
% user is prompted if they want to change the label.  If this is desired
% the current frame label is changed and all proceeding frames are
% unlabeled.
if((strcmp(handles.state,'waitForLabel') || strcmp(handles.state,'automation')) && ...
          handles.currFrame <= handles.currentLabel)
  if(handles.currFrame == handles.currentLabel)
       % Label The Unlabeled Frames
       handles.label(handles.currFrame) = labelNum;
       handles = performTracking(handles,handles.currFrame + 1);
       handles.currFrame = handles.currFrame + 1;
       handles.currentLabel = handles.currFrame;
       %displayFun(handles);
  else
      % Label & Clear
      btn = questdlg('This frame has already been labeled.  By labeling it, all labels set after it will be lost.  Are you sure?',...
          'Change already set label?','No');
      if(strcmp(btn,'Yes'))
          handles.label(handles.currFrame) = labelNum;
          handles = clearAfter(handles,handles.currFrame);
          handles = performTracking(handles,handles.currFrame + 1);
          handles.currFrame = handles.currFrame + 1;
          handles.currentLabel = handles.currFrame;
          %displayFun(handles);
      end 
  end  
end

function labelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to labelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of labelEdit as text
%        str2double(get(hObject,'String')) returns contents of labelEdit as a double


% --- Executes during object creation, after setting all properties.
function labelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function previousFrameLabelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to previousFrameLabelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of previousFrameLabelEdit as text
%        str2double(get(hObject,'String')) returns contents of previousFrameLabelEdit as a double


% --- Executes during object creation, after setting all properties.
function previousFrameLabelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to previousFrameLabelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function algorithmTypeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to algorithmTypeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of algorithmTypeEdit as text
%        str2double(get(hObject,'String')) returns contents of algorithmTypeEdit as a double


% --- Executes during object creation, after setting all properties.
function algorithmTypeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to algorithmTypeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function fileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to fileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function currVidDir = getCurrVidDir(handles)
% GETCURRVIDDIR Gets the Current Video's Directory
% 
% CURRVIDDIR = getCurrVidDir(HANDLES) Checks if a video file is associated
% with the HANDLES structure.  If so returns directory name containing
% video in CURRVIDDIR, otherwise returns 0.

% Check if vidName exists
if(isfield(handles,'vidName'))
    % Check that it is not empty
    if(~isempty(handles.vidName))
        % Get Directory and Return
        currVidDir = fileparts(handles.vidName);
        return;
    end
end

% Otherwise Return 0
currVidDir = 0;

function handles = initNonVideoStructures(handles)
% INITNONVIDEOSTRUCTURES Initializes all non-video data structures and vars
%
% HANDLES = initNonVideoStructures(handles)  Initializes tracking
% structures : label(NaNs), mark(NaNs), inst(NaNs), and algoInfo(NaNs).
% Also initializes template matching structures neighborhoodSize(2),
% pc(empty), and templateFrame(empty).  Sets the axis to display frames at
% the native size.  Initializes the currentLabel variable (empty).  Sets
% system state to 'waitForTemplate'.

nFrames = handles.numFrames;

% Analysis File
handles.analysisFileName = [];
set(handles.analysisFileEdit,'String','');

% Set Axis Dimensions to that of Image
set(handles.dispAxes,'Units','pixels');
pos = get(handles.dispAxes,'Position');
pos(3:4) = handles.imSize([2,1]);
set(handles.dispAxes,'Position',pos);

% Init Structs
handles.label = nan*ones(nFrames,1);
handles.mark = struct('tCorn',nan*ones(nFrames,2),'wStats',[],...
    'nCorn',nan*ones(nFrames,2),'subT',nan*ones(nFrames,2));
handles.inst = struct('rho',nan*ones(nFrames,2),'theta',...
    nan*ones(nFrames,2),'trackPt',nan*ones(nFrames,2));
handles.algoInfo = struct('nccScore',nan*ones(nFrames,1),'fitMSE',...,
    nan*ones(nFrames,1),'votesLeft',nan*ones(nFrames,1),'votesRight',...
    nan*ones(nFrames,1),'leftNumInliers',nan*ones(nFrames,1),...
    'leftRefitMSE',nan*ones(nFrames,1),'rightNumInliers',...
    nan*ones(nFrames,1),'rightRefitMSE',nan*ones(nFrames,1));

% Template Related
handles.neighborhoodSize = 2;
handles.pc = []; % Precomputation Struct Used for Template Tracking
handles.templateFrame = [];

% Analysis Progress Variables
handles.currentLabel = [];

% Set State
handles.state = 'waitForTemplate';
set(handles.statusEdit,'String','Waiting For Initialize');

% -------------------------------------------------------------------
function newTrackingMenu_Callback(hObject, eventdata, handles)
% hObject    handle to newTrackingMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Prompt User for Video File Location
currVidDir = getCurrVidDir(handles);
if(currVidDir ~= 0)
    % Open File Browser in Current Video Directory
    [vidName,vidPath] = uigetfile('*.avi','Load Video File',currVidDir);
else
    [vidName,vidPath] = uigetfile('*.avi','Load Video File');
end

if(vidName == 0); return; end;

[backName,backPath] = uigetfile('*.avi','Load Background File',vidPath);
if(backName == 0); return; end;

% Load Video & Store Video Parameters
handles.vidName = fullfile(vidPath,vidName);
handles.backName = fullfile(backPath,backName);
handles.vidObj = VideoReader(handles.vidName);
handles.backObj = VideoReader(handles.backName);
handles.imSize = [handles.vidObj.Height,handles.vidObj.Width];
nFrames = handles.vidObj.NumberOfFrames; handles.numFrames = nFrames;

% Store Backframe
handles.backIm = rgb2gray(handles.backObj.read(1));

% Display Video File Info
set(handles.videoFileEdit,'String',handles.vidName);

% Display Variables
handles.currFrame = 1;

% Set Total Frame Edit
set(handles.totalFramesEdit,'String',num2str(nFrames));

% Initialize all non-video data structures.  These structures are used to
% store tracking data and the tracking file name
handles = initNonVideoStructures(handles);
guidata(hObject,handles);

% Show First Frame
displayFun(handles);

function currTrackFileDir = getCurrTrackFileDir(handles)
% GETCURRTRACKFILEDIR Gets the Current Track File's Directory
% 
% CURRTRACKFILEDIR = getCurrTrackFileDir(HANDLES) Checks if a track file is
% associated with the HANDLES structure.  If so returns directory name
% containing track file in CURRTRACKFILEDIR, otherwise returns 0.

% Check if analysisFileName exists
if(isfield(handles,'analysisFileName'))
    % Check that it is not empty
    if(~isempty(handles.analysisFileName))
        % Get Directory and Return
        currTrackFileDir = fileparts(handles.analysisFileName);
        return;
    end
end

% Otherwise Return 0
currTrackFileDir = 0;



% --------------------------------------------------------------------
function loadTrackingMenu_Callback(hObject, eventdata, handles)
% hObject    handle to loadTrackingMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currTrackFileDir = getCurrTrackFileDir(handles);
if(currTrackFileDir ~= 0)
    % Open File Browser in Current Track File Directory
    [fName,pName] = uigetfile('*.mat','Select Tracking File',...                                                    
                                                    currTrackFileDir);
else
    % No Tracking File Avaialble.  Check Video File
    currVidDir = getCurrVidDir(handles);
    if(currVidDir ~= 0)
        [fName,pName] = uigetfile('*.mat','Select Tracking File',...
                                        currVidDir);
    else
        [fName,pName] = uigetfile('*.mat','Select Tracking File');
    end
end

if(fName == 0)
    return;
end

load(fullfile(pName,fName),'saveStruct');

saveFields = {'mark','pc','neighborhoodSize','templateFrame','inst',...
'algoInfo','vidName','backName','imSize','numFrames',...
'backIm','state','currentLabel','label','currFrame','analysisFileName'};

% Load all Save Fields
for k = 1:numel(saveFields)
    handles.(saveFields{k}) = saveStruct.(saveFields{k});
end

% Check if tracking file name or path has changed

% Setup Video Objects
bExist = exist(handles.vidName,'file');
if bExist == 0
    [vidName,vidPath] = uigetfile('*.avi','Load Video File');
    handles.vidName = fullfile(vidPath,vidName);
end
handles.vidObj = VideoReader(handles.vidName);

bExist = exist(handles.backName,'file');
if bExist == 0
    [backName,backPath] = uigetfile('*.avi','Load Video File');
    handles.backName = fullfile(backPath,backName);
end
handles.backObj = VideoReader(handles.backName);

displayFun(handles);

% Set Axis Dimensions to that of Image
set(handles.dispAxes,'Units','pixels');
pos = get(handles.dispAxes,'Position');
pos(3:4) = handles.imSize([2,1]);
set(handles.dispAxes,'Position',pos);

% Display Video File & Analysis File Names
set(handles.videoFileEdit,'String',handles.vidName);
set(handles.analysisFileEdit,'String',handles.analysisFileName);

% Set Label to waitForLabel
handles.state = 'waitForLabel';
set(handles.statusEdit,'String','Waiting for Label');

guidata(hObject,handles);

% --------------------------------------------------------------------
function saveTrackingMenu_Callback(hObject, eventdata, handles)
% hObject    handle to saveTrackingMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Do Nothing if Idle
if(strcmp(handles.state,'idle'))
    return;
end

if(isempty(handles.analysisFileName))
    currVidDir = getCurrVidDir(handles);
    if(currVidDir ~= 0)
        % Open File Browser in current video directory
        [fName,pName] = uiputfile('*.mat','Save Tracking File',currVidDir);
    else
        [fName,pName] = uiputfile('*.mat','Save Analysis File');
    end
    
    % Check for User Cancel
    if(fName == 0)
        return;
    end
    
    handles.analysisFileName = fullfile(pName,fName);
end

saveFields = {'mark','pc','neighborhoodSize','templateFrame','inst',...
'algoInfo','vidName','backName','imSize','numFrames',...
'backIm','state','currentLabel','label','currFrame','analysisFileName'};

saveStruct = handles;
allFields = fieldnames(saveStruct);

% Remove all non save fields
for k = 1:numel(allFields)
    cName = allFields{k};
    if(~any(strcmp(cName,saveFields)))
        saveStruct = rmfield(saveStruct,cName);
    end
end

save(handles.analysisFileName, 'saveStruct');
set(handles.analysisFileEdit,'String',handles.analysisFileName);
guidata(hObject,handles);

% --------------------------------------------------------------------
function saveTrackingAsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to saveTrackingAsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Do Nothing if Idle
if(strcmp(handles.state,'idle'))
    return;
end

currTrackFileDir = getCurrTrackFileDir(handles);
if(currTrackFileDir ~= 0)
    [fName,pName] = uiputfile('*.mat','Save Tracking File',...
                                                    currTrackFileDir);
else
    % No Tracking File Available
    currVidDir = getCurrVidDir(handles);
    if(currVidDir ~= 0)
        % Open File Browswer in current video directory
        [fName,pName] = uiputfile('*.mat','Save Tracking File',currVidDir);
    else
        [fName,pName] = uiputfile('*.mat','Save Tracking File');
    end
end

% Check for User Cancel
if(fName == 0)
    return;
end

handles.analysisFileName = fullfile(pName,fName);

saveFields = {'mark','pc','neighborhoodSize','templateFrame','inst',...
'algoInfo','vidName','backName','imSize','numFrames',...
'backIm','state','currentLabel','label','currFrame','analysisFileName'};

saveStruct = handles;
allFields = fieldnames(saveStruct);

% Remove all non save fields
for k = 1:numel(allFields)
    cName = allFields{k};
    if(~any(strcmp(cName,saveFields)))
        saveStruct = rmfield(saveStruct,cName);
    end
end

save(handles.analysisFileName, 'saveStruct');
set(handles.analysisFileEdit,'String',handles.analysisFileName);
guidata(hObject,handles);

% --------------------------------------------------------------------
function exitMenu_Callback(hObject, eventdata, handles)
% hObject    handle to exitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% If a tracking is open double check that the user want to quit
if(strcmp(handles.state,'idle'))
    delete(handles.UGT);
else
    user_response = questdlg('Are you sure you want to quit?',...
        'Exit Confirmation','Yes','No','No');
    switch(user_response)
        case 'Yes'
            delete(handles.UGT);
    end
end

function videoFileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to videoFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of videoFileEdit as text
%        str2double(get(hObject,'String')) returns contents of videoFileEdit as a double

% --- Executes during object creation, after setting all properties.
function videoFileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to videoFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function analysisFileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to analysisFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of analysisFileEdit as text
%        str2double(get(hObject,'String')) returns contents of analysisFileEdit as a double


% --- Executes during object creation, after setting all properties.
function analysisFileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analysisFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in currentLabelButton.
function currentLabelButton_Callback(hObject, eventdata, handles)
% hObject    handle to currentLabelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Do Nothing if not waiting for label
if(~strcmp(handles.state,'waitForLabel'))
    return;
end

% Go to frame that needs to be labeled
handles.currFrame = handles.currentLabel;

% Display Frame
displayFun(handles);
guidata(hObject,handles);

function displayFun(handles)
% DISPLAYFUN Displays current frame and overlays tracking data
%
% displayFun(HANDLES) Displays the current frame and overlays the template
% window and boundary line estimates if available.

currFrame = handles.currFrame;

% Draw Marker Mask %
tCurr = handles.mark.tCorn(currFrame,:);
nCurr = handles.mark.nCorn(currFrame,:);
imSize = handles.imSize;

if(~any(isnan(tCurr)))
    tSize = handles.mark.wStats.wSize;
    nSize = handles.mark.wStats.nx-1;
    markMask = trackMask(tCurr,nCurr,tSize,nSize,imSize);
else
    markMask = false(imSize);
end

% Draw Boundary Line Mask %
rho = handles.inst.rho(currFrame,:);
theta = handles.inst.theta(currFrame,:);
trackPt = handles.inst.trackPt(currFrame,:);
if(~any(isnan([rho,theta,trackPt])))
    lineMask = drawLineMask(imSize,[rho,mean(rho)],[theta,mean(theta)]);
    lineMask(round(tCurr(2)):end,:) = 0;
    lineMask((-5:5) + round(trackPt(2)),(-5:5) + round(trackPt(1))) = 1;
else
   lineMask =  false(imSize);
end

% Draw Full Image
frameIm = rgb2gray(handles.vidObj.read(currFrame));
fullIm = genOverlayIm(frameIm,(lineMask | markMask));
imshow(fullIm,'InitialMagnification',100,'Parent',handles.dispAxes);

% Update GUI Objects
% Algorithm Edit
if(currFrame == handles.templateFrame)
    set(handles.algorithmTypeEdit,'String','Template Frame');
elseif(currFrame == 1)
    set(handles.algorithmTypeEdit,'String','');
else
    set(handles.algorithmTypeEdit,'String',...
        getAlgoTypeString(handles.label(currFrame-1)));
end

% Label & Previous Frame Label
set(handles.labelEdit,'String',getLabelString(handles.label(currFrame)));
if(currFrame ~= 1)
    set(handles.previousFrameLabelEdit,'String',...
        getLabelString(handles.label(currFrame-1)));
else
    set(handles.previousFrameLabelEdit,'String','');
end

% Current Frame
set(handles.currFrameEdit,'String',num2str(currFrame));
set(handles.totalFramesEdit,'String',num2str(handles.numFrames));

function algoString = getAlgoTypeString(prevFrameLabel)
% GETALGOTYPESTRING Determines Algo Type based on previous frame label
%
% algoString = getAlgoTypeString(PREVFRAMELABEL) Converts numeric previous
% frame label to the algorithm type used for template & boundary line
% detection in the current frame.

switch(prevFrameLabel)
    case 0
        algoString = 'Tracking';
    case {1,2,3}
        algoString = 'Detection';
    otherwise
        algoString = '';
end

function labelString = getLabelString(label)
% GETLABELSTRING Gets the label string corresponding to numeric label
%
% labelString = getLabelString(LABEL)  Converts numeric label to its
% corresponding string descriptor.

switch(label)
    case 0
        labelString = 'Correct';
    case 1
        labelString = 'Error Blur';
    case 2
        labelString = 'Error Occlusion';
    case 3
        labelString = 'Error Out of Frame';
    case 4
        labelString = 'Error Other';
    otherwise
        labelString = '';
end

function handles = performTracking(handles,frameNum)
% PERFORMTRACKING Runs Instrument Tracking/Detection on single frame
%
% HANDLES = performTracking(HANDLES,FRAMENUM) Instrument tracking/detection
% algorithm.  The previous frame label determines whether detection or
% tracking is performed.  All instrument tracking data structures in
% returned HANDLES are updated.

if(isnan(handles.label(frameNum-1)))
    return;
end

frameIm = rgb2gray(handles.vidObj.read(frameNum));

fail_bLine = 0;
switch(handles.label(frameNum-1))
    case 0
        % Perform Tracking
        % Template
        tPrev = handles.mark.tCorn(frameNum-1,:);
        if(frameNum > 2 && handles.label(frameNum - 2) == 0)
            tPrev2 = handles.mark.tCorn(frameNum-2,:);
        else
            tPrev2 = tPrev;
        end
        trackPtPrev = handles.inst.trackPt(frameNum-1,:);
        thetaPrev = handles.inst.theta(frameNum-1,:);
        nPrev = handles.mark.nCorn(frameNum-1,:);
        nSize = handles.mark.wStats.nx-1;
        
        [tCurr,nCurr,maxScore,fitMSE] = track_temp(frameIm,handles.pc,tPrev,...
                                            tPrev2,nPrev,nSize);
                                        
        % Boundary Line Estimates
        try
            [r,t,dbl_ais] = track_bLines(frameIm,handles.backIm,...
                                    round(tCurr),tPrev,trackPtPrev,thetaPrev);
        catch ME
            fail_bLine = 1;
        end
    otherwise
        % Perform Detection
        % Template Detection
        [tCurr,nCurr,maxScore,fitMSE] = detection_temp(frameIm,...
            handles.pc,handles.mark.wStats);
        
        % Boundary Lines
        % Use Boundary Line Detector Based on Template Mask Orientation
        
        % Determine Last Correct Labeling
        % Use Last Correct to Seed Detection
        %lastCorr = find(handles.label(1:frameNum-1) == 0,1,'last');
        %tCorn_lastCorr = handles.mark.tCorn(lastCorr,:);
        %trackPt_lastCorr = handles.inst.trackPt(lastCorr,:);
        %thetaPrev = mean(handles.inst.theta(lastCorr,:));
        %trackPtDelta = trackPt_lastCorr - tCorn_lastCorr;
        try
            [r,t,dbl_ais] = detection_bLines(frameIm,handles.backIm,...
                round(tCurr),handles.mark.wStats,handles.pc);
            %[r,t,dbl_ais] = track_bLines(frameIm,handles.backIm,...
            %    round(tCurr),tCurr,tCurr + trackPtDelta,thetaPrev);
        catch ME
            fail_bLine = 1;
        end
end

% Store Results
handles.mark.tCorn(frameNum,:) = round(tCurr);
handles.mark.subT(frameNum,:) = tCurr;
handles.mark.nCorn(frameNum,:) = round(nCurr);

% Template Related Algo Info
handles.algoInfo.nccScore(frameNum,1) = maxScore;
handles.algoInfo.fitMSE(frameNum,1) = fitMSE;

if(fail_bLine == 0)
% Boundary Line Estimates
handles.inst.rho(frameNum,:) = r;
handles.inst.theta(frameNum,:) = t;
halfSize = round((handles.mark.wStats.wSize-1)/2);
handles.inst.trackPt(frameNum,:) = computeTrackPt(r,t,tCurr,halfSize);

% Boundary Line Algo Info
handles.algoInfo.votesLeft(frameNum,:) = dbl_ais.votesLeft;
handles.algoInfo.votesRight(frameNum,:) = dbl_ais.votesRight;
handles.algoInfo.leftNumInliers(frameNum,:) = dbl_ais.leftNumInliers;
handles.algoInfo.leftRefitMSE(frameNum,:) = dbl_ais.leftRefitMSE;
handles.algoInfo.rightNumInliers(frameNum,:) = dbl_ais.rightNumInliers;
handles.algoInfo.rightRefitMSE(frameNum,:) = dbl_ais.rightRefitMSE;
end

function handles = clearAfter(handles,afterFrame)
% CLEARAFTER Clears all stored data after a specific frame
%
% HANDLES = clearAfter(HANDLES,AFTERFRAME) Clears labels and stored
% tracking data in all frames after AFTERFRAME.

f = afterFrame + 1;

handles.label(f:end,:) = nan;

% Store Results
handles.mark.tCorn(f:end,:) = nan;
handles.mark.subT(f:end,:) = nan;
handles.mark.nCorn(f:end,:) = nan;

% Template Related Algo Info
handles.algoInfo.nccScore(f:end,1) = nan;
handles.algoInfo.fitMSE(f:end,1) = nan;

% Boundary Line Estimates
handles.inst.rho(f:end,:) = nan;
handles.inst.theta(f:end,:) = nan;
%halfSize = round((handles.mark.wStats.wSize-1)/2);
handles.inst.trackPt(f:end,:) = nan;

% Boundary Line Algo Info
handles.algoInfo.votesLeft(f:end,:) = nan;
handles.algoInfo.votesRight(f:end,:) = nan;
handles.algoInfo.leftNumInliers(f:end,:) = nan;
handles.algoInfo.leftRefitMSE(f:end,:) = nan;
handles.algoInfo.rightNumInliers(f:end,:) = nan;
handles.algoInfo.rightRefitMSE(f:end,:) = nan;

% --- Executes on key press with focus on UGT or any of its controls.
function UGT_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to UGT (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% Allow Prev Frame & Next Frame to be controlled by 'o' and 'p'
switch(eventdata.Key)
    case 'p'
        if(handles.currFrame < handles.numFrames)
            handles.currFrame = handles.currFrame + 1;
            displayFun(handles);
            guidata(hObject,handles);
        end
    case 'o'
        if(handles.currFrame > 1)
            handles.currFrame = handles.currFrame - 1;
            displayFun(handles);
            guidata(hObject,handles);
        end
end

function statusEdit_Callback(hObject, eventdata, handles)
% hObject    handle to statusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of statusEdit as text
%        str2double(get(hObject,'String')) returns contents of statusEdit as a double


% --- Executes during object creation, after setting all properties.
function statusEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statusEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% autoPilotButton_Callback
% Runs the Auto-Pilot automation routine.
%
% The auto-pilot routine is ran from the current frame waiting for a label
% to a user supplied endFrame.  It is assumed the first frame is tracked
% correctly.  For each proceeding frame, tracking is performed.  Algorithm
% confidence parameters are calculated from the newly found tracking
% parameters.  The confidence parameters are used to classify algorithm
% confidence.  If it is not confident, the user is prompted for
% intervention.  If the algorithm is confident, it continues this process
% until the end frame.  Auto-Pilot can be stopped by pressing the "Stop
% Auto" button.

% --- Executes on button press in autoPilotButton.
function autoPilotButton_Callback(hObject, eventdata, handles)
% hObject    handle to autoPilotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Only run if template has been initialized
if(~strcmp(handles.state,'waitForLabel'))
    return;
end

% Get End Frame
startFrame = handles.currentLabel;
handles.currFrame = startFrame;
displayFun(handles);
dlgString = ['Start Frame : ' num2str(startFrame) ...
             '.  Please Select End Frame'];
         
 if(isfield(handles,'prevAutoPilotStop'))
     answer = inputdlg(dlgString,'Auto-Pilot',1,{num2str(handles.prevAutoPilotStop)});
 else
     answer = inputdlg(dlgString,'Auto-Pilot');
 end

if(isempty(answer))
    % User Cancelled
    return;
end

% Check Valid Answer
stopFrame = round(str2num(answer{1}));
if(stopFrame > handles.currentLabel && stopFrame <= handles.numFrames)
   
    % Set Program State to Automation
    handles.state = 'automation';
    
    % Set Status Indicator
    set(handles.statusEdit,'String','Auto Pilot');
    handles.prevAutoPilotStop = stopFrame;
    
    % Set Automation Indicators
    set(handles.automationTypeEdit','String','Auto Pilot');
    set(handles.automationStartEdit,'String',num2str(startFrame));
    set(handles.automationStopEdit,'String',num2str(stopFrame));
    
    % Label the First Frame Correct
    handles = labelCurr(handles,0);
    set(handles.currFrameEdit,'String',num2str(handles.currFrame));
    
    % Save guidata
    guidata(hObject,handles);
    
    % Get guidata to be used in loop
    % Note a local handles copy (h_loop) is used.  The handles structure is
    % updated at the end of this function.  This is done to avoid potential
    % problems with the handles structure updating when the "Stop Auto"
    % button is pressed.
    h_loop = guidata(hObject);
    
    % Auto-Pilot Loop
    for k = h_loop.currFrame:stopFrame
        % Form Feature Vector for Confidence Classifier
        % fVec = [delta_NCC,delta_width,delta_orient,delta_inliers,...
                                                    %  inliersL,inliersR]
        %{
        f = [k-1,k];
        fVec = [abs(diff(h_loop.algoInfo.nccScore(f),1,1)),...
            abs(diff(diff(h_loop.inst.rho(f,:),1,2),1,1)),...
            abs(diff(mean(h_loop.inst.theta(f,:),2),1,1)),...
            diff(h_loop.algoInfo.leftNumInliers(f) + h_loop.algoInfo.rightNumInliers(f),1,1),...
            h_loop.algoInfo.leftNumInliers(k),...
            h_loop.algoInfo.rightNumInliers(k)];
        
        % Get Result of Confidence Classification
        result = threshClassifyFun(fVec,h_loop.threshVec,h_loop.typeVec);
        %}
                                                    
        % Check Algorithm Confidence Parameters
        result = checkAlgoConf(h_loop,k);
                                                    
        if(any(result))
            % Algorithm is Underconfident %
            % Display feature vector thresholding in command window
            disp(result);
            
            % Update Display with current frame
            displayFun(h_loop);
            
            % Display Dialog Asking for labelling
            figPos = get(h_loop.UGT,'Position');
            btn = lowAlgoConfDialog(figPos(1),figPos(2));
            
            if(strcmp(btn,'No'))
                % User Wants to Stop
                break;
            end
        end
        
        % Algorithm is Confident %
        % Label frame as correct
        h_loop = labelCurr(h_loop,0);
        
        set(h_loop.currFrameEdit,'String',num2str(h_loop.currFrame));
        
        % Stop Auto-Pilot if user presses Stop Automation %
        % Get Latest Handles
        h_latest = guidata(hObject);
        
        % Check Stop Flag
        if(h_latest.autoStopFlag == 1)
            % Stop Pressed.  Break
            h_loop.autoStopFlag = 0;
            break;
        end
    end
    
    % Display a Messge Box if auto-pilot ran until stopframe
    if(k == stopFrame)
        % Display Message Box Indicating Completion
        msgbox('Auto-Pilot Complete!','Auto-Pilot Complete!','modal');
    end
    
    % Set State to Wait For Label
    h_loop.state = 'waitForLabel';
    
    handles = h_loop;
    
    % Update Display
    displayFun(handles);
    
    % Update Automation Indicators
    set(handles.automationTypeEdit,'String','');
    set(handles.automationStartEdit,'String','');
    set(handles.automationStopEdit,'String','');
    
    % Update Status
    set(handles.statusEdit,'String','Waiting for Label');
    
else
    msgbox('Incorrect End Frame Entered','Incorrect End Frame','error','modal');
end
guidata(hObject,handles);
    
% --- Executes on button press in haveFaithButton.
function haveFaithButton_Callback(hObject, eventdata, handles)
% hObject    handle to haveFaithButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Only run if template has been initialized
if(~strcmp(handles.state,'waitForLabel'))
    return;
end

startFrame = handles.currentLabel;
dlgString = ['Start Frame : ' num2str(startFrame) ...
             '.  Please Select End Frame'];
answer = inputdlg(dlgString,'Have Faith');

if(isempty(answer))
    % User Cancelled
    return;
end

% Check Valid Answer
stopFrame = round(str2num(answer{1}));
if(stopFrame > handles.currentLabel && stopFrame <= handles.numFrames)
    
    % Set Program State to Automation
    handles.state = 'automation';
    guidata(hObject,handles);
    
    % Set Status Indicator
    set(handles.statusEdit,'String','Have Faith');
    
    % Set Automation Indicators
    set(handles.automationTypeEdit','String','Have Faith');
    set(handles.automationStartEdit,'String',num2str(startFrame));
    set(handles.automationStopEdit,'String',num2str(stopFrame));
    
    % Save guidata
    guidata(hObject,handles);
    
    % Get guidata to be used in loop
    % Note a local handles copy (h_loop) is used.  The handles structure is
    % updated at the end of this function.  This is done to avoid potential
    % problems with the handles structure updating when the "Stop Auto"
    % button is pressed.
    h_loop = guidata(hObject);
    
    % Run Have Faith %
    h_loop.currFrame = h_loop.currentLabel;
    displayFun(h_loop);
    
    for k = h_loop.currFrame:stopFrame
        h_loop = labelCurr(h_loop,0);
        set(h_loop.currFrameEdit,'String',num2str(h_loop.currFrame));
        
        % Stop Have Faith if user presses Stop Automation %
        % Get Latest Handles
        h_latest = guidata(hObject);
        
        % Check Stop Flag
        if(h_latest.autoStopFlag == 1)
            % Stop Pressed.  Break
            h_loop.autoStopFlag = 0;
            break;
        end    
    end
    
    % Display a Messge Box if Have Faith ran until stopframe
    if(k == stopFrame)
        % Display Message Box Indicating Completion
        msgbox('Have Faith Complete!','Have Faith Complete!','modal');
    end
    
    % Set State to Wait For Label
    h_loop.state = 'waitForLabel';
    
    % This is the structure that will be used to update guidata
    handles = h_loop;
    
    % Update Display
    displayFun(handles);
    
    % Update Automation Indicators
    set(handles.automationTypeEdit,'String','');
    set(handles.automationStartEdit,'String','');
    set(handles.automationStopEdit,'String','');
    
    % Update Status
    set(handles.statusEdit,'String','Waiting for Label');

else
    msgbox('Incorrect End Frame Entered','Incorrect End Frame','error','modal'); 
end
% Update guidata
guidata(hObject,handles);

% --- Executes on button press in multiErrorButton.
function multiErrorButton_Callback(hObject, eventdata, handles)
% hObject    handle to multiErrorButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Only run if template has been initialized
if(~strcmp(handles.state,'waitForLabel'))
    return;
end

startFrame = handles.currentLabel;

% Get End Frame and Error Type from user using modal dialog
multiErrorInfo = multiErrorDlg('StartFrame',startFrame);
stopFrame = multiErrorInfo(1); 
errorType = multiErrorInfo(2);

% Check for user cancel
if(all(multiErrorInfo == 0))
    return;
end

% Check Valid stopFrame
if(stopFrame > handles.currentLabel && stopFrame <= handles.numFrames)
    
    % Note : Multi Error ignores the Stop Auto Button.  This is because it
    % labels all frames but the last without running the detector.
    % Therefore it should not take long to run.
    
    % Set Program State to Automation
    handles.state = 'automation';
    guidata(hObject,handles);
    
    % Set Status Indicator
    set(handles.statusEdit,'String','Multi Error');
    
    % Set Automation Indicators
    set(handles.automationTypeEdit','String','Multi Error');
    set(handles.automationStartEdit,'String',num2str(startFrame));
    set(handles.automationStopEdit,'String',num2str(stopFrame));
    
    % Save guidata
    guidata(hObject,handles);
    
    % Get guidata to be used in loop
    % Note, a local handles copy (h_loop) is used.  The handles structure
    % is updated at the end of this function.  This is done to avoid
    % potential problems with the handles structure updating when the "Stop
    % Auto" button is pressed.
    h_loop = guidata(hObject);
    
    % Run Multi Error %
    h_loop.currFrame = h_loop.currentLabel;
    displayFun(h_loop);
    
    % Label all but stopFrame without running the detector %
    if(stopFrame >= startFrame + 1)
        h_loop.label(startFrame:(stopFrame-1)) = errorType;
        h_loop.currentLabel = stopFrame;
        h_loop.currFrame = stopFrame;
    end
     
    % Label stopFrame with labelCurr so detector is run on next frame
    h_loop = labelCurr(h_loop,errorType);
    set(h_loop.currFrameEdit,'String',num2str(h_loop.currFrame));
        
    % Set State to Wait For Label
    h_loop.state = 'waitForLabel';
    
    % This is the structure that will be used to update guidata
    handles = h_loop;
    
    % Update Display
    displayFun(handles);
    
    % Update Automation Indicators
    set(handles.automationTypeEdit,'String','');
    set(handles.automationStartEdit,'String','');
    set(handles.automationStopEdit,'String','');
    
    % Display Message Box Indicating Completion
    msgbox('Multi Error Complete!','Multi Error Complete!','modal');
    
    % Update Status
    set(handles.statusEdit,'String','Waiting for Label');
%}
    
else
    msgbox('Incorrect End Frame Entered','Incorrect End Frame','error','modal'); 
end
% Update guidata
guidata(hObject,handles);


% --- Executes on button press in viewTrackingDataButton.
function viewTrackingDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to viewTrackingDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in stopAutoButton.
function stopAutoButton_Callback(hObject, eventdata, handles)
% hObject    handle to stopAutoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% If currently in an automation mode, set stop flag in handles
% structure.
if(strcmp(handles.state,'automation'))
    % Update autoStopFlag
    handles.autoStopFlag = 1;
    
    % Save Change to handles
    guidata(hObject,handles);
end

function automationTypeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to automationTypeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of automationTypeEdit as text
%        str2double(get(hObject,'String')) returns contents of automationTypeEdit as a double


% --- Executes during object creation, after setting all properties.
function automationTypeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to automationTypeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function automationStartEdit_Callback(hObject, eventdata, handles)
% hObject    handle to automationStartEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of automationStartEdit as text
%        str2double(get(hObject,'String')) returns contents of automationStartEdit as a double


% --- Executes during object creation, after setting all properties.
function automationStartEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to automationStartEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function automationStopEdit_Callback(hObject, eventdata, handles)
% hObject    handle to automationStopEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of automationStopEdit as text
%        str2double(get(hObject,'String')) returns contents of automationStopEdit as a double


% --- Executes during object creation, after setting all properties.
function automationStopEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to automationStopEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function utilityMenu_Callback(hObject, eventdata, handles)
% hObject    handle to utilityMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function viewTrackingDataMenu_Callback(hObject, eventdata, handles)
% hObject    handle to viewTrackingDataMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~strcmp(handles.state,'waitForLabel'))
    return;
end

viewInstTrackData(handles);
