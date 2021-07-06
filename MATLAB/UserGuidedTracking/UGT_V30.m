function varargout = UGT_V30(varargin)
% UGT_V30 MATLAB code for UGT_V30.fig
%      UGT_V30, by itself, creates a new UGT_V30 or raises the existing
%      singleton*.
%
%      H = UGT_V30 returns the handle to a new UGT_V30 or the handle to
%      the existing singleton*.
%
%      UGT_V30('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UGT_V30.M with the given input arguments.
%
%      UGT_V30('Property','Value',...) creates a new UGT_V30 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UGT_V30_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UGT_V30_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singlet on)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UGT_V30

% Last Modified by GUIDE v2.5 18-Feb-2016 14:02:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UGT_V30_OpeningFcn, ...
                   'gui_OutputFcn',  @UGT_V30_OutputFcn, ...
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


% --- Executes just before UGT_V30 is made visible.
function UGT_V30_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UGT_V30 (see VARARGIN)

% Choose default command line output for UGT_V30
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
% Setup Auto-Pilot
% threshVec = [delta_NCC,delta_width,delta_orient,delta_inliersDistL,...
%              delta_inliersDistR,inliersL,inliersR]
% Threshold Parameters Set 2011-11-25
handles.threshVec = [0.0332,1.72,0.00987,10,10,10,10];
handles.threshVec_ini = [0.0332,1.72,0.00987,10,10,10,10];
handles.threshVec_org = [0.0332,1.72,0.00987,10,10,10,10];
handles.typeVec = [0,0,0,0,0,1,1];

handles = iniProgressBar(hObject,handles);

set(handles.checkboxCeasing,'Value',0);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UGT_V30 wait for user response (see UIRESUME)
% uiwait(handles.UGT_V30);


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
           
function [result,fVec] = checkAlgoConf(handles,frameNum)
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

% paint the bad element in red
paintAbnormalEdit(handles,result);
% display the  Confidence Parameters to the bottom edit line
bTop = 0;
displayThreshholdVecToEdits(handles,bTop,2,fVec);

% --- Outputs from this function are returned to the command line.
function varargout = UGT_V30_OutputFcn(hObject, eventdata, handles) 
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

if ~isfield(handles,'batchRun')
    msgbox('No task !','Error');
    return;
end

% if there is not task, return
taskSum = handles.batchRun.totalTaskNum;
if taskSum == 0
    msgbox('No task !','Error');
    return;
end
% get current task
videoTaskList = handles.batchRun.videoTaskList;
videoTask = videoTaskList{handles.currentTaskNo};

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
%     % Clear Analysis File Name
%     handles.analysisFileName = [];
%     set(handles.analysisFileEdit,'String','');
    
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
h_weightWindow = figure('Name','Define Interest Region');
imshow(temp);
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

% Record the template to videoTask
videoTask.pc = handles.pc;
videoTask.mark = handles.mark;
videoTask.inst = handles.inst;
videoTask.algoInfo = handles.algoInfo;
videoTask.templateFrame = handles.templateFrame;

videoTaskList{handles.currentTaskNo} = videoTask;
handles.batchRun.videoTaskList = videoTaskList;

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
       if handles.currFrame ~= handles.numFrames
           handles = performTracking(handles,handles.currFrame + 1);
           handles.currFrame = handles.currFrame + 1;
       end       
       handles.currentLabel = handles.currFrame;
       %displayFun(handles);
  else
      % Label & Clear
%       btn = questdlg('This frame has already been labeled.  By labeling it, all labels set after it will be lost.  Are you sure?',...
%           'Change already set label?','No');
        btn = 'Yes';
        if(strcmp(btn,'Yes'))
            handles.label(handles.currFrame) = labelNum;
            handles = clearAfter(handles,handles.currFrame);
            if handles.currFrame ~= handles.numFrames
                handles = performTracking(handles,handles.currFrame + 1);
                handles.currFrame = handles.currFrame + 1;
            end
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
handles.vidObj = VideoReader(handles.vidName);
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

% --- Executes on key press with focus on UGT_V30 or any of its controls.
function UGT_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to UGT_V30 (see GCBO)
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



function batchTxtFileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to batchTxtFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of batchTxtFileEdit as text
%        str2double(get(hObject,'String')) returns contents of batchTxtFileEdit as a double


% --- Executes during object creation, after setting all properties.
function batchTxtFileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to batchTxtFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in batchTxtFileSelBtn.
function batchTxtFileSelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to batchTxtFileSelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This Button select batch txt files
[fName,pName] = uigetfile({'*.txt','(*.txt) Batch Text Files'}...
            ,'Load left Threshold Vector');
% Check for User Cancel
if(fName == 0)
    return;
end

handles.batchRun.currentTaskNo = 0;
set(handles.currentTaskEdit,'String',num2str(handles.batchRun.currentTaskNo));
handles.batchRun.totalTaskNum = 0;
set(handles.currentTaskEdit,'String',num2str(handles.batchRun.totalTaskNum));

batchTxtFileName = fullfile(pName,fName);

% read in the data in the txt file
[cellFileContent,numFileLine] = readTxtFile(batchTxtFileName);
set(handles.batchTxtFileEdit,'String',batchTxtFileName);

% begin to check if the task file is valid.
nmax = length(cellFileContent);
validFileSum = 0;
for n = 1 : nmax
    listFilename = cellFileContent{n};
    if ( isa(listFilename,'char')== 0)
        continue;
    end
    [pName,fName,extName] = fileparts(listFilename);
    extName = lower(extName);
    if (strcmp(extName,'.avi')~=1)
        continue;
    end
    
    hFile = fopen(listFilename,'r');
    % if exist(taskFilename,'file')~=0
    if hFile ~= 0 && hFile ~= -1
        validFileSum = validFileSum + 1;
        validFileList{validFileSum} = listFilename;
        fclose (hFile);
    end    
end

nmax = floor(validFileSum / 2);
if nmax == 0
    msgbox('No valid task!','No task;')
    set(handles.totalTaskEdit,'String','0');
    set(handles.currentTaskEdit,'String','0');
    return;
end

taskSum = nmax;
for n = 1 : nmax
    videoTask.vidName = validFileList{2*n-1};
    % take the 2nd file name as the background video file;
    videoTask.backName = validFileList{2*n}; 
    displayTaskList{n} = videoTask.vidName;
    
    [pName,fName,extName] = fileparts(videoTask.vidName);
    extName = '.mat';
    videoTask.analysisFileName = fullfile(pName,[fName '_analysis_task' num2str(n) extName]);        
    
    set(handles.outputFileEdit,'String',videoTask.analysisFileName);
    
    videoTaskList{n} = videoTask;
end
% set display
set(handles.taskListbox,'String',displayTaskList);
set(handles.totalTaskEdit,'String',num2str(taskSum));
set(handles.currentTaskEdit,'String','0');

% save data
handles.batchRun.videoTaskList = videoTaskList;
handles.batchRun.totalTaskNum = taskSum;
handles.batchRun.currentTaskNo = 1;

handles.currentTaskNo = 1;
[handles] = loadVedioTaskData(hObject, eventdata, handles, 1);
guidata(hObject,handles);




% --- Executes on selection change in taskListbox.
function taskListbox_Callback(hObject, eventdata, handles)
% hObject    handle to taskListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns taskListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from taskListbox
idxClickSel = get(hObject,'value');
handles.currentTaskNo = idxClickSel;
set(handles.currentTaskEdit,'String',num2str(idxClickSel));
[handles] = loadVedioTaskData(hObject, eventdata, handles, idxClickSel);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function taskListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to taskListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function currentTaskEdit_Callback(hObject, eventdata, handles)
% hObject    handle to currentTaskEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentTaskEdit as text
%        str2double(get(hObject,'String')) returns contents of currentTaskEdit as a double


% --- Executes during object creation, after setting all properties.
function currentTaskEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentTaskEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function totalTaskEdit_Callback(hObject, eventdata, handles)
% hObject    handle to totalTaskEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of totalTaskEdit as text
%        str2double(get(hObject,'String')) returns contents of totalTaskEdit as a double


% --- Executes during object creation, after setting all properties.
function totalTaskEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to totalTaskEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkReadyBtn.
function checkReadyBtn_Callback(hObject, eventdata, handles)
% hObject    handle to checkReadyBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This function build the task with all the input data needed
% including:
% 1. background picture.
% 2. template data
% 3. threshold setting
if ~isfield(handles,'batchRun')
    return;
end

taskSum = handles.batchRun.totalTaskNum;
totalFrameNum = 0;
bHasTemplete = zeros(1,taskSum);
bHasThreshVec = zeros(1,taskSum);
for n = 1 : taskSum
    handles.batchRun.currentTaskNo = n;
    set(handles.currentTaskEdit,'String',num2str(n));
       
  % Load Video & Store Video Parameters
    handles = loadVedioTaskData(hObject, eventdata, handles, n); 
    videoTaskList = handles.batchRun.videoTaskList;
    videoTask = videoTaskList{n};
    % Calculate total frame number
    totalFrameNum = totalFrameNum + handles.numFrames
   
    % Clear Data Structures
    % handles = initNonVideoStructures(handles);
    % guidata(hObject,handles);
    displayFun(handles);
    
    % currFrame = handles.currFrame;

    % Check if there is a template
    if(isfield(videoTask,'mark'))
        bHasTemplete(n) = 1; 
    end
    
    if(isfield(videoTask,'theThresholdVector'))
        bHasThreshVec(n) = 1;
    end
%     % Begin to prepare template
%     % Interactively Get Template Window
%     defStartSize = 99; % Default Window Starting Size
%     frameIm = rgb2gray(handles.vidObj.read(handles.currFrame));
%     [centerPt,tSize] = squareDraw(frameIm,defStartSize);
% 
%     % Check for user cancel
%     if(isempty(centerPt))
%         errordlg('Cancelled Template Initializiation Process','Cancelled Init','modal');
%         return;
%     end
% 
%     % Upper Left Corner of Template
%     tempWindow = centerPt - ((tSize-1)/2);
% 
%     nSize = round(tSize * handles.neighborhoodSize)-1;
%     wStats = struct('wSize',tSize,'nx',nSize+1,'ny',nSize+1);
%     nCurr = tempWindow - repmat(round((nSize - tSize)/2),1,2);
%     temp = frameIm((0:tSize-1)+tempWindow(2),(0:tSize-1)+tempWindow(1),:);
% 
%     % Interactively get Template Weighting Mask
%     h_weightWindow = figure('Name','Define Interest Region');imshow(temp);
%     h_poly = impoly(gca,'Closed',1);
%     wait(h_poly);
% 
%     % If User closed window do not continue template initialization process
%     if(isempty(h_poly))
%         if(ishghandle(h_weightWindow))
%             close(h_weightWindow);
%         end
%         errordlg('Cancelled Template Initializiation Process','Cancelled Init','modal');
%         return;
%     end
% 
%     % Store Weighting Mask and Close Window
%     weights = double(h_poly.createMask());
%     close(h_weightWindow);
% 
%     % Template Related Data Structures
%     handles.pc = wMatchBuildPcStruct(double(temp),weights);
%     handles.mark.tCorn(currFrame,:) = tempWindow;
%     handles.mark.subT(currFrame,:) = tempWindow;
%     handles.mark.wStats = wStats;
%     handles.mark.nCorn(currFrame,:) = nCurr;
%     handles.templateFrame = currFrame;
%     
%     % Copy to videoTask
%     videoTask.pc = handles.pc;
%     videoTask.mark = handles.mark;
%     % templateFrame = handles.templateFrame;
% 
%     % Template Related Algo Info
%     handles.algoInfo.nccScore(currFrame,1) = 1;
%     handles.algoInfo.fitMSE(currFrame,1) = 0;
% 
%     % Estimate Boundary Lines
%     [r,t,dbl_ais] = detection_bLines(frameIm,handles.backIm,tempWindow,wStats,handles.pc);
%     handles.inst.rho(currFrame,:) = r;
%     handles.inst.theta(currFrame,:) = t;
%     handles.inst.trackPt(currFrame,:) = computeTrackPt(r,t,tempWindow,...
%                                                         round((tSize-1)/2));
%     % Boundary Line Algo Info
%     handles.algoInfo.votesLeft(currFrame,:) = dbl_ais.votesLeft;
%     handles.algoInfo.votesRight(currFrame,:) = dbl_ais.votesRight;
%     handles.algoInfo.leftNumInliers(currFrame,:) = dbl_ais.leftNumInliers;
%     handles.algoInfo.leftRefitMSE(currFrame,:) = dbl_ais.leftRefitMSE;
%     handles.algoInfo.rightNumInliers(currFrame,:) = dbl_ais.rightNumInliers;
%     handles.algoInfo.rightRefitMSE(currFrame,:) = dbl_ais.rightRefitMSE;
%     
%     % Set State to waitForLabel
%     handles.currentLabel = currFrame;
%     handles.state = 'waitForLabel';
%     set(handles.statusEdit,'String','Waiting for Label');
%     
%     % 2nd Copy to videoTask
%     videoTask.inst = handles.inst;
%     videoTask.algoInfo = handles.algoInfo;
%     videoTask.templateFrame = handles.templateFrame;
%     
%     displayFun(handles);
%     guidata(hObject,handles);
%     % Finish preparing template    
%     
%     % strMsg = ['Load Background File for Task ' num2str(n)];
% 
%     videoTaskList{n} = videoTask;
%     handles.batchRun.videoTaskList = videoTaskList;
end
% template parameter
    % handles.inst
    % handles.algoInfo
    % handles.pc
    % handles.mark
    % handles.templateFrame
% Update elements
handles.batchRun.totalFrameNum = totalFrameNum;
% handles.batchRun.videoTaskList = videoTaskList;
% Save data
guidata(hObject,handles);

bTempleteOK = 0;
strNum = '';
if sum(bHasTemplete)==taskSum
    bTempleteOK = 1;    
else
    bTempleteOK = 0;
    strNum = '';
    for n = 1 : taskSum
        if bHasTemplete(n)==0
            strNum = [strNum num2str(n) ', '];
        end
    end
    if length(strNum)>=2
        strNum(end-1:end)=[]; 
    end
end

bThreshVecOK = 0;
strNum1 = '';
if sum(bHasThreshVec)==taskSum
    bThreshVecOK = 1;
else
    bThreshVecOK = 0;
    strNum1 = '';
    for n = 1 : taskSum
        if bHasThreshVec(n)==0
            strNum1 = [strNum1 num2str(n) ', '];
        end
    end
    if length(strNum1)>=2
        strNum1(end-1:end)=[];
    end
end

if bTempleteOK == 1 && bThreshVecOK==1
    msgbox('All tasks are ready!','Ready');
    return;
elseif bTempleteOK == 1 && bThreshVecOK==0
    strmsg = {'Templete setting is OK.'...
        ,['No Threshold Vector is set for No. ' strNum1 ' task ! Please initialize!']};
elseif bTempleteOK == 0 && bThreshVecOK==1 
    strmsg = {['No templete is set for No. ' strNum ' task ! Please initialize!']...
        ,'Thresh vector setting is OK! '};
else
    strmsg = {['No templete is set for No. ' strNum ' task ! Please initialize!']...
        ,['No Threshold Vector is set for No. ' strNum1 ' task ! Please initialize!']};
    
end

msgbox(strmsg,'Not Ready');


%--------------------------------------------------------------------------
function usedThreshVecEdit01_Callback(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of usedThreshVecEdit01 as text
%        str2double(get(hObject,'String')) returns contents of usedThreshVecEdit01 as a double


% --- Executes during object creation, after setting all properties.
function usedThreshVecEdit01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function usedThreshVecEdit02_Callback(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of usedThreshVecEdit02 as text
%        str2double(get(hObject,'String')) returns contents of usedThreshVecEdit02 as a double


% --- Executes during object creation, after setting all properties.
function usedThreshVecEdit02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function usedThreshVecEdit03_Callback(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of usedThreshVecEdit03 as text
%        str2double(get(hObject,'String')) returns contents of usedThreshVecEdit03 as a double


% --- Executes during object creation, after setting all properties.
function usedThreshVecEdit03_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function usedThreshVecEdit04_Callback(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit04 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of usedThreshVecEdit04 as text
%        str2double(get(hObject,'String')) returns contents of usedThreshVecEdit04 as a double


% --- Executes during object creation, after setting all properties.
function usedThreshVecEdit04_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit04 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function usedThreshVecEdit05_Callback(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit05 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of usedThreshVecEdit05 as text
%        str2double(get(hObject,'String')) returns contents of usedThreshVecEdit05 as a double


% --- Executes during object creation, after setting all properties.
function usedThreshVecEdit05_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit05 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function usedThreshVecEdit06_Callback(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit06 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of usedThreshVecEdit06 as text
%        str2double(get(hObject,'String')) returns contents of usedThreshVecEdit06 as a double


% --- Executes during object creation, after setting all properties.
function usedThreshVecEdit06_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit06 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function usedThreshVecEdit07_Callback(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit07 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of usedThreshVecEdit07 as text
%        str2double(get(hObject,'String')) returns contents of usedThreshVecEdit07 as a double


% --- Executes during object creation, after setting all properties.
function usedThreshVecEdit07_CreateFcn(hObject, eventdata, handles)
% hObject    handle to usedThreshVecEdit07 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function runingParaEdit01_Callback(hObject, eventdata, handles)
% hObject    handle to runingParaEdit01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of runingParaEdit01 as text
%        str2double(get(hObject,'String')) returns contents of runingParaEdit01 as a double


% --- Executes during object creation, after setting all properties.
function runingParaEdit01_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runingParaEdit01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function runingParaEdit02_Callback(hObject, eventdata, handles)
% hObject    handle to runingParaEdit02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of runingParaEdit02 as text
%        str2double(get(hObject,'String')) returns contents of runingParaEdit02 as a double


% --- Executes during object creation, after setting all properties.
function runingParaEdit02_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runingParaEdit02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function runingParaEdit03_Callback(hObject, eventdata, handles)
% hObject    handle to runingParaEdit03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of runingParaEdit03 as text
%        str2double(get(hObject,'String')) returns contents of runingParaEdit03 as a double


% --- Executes during object creation, after setting all properties.
function runingParaEdit03_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runingParaEdit03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function runingParaEdit04_Callback(hObject, eventdata, handles)
% hObject    handle to runingParaEdit04 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of runingParaEdit04 as text
%        str2double(get(hObject,'String')) returns contents of runingParaEdit04 as a double


% --- Executes during object creation, after setting all properties.
function runingParaEdit04_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runingParaEdit04 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function runingParaEdit05_Callback(hObject, eventdata, handles)
% hObject    handle to runingParaEdit05 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of runingParaEdit05 as text
%        str2double(get(hObject,'String')) returns contents of runingParaEdit05 as a double


% --- Executes during object creation, after setting all properties.
function runingParaEdit05_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runingParaEdit05 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function runingParaEdit06_Callback(hObject, eventdata, handles)
% hObject    handle to runingParaEdit06 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of runingParaEdit06 as text
%        str2double(get(hObject,'String')) returns contents of runingParaEdit06 as a double


% --- Executes during object creation, after setting all properties.
function runingParaEdit06_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runingParaEdit06 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function runingParaEdit07_Callback(hObject, eventdata, handles)
% hObject    handle to runingParaEdit07 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of runingParaEdit07 as text
%        str2double(get(hObject,'String')) returns contents of runingParaEdit07 as a double


% --- Executes during object creation, after setting all properties.
function runingParaEdit07_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runingParaEdit07 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in setThreshVecBtn01.
function setThreshVecBtn01_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setRuningEditToUsedEdit(handles,1);

% --- Executes on button press in setThreshVecBtn02.
function setThreshVecBtn02_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setRuningEditToUsedEdit(handles,2);

% --- Executes on button press in setThreshVecBtn03.
function setThreshVecBtn03_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setRuningEditToUsedEdit(handles,3);

% --- Executes on button press in setThreshVecBtn04.
function setThreshVecBtn04_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn04 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setRuningEditToUsedEdit(handles,4);

% --- Executes on button press in setThreshVecBtn05.
function setThreshVecBtn05_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn05 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setRuningEditToUsedEdit(handles,5);

% --- Executes on button press in setThreshVecBtn06.
function setThreshVecBtn06_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn06 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setRuningEditToUsedEdit(handles,6);

% --- Executes on button press in setThreshVecBtn07.
function setThreshVecBtn07_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn07 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setRuningEditToUsedEdit(handles,7);

% --- Executes on button press in applyThresholdToVedioBtn.
function applyThresholdToVedioBtn_Callback(hObject, eventdata, handles)
% hObject    handle to applyThresholdToVedioBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
theVec = updateThreshVecFromEdits(handles);

handles.threshVec = theVec;
handles.threshVec_ini = handles.threshVec;

guidata(hObject,handles);

bTop = 1; nVecIdx = 0;
displayThreshholdVecToEdits(handles,bTop,nVecIdx) ;
bTop = 0; nVecIdx = 1;
displayThreshholdVecToEdits(handles,bTop,nVecIdx) ;

%--------------------------------------------------------------------------
% --- Executes on button press in saveThresholdToFileBtn.
function saveThresholdToFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to saveThresholdToFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currVidDir = getCurrVidDir(handles);
if(currVidDir ~= 0)
    % Open File Browser in current video directory
    [fName,pName] = uiputfile({'*.mat','(*.mat) Mat Files'; '*.txt','(*.txt) Text Files'}...
        ,'Save Left Threshold Vector',currVidDir);
else
    [fName,pName] = uiputfile({'*.mat','(*.mat) Mat Files'; '*.txt','(*.txt) Text Files'}...
        ,'Save Left Threshold Vector');
end

% Check for User Cancel
if(fName == 0)
    return;
end
thresholdVectorFileName = fullfile(pName,fName);
saveVec = handles.threshVec_ini;

[pName,fName,extName] = fileparts(thresholdVectorFileName);
extName = lower(extName);
if (strcmp(extName,'.mat')==1)
    % mat file save
    save(thresholdVectorFileName, 'saveVec');
else 
    %txt file save
    fID = fopen(thresholdVectorFileName,'wt');
    fprintf(fID,'%f\t',saveVec);
    fclose(fID);        
end    

%--------------------------------------------------------------------------
% --- Executes on button press in loadThresholdFromFileBtn.
function loadThresholdFromFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to loadThresholdFromFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currVidDir = getCurrVidDir(handles);
if(currVidDir ~= 0)
    % Open File Browser in current video directory
    [fName,pName] = uigetfile({'*.mat','(*.mat) Mat Files'; '*.txt','(*.txt) Text Files'}...
        ,'Load left Threshold Vector',currVidDir);
else
    [fName,pName] = uigetfile({'*.mat','(*.mat) Mat Files'; '*.txt','(*.txt) Text Files'}...
        ,'Load left Threshold Vector');
end

% Check for User Cancel
if(fName == 0)
    return;
end
thresholdVectorFileName = fullfile(pName,fName);    
    
[pName,fName,extName] = fileparts(thresholdVectorFileName);
extName = lower(extName);
if (strcmp(extName,'.mat')==1)
    % mat file load
    clear saveVec;
    load(thresholdVectorFileName);
    if ~exist('saveVec')
        errordlg('Wrong file is loaded. Please reload all data files !!','File loading error'); 
        return;
    end
else 
    %txt file load
    saveVec = load(thresholdVectorFileName); 
    [ro,co] = size(saveVec);
    if ro~= 1 || co~= 7
        errordlg('Vector size should be 1*7','File loading error'); 
        return;
    end
end
handles.threshVec_ini = saveVec;
bLeft = 1; bTop = 1; nVecIdx = 0;
displayThreshholdVecToEdits(handles,bTop,nVecIdx) ;
bTop = 0;
displayThreshholdVecToEdits(handles,bTop,nVecIdx) ;

guidata(hObject,handles);
%--------------------------------------------------------------------------
function setRuningEditToUsedEdit(handles,idx)

switch idx
    case 1
        theStr = get(handles.runingParaEdit01,'String');
        if isempty(str2num(theStr))
            warning('Not a number.');
            return;
        end
        set(handles.usedThreshVecEdit01,'String',theStr);
    case 2
        theStr = get(handles.runingParaEdit02,'String');
        if isempty(str2num(theStr))
            warning('Not a number.');
            return;
        end
        set(handles.usedThreshVecEdit02,'String',theStr);
    case 3
        theStr = get(handles.runingParaEdit03,'String');
        if isempty(str2num(theStr))
            warning('Not a number.');
            return;
        end
        set(handles.usedThreshVecEdit03,'String',theStr);
    case 4
        theStr = get(handles.runingParaEdit04,'String');
        if isempty(str2num(theStr))
            warning('Not a number.');
            return;
        end
        set(handles.usedThreshVecEdit04,'String',theStr);
    case 5
        theStr = get(handles.runingParaEdit05,'String');
        if isempty(str2num(theStr))
            warning('Not a number.');
            return;
        end
        set(handles.usedThreshVecEdit05,'String',theStr);
    case 6
        theStr = get(handles.runingParaEdit06,'String');
        if isempty(str2num(theStr))
            warning('Not a number.');
            return;
        end
        set(handles.usedThreshVecEdit06,'String',theStr);
    case 7
        theStr = get(handles.runingParaEdit07,'String');
        if isempty(str2num(theStr))
            warning('Not a number.');
            return;
        end
        set(handles.usedThreshVecEdit07,'String',theStr);
    otherwise
        warning('Not valid value.');
end
%--------------------------------------------------------------------------
% This function execute 1 video analysis task
function [handles] = executeTask(hObject, eventdata, handles, taskNo)
% hObject    handle to loadThresholdFromFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% taskNo     input task No. to locate task index

if ~isfield(handles,'batchRun')
    msgbox('No task !','Error');
    return;
end

taskSum = handles.batchRun.totalTaskNum;
if taskSum == 0
    msgbox('No task !','Error');
    return;
end
videoTaskList = handles.batchRun.videoTaskList;
% totalFrameNum = handles.batchRun.totalFrameNum;

theLabelCorrect = 0;
theLabelError = 4;

if taskNo > taskSum
    return;
end

% get the videoTask
videoTask = videoTaskList{taskNo};

startFrame = 1;
stopFrame = videoTask.numFrames ;
% stopFrame = 30;
% startFrame = 11000;
% if stopFrame > 500
%     stopFrame = 500;
% end
% stopFrame = 15;
% handles.batchRun.totalFrameNum = 30;

% start of progress bar setting
setEvaluateNum(0,videoTask.numFrames,...
                handles.hPatchCurrentTaskPercentage,...
                handles.textPercentageCurrTask);
            
setEvaluateNum(handles.batchRun.nCompletedFrame,handles.batchRun.totalFrameNum,...
                handles.hPatchTotalTaskPercentage,...
                handles.textPercentageTotalTask);
% end of progress bar setting 

% Now Load data from vedioTask
handles = loadVedioTaskData(hObject, eventdata, handles, taskNo);

% Load template data
handles.mark = videoTask.mark;
handles.inst = videoTask.inst;
handles.algoInfo = videoTask.algoInfo;
handles.pc = videoTask.pc;

handles.threshVec_ini = videoTask.theThresholdVector;
% Save guidata
guidata(hObject,handles);

% Added on Jan 25 2016, by Jonathan
% This part od program set the checking range of a task. If rangeInfo has
% been set, use the set info to findout the index and record to variance 
% 'checkingRange'. If there is not rangeInfo, 'checkingRange' will be set
% by the whole video frame index;
checkingRange = [];
% Orgnize range info, only being done when videoTask has rangeInfo
if (isfield(videoTask,'rangeInfo'))
    [rangeInfoRow,rangeInfoCol] = size(videoTask.rangeInfo);
    for n = 1 :rangeInfoRow
        infoRow = videoTask.rangeInfo(n,:)
        checkingRange = [checkingRange,infoRow(2):infoRow(3)];
    end
else
    checkingRange = [startFrame:stopFrame];
end
% reset the startFrame to first frame of the checkingRange
startFrame = checkingRange(1);
% reset the stopFrame
stopFrame = checkingRange(end);
% End adding on Jan 25 2016,


if(~strcmp(handles.state,'waitForLabel'))
    return;
end

% Set Program State to Automation
% handles.state = 'automation';

% Set Status Indicator
set(handles.statusEdit,'String','Auto Pilot');
handles.prevAutoPilotStop = stopFrame;

% Set Automation Indicators
set(handles.automationTypeEdit','String','Auto Pilot');

% Move to Frame No. startFrame
handles.currFrame = startFrame;
handles.currentLabel = handles.currFrame;
set(handles.currFrameEdit,'String',num2str(startFrame));
displayFun(handles);

% Label the First Frame Correct
% if ~any(checkingRange(:)==1)
%     handles = labelCurr(handles,theLabelError);
% else
%     handles = labelCurr(handles,theLabelCorrect);
% end
handles = labelCurr(handles,theLabelCorrect);

setEvaluateNum(1,videoTask.numFrames,...
                handles.hPatchCurrentTaskPercentage,...
                handles.textPercentageCurrTask);
            
handles.batchRun.nCompletedFrame = handles.batchRun.nCompletedFrame + 1;
setEvaluateNum(handles.batchRun.nCompletedFrame,handles.batchRun.totalFrameNum,...
                handles.hPatchTotalTaskPercentage,...
                handles.textPercentageTotalTask);
            
set(handles.currFrameEdit,'String',num2str(handles.currFrame));
set(handles.currFrameEdit,'String',num2str(startFrame));
% Save guidata
guidata(hObject,handles);

% Label all frames which are out of checking range error;
% find out 'outRange'
outRange = [1:videoTask.numFrames];
outRange(checkingRange) = [];
% Label them Error
if ~isempty(outRange)
    handles.label(outRange) = theLabelError;
    
    % start of progressbar setting
    setEvaluateNum(1+length(outRange),videoTask.numFrames,...
                handles.hPatchCurrentTaskPercentage,...
                handles.textPercentageCurrTask);
            
    handles.batchRun.nCompletedFrame = handles.batchRun.nCompletedFrame + length(outRange);
    setEvaluateNum(handles.batchRun.nCompletedFrame,handles.batchRun.totalFrameNum,...
                    handles.hPatchTotalTaskPercentage,...
                    handles.textPercentageTotalTask);
    % end of progressbar setting

    set(handles.currFrameEdit,'String',num2str(handles.currFrame));
    % Save guidata
    guidata(hObject,handles);
end

% Get guidata to be used in loop
% Note a local handles copy (h_loop) is used.  The handles structure is
% updated at the end of this function.  This is done to avoid potential
% problems with the handles structure updating when the "Stop Auto"
% button is pressed.
    
% set 'autoStopFlag' to zero 
handles.autoStopFlag = 0;

h_loop = guidata(hObject);   
    
bTop = 1; nVecIdx = 0;
displayThreshholdVecToEdits(handles,bTop,nVecIdx);
handles.threshVec = handles.threshVec_ini;
  
% Auto-Pilot Loop
kmax = length(checkingRange);
for k = 2 : kmax    
    % Auto pilot stop check   
    if(handles.autoStopFlag == 1)
         set(handles.automationStartEdit,'String',num2str(checkingRange(k)));
         set(handles.currFrameEdit,'String',num2str(checkingRange(k)));
         set(handles.statusEdit,'String','waitForLabel');
         handles.autoStopFlag = 0;
         break;
    end
    
    if get(handles.checkboxCeasing,'Value')~=0
        disp('Procesure is ceasing!');
        return;
    end
    
    % Added on Jan 25 2016, by Jonathan
    % Check whether 'k' is in errRange
%     if ~any(checkingRange(:)==k)
%         displayFun(h_loop);
%         % Mark frame error
%         h_loop = labelCurr(h_loop,theLabelError);
%         % jump to next frame directly
%         continue;
%     end
    % End adding on Jan 25 2016,
    
    % if multiple range is met, move to the startFrame and label it correct
    if (checkingRange(k)-checkingRange(k-1)~=1)
        h_loop.currFrame = checkingRange(k);
        h_loop.currentLabel = h_loop.currFrame;
        displayFun(h_loop);
        h_loop = labelCurr(h_loop,theLabelCorrect);
        continue;
    end
    % Check Algorithm Confidence Parameters
    [result,fVec] = checkAlgoConf(h_loop,checkingRange(k));

    if(any(result))
        % Algorithm is Underconfident %
        % Display feature vector thresholding in command window
        % disp(result);
        % Update Display with current frame
        displayFun(h_loop);

        % Mark frame error
        h_loop = labelCurr(h_loop,theLabelError);
    else   
        displayFun(h_loop);
        try 
            h_loop = labelCurr(h_loop,theLabelCorrect);
        catch e
            h_loop = labelCurr(h_loop,theLabelError);
        end
    end
        
    bTop = 0; nVecIdx = 1;
    displayThreshholdVecToEdits(h_loop,bTop,nVecIdx);

    % Algorithm is Confident %
    % Label frame as correct

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
    
    setEvaluateNum(k + length(outRange),videoTask.numFrames,...
                handles.hPatchCurrentTaskPercentage,...
                handles.textPercentageCurrTask);
    handles.batchRun.nCompletedFrame = handles.batchRun.nCompletedFrame + 1;
    setEvaluateNum(handles.batchRun.nCompletedFrame,handles.batchRun.totalFrameNum,...
                handles.hPatchTotalTaskPercentage,...
                handles.textPercentageTotalTask);
    guidata(hObject,handles);
end % end of the loop for k

% Display a Messge Box if auto-pilot ran until stopframe
% if(k ~= stopFrame)
%     % Display Message Box Indicating Completion
%     msgbox('Video Auto-Pilot is not Completed!','Auto-Pilot not Complete!','modal');
%     % Move to Frame No. stopFrame
%     h_loop.currFrame = stopFrame + 1;
% end

% Set State to Wait For Label
h_loop.state = 'waitForLabel'; 

if k == stopFrame
    h_loop.currFrame = h_loop.currFrame - 1;
end
set(h_loop.currFrameEdit,'String',num2str(h_loop.currFrame));

handles = h_loop;

% Update Display
displayFun(handles);
displayThreshholdVecToEdits(handles,bTop,nVecIdx);

% Update Status
set(handles.statusEdit,'String','Waiting for Label');
handles.currFrame = 1;

guidata(hObject,handles);

% Begin to save data
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
set(handles.outputFileEdit,'String',handles.analysisFileName);
guidata(hObject,handles);

display('The video analysis data has been saved');
    

%--------------------------------------------------------------------------
function displayThreshholdVecToEdits(handles,bTop,nVecIdx,vecIn)
% handles    structure with handles and user data (see GUIDATA)
% bTop  :  0: display the Threshold vector to bottom edit line; 
%          1: display the Threshold vector to top edit line;
% nVecIdx: 0: display threshVec_ini 
%          1: display threshVec 
%          2: display input vec
% vecIn :  the input vector to display
if nargin<3
    warning('Threshold edits display error. Input is not enough.');
    return;
end
if nargin==3
    if nVecIdx > 1
        warning('Threshold edits display error. Input is not enough.');
        return;
    end
end

if nVecIdx==0
    theVec = handles.threshVec_ini;
end

if nVecIdx==1 
    theVec = handles.threshVec;
end

if nVecIdx==2
    theVec = vecIn;
end

if bTop ~= 0
    set(handles.usedThreshVecEdit01,'String',num2str(theVec(1)));
    set(handles.usedThreshVecEdit02,'String',num2str(theVec(2)));
    set(handles.usedThreshVecEdit03,'String',num2str(theVec(3)));
    set(handles.usedThreshVecEdit04,'String',num2str(theVec(4)));
    set(handles.usedThreshVecEdit05,'String',num2str(theVec(5)));
    set(handles.usedThreshVecEdit06,'String',num2str(theVec(6)));
    set(handles.usedThreshVecEdit07,'String',num2str(theVec(7)));
else
    set(handles.runingParaEdit01,'String',num2str(theVec(1)));
    set(handles.runingParaEdit02,'String',num2str(theVec(2)));
    set(handles.runingParaEdit03,'String',num2str(theVec(3)));
    set(handles.runingParaEdit04,'String',num2str(theVec(4)));
    set(handles.runingParaEdit05,'String',num2str(theVec(5)));
    set(handles.runingParaEdit06,'String',num2str(theVec(6)));
    set(handles.runingParaEdit07,'String',num2str(theVec(7)));
end

%--------------------------------------------------------------------------
% --- Executes on button press in executeBtn.
function executeBtn_Callback(hObject, eventdata, handles)
% hObject    handle to executeBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
taskSum = handles.batchRun.totalTaskNum;
if taskSum < 1
    return;
end

% if ~isfield(handles.batchRun,'totalFrameNum')
    totalFrameNum = 0;
    for n = 1 : taskSum
        nFrames = handles.batchRun.videoTaskList{n}.numFrames;
        totalFrameNum = totalFrameNum + nFrames;
    end
    handles.batchRun.totalFrameNum = totalFrameNum;
    guidata(hObject,handles);
% end

handles.batchRun.nCompletedFrame = 0;

set(handles.checkboxCeasing,'Value',0);
% execute task one by one
for n = 1 : taskSum    
    t1 = clock;
    handles.batchRun.idxCurrentTask = n;
    set(handles.currentTaskEdit,'String',num2str(handles.batchRun.idxCurrentTask));       
    
    % reset the labels;
    handles.label = nan*ones(nFrames,1);
    
    try
        % use function 'executeTask' to execute single task
        executeTask(hObject, eventdata, handles, n);
    catch e
        display(['Task No.' ,num2str(n),' error!']);
        continue;
    end
    
    if get(handles.checkboxCeasing,'Value')~=0
        disp('Procesure was ceased by user!');
        disp([num2str(n-1) ' task(s) had been done.' ]);
        set(handles.checkboxCeasing,'Value',0);
        break;
    end
    
    handles.batchRun.nCompletedFrame = handles.batchRun.nCompletedFrame + nFrames;
    t2 = clock;
    durationInSecond = etime(t2,t1);
    msgStr = ['No.' ,num2str(n), ' task spent ' ,num2str(durationInSecond) ,' seconds.\n Every frame needs ',num2str(durationInSecond/nFrames), ' seconds.\n'];
    display(msgStr);
                  
    guidata(hObject,handles);
end
guidata(hObject,handles);

%--------------------------------------------------------------------------
function paintAbnormalEdit(handles,vecRedOn1)
theRed = [1,0,0];
theBlack = [0,0,0];
hEdit(1) = handles.runingParaEdit01;
hEdit(2) = handles.runingParaEdit02;
hEdit(3) = handles.runingParaEdit03;
hEdit(4) = handles.runingParaEdit04;
hEdit(5) = handles.runingParaEdit05;
hEdit(6) = handles.runingParaEdit06;
hEdit(7) = handles.runingParaEdit07;
for n=1:7
    if vecRedOn1(n)==1
        set(hEdit(n),'ForegroundColor',theRed);
    else
        set(hEdit(n),'ForegroundColor',theBlack);
    end
end
%--------------------------------------------------------------------------
function [theVec] = updateThreshVecFromEdits(handles)

theVec(1) = str2num(get(handles.usedThreshVecEdit01,'String'));
theVec(2) = str2num(get(handles.usedThreshVecEdit02,'String'));
theVec(3) = str2num(get(handles.usedThreshVecEdit03,'String'));
theVec(4) = str2num(get(handles.usedThreshVecEdit04,'String'));
theVec(5) = str2num(get(handles.usedThreshVecEdit05,'String'));
theVec(6) = str2num(get(handles.usedThreshVecEdit06,'String'));
theVec(7) = str2num(get(handles.usedThreshVecEdit07,'String'));
%--------------------------------------------------------------------------
function [handles] = loadVedioTaskData(hObject, eventdata, handles, taskNo)
% hObject    handle to executeBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function Load Video & Store Video Parameters according to taskNo
% Also template information
% If there is not template information , it will do a "initNonVideoStructures" 
% which means "Clear Data Structures"
% Then display;

if ~isfield(handles,'batchRun')
    msgbox('No task !','Error');
    return;
end

taskSum = handles.batchRun.totalTaskNum;
if taskSum == 0
    msgbox('No task !','Error');
    return;
end
videoTaskList = handles.batchRun.videoTaskList;
videoTask = videoTaskList{taskNo};

handles.vidName = videoTask.vidName;
handles.backName = videoTask.backName;
handles.vidObj = VideoReader(handles.vidName);
handles.backObj = VideoReader(handles.backName);
handles.imSize = [handles.vidObj.Height,handles.vidObj.Width];
nFrames = handles.vidObj.NumberOfFrames; 
handles.numFrames = nFrames;
videoTask.numFrames = nFrames;

% Store Backframe
handles.backIm = rgb2gray(handles.backObj.read(1));
% Display Video File Info
set(handles.videoFileEdit,'String',handles.vidName);
% Display Variables
handles.currFrame = 1;
handles.currentLabel = handles.currFrame;

% Set Total Frame Edit
set(handles.totalFramesEdit,'String',num2str(nFrames));

if isfield(videoTask,'analysisFileName')
    set(handles.outputFileEdit,'String',videoTask.analysisFileName);        
end

% Update the videoTask
videoTaskList{taskNo} = videoTask;
handles.batchRun.videoTaskList = videoTaskList;

handles.analysisFileName = videoTask.analysisFileName;
set(handles.analysisFileEdit,'String',videoTask.analysisFileName);

if isfield(videoTask,'mark')
    % Load template data
    handles.mark = videoTask.mark;
    handles.inst = videoTask.inst;
    handles.algoInfo = videoTask.algoInfo; 
    handles.pc = videoTask.pc;
    handles.templateFrame = videoTask.templateFrame;
    
    handles.state = 'waitForLabel';
    set(handles.statusEdit,'String','Wait for Label');
else
    handles = initNonVideoStructures(handles);
end

% begin to load threshold vector 
if(isfield(videoTask,'theThresholdVector'))
    handles.threshVec_ini = videoTask.theThresholdVector;
    guidata(hObject,handles);

    bTop = 1; nVecIdx = 0;
    displayThreshholdVecToEdits(handles,bTop,nVecIdx) ;
    bTop = 0; nVecIdx = 1;
    displayThreshholdVecToEdits(handles,bTop,nVecIdx) ;
end

% end of loading threshold vector

guidata(hObject,handles);
displayFun(handles);

%--------------------------------------------------------------------------
% --- Executes on button press in testTemplateBtn.
function testTemplateBtn_Callback(hObject, eventdata, handles)
% hObject    handle to testTemplateBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% If idle or running automation do nothing
if(any(strcmp(handles.state,{'idle','automation'})))
    return;
end

if ~isfield(handles,'batchRun')
    msgbox('No task !','Error');
    return;
end

taskSum = handles.batchRun.totalTaskNum;
if taskSum == 0
    msgbox('No task !','Error');
    return;
end

% Backup the handles;
handles_Backup = handles;

theLabelCorrect = 0;
theLabelError = 4;

% Now get current task
handles = loadVedioTaskData(hObject, eventdata, handles, handles.currentTaskNo);

videoTaskList = handles.batchRun.videoTaskList;
videoTask = videoTaskList{handles.currentTaskNo};

if ~isfield(videoTask,'mark')
    msgbox('No template in current task!','error');
    return;
end
% Load template data
handles.mark = videoTask.mark;
handles.inst = videoTask.inst;
handles.algoInfo = videoTask.algoInfo;
handles.pc = videoTask.pc;

if isfield('videoTask','theThresholdVector')
    handles.threshVec_ini = videoTask.theThresholdVector;
end
% Save guidata
guidata(hObject,handles);


% Set the test range
startFrame = videoTask.templateFrame;
if videoTask.numFrames - startFrame < 30
    stopFrame =  videoTask.numFrames;
else
    stopFrame = startFrame + 25;
end

if(~strcmp(handles.state,'waitForLabel'))
    return;
end

% Set Program State to Automation
% handles.state = 'automation';

% Set Status Indicator
set(handles.statusEdit,'String','Testing template');
handles.prevAutoPilotStop = stopFrame;

% Move to Frame No. startFrame
handles.currFrame = startFrame;
% Added by Jonathan Feb 19th, 2016
handles.currentLabel = startFrame;
% End of adding Feb 19th, 2016
set(handles.currFrameEdit,'String',num2str(startFrame));
% Save guidata
handles.label = nan*ones(videoTask.numFrames,1);
% guidata(hObject,handles);
displayFun(handles);
% Label the First Frame Correct
handles = labelCurr(handles,theLabelCorrect);
set(handles.currFrameEdit,'String',num2str(startFrame));

% Save guidata
guidata(hObject,handles);

% Get guidata to be used in loop
% Note a local handles copy (h_loop) is used.  The handles structure is
% updated at the end of this function.  This is done to avoid potential
% problems with the handles structure updating when the "Stop Auto"
% button is pressed.

h_loop = guidata(hObject);   
    
bTop = 1; nVecIdx = 0;
displayThreshholdVecToEdits(handles,bTop,nVecIdx);
handles.threshVec = handles.threshVec_ini;
  
% Auto-Pilot Loop
for k = startFrame+1:stopFrame
    % Auto pilot stop check
    if(handles.autoStopFlag == 1)
         set(handles.automationStartEdit,'String',num2str(k));
         set(handles.currFrameEdit,'String',num2str(k));
         set(handles.statusEdit,'String','waitForLabel');
         handles.autoStopFlag = 0;
         break;
    end
    
    % Check Algorithm Confidence Parameters

    [result,fVec] = checkAlgoConf(h_loop,k);

    if(any(result))
        % Algorithm is Underconfident %
        % Display feature vector thresholding in command window
        disp(result);
        % Update Display with current frame
        displayFun(h_loop);

        % Mark frame error
        h_loop = labelCurr(h_loop,theLabelError);
    else   
        displayFun(h_loop);
        h_loop = labelCurr(h_loop,theLabelCorrect);
    end
        
    bTop = 0; nVecIdx = 1;
    displayThreshholdVecToEdits(h_loop,bTop,nVecIdx);

    % Algorithm is Confident %
    % Label frame as correct
    

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
end % end of the loop for

% Display a Messge Box if auto-pilot ran until stopframe
% if(k ~= stopFrame)
%     % Display Message Box Indicating Completion
%     msgbox('Video Auto-Pilot is not Completed!','Auto-Pilot not Complete!','modal');
%     % Move to Frame No. stopFrame
%     h_loop.currFrame = stopFrame + 1;
% end

% Set State to Wait For Label
h_loop.state = 'waitForLabel'; 

if k == stopFrame
    h_loop.currFrame = h_loop.currFrame - 1;
end
set(h_loop.currFrameEdit,'String',num2str(h_loop.currFrame));

handles = h_loop;

% Update Display
displayFun(handles);
displayThreshholdVecToEdits(handles,bTop,nVecIdx);

% Update Status
handles = handles_Backup;
handles = loadVedioTaskData(hObject, eventdata, handles, handles.currentTaskNo);

guidata(hObject,handles);

display('The template test completed. If you need adjustment, push "Initialize Template" button, please.');
msgbox('The template test completed. If you need adjustment, push "Initialize Template" button, please.','Test completed.');


%--------------------------------------------------------------------------
function [handles] = iniProgressBar(hObject,handles)
% hObject    handle to testTemplateBtn (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% This function initialize progress bar "TotalTaskPercentageAxes" and 
% "CurrentTaskPercentageAxes"

set(handles.CurrentTaskPercentageAxes,...
    'XLim', [0 1000], ...
    'YLim', [0 1],...
    'Box', 'on', ...
    'ytick', [], ...
    'xtick', [0,250,500,750]);
set(handles.CurrentTaskPercentageAxes,'XTickLabel',...
    {' ';'25%';'50%';'75%'},'FontSize',8);
axes(handles.CurrentTaskPercentageAxes);

hPatchCurrentTaskPercentage = patch( ...
                           'XData', [0 1 1 0], ...
                           'YData', [0 0 1 1] );
set(handles.textPercentageCurrTask,'string','0%');
handles.hPatchCurrentTaskPercentage = hPatchCurrentTaskPercentage;
                       
set(handles.TotalTaskPercentageAxes,...
    'XLim', [0 1000], ...
    'YLim', [0 1],...
    'Box', 'on', ...
    'ytick', [], ...
    'xtick', [0,250,500,750]);
set(handles.TotalTaskPercentageAxes,'XTickLabel',...
    {' ';'25%';'50%';'75%'},'FontSize',8);
axes(handles.TotalTaskPercentageAxes);
hPatchTotalTaskPercentage = patch( ...
                           'XData', [0 1 1 0], ...
                           'YData', [0 0 1 1] );
set(handles.textPercentageTotalTask,'string','0%');
handles.hPatchTotalTaskPercentage = hPatchTotalTaskPercentage;

guidata(hObject,handles);
%--------------------------------------------------------------------------
function setEvaluateNum(currNum,totalNum,hNumPatch,hTxt)
% Calculate Percentage
x = currNum / totalNum * 1000;
theNum = round(x);
displayStr = [num2str(theNum / 10),'%'];

thiscolor = [];
Color_0_25 = [0.5437    0.9848    0.7157];
Color_25_49 =[ 0.1190    0.4984    0.95970];
Color_50_74 = [0.6074    0.1917    0.7384];
Color_75_100 = [0.8584    0.1253    0.0897];

if theNum <= 250
    thiscolor = Color_0_25;
elseif theNum <500 && theNum > 250
    thiscolor = Color_25_49;
elseif theNum <750 && theNum >= 500
    thiscolor = Color_50_74;
else
    thiscolor = Color_75_100;
end

% rectangle('Position',[0,0,theNum,1],'facecolor',thiscolor);
set(hNumPatch,...
    'XData', [0 theNum theNum 0], ...
    'YData', [0 0 1 1] );
set(hNumPatch, 'FaceColor', thiscolor);
set(hTxt,'String',displayStr);

%--------------------------------------------------------------------------
function [cellFileContent,numFileLine] = readTxtFile(theTxtFileName)
% this fuction read in lines from a txt file
% input :
% theTxtFileName: the name of the txt file, including the path
% output:
% cellFileContent: the Content of the txt file, erery line is a member if
%                   this cell.
% numFileLine: the number of the lines in this file
fid1=fopen(theTxtFileName,'rt');
if fid1 == -1
    display('Cannot open the txtFile!');
    cellFileContent = [];
    numFileLine = -1;
    return;
end
% begin to read the file
numFileLine = 0;
while ~feof(fid1)
    aline=fgetl(fid1);
    numFileLine = numFileLine + 1;
    cellFileContent{numFileLine} = aline;
end
fclose(fid1);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% --- Executes on selection change in comboBoxPresetThresholdVector.
function comboBoxPresetThresholdVector_Callback(hObject, eventdata, handles)
% hObject    handle to comboBoxPresetThresholdVector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns comboBoxPresetThresholdVector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from comboBoxPresetThresholdVector

bUseRecommneded = get(handles.checkboxUseRecommended,'Value');
if bUseRecommneded == 0;
    return;
end

% if bUsegetRecommneded==1, % if use Recommended threshold vector
% 1. check the selected item of the comboBox
selectedThreshVecNo = get(handles.comboBoxPresetThresholdVector,'Value');
% 2. load the pre-set threshold vector
selRecommendThreshVec = getPresetThreshVec(selectedThreshVecNo);

% apply the threshold vector and display
handles.threshVec_ini = selRecommendThreshVec;
bTop = 1; nVecIdx = 0;
displayThreshholdVecToEdits(handles,bTop,nVecIdx) ;
bTop = 0;
displayThreshholdVecToEdits(handles,bTop,nVecIdx) ;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function comboBoxPresetThresholdVector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to comboBoxPresetThresholdVector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in checkboxUseRecommended.
function checkboxUseRecommended_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxUseRecommended (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxUseRecommended
bUseRecommneded = get(handles.checkboxUseRecommended,'Value');

if bUseRecommneded==1
    % disable Apply button and load button
    set(handles.applyThresholdToVedioBtn,'Enable','off');
    set(handles.loadThresholdFromFileBtn,'Enable','off');
    % if use Recommended threshold vector
    % 1. check the selected item of the comboBox
    selectedThreshVecNo = get(handles.comboBoxPresetThresholdVector,'Value');
    % 2. load the pre-set threshold vector
    selRecommendThreshVec = getPresetThreshVec(selectedThreshVecNo);
    
    % apply the threshold vector and display
    handles.threshVec_ini = selRecommendThreshVec;
    
    bTop = 1; nVecIdx = 0;
    displayThreshholdVecToEdits(handles,bTop,nVecIdx) ;
    bTop = 0;
    displayThreshholdVecToEdits(handles,bTop,nVecIdx) ;

    guidata(hObject,handles);
else
    set(handles.applyThresholdToVedioBtn,'Enable','on');
    set(handles.loadThresholdFromFileBtn,'Enable','on');
end
%--------------------------------------------------------------------------
function [selRecommendThreshVec] = getPresetThreshVec(idx)
knifeThreshVec      = [1,1,1,1,1,1,1];
scissorsThreshVec   = [2,2,2,2,2,2,2];
forcepsThreshVec    = [3,3,3,3,3,3,3];

if idx == 1
    selRecommendThreshVec = knifeThreshVec;
elseif idx == 2
    selRecommendThreshVec = scissorsThreshVec;
else % idx == 3
    selRecommendThreshVec = forcepsThreshVec;
end

%--------------------------------------------------------------------------
% --- Executes on button press in setThreshVecToTaskBtn.
function setThreshVecToTaskBtn_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecToTaskBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This function set the threshold vector to current Task
% If idle or running automation do nothing
if(any(strcmp(handles.state,{'idle','automation'})))
    return;
end

if ~isfield(handles,'batchRun')
    msgbox('No task !','Error');
    return;
end

% if there is not task, return
taskSum = handles.batchRun.totalTaskNum;
if taskSum == 0
    msgbox('No task !','Error');
    return;
end
% get current task
videoTaskList = handles.batchRun.videoTaskList;
if ~isfield(handles,'currentTaskNo')
    handles.currentTaskNo = 1;
end
videoTask = videoTaskList{handles.currentTaskNo};

% If waiting for label, verify with user that re-initializing the template
% will result in a loss of tracking data.
if(isfield(videoTask,'theThresholdVector'))
    % Prompt User
    msgAns = questdlg('A Threshold Vector has been set for this video. Reset it?',...
    'Reset Threshold Vector','Yes','No','No');
    
    % Do nothing if user does not want to reinitialize template
    if(~strcmp('Yes',msgAns))
        return;
    end
 
end

videoTask.theThresholdVector = handles.threshVec_ini;
videoTaskList{handles.currentTaskNo} = videoTask;
handles.batchRun.videoTaskList = videoTaskList;

guidata(hObject,handles);


% --- Executes on button press in addTaskRangBtn.
function addTaskRangBtn_Callback(hObject, eventdata, handles)
% This fuction set up the frame range of every task. One Task can have
% several range. Frames in these range will be taken as 'Error'. If there
% is no range setting, the program will check every frame of the video.
% This function is optional. 
% hObject    handle to addTaskRangBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% If idle or running automation do nothing
if(any(strcmp(handles.state,{'idle','automation'})))
    return;
end

if ~isfield(handles,'batchRun')
    msgbox('No task !','Error');
    return;
end

% if there is not a task, return
taskSum = handles.batchRun.totalTaskNum;
if taskSum == 0
    msgbox('No task !','Error');
    return;
end
% get current task
videoTaskList = handles.batchRun.videoTaskList;
if ~isfield(handles,'currentTaskNo')
    handles.currentTaskNo = 1;
end
videoTask = videoTaskList{handles.currentTaskNo};

% Get range info from user using modal dialog
% rangeInfo = [1,50,75;2,78,103];
% prepare rangeInfo
if (~isfield(videoTask, 'rangeInfo'))
    rangeInfo = 0;
else
    rangeInfo = videoTask.rangeInfo;
end

% Now jumps out the range setting window to set the range
rangeSetting = taskRangeSettingDlg('rangeInfo',rangeInfo,'numberOfFrames',videoTask.numFrames);
bCancel = rangeSetting{1};
% Check for user cancel
if(bCancel)
    return;
end
% if not cancel
rangeInfo = rangeSetting{2};
videoTask.rangeInfo = rangeInfo;
videoTaskList{handles.currentTaskNo} = videoTask;
handles.batchRun.videoTaskList = videoTaskList;

guidata(hObject,handles);


% --- Executes on button press in ceaseBtn.
function ceaseBtn_Callback(hObject, eventdata, handles)
% hObject    handle to ceaseBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.checkboxCeasing,'Value',1);


% --- Executes on button press in checkboxCeasing.
function checkboxCeasing_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCeasing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxCeasing



function outputFileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to outputFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputFileEdit as text
%        str2double(get(hObject,'String')) returns contents of outputFileEdit as a double


% --- Executes during object creation, after setting all properties.
function outputFileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in outputFileSelBtn.
function outputFileSelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to outputFileSelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% This Button select batch txt files
[fName,pName] = uiputfile({'*.mat','(*.mat) Analysis Files'}...
            ,'Change the output file name');
% Check for User Cancel
if(fName == 0)
    return;
end

analysisFileName = fullfile(pName,fName);

% get current task
videoTaskList = handles.batchRun.videoTaskList;
if ~isfield(handles,'currentTaskNo')
    handles.currentTaskNo = 1;
end
videoTask = videoTaskList{handles.currentTaskNo};

videoTask.analysisFileName = analysisFileName;
set(handles.outputFileEdit,'String',analysisFileName);

videoTaskList{handles.currentTaskNo} = videoTask;
handles.batchRun.videoTaskList = videoTaskList;

guidata(hObject,handles);
