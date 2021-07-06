function varargout = UGT_V20(varargin)
% UGT_V20 MATLAB code for UGT_V20.fig
%      UGT_V20, by itself, creates a new UGT_V20 or raises the existing
%      singleton*.
%
%      H = UGT_V20 returns the handle to a new UGT_V20 or the handle to
%      the existing singleton*.
%
%      UGT_V20('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UGT_V20.M with the given input arguments.
%
%      UGT_V20('Property','Value',...) creates a new UGT_V20 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UGT_V20_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UGT_V20_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singlet on)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UGT_V20

% Last Modified by GUIDE v2.5 16-Nov-2015 16:11:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UGT_V20_OpeningFcn, ...
                   'gui_OutputFcn',  @UGT_V20_OutputFcn, ...
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


% --- Executes just before UGT_V20 is made visible.
function UGT_V20_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UGT_V20 (see VARARGIN)

% Choose default command line output for UGT_V20
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

% right part
handles.state_Right = 'idle';
set(handles.statusEdit_Right,'String','Waiting for File');

% Automation Stop Flag
% Used to stop automation prior to reaching the stop frame
handles.autoStopFlag = 0;

% Setup Auto-Pilot
% threshVec = [delta_NCC,delta_width,delta_orient,delta_inliersDistL,...
%              delta_inliersDistR,inliersL,inliersR]
% Threshold Parameters Set 2011-11-25
handles.threshVec = [0.0332,1.72,0.00987,10,10,10,10];
handles.threshVec_Right = [0.0332,1.72,0.00987,10,10,10,10];
handles.threshVec_ini = [0.0332,1.72,0.00987,10,10,10,10];
handles.threshVec_ini_Right = [0.0332,1.72,0.00987,10,10,10,10];
handles.threshVec_org = [0.0332,1.72,0.00987,10,10,10,10];
handles.typeVec = [0,0,0,0,0,1,1];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UGT_V20 wait for user response (see UIRESUME)
% uiwait(handles.UGT_V20);


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
% NCCy
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

% this function Checks Algorithm Confidence Parameters 
function [result, fVec] = checkAlgoConf(handles,frameNum,bLeft)
% CHECKALGOCONF Checks Algorithm Confidence Parameters
%
% RESULT = CHECKALGOCONF(HANDLES,FRAMENUM) Checks algorithm tracking
% confidence at FRAMENUM.  HANDLES is the gui HANDLES structure containing
% mark, inst, algoInfo data structures with tracking data and threshVec and
% typeVec vectors.  RESULT is a binary vector, value = 0 means that the
% algorithm is confident about the parameter.  A value = 1 means the
% algorithm underconfident about the parameter.

% Get Confidence Parameter Vector

if nargin<=2
    bLeft = 1;
end

if bLeft==1    
    % Get Confidence Parameters of left 
    fVec = getConfParams(handles.mark,handles.inst,handles.algoInfo,frameNum);
    % Get Result of Confidence Classification
    result = threshClassifyFun(fVec,handles.threshVec,handles.typeVec);
else
    % get Confidence Parameters of right
    fVec = getConfParams(handles.mark_Right,handles.inst_Right,handles.algoInfo_Right,frameNum);    
    % Get Result of Confidence Classification
    result = threshClassifyFun(fVec,handles.threshVec_Right,handles.typeVec);
end

% paint the bad element in red
paintAbnormalEdit(handles,result);
% display the  Confidence Parameters to the bottom edit line
bTop = 0;
displayThreshholdVecToEdits(handles,bLeft,bTop,2,fVec);

% --- Outputs from this function are returned to the command line.
function varargout = UGT_V20_OutputFcn(hObject, eventdata, handles) 
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
iModeOperation = getSelObject(handles);
if  iModeOperation == 1 % Right only mode
    return;
end

if(strcmp(handles.state,{'automation','idle'}))
    % Do Nothing if tracking file has not been started or
    % automation is running
    return;
end
if  iModeOperation == 2 % Both & Synchro mode
    if(strcmp(handles.state_Right,{'automation','idle'}))
        % Do Nothing if tracking file has not been started or
        % automation is running
        return;
    end	
end

if(handles.currFrame > 1)
    handles.currFrame = handles.currFrame - 1;
    displayFun(handles);
    guidata(hObject,handles);
end

if  iModeOperation == 2 % Both & Synchro mode
    if(handles.currFrame_Right > 1)
        handles.currFrame_Right = handles.currFrame_Right - 1;
        displayFun(handles,0);
        guidata(hObject,handles);
    end
end

% --- Executes on button press in nextFrameButton.
function nextFrameButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextFrameButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iModeOperation = getSelObject(handles);
if  iModeOperation == 1 % Right only mode
    return;
end

if(strcmp(handles.state,{'automation','idle'}))
    % Do Nothing if tracking file has not been started or
    % automation is running
    return;
end

if  iModeOperation == 2 % Both & Synchro mode
    if(strcmp(handles.state_Right,{'automation','idle'}))
        % Do Nothing if tracking file has not been started or
        % automation is running
        return;
    end
end

if(handles.currFrame < handles.numFrames)
    handles.currFrame = handles.currFrame + 1;
    displayFun(handles);
    guidata(hObject,handles);
end

if  iModeOperation == 2 % Both & Synchro mode
    if(handles.currFrame_Right < handles.numFrames_Right)
        handles.currFrame_Right = handles.currFrame_Right + 1;
        displayFun(handles,0);
        guidata(hObject,handles);
    end
end


% -------------------------------------------------------------------------
% --- Executes on button press in initializeTemplateButton.
function initializeTemplateButton_Callback(hObject, eventdata, handles)
% hObject    handle to initializeTemplateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% If idle or running automation do nothing
if(any(strcmp(handles.state,{'idle','automation'})))
    return;
end

iModeOperation = getSelObject(handles);
if  iModeOperation == 0 || iModeOperation == 2% Left only mode or Both & Synchro mode
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
end

% -------------------------------------------------------------------------

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

%--------------------------------------------------------------------------
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


%--------------------------------------------------------------------------
function totalFramesEdit_Callback(hObject, eventdata, handles)
% hObject    handle to totalFramesEdit_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of totalFramesEdit_R as text
%        str2double(get(hObject,'String')) returns contents of totalFramesEdit_R as a double

%--------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function totalFramesEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to totalFramesEdit_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% --- Executes on button press in labelCorrectButton.
function labelCorrectButton_Callback(hObject, eventdata, handles)
% hObject    handle to labelCorrectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% get Operation Model, Left, Right or Both
iModelOperation = getSelObject(handles);
if iModelOperation == 0 || iModelOperation == 2
    if(~strcmp(handles.state,'waitForLabel'))
        return;
    end
    
    handles = labelCurr(handles,0);
    displayFun(handles);
    guidata(hObject,handles);
end
if iModelOperation == 1 || iModelOperation == 2
    if(~strcmp(handles.state_Right,'waitForLabel'))
        return;
    end
    
    handles = labelCurr(handles,0,0);
    displayFun(handles,0);
    guidata(hObject,handles);
end

%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
% --- Executes on button press in labelOtherButton.
function labelOtherButton_Callback(hObject, eventdata, handles)
% hObject    handle to labelOtherButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iModelOperation = getSelObject(handles);
if iModelOperation == 0 || iModelOperation == 2
    if(~strcmp(handles.state,'waitForLabel'))
        return;
    end
    
    handles = labelCurr(handles,4);
    displayFun(handles);
    guidata(hObject,handles);
end
if iModelOperation == 1 || iModelOperation == 2
    if(~strcmp(handles.state_Right,'waitForLabel'))
        return;
    end
    
    handles = labelCurr(handles,4,0);
    displayFun(handles,0);
    guidata(hObject,handles);
end

%--------------------------------------------------------------------------
function handles = labelCurr(handles,labelNum,bLeft)
% LABELCURR Labels the Current Frame & Runs tracking on next frame
%
% HANDLES = labelCurr(HANDLES,LABELNUM)  Labels the current frame as
% integer LABELNUM -> {0:Correct,1:Blur,2:Occlusion,3:OtherError}.  After
% labeling the tracking is run on the next frame and the current frame and
% label are updated.  If the current frame is not the current label the
% user is prompted if they want to change the label.  If this is desired
% the current frame label is changed and all proceeding frames are
% unlabeled.
if nargin<=2
    bLeft = 1;
end

if bLeft==1
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
else % Right
    if((strcmp(handles.state_Right,'waitForLabel') || strcmp(handles.state_Right,'automation')) && ...
          handles.currFrame_Right <= handles.currentLabel_Right)
      if(handles.currFrame_Right == handles.currentLabel_Right)
           % Label The Unlabeled Frames
           handles.label_Right(handles.currFrame_Right) = labelNum;
           handles = performTracking(handles,handles.currFrame_Right + 1,0);
           handles.currFrame_Right = handles.currFrame_Right + 1;
           handles.currentLabel_Right = handles.currFrame_Right;
           %displayFun(handles);
      else
          % Label & Clear
          btn = questdlg('This frame in right has already been labeled.  By labeling it, all labels set after it will be lost.  Are you sure?',...
              'Change already set label?','No');
          if(strcmp(btn,'Yes'))
              handles.label_Right(handles.currFrame_Right) = labelNum;
              handles = clearAfter(handles,handles.currFrame_Right,0);
              handles = performTracking(handles,handles.currFrame_Right + 1,0);
              handles.currFrame_Right = handles.currFrame_Right + 1;
              handles.currentLabel_Right = handles.currFrame_Right;
              %displayFun(handles);
          end 
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


function currVidDir = getCurrVidDir(handles,bLeft)
% GETCURRVIDDIR Gets the Current Video's Directory
% 
% CURRVIDDIR = getCurrVidDir(HANDLES) Checks if a video file is associated
% with the HANDLES structure.  If so returns directory name containing
% video in CURRVIDDIR, otherwise returns 0.
if nargin<2
    bLeft = 1;
end

if bLeft == 1
    % Check if vidName exists
    if(isfield(handles,'vidName'))
        % Check that it is not empty
        if(~isempty(handles.vidName))
            % Get Directory and Return
            currVidDir = fileparts(handles.vidName);
            return;
        end
    end
end

if bLeft == 0
    % Check if vidName exists
    if(isfield(handles,'vidName_Right'))
        % Check that it is not empty
        if(~isempty(handles.vidName_Right))
            % Get Directory and Return
            currVidDir = fileparts(handles.vidName_Right);
            return;
        end
    end
end

% Otherwise Return 0
currVidDir = 0;

function handles = initNonVideoStructures(handles,bLeft)
% INITNONVIDEOSTRUCTURES Initializes all non-video data structures and vars
%
% HANDLES = initNonVideoStructures(handles)  Initializes tracking
% structures : label(NaNs), mark(NaNs), inst(NaNs), and algoInfo(NaNs).
% Also initializes template matching structures neighborhoodSize(2),
% pc(empty), and templateFrame(empty).  Sets the axis to display frames at
% the native size.  Initializes the currentLabel variable (empty).  Sets
% system state to 'waitForTemplate'.
if nargin==1
    bLeft = 1;
end

if bLeft == 1
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
    % Left init completed
else
    nFrames = handles.numFrames_Right;

    % Analysis File
    handles.analysisFileName_Right = [];
    set(handles.RightAnalysisFileEdit,'String','');

    % Set Axis Dimensions to that of Image
    set(handles.dispAxesRight,'Units','pixels');
    pos = get(handles.dispAxesRight,'Position');
    pos(3:4) = handles.imSize_Right([2,1]);
    set(handles.dispAxesRight,'Position',pos);

    % Init Structs
    handles.label_Right = nan*ones(nFrames,1);
    handles.mark_Right = struct('tCorn',nan*ones(nFrames,2),'wStats',[],...
        'nCorn',nan*ones(nFrames,2),'subT',nan*ones(nFrames,2));
    handles.inst_Right = struct('rho',nan*ones(nFrames,2),'theta',...
        nan*ones(nFrames,2),'trackPt',nan*ones(nFrames,2));
    handles.algoInfo_Right = struct('nccScore',nan*ones(nFrames,1),'fitMSE',...,
        nan*ones(nFrames,1),'votesLeft',nan*ones(nFrames,1),'votesRight',...
        nan*ones(nFrames,1),'leftNumInliers',nan*ones(nFrames,1),...
        'leftRefitMSE',nan*ones(nFrames,1),'rightNumInliers',...
        nan*ones(nFrames,1),'rightRefitMSE',nan*ones(nFrames,1));

    % Template Related
    handles.neighborhoodSize_Right = 2;
    handles.pc_Right = []; % Precomputation Struct Used for Template Tracking
    handles.templateFrame_Right = [];

    % Analysis Progress Variables
    handles.currentLabel_Right = [];

    % Set State
    handles.state_Right = 'waitForTemplate';
    set(handles.statusEdit_Right,'String','Waiting For Initialize');

    % set up display para for right video
    % bHasRight = get(handles.ConnectedWithRightCheckbox,'Value');
end


% -------------------------------------------------------------------
% This function reply the Menu item newTracking click
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
set(handles.totalFramesEdit_R,'String',num2str(nFrames));

% Initialize all non-video data structures.  These structures are used to
% store tracking data and the tracking file name
handles = initNonVideoStructures(handles);
guidata(hObject,handles);

% Show First Frame
displayFun(handles);

function currTrackFileDir = getCurrTrackFileDir(handles,bLeft)
% GETCURRTRACKFILEDIR Gets the Current Track File's Directory
% 
% CURRTRACKFILEDIR = getCurrTrackFileDir(HANDLES) Checks if a track file is
% associated with the HANDLES structure.  If so returns directory name
% containing track file in CURRTRACKFILEDIR, otherwise returns 0.

if nargin<2
    bLeft = 1;
end

if bLeft == 1
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
end

if bLeft == 0
    % Check if analysisFileName exists
    if(isfield(handles,'analysisFileName_Right'))
        % Check that it is not empty
        if(~isempty(handles.analysisFileName_Right))
            % Get Directory and Return
            currTrackFileDir = fileparts(handles.analysisFileName_Right);
            return;
        end
    end

    % Otherwise Return 0
    currTrackFileDir = 0;
end



% --------------------------------------------------------------------
function loadTrackingMenu_Callback(hObject, eventdata, handles)
% hObject    handle to loadTrackingMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% currTrackFileDir = getCurrTrackFileDir(handles);
% if(currTrackFileDir ~= 0)
%     % Open File Browser in Current Track File Directory
%     [fName,pName] = uigetfile('*.mat','Select Tracking File',...                                                    
%                                                     currTrackFileDir);
% else
%     % No Tracking File Avaialble.  Check Video File
%     currVidDir = getCurrVidDir(handles);
%     if(currVidDir ~= 0)
%         [fName,pName] = uigetfile('*.mat','Select Tracking File',...
%                                         currVidDir);
%     else
%         [fName,pName] = uigetfile('*.mat','Select Tracking File');
%     end
% end
% 
% if(fName == 0)
%     return;
% end
% 
% load(fullfile(pName,fName),'saveStruct');
% 
% saveFields = {'mark','pc','neighborhoodSize','templateFrame','inst',...
% 'algoInfo','vidName','backName','imSize','numFrames',...
% 'backIm','state','currentLabel','label','currFrame','analysisFileName'};
% 
% % Load all Save Fields
% for k = 1:numel(saveFields)
%     handles.(saveFields{k}) = saveStruct.(saveFields{k});
% end
% 
% % Check if tracking file name or path has changed
% 
% 
% 
% 
% % Setup Video Objects
% handles.vidObj = VideoReader(handles.vidName);
% handles.backObj = VideoReader(handles.backName);
% displayFun(handles);
% 
% % Set Axis Dimensions to that of Image
% set(handles.dispAxes,'Units','pixels');
% pos = get(handles.dispAxes,'Position');
% pos(3:4) = handles.imSize([2,1]);
% set(handles.dispAxes,'Position',pos);
% 
% % Display Video File & Analysis File Names
% set(handles.videoFileEdit,'String',handles.vidName);
% set(handles.analysisFileEdit,'String',handles.analysisFileName);
% 
% % Set Label to waitForLabel
% handles.state = 'waitForLabel';
% set(handles.statusEdit,'String','Waiting for Label');
% 
% guidata(hObject,handles);
loadTrackFileBtn_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function saveTrackingMenu_Callback(hObject, eventdata, handles)
% hObject    handle to saveTrackingMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Do Nothing if Idle
% if(strcmp(handles.state,'idle'))
%     return;
% end
% 
% if(isempty(handles.analysisFileName))
%     currVidDir = getCurrVidDir(handles);
%     if(currVidDir ~= 0)
%         % Open File Browser in current video directory
%         [fName,pName] = uiputfile('*.mat','Save Tracking File',currVidDir);
%     else
%         [fName,pName] = uiputfile('*.mat','Save Analysis File');
%     end
%     
%     % Check for User Cancel
%     if(fName == 0)
%         return;
%     end
%     
%     handles.analysisFileName = fullfile(pName,fName);
% end
% 
% saveFields = {'mark','pc','neighborhoodSize','templateFrame','inst',...
% 'algoInfo','vidName','backName','imSize','numFrames',...
% 'backIm','state','currentLabel','label','currFrame','analysisFileName'};
% 
% saveStruct = handles;
% allFields = fieldnames(saveStruct);
% 
% % Remove all non save fields
% for k = 1:numel(allFields)
%     cName = allFields{k};
%     if(~any(strcmp(cName,saveFields)))
%         saveStruct = rmfield(saveStruct,cName);
%     end
% end
% 
% save(handles.analysisFileName, 'saveStruct');
% set(handles.analysisFileEdit,'String',handles.analysisFileName);
% guidata(hObject,handles);
LTrackFileSelBtn_Callback(hObject, eventdata, handles);

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

function displayFun(handles,bLeft)
% DISPLAYFUN Displays current frame and overlays tracking data
%
% displayFun(HANDLES) Displays the current frame and overlays the template
% window and boundary line estimates if available.
if nargin<2
    bLeft = 1;
end
if bLeft==1
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
else
    currFrame = handles.currFrame_Right;

    % Draw Marker Mask %
    tCurr = handles.mark_Right.tCorn(currFrame,:);
    nCurr = handles.mark_Right.nCorn(currFrame,:);
    imSize = handles.imSize_Right;

    if(~any(isnan(tCurr)))
        tSize = handles.mark_Right.wStats.wSize;
        nSize = handles.mark_Right.wStats.nx-1;
        markMask = trackMask(tCurr,nCurr,tSize,nSize,imSize);
    else
        markMask = false(imSize);
    end

    % Draw Boundary Line Mask %
    rho = handles.inst_Right.rho(currFrame,:);
    theta = handles.inst_Right.theta(currFrame,:);
    trackPt = handles.inst_Right.trackPt(currFrame,:);
    if(~any(isnan([rho,theta,trackPt])))
        lineMask = drawLineMask(imSize,[rho,mean(rho)],[theta,mean(theta)]);
        lineMask(round(tCurr(2)):end,:) = 0;
        lineMask((-5:5) + round(trackPt(2)),(-5:5) + round(trackPt(1))) = 1;
    else
       lineMask =  false(imSize);
    end

    % Draw Full Image
    frameIm = rgb2gray(handles.vidObj_Right.read(currFrame));
    fullIm = genOverlayIm(frameIm,(lineMask | markMask));
    imshow(fullIm,'InitialMagnification',100,'Parent',handles.dispAxesRight);

    % Update GUI Objects
    % Algorithm Edit
%     if(currFrame == handles.templateFrame_Right)
%         set(handles.algorithmTypeEdit,'String','Template Frame');
%     elseif(currFrame == 1)
%         set(handles.algorithmTypeEdit,'String','');
%     else
%         set(handles.algorithmTypeEdit,'String',...
%             getAlgoTypeString(handles.label_Right(currFrame-1)));
%     end

    % Label & Previous Frame Label
    set(handles.labelEdit_Right,'String',getLabelString(handles.label_Right(currFrame)));
    if(currFrame ~= 1)
        set(handles.previousFrameLabelEdit_Right,'String',...
            getLabelString(handles.label_Right(currFrame-1)));
    else
        set(handles.previousFrameLabelEdit_Right,'String','');
    end

    % Current Frame
    set(handles.currFrameEdit_R,'String',num2str(currFrame));
    set(handles.totalFramesEdit_R,'String',num2str(handles.numFrames_Right));
end


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

function handles = performTracking(handles,frameNum,bLeft)
% PERFORMTRACKING Runs Instrument Tracking/Detection on single frame
%
% HANDLES = performTracking(HANDLES,FRAMENUM) Instrument tracking/detection
% algorithm.  The previous frame label determines whether detection or
% tracking is performed.  All instrument tracking data structures in
% returned HANDLES are updated.
if nargin<=2
    bLeft = 1;
end

if bLeft==1
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
end

if bLeft==0
    if(isnan(handles.label_Right(frameNum-1)))
        return;
    end

    frameIm = rgb2gray(handles.vidObj_Right.read(frameNum));

    fail_bLine = 0;
    
    switch(handles.label_Right(frameNum-1))
        case 0
            % Perform Tracking
            % Template
            tPrev = handles.mark_Right.tCorn(frameNum-1,:);
            if(frameNum > 2 && handles.label_Right(frameNum - 2) == 0)
                tPrev2 = handles.mark_Right.tCorn(frameNum-2,:);
            else
                tPrev2 = tPrev;
            end
            trackPtPrev = handles.inst_Right.trackPt(frameNum-1,:);
            thetaPrev = handles.inst_Right.theta(frameNum-1,:);
            nPrev = handles.mark_Right.nCorn(frameNum-1,:);
            nSize = handles.mark_Right.wStats.nx-1;

            [tCurr,nCurr,maxScore,fitMSE] = track_temp(frameIm,handles.pc_Right,tPrev,...
                                                tPrev2,nPrev,nSize);

            % Boundary Line Estimates
            try
                [r,t,dbl_ais] = track_bLines(frameIm,handles.backIm_Right,...
                                        round(tCurr),tPrev,trackPtPrev,thetaPrev);
            catch ME
                fail_bLine = 1;
            end
        otherwise
            % Perform Detection
            % Template Detection
            [tCurr,nCurr,maxScore,fitMSE] = detection_temp(frameIm,...
                handles.pc_Right,handles.mark_Right.wStats);

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
                [r,t,dbl_ais] = detection_bLines(frameIm,handles.backIm_Right,...
                    round(tCurr),handles.mark_Right.wStats,handles.pc_Right);
                %[r,t,dbl_ais] = track_bLines(frameIm,handles.backIm,...
                %    round(tCurr),tCurr,tCurr + trackPtDelta,thetaPrev);
            catch ME
                fail_bLine = 1;
            end
    end
    % Store Results
    handles.mark_Right.tCorn(frameNum,:) = round(tCurr);
    handles.mark_Right.subT(frameNum,:) = tCurr;
    handles.mark_Right.nCorn(frameNum,:) = round(nCurr);

    % Template Related Algo Info
    handles.algoInfo_Right.nccScore(frameNum,1) = maxScore;
    handles.algoInfo_Right.fitMSE(frameNum,1) = fitMSE;

    if(fail_bLine == 0)
        % Boundary Line Estimates
        handles.inst_Right.rho(frameNum,:) = r;
        handles.inst_Right.theta(frameNum,:) = t;
        halfSize = round((handles.mark_Right.wStats.wSize-1)/2);
        handles.inst_Right.trackPt(frameNum,:) = computeTrackPt(r,t,tCurr,halfSize);

        % Boundary Line Algo Info
        handles.algoInfo_Right.votesLeft(frameNum,:) = dbl_ais.votesLeft;
        handles.algoInfo_Right.votesRight(frameNum,:) = dbl_ais.votesRight;
        handles.algoInfo_Right.leftNumInliers(frameNum,:) = dbl_ais.leftNumInliers;
        handles.algoInfo_Right.leftRefitMSE(frameNum,:) = dbl_ais.leftRefitMSE;
        handles.algoInfo_Right.rightNumInliers(frameNum,:) = dbl_ais.rightNumInliers;
        handles.algoInfo_Right.rightRefitMSE(frameNum,:) = dbl_ais.rightRefitMSE;
    end
end



function handles = clearAfter(handles,afterFrame,bLeft)
% CLEARAFTER Clears all stored data after a specific frame
%
% HANDLES = clearAfter(HANDLES,AFTERFRAME) Clears labels and stored
% tracking data in all frames after AFTERFRAME.
if nargin<=2
    bLeft = 1;
end

if bLeft==1
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
else
    f = afterFrame + 1;

    handles.label_Right(f:end,:) = nan;

    % Store Results
    handles.mark_Right.tCorn(f:end,:) = nan;
    handles.mark_Right.subT(f:end,:) = nan;
    handles.mark_Right.nCorn(f:end,:) = nan;

    % Template Related Algo Info
    handles.algoInfo_Right.nccScore(f:end,1) = nan;
    handles.algoInfo_Right.fitMSE(f:end,1) = nan;

    % Boundary Line Estimates
    handles.inst_Right.rho(f:end,:) = nan;
    handles.inst_Right.theta(f:end,:) = nan;
    %halfSize = round((handles.mark.wStats.wSize-1)/2);
    handles.inst_Right.trackPt(f:end,:) = nan;

    % Boundary Line Algo Info
    handles.algoInfo_Right.votesLeft(f:end,:) = nan;
    handles.algoInfo_Right.votesRight(f:end,:) = nan;
    handles.algoInfo_Right.leftNumInliers(f:end,:) = nan;
    handles.algoInfo_Right.leftRefitMSE(f:end,:) = nan;
    handles.algoInfo_Right.rightNumInliers(f:end,:) = nan;
    handles.algoInfo_Right.rightRefitMSE(f:end,:) = nan;
end



% --- Executes on key press with focus on UGT_V20 or any of its controls.
function UGT_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to UGT_V20 (see GCBO)
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
% Get Start Frame from 'automationStartEdit'

startFrame = getNumfromEdit(handles.automationStartEdit,'Automation Frame From Edit');
% Get End Frame from 'automationStopEdit'
stopFrame = getNumfromEdit(handles.automationStopEdit,'Automation End Frame Edit');

iModelOperation = getSelObject(handles);

% Check startFrame & stopFrame
if iModelOperation == 0 || iModelOperation == 2
    if startFrame<1 || startFrame>handles.numFrames || startFrame>stopFrame ...
            || stopFrame<1 || stopFrame>handles.numFrames 
        msgbox('Incorrect End Frame Entered','Incorrect End Frame','error','modal'); 
        return;
    end
else 
    if iModelOperation == 0 || iModelOperation == 2
        if startFrame<1 || startFrame>handles.numFrames_Right || startFrame>stopFrame ...
                || stopFrame<1 || stopFrame>handles.numFrames_Right 
            msgbox('Incorrect End Frame Entered','Incorrect End Frame','error','modal'); 
            return;
        end
    end
end

% if Left Only
if iModelOperation == 0 
    if(~strcmp(handles.state,'waitForLabel'))
        return;
    end
    
    if iModelOperation == 2
        if(~strcmp(handles.state_Right,'waitForLabel'))
            return;
        end
    end
    
    set(handles.radiobuttonRightThreshold,'Value',0);
    set(handles.radiobuttonLeftThreshold,'Value',1);

    % Get End Frame
%     startFrame = handles.currentLabel;
%     handles.currFrame = startFrame;
%     displayFun(handles);
%     dlgString = ['Start Frame : ' num2str(startFrame) ...
%                  '.  Please Select End Frame'];
             
%     % Auto pilot stop check
%     if(isfield(handles,'prevAutoPilotStop'))
%          answer = inputdlg(dlgString,'Auto-Pilot',1,{num2str(handles.prevAutoPilotStop)});
%     else
%          answer = inputdlg(dlgString,'Auto-Pilot');
%     end

%     if(isempty(answer))
%         % User Cancelled
%         return;
%     end
    
    % Check Valid Answer
    %stopFrame = round(str2num(answer{1}));
    %if(stopFrame > handles.currentLabel && stopFrame <= handles.numFrames)
   
    % Set Program State to Automation
    handles.state = 'automation';

    % Set Status Indicator
    set(handles.statusEdit,'String','Auto Pilot');
    handles.prevAutoPilotStop = stopFrame;

    % Set Automation Indicators
    set(handles.automationTypeEdit','String','Auto Pilot');
%     set(handles.automationStartEdit,'String',num2str(startFrame));
%     set(handles.automationStopEdit,'String',num2str(stopFrame));

    % Move to Frame No. startFrame
    handles.currFrame = startFrame;
    set(handles.currFrameEdit,'String',num2str(startFrame));
    displayFun(handles);
    % Label the First Frame Correct
    handles = labelCurr(handles,0);
    %set(handles.currFrameEdit,'String',num2str(handles.currFrame));
    set(handles.currFrameEdit,'String',num2str(startFrame));

    % Save guidata
    guidata(hObject,handles);

    % Get guidata to be used in loop
    % Note a local handles copy (h_loop) is used.  The handles structure is
    % updated at the end of this function.  This is done to avoid potential
    % problems with the handles structure updating when the "Stop Auto"
    % button is pressed.
    
    % set 'autoStopFlag' to zero 
    handles.autoStopFlag = 0;
    
    h_loop = guidata(hObject);   
    
    bLeft = 1; bTop = 1; nVecIdx = 0;
    displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx);
    handles.threshVec = handles.threshVec_ini;
    
    % Auto-Pilot Loop
    for k = startFrame+1:stopFrame
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

        %if(any(result))
        if(any(result))
            % Algorithm is Underconfident %
            % Display feature vector thresholding in command window
            disp(result);

            % Update Display with current frame
            displayFun(h_loop);

            % Display Dialog Asking for labelling
            figPos = get(h_loop.UGT,'Position');
            % btn = 'Yes';
            btn = lowAlgoConfDialog(figPos(1),figPos(2));
           
            if(strcmp(btn,'No'))
                % User Wants to Stop
                break;
            end
            
            if get(handles.autoThresholdCheckBox,'Value')==1
                bLeft = 1; 
                h_loop = setLowerThreshold(h_loop,fVec,bLeft);
            end
        end
        
        bLeft = 1; bTop = 1; nVecIdx = 1;
        displayThreshholdVecToEdits(h_loop,bLeft,bTop,nVecIdx);

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
        msgbox('Left Video Auto-Pilot Complete!','Auto-Pilot Complete!','modal');
        % Move to Frame No. stopFrame
        h_loop.currFrame = stopFrame + 1;
    else
        % Display Message Box Indicating Break
        msgbox('Left Video Auto-Pilot was stopped by user!','Auto-Pilot Complete!','modal');
        h_loop.currFrame = k;
    end

    % Set State to Wait For Label
    h_loop.state = 'waitForLabel';    
    
    set(h_loop.currFrameEdit,'String',num2str(h_loop.currFrame));

    handles = h_loop;

    % Update Display
    displayFun(handles);
    displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx);

    % Update Automation Indicators
%     set(handles.automationTypeEdit,'String','');
%     set(handles.automationStartEdit,'String','');
%     set(handles.automationStopEdit,'String','');

    % Update Status
    set(handles.statusEdit,'String','Waiting for Label');

%     else
%         msgbox('Incorrect End Frame Entered','Incorrect End Frame','error','modal');
%     end
    guidata(hObject,handles);
    
end

% if Right Only 
if iModelOperation == 1
    if(~strcmp(handles.state_Right,'waitForLabel'))
        return;
    end

    set(handles.radiobuttonLeftThreshold,'Value',0);
    set(handles.radiobuttonRightThreshold,'Value',1);
    
    handles.state_Right = 'automation';

    % Set Status Indicator
    set(handles.statusEdit_Right,'String','Auto Pilot');
    handles.prevAutoPilotStop_Right = stopFrame;

    % Set Automation Indicators
    % set(handles.automationTypeEdit,'String','Auto Pilot');
%     set(handles.automationStartEdit,'String',num2str(startFrame));
%     set(handles.automationStopEdit,'String',num2str(stopFrame));

    % Move to Frame No. startFrame
    handles.currFrame_Right = startFrame;
    set(handles.currFrameEdit_R,'String',num2str(startFrame));
    displayFun(handles,0);
    
    % Label the First Frame Correct
    handles = labelCurr(handles,0,0);
    set(handles.currFrameEdit_R,'String',num2str(startFrame));
    
    % Save guidata
    guidata(hObject,handles);

    % Get guidata to be used in loop
    % Note a local handles copy (h_loop) is used.  The handles structure is
    % updated at the end of this function.  This is done to avoid potential
    % problems with the handles structure updating when the "Stop Auto"
    % button is pressed.
    
    % set 'autoStopFlag' to zero 
    handles.autoStopFlag = 0;
    
    h_loop = guidata(hObject); 
    bLeft = 0; bTop = 1; nVecIdx = 1;
    displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx);
    handles.threshVec_Right = handles.threshVec_ini_Right;
    
    % Auto-Pilot Loop
    for k = startFrame+1:stopFrame
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
        
        % Auto pilot stop check
        if(handles.autoStopFlag == 1)
             set(handles.automationStartEdit,'String',num2str(k));
             set(handles.currFrameEdit_R,'String',num2str(k));
             set(handles.statusEdit_Right,'String','waitForLabel');
             break;
        end

        % Check Algorithm Confidence Parameters
        [result,fVec] = checkAlgoConf(h_loop,k,0);

        if(any(result))
            % Algorithm is Underconfident %
            % Display feature vector thresholding in command window
            disp(result);

            % Update Display with current frame
            displayFun(h_loop,0);

            % Display Dialog Asking for labelling
            figPos = get(h_loop.UGT,'Position');
            % btn = 'Yes';
            btn = lowAlgoConfDialog(figPos(1),figPos(2));

            if(strcmp(btn,'No'))
                % User Wants to Stop
                break;
            end
            
            if get(handles.autoThresholdCheckBox,'Value')==1
                bLeft = 0; 
                h_loop = setLowerThreshold(h_loop,fVec,bLeft);               
            end
        end
        
        bLeft = 0; bTop = 1; nVecIdx = 1;
        displayThreshholdVecToEdits(h_loop,bLeft,bTop,nVecIdx);

        % Algorithm is Confident %
        % Label frame as correct
        h_loop = labelCurr(h_loop,0,0);

        set(h_loop.currFrameEdit_R,'String',num2str(h_loop.currFrame_Right));

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
        msgbox('Right Video Auto-Pilot Complete!','Auto-Pilot Complete!','modal');
        h_loop.currFrame_Right = stopFrame + 1;
    else
        % Display Message Box Indicating Break
        msgbox('Right Video Auto-Pilot was stopped by user!','Auto-Pilot Complete!','modal');
        h_loop.currFrame_Right = k;
    end

    % Set State to Wait For Label
    h_loop.state_Right = 'waitForLabel';
    
    % Move to Frame No. stopFrame
    set(h_loop.currFrameEdit_R,'String',num2str(h_loop.currFrame_Right));

    handles = h_loop;

    % Update Display
    displayFun(handles,0);

    % Update Automation Indicators
%     set(handles.automationTypeEdit,'String','');
%     set(handles.automationStartEdit,'String','');
%     set(handles.automationStopEdit,'String','');

    % Update Status
    set(handles.statusEdit_Right,'String','Waiting for Label');

%     else
%         msgbox('Incorrect End Frame Entered','Incorrect End Frame','error','modal');
%     end
    guidata(hObject,handles);
    
end
         
% if Both
if iModelOperation == 2
    if(~strcmp(handles.state,'waitForLabel')||~strcmp(handles.state_Right,'waitForLabel'))
        return;
    end
    
    set(handles.radiobuttonLeftThreshold,'Value',1);
    set(handles.radiobuttonRightThreshold,'Value',0);

    handles.state = 'automation';
    handles.state_Right = 'automation';

    % Set Status Indicator
    set(handles.statusEdit,'String','Auto Pilot');
    handles.prevAutoPilotStop = stopFrame;
    set(handles.statusEdit_Right,'String','Auto Pilot');
    handles.prevAutoPilotStop_Right = stopFrame;

    % Move to Frame No. startFrame
    handles.currFrame = startFrame;
    set(handles.currFrameEdit,'String',num2str(startFrame));
    displayFun(handles);
    handles.currFrame_Right = startFrame;
    set(handles.currFrameEdit_R,'String',num2str(startFrame));
    displayFun(handles,0);
    
    % Label the First Frame Correct
    handles = labelCurr(handles,0);
    set(handles.currFrameEdit,'String',num2str(startFrame));
    handles = labelCurr(handles,0,0);
    set(handles.currFrameEdit_R,'String',num2str(startFrame));
    
    % Save guidata
    guidata(hObject,handles);

    % Get guidata to be used in loop
    % Note a local handles copy (h_loop) is used.  The handles structure is
    % updated at the end of this function.  This is done to avoid potential
    % problems with the handles structure updating when the "Stop Auto"
    % button is pressed.
    
    % set 'autoStopFlag' to zero 
    handles.autoStopFlag = 0;    
    h_loop = guidata(hObject);    
    
    bLeft = 1; bTop = 1; nVecIdx = 0;
    displayThreshholdVecToEdits(h_loop,bLeft,bTop,nVecIdx);
    h_loop.threshVec = h_loop.threshVec_ini;
    h_loop.threshVec_Right = h_loop.threshVec_ini_Right;
    
    % Auto-Pilot Loop
    for k = startFrame+1:stopFrame        
        
        % Auto pilot stop check
        if(handles.autoStopFlag == 1)
             set(h_loop.automationStartEdit,'String',num2str(k));
             set(h_loop.currFrameEdit,'String',num2str(k));
             set(h_loop.statusEdit,'String','waitForLabel');
             set(h_loop.currFrameEdit_R,'String',num2str(k));
             set(h_loop.statusEdit_Right,'String','waitForLabel');
             break;
        end

        set(handles.radiobuttonLeftThreshold,'Value',1);
        set(handles.radiobuttonRightThreshold,'Value',0);
        
        % display left ini
        bLeft = 1; bTop = 1; nVecIdx = 0;
        displayThreshholdVecToEdits(h_loop,bLeft,bTop,nVecIdx);
        
        % Check Algorithm Confidence Parameters
        [result,fVec] = checkAlgoConf(h_loop,k);
        if(any(result))
            % Algorithm is Underconfident %
            % Display feature vector thresholding in command window
            disp(result);

            % Update Display with current frame
            displayFun(h_loop);

            % Display Dialog Asking for labelling
            figPos = get(h_loop.UGT,'Position');
            % btn = 'Yes';
            btn = lowAlgoConfDialog(figPos(1),figPos(2));

            if(strcmp(btn,'No'))
                % User Wants to Stop
                break;
            end
            
            if get(handles.autoThresholdCheckBox,'Value')==1
                bLeft = 1; 
                h_loop = setLowerThreshold(h_loop,fVec,bLeft);            
            end
        end
        
        % display left relsult
        bLeft = 1; bTop = 1; nVecIdx = 0;
        displayThreshholdVecToEdits(h_loop,bLeft,bTop,nVecIdx);
        
        h_loop = labelCurr(h_loop,0);       
        set(h_loop.currFrameEdit,'String',num2str(h_loop.currFrame));        
        
        set(handles.radiobuttonLeftThreshold,'Value',0);
        set(handles.radiobuttonRightThreshold,'Value',1);
    
        % display right ini
        bLeft = 0; bTop = 1; nVecIdx = 0;
        displayThreshholdVecToEdits(h_loop,bLeft,bTop,nVecIdx);
        
        % Stop Auto-Pilot if user presses Stop Automation %
        
        [result_R,fVec_R] = checkAlgoConf(h_loop,k,0);

        if(any(result_R))
            % Algorithm is Underconfident %
            % Display feature vector thresholding in command window
            disp(result_R);

            % Update Display with current frame
            displayFun(h_loop,0);

            % Display Dialog Asking for labelling
            figPos = get(h_loop.UGT,'Position');
            % btn = 'Yes';
            btn = lowAlgoConfDialog(figPos(1),figPos(2));

            if(strcmp(btn,'No'))
                % User Wants to Stop
                break;
            end
            
            if get(handles.autoThresholdCheckBox,'Value')==1
                bLeft = 0; 
                h_loop = setLowerThreshold(h_loop,fVec_R,bLeft);             
            end
        end
        
        % display right relsult
        bLeft = 0; bTop = 1; nVecIdx = 0;
        displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx);

        % Algorithm is Confident %
        % Label frame as correct
        h_loop = labelCurr(h_loop,0,0);
        set(h_loop.currFrameEdit_R,'String',num2str(h_loop.currFrame_Right));

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
        msgbox('Both Video Auto-Pilot Complete!','Auto-Pilot Complete!','modal');
        h_loop.currFrame = stopFrame + 1;
        h_loop.currFrame_Right = stopFrame + 1;
    else
        % Display Message Box Indicating Break
        msgbox('Both Video Auto-Pilot was stopped by user!','Auto-Pilot Complete!','modal');
        h_loop.currFrame = k;
        h_loop.currFrame_Right = k;
    end

    % Set State to Wait For Label
    h_loop.state = 'waitForLabel';
    h_loop.state_Right = 'waitForLabel';
    
    % Move to Frame No. stopFrame
    set(h_loop.currFrameEdit,'String',num2str(h_loop.currFrame));
    set(h_loop.currFrameEdit_R,'String',num2str(h_loop.currFrame_Right));

    handles = h_loop;

    % Update Display
    displayFun(handles);
    displayFun(handles,0);
    bLeft = 1; bTop = 1;nVecIdx = 1;
    displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx);

    % Update Status
    set(handles.statusEdit_Right,'String','Waiting for Label');
    guidata(hObject,handles);
    
end 


    
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
theSelLabel = get(handles.labelSelComboEdit,'Value');
if theSelLabel == 1 % Correct
    theLabelToSet = 0;
end
if theSelLabel== 2 % Error
    theLabelToSet = 4;
end

% Get Start Frame from 'automationStartEdit'
startFrame = getNumfromEdit(handles.automationStartEdit,'Automation From Frame Edit');
% Get End Frame from 'automationStopEdit'
stopFrame = getNumfromEdit(handles.automationStopEdit,'Automation End Frame Edit');
errorType = theLabelToSet;

% Check startFrame & stopFrame
if startFrame<1 || startFrame>handles.numFrames || startFrame>stopFrame 
    msgbox('Incorrect Start Frame Entered','Incorrect Start Frame','error','modal'); 
    return;
end

if startFrame>stopFrame || stopFrame<1 || stopFrame>handles.numFrames 
    msgbox('Incorrect End Frame Entered','Incorrect End Frame','error','modal'); 
    return;
end

iModelOperation = getSelObject(handles);
if iModelOperation == 0
    if(~strcmp(handles.state,'waitForLabel'))
        return;
    end
    if iModelOperation == 2
        if(~strcmp(handles.state_Right,'waitForLabel'))
            return;
        end
    end
    
    % startFrame = handles.currentLabel;
    

    % Get End Frame and Error Type from user using modal dialog
%     multiErrorInfo = multiErrorDlg('StartFrame',startFrame);
%     stopFrame = multiErrorInfo(1); 
%     errorType = multiErrorInfo(2);
    
    

    % Check for user cancel
%     if(all(multiErrorInfo == 0))
%         return;
%     end

    % Check Valid stopFrame
     % if(stopFrame > handles.currentLabel && stopFrame <= handles.numFrames)

    % Note : Multi Error ignores the Stop Auto Button.  This is because it
    % labels all frames but the last without running the detector.
    % Therefore it should not take long to run.

    % Set Program State to Automation
    handles.state = 'automation';
    guidata(hObject,handles);

    % Set Status Indicator
    set(handles.statusEdit,'String','Batch Labeling (automation)');

    % Set Automation Indicators
%     set(handles.automationTypeEdit','String','waitForLabel');
%     set(handles.automationStartEdit,'String',num2str(startFrame));
%     set(handles.automationStopEdit,'String',num2str(stopFrame));

    % Save guidata
    guidata(hObject,handles);

    % Get guidata to be used in loop
    % Note, a local handles copy (h_loop) is used.  The handles structure
    % is updated at the end of this function.  This is done to avoid
    % potential problems with the handles structure updating when the "Stop
    % Auto" button is pressed.
    h_loop = guidata(hObject);

    % Run Multi Error %
    h_loop.currFrame = startFrame;
    displayFun(h_loop);

    if errorType == 0
        for k = startFrame:stopFrame
            h_loop = labelCurr(h_loop,errorType);
            set(h_loop.currFrameEdit,'String',num2str(k));

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
    else
        % Label all but stopFrame without running the detector %
        if(stopFrame >= startFrame + 1)
            h_loop.label(startFrame:(stopFrame-1)) = errorType;
            h_loop.currentLabel = stopFrame;
            h_loop.currFrame = stopFrame;
        end
        % Label stopFrame with labelCurr so detector is run on next frame
        h_loop = labelCurr(h_loop,errorType);
        set(h_loop.currFrameEdit,'String',num2str(h_loop.currFrame));
    end
    
    

    % Set State to Wait For Label
    h_loop.state = 'waitForLabel';

    % This is the structure that will be used to update guidata
    handles = h_loop;

    % Update Display
    displayFun(handles);

    % Update Automation Indicators
%     set(handles.automationTypeEdit,'String','');
%     set(handles.automationStartEdit,'String','');
%     set(handles.automationStopEdit,'String','');

    % Display Message Box Indicating Completion
    msgbox('Left Video Batch Label Complete!','Left Video Batch Label Complete!','modal');

    % Update Status
    set(handles.statusEdit,'String','Waiting for Label');
    %}

%     else
%         msgbox('Incorrect End Frame Entered','Incorrect End Frame','error','modal'); 
%     end
    % Update guidata
    guidata(hObject,handles);
end

if iModelOperation == 1
    if(~strcmp(handles.state_Right,'waitForLabel'))
        return;
    end
    
    handles.state_Right = 'automation';
    guidata(hObject,handles);

    % Set Status Indicator
    set(handles.statusEdit_Right,'String','Batch Labeling (Automation)');

    % Set Automation Indicators
%     set(handles.automationTypeEdit','String','Multi Error');
%     set(handles.automationStartEdit,'String',num2str(startFrame));
%     set(handles.automationStopEdit,'String',num2str(stopFrame));

    % Save guidata
    guidata(hObject,handles);

    % Get guidata to be used in loop
    % Note, a local handles copy (h_loop) is used.  The handles structure
    % is updated at the end of this function.  This is done to avoid
    % potential problems with the handles structure updating when the "Stop
    % Auto" button is pressed.
    h_loop = guidata(hObject);

    % Run Multi Error %
    h_loop.currFrame_Right = startFrame;
    displayFun(h_loop,0);

    if errorType == 0;
        for k = startFrame:stopFrame
            h_loop = labelCurr(h_loop,errorType,0);
            set(h_loop.currFrameEdit_R,'String',num2str(k));

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
    else
        % Label all but stopFrame without running the detector %
        if(stopFrame >= startFrame + 1)
            h_loop.label_Right(startFrame:(stopFrame-1)) = errorType;
            h_loop.currentLabel_Right = stopFrame;
            h_loop.currFrame_Right = stopFrame;
        end
        % Label stopFrame with labelCurr so detector is run on next frame
        h_loop = labelCurr(h_loop,errorType,0);
        set(h_loop.currFrameEdit_R,'String',num2str(h_loop.currFrame_Right));
    end

    % Set State to Wait For Label
    h_loop.state_Right = 'waitForLabel';

    % This is the structure that will be used to update guidata
    handles = h_loop;

    % Update Display
    displayFun(handles,0);

    % Update Automation Indicators
%     set(handles.automationTypeEdit,'String','');
%     set(handles.automationStartEdit,'String','');
%     set(handles.automationStopEdit,'String','');

    % Display Message Box Indicating Completion
    msgbox('Right Video Batch Label Complete!','Right Video Batch Label Complete!','modal');

    % Update Status
    set(handles.statusEdit_Right,'String','Waiting for Label');
    %}

    guidata(hObject,handles);

end

if iModelOperation == 2
    if(~strcmp(handles.state,'waitForLabel')||~strcmp(handles.state_Right,'waitForLabel'))
        return;
    end
   
    handles.state = 'automation';
    handles.state_Right = 'automation';
    guidata(hObject,handles);

    % Set Status Indicator
    set(handles.statusEdit,'String','Batch Labeling (Automation)');
    set(handles.statusEdit_Right,'String','Batch Labeling (Automation)');

    % Save guidata
    guidata(hObject,handles);

    % Get guidata to be used in loop
    % Note, a local handles copy (h_loop) is used.  The handles structure
    % is updated at the end of this function.  This is done to avoid
    % potential problems with the handles structure updating when the "Stop
    % Auto" button is pressed.
    h_loop = guidata(hObject);

    % Run Multi Error %
    h_loop.currFrame = startFrame;
    displayFun(h_loop);
    h_loop.currFrame_Right = startFrame;
    displayFun(h_loop,0);

    if errorType == 0
        for k = startFrame:stopFrame
            h_loop = labelCurr(h_loop,errorType);
            set(h_loop.currFrameEdit,'String',num2str(k));
            h_loop = labelCurr(h_loop,errorType,0);
            set(h_loop.currFrameEdit_R,'String',num2str(k));

            % Stop Have Faith if user presses Stop Automation %
            % Get Latest Handles
            h_latest = guidata(hObject);

            % Check Stop Flag
            if(h_latest.autoStopFlag == 1)
                % Stop Pressed.  Break
                h_loop.autoStopFlag = 0;
                break;
            end 

            displayFun(h_loop);
            displayFun(h_loop,0);
        end
    else
        % Label all but stopFrame without running the detector %
        if(stopFrame >= startFrame + 1)
            h_loop.label(startFrame:(stopFrame-1)) = errorType;
            h_loop.currentLabel = stopFrame;
            h_loop.currFrame = stopFrame;
            h_loop.label_Right(startFrame:(stopFrame-1)) = errorType;
            h_loop.currentLabel_Right = stopFrame;
            h_loop.currFrame_Right = stopFrame;
        end
    % Label stopFrame with labelCurr so detector is run on next frame
        h_loop = labelCurr(h_loop,errorType);
        set(h_loop.currFrameEdit,'String',num2str(h_loop.currFrame));
        h_loop = labelCurr(h_loop,errorType,0);
        set(h_loop.currFrameEdit_R,'String',num2str(h_loop.currFrame_Right));
    end

    % Set State to Wait For Label
    h_loop.state = 'waitForLabel';
    h_loop.state_Right = 'waitForLabel';

    % This is the structure that will be used to update guidata
    handles = h_loop;

    % Update Display
    displayFun(handles);
    displayFun(handles,0);

    % Update Automation Indicators
%     set(handles.automationTypeEdit,'String','');
%     set(handles.automationStartEdit,'String','');
%     set(handles.automationStopEdit,'String','');

    % Display Message Box Indicating Completion
    msgbox('Both Video Batch Label Complete!','Right Video Batch Label Complete!','modal');

    % Update Status
    set(handles.statusEdit,'String','Waiting for Label');
    set(handles.statusEdit_Right,'String','Waiting for Label');
    %}

    guidata(hObject,handles);

end


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



function rightVideoFileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to rightVideoFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rightVideoFileEdit as text
%        str2double(get(hObject,'String')) returns contents of rightVideoFileEdit as a double


% --- Executes during object creation, after setting all properties.
function rightVideoFileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rightVideoFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RightAnalysisFileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to RightAnalysisFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RightAnalysisFileEdit as text
%        str2double(get(hObject,'String')) returns contents of RightAnalysisFileEdit as a double


% --- Executes during object creation, after setting all properties.
function RightAnalysisFileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RightAnalysisFileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ConnectedWithRightCheckbox.
function ConnectedWithRightCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to ConnectedWithRightCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ConnectedWithRightCheckbox


%--------------------------------------------------------------------------
% --- Executes on button press in PrcsRightButton.
function PrcsRightButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrcsRightButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Prompt User for Video File Location

% idx = getSelObject(handles)
% theNum = getNumfromEdit(handles.automationStartEdit)
% displayThreshholdVecToEdits(handles,0,2,[11:17]);
% updateThreshVecFromEdits(handles)
get(handles.autoThresholdCheckBox,'Value')
return;

currVidDir_R = getCurrVidDir(handles);
if(currVidDir_R ~= 0)
    % Open File Browser in Current Video Directory
    [vidName,vidPath] = uigetfile('*.avi','Load Video File',currVidDir);
else
    [vidName,vidPath] = uigetfile('*.avi','Load Video File');
end

if(vidName == 0); 
    return; 
end;

[backName,backPath] = uigetfile('*.avi','Load Background File',vidPath);
if(backName == 0); 
    return; 
end;

% Load Video & Store Video Parameters
% different to left, they are for right video
handles.vidName_Right   = fullfile(vidPath,vidName);
handles.backName_Right  = fullfile(backPath,backName);
handles.vidObj_Right    = VideoReader(handles.vidName_Right);
handles.backObj_Right   = VideoReader(handles.backName_Right);
handles.imSize_Right    = [handles.vidObj_Right.Height,handles.vidObj_Right.Width];
nFrames = handles.vidObj_Right.NumberOfFrames;
handles.numFrames_Right = nFrames;

% Store Backframe
handles.backIm_Right = rgb2gray(handles.backObj_Right.read(1));

% Display Video File Info
set(handles.rightVideoFileEdit,'String',handles.vidName_Right);

% Display Variables
handles.currFrame_Right = 1;

% Set Total Frame Edit
set(handles.totalFramesEdit_R,'String',num2str(nFrames));

% Initialize all non-video data structures.  These structures are used to
% store tracking data and the tracking file name
handles = initNonVideoStructures(handles);
guidata(hObject,handles);

% Show First Frame
displayFun(handles);


% --- Executes on button press in previousFrameButton_R.
function previousFrameButton_R_Callback(hObject, eventdata, handles)
% hObject    handle to previousFrameButton_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iModeOperation = getSelObject(handles);
if  iModeOperation == 0 % Left only mode
    return;
end

if(strcmp(handles.state_Right,{'automation','idle'}))
    % Do Nothing if tracking file has not been started or
    % automation is running
    return;
end
if  iModeOperation == 2 % Both & Synchro mode
    if(strcmp(handles.state,{'automation','idle'}))
        % Do Nothing if tracking file has not been started or
        % automation is running
        return;
    end
end

if(handles.currFrame_Right > 1)
    handles.currFrame_Right = handles.currFrame_Right - 1;
    displayFun(handles,0);
    guidata(hObject,handles);
end

if  iModeOperation == 2 % Both & Synchro mode
    if(handles.currFrame > 1)
        handles.currFrame = handles.currFrame - 1;
        displayFun(handles);
        guidata(hObject,handles);
    end
end


% --- Executes on button press in nextFrameButton_R.
function nextFrameButton_R_Callback(hObject, eventdata, handles)
% hObject    handle to nextFrameButton_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iModeOperation = getSelObject(handles);
if  iModeOperation == 0 % Left only mode
    return;
end

if(strcmp(handles.state_Right,{'automation','idle'}))
    % Do Nothing if tracking file has not been started or
    % automation is running
    return;
end
if  iModeOperation == 2 % Both & Synchro mode
    if(strcmp(handles.state,{'automation','idle'}))
        % Do Nothing if tracking file has not been started or
        % automation is running
        return;
    end
end

if(handles.currFrame_Right < handles.numFrames_Right)
    handles.currFrame_Right = handles.currFrame_Right + 1;
    displayFun(handles,0);
    guidata(hObject,handles);
end

if  iModeOperation == 2 % Both & Synchro mode
    if(handles.currFrame < handles.numFrames)
        handles.currFrame = handles.currFrame + 1;
        displayFun(handles);
        guidata(hObject,handles);
    end
end


function currFrameEdit_R_Callback(hObject, eventdata, handles)
% hObject    handle to currFrameEdit_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currFrameEdit_R as text
%        str2double(get(hObject,'String')) returns contents of currFrameEdit_R as a double
% Do Nothing if idle or automation
if(strcmp(handles.state_Right,{'idle','automation'}))
    return;
end

frameNum = round(str2double(get(hObject,'String')));
if(frameNum >= 1 && frameNum <= handles.numFrames_Right)
    handles.currFrame_Right = frameNum;
    set(handles.currFrameEdit_R,'String',num2str(frameNum));
    displayFun(handles,0);
else
    set(handles.currFrameEdit_R,'String',num2str(handles.currFrame_Right));
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function currFrameEdit_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currFrameEdit_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in currentLabelButton_R.
function currentLabelButton_R_Callback(hObject, eventdata, handles)
% hObject    handle to currentLabelButton_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~strcmp(handles.state_Right,'waitForLabel'))
    return;
end

% Go to frame that needs to be labeled
handles.currFrame_Right = handles.currentLabel_Right;

% Display Frame
displayFun(handles,0);
guidata(hObject,handles);

function totalFramesEdit_R_Callback(hObject, eventdata, handles)
% hObject    handle to totalFramesEdit_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of totalFramesEdit_R as text
%        str2double(get(hObject,'String')) returns contents of totalFramesEdit_R as a double


% --- Executes during object creation, after setting all properties.
function totalFramesEdit_R_CreateFcn(hObject, eventdata, handles)
% hObject    handle to totalFramesEdit_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in CommentListbox.
function CommentListbox_Callback(hObject, eventdata, handles)
% hObject    handle to CommentListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CommentListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CommentListbox


% --- Executes during object creation, after setting all properties.
function CommentListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CommentListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- this function get the status of the "Object Select " button group,
% return the selected group.
function [idx] = getSelObject(handles)
% handles    structure with handles and user data (see GUIDATA)
if (get(handles.radiobuttonRightOnly,'Value'))==1
    idx = 1;
elseif (get(handles.radiobuttonBoth,'Value'))==1
    idx = 2;
else 
    idx = 0;
end


% --- Executes on button press in LVideoSelBtn.
function LVideoSelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to LVideoSelBtn (see GCBO)
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


% --- Executes on button press in RVideoSelBtn.
function RVideoSelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to RVideoSelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Prompt User for Video File Location
% idx = getSelObject(handles)
% return;

currVidDir_R = getCurrVidDir(handles);
if(currVidDir_R ~= 0)
    % Open File Browser in Current Video Directory
    [vidName,vidPath] = uigetfile('*.avi','Load Right Video File',currVidDir_R);
else
    [vidName,vidPath] = uigetfile('*.avi','Load Right Video File');
end

if(vidName == 0); 
    return; 
end;

[backName,backPath] = uigetfile('*.avi','Load Background File',vidPath);
if(backName == 0); 
    return; 
end;

% Load Video & Store Video Parameters
% different to left, they are for right video
handles.vidName_Right   = fullfile(vidPath,vidName);
handles.backName_Right  = fullfile(backPath,backName);
handles.vidObj_Right    = VideoReader(handles.vidName_Right);
handles.backObj_Right   = VideoReader(handles.backName_Right);
handles.imSize_Right    = [handles.vidObj_Right.Height,handles.vidObj_Right.Width];
nFrames = handles.vidObj_Right.NumberOfFrames;
handles.numFrames_Right = nFrames;

% Store Backframe
handles.backIm_Right = rgb2gray(handles.backObj_Right.read(1));

% Display Video File Info
set(handles.rightVideoFileEdit,'String',handles.vidName_Right);

% Display Variables
handles.currFrame_Right = 1;

% Set Total Frame Edit
set(handles.totalFramesEdit_R,'String',num2str(nFrames));

% Initialize all non-video data structures.  These structures are used to
% store tracking data and the tracking file name
handles = initNonVideoStructures(handles,0);
guidata(hObject,handles);

% Show First Frame
displayFun(handles,0);



% --- Executes on button press in LTrackFileSelBtn.
function LTrackFileSelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to LTrackFileSelBtn (see GCBO)
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

msgbox('Left video data has been saved','Data is Saved');


% --- Executes on button press in RTrackFileSelBtn.
function RTrackFileSelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to RTrackFileSelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(strcmp(handles.state_Right,'idle'))
    return;
end

if(isempty(handles.analysisFileName_Right))
    currVidDir_R = getCurrVidDir(handles);
    if(currVidDir_R ~= 0)
        % Open File Browser in current video directory
        [fName,pName] = uiputfile('*.mat','Save Tracking File',currVidDir_R);
    else
        [fName,pName] = uiputfile('*.mat','Save Analysis File');
    end
    
    % Check for User Cancel

    if(fName == 0)
        return;
    end
    
    handles.analysisFileName_Right = fullfile(pName,fName);
end

% The following fields had been set in function 'iniNonVideoStructure'
handlesForSave.mark = handles.mark_Right;
handlesForSave.pc = handles.pc_Right;
handlesForSave.neighborhoodSize = handles.neighborhoodSize_Right;
handlesForSave.templateFrame = handles.templateFrame_Right;
handlesForSave.inst = handles.inst_Right;

%
handlesForSave.algoInfo = handles.algoInfo;
% The following fields had been set in function 'Select Right Video Button' 
handlesForSave.vidName = handles.vidName_Right;
handlesForSave.backName = handles.backName_Right;
handlesForSave.imSize = handles.imSize_Right;
handlesForSave.numFrames = handles.numFrames_Right;
handlesForSave.backIm = handles.backIm_Right;
% state_Right has been initialized in Open Fcn
handlesForSave.state = handles.state_Right;
handlesForSave.currentLabel = handles.currentLabel_Right;
handlesForSave.label = handles.label_Right;
handlesForSave.currFrame = handles.currFrame_Right;
handlesForSave.analysisFileName = handles.analysisFileName_Right;
% saveFileds must modified
saveFields = {'mark','pc','neighborhoodSize','templateFrame','inst',...
'algoInfo','vidName','backName','imSize','numFrames',...
'backIm','state','currentLabel','label','currFrame','analysisFileName'};

saveStruct = handlesForSave;
allFields = fieldnames(saveStruct);

% Remove all non save fields
for k = 1:numel(allFields)
    cName = allFields{k};
    if(~any(strcmp(cName,saveFields)))
        saveStruct = rmfield(saveStruct,cName);
    end
end

save(handles.analysisFileName_Right, 'saveStruct');
set(handles.RightAnalysisFileEdit,'String',handles.analysisFileName_Right);
guidata(hObject,handles);

msgbox('Right video data has been saved','Data is Saved');



function labelEdit_Right_Callback(hObject, eventdata, handles)
% hObject    handle to labelEdit_Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of labelEdit_Right as text
%        str2double(get(hObject,'String')) returns contents of labelEdit_Right as a double


% --- Executes during object creation, after setting all properties.
function labelEdit_Right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labelEdit_Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function previousFrameLabelEdit_Right_Callback(hObject, eventdata, handles)
% hObject    handle to previousFrameLabelEdit_Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of previousFrameLabelEdit_Right as text
%        str2double(get(hObject,'String')) returns contents of previousFrameLabelEdit_Right as a double


% --- Executes during object creation, after setting all properties.
function previousFrameLabelEdit_Right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to previousFrameLabelEdit_Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function statusEdit_Right_Callback(hObject, eventdata, handles)
% hObject    handle to statusEdit_Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of statusEdit_Right as text
%        str2double(get(hObject,'String')) returns contents of statusEdit_Right as a double


% --- Executes during object creation, after setting all properties.
function statusEdit_Right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statusEdit_Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in initializeTemplateButton_Right.
function initializeTemplateButton_Right_Callback(hObject, eventdata, handles)
% hObject    handle to initializeTemplateButton_Right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(any(strcmp(handles.state_Right,{'idle','automation'})))
    return;
end

iModeOperation = getSelObject(handles);
if  iModeOperation == 1 || iModeOperation == 2% Right only mode or Both & Synchro mode
    % If waiting for label, verify with user that re-initializing the template
    % will result in a loss of tracking data.
    if(strcmp(handles.state_Right,'waitForLabel'))
        % Prompt User
        msgAns = questdlg('Re-initializing the right template will close the current right tracking file.  Any unsaved progress will be lost.  Continue?',...
        'Re-initialize Template','Yes','No','No');

        % Do nothing if user does not want to reinitialize template
        if(~strcmp('Yes',msgAns))
            return;
        end

        % Clear Right Analysis File Name
        handles.analysisFileName_Right = [];
        set(handles.RightAnalysisFileEdit,'String','');

        % Clear Data Structures
        handles = initNonVideoStructures(handles,0);
        guidata(hObject,handles);
        displayFun(handles,0);

    end

    currFrame = handles.currFrame_Right;

    % Interactively Get Template Window
    defStartSize = 99; % Default Window Starting Size
    frameIm = rgb2gray(handles.vidObj_Right.read(handles.currFrame_Right));
    [centerPt,tSize] = squareDraw(frameIm,defStartSize);

    % Check for user cancel
    if(isempty(centerPt))
        errordlg('Cancelled Template Initializiation Process','Cancelled Init','modal');
        return;
    end

    % Upper Left Corner of Template
    tempWindow = centerPt - ((tSize-1)/2);

    nSize = round(tSize * handles.neighborhoodSize_Right)-1;
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
    handles.pc_Right = wMatchBuildPcStruct(double(temp),weights);
    handles.mark_Right.tCorn(currFrame,:) = tempWindow;
    handles.mark_Right.subT(currFrame,:) = tempWindow;
    handles.mark_Right.wStats = wStats;
    handles.mark_Right.nCorn(currFrame,:) = nCurr;
    handles.templateFrame_Right = currFrame;

    % Template Related Algo Info
    handles.algoInfo_Right.nccScore(currFrame,1) = 1;
    handles.algoInfo_Right.fitMSE(currFrame,1) = 0;

    % Estimate Boundary Lines
    [r,t,dbl_ais] = detection_bLines(frameIm,handles.backIm_Right,tempWindow,wStats,handles.pc_Right);
    handles.inst_Right.rho(currFrame,:) = r;
    handles.inst_Right.theta(currFrame,:) = t;
    handles.inst_Right.trackPt(currFrame,:) = computeTrackPt(r,t,tempWindow,...
                                                        round((tSize-1)/2));
    % Boundary Line Algo Info
    handles.algoInfo_Right.votesLeft(currFrame,:) = dbl_ais.votesLeft;
    handles.algoInfo_Right.votesRight(currFrame,:) = dbl_ais.votesRight;
    handles.algoInfo_Right.leftNumInliers(currFrame,:) = dbl_ais.leftNumInliers;
    handles.algoInfo_Right.leftRefitMSE(currFrame,:) = dbl_ais.leftRefitMSE;
    handles.algoInfo_Right.rightNumInliers(currFrame,:) = dbl_ais.rightNumInliers;
    handles.algoInfo_Right.rightRefitMSE(currFrame,:) = dbl_ais.rightRefitMSE;

    % Set State to waitForLabel
    handles.currentLabel_Right = currFrame;
    handles.state_Right = 'waitForLabel';
    set(handles.statusEdit_Right,'String','Waiting for Label');

    displayFun(handles,0);
    guidata(hObject,handles);
end


% --- Executes on selection change in labelSelComboEdit.
function labelSelComboEdit_Callback(hObject, eventdata, handles)
% hObject    handle to labelSelComboEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns labelSelComboEdit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from labelSelComboEdit


% --- Executes during object creation, after setting all properties.
function labelSelComboEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labelSelComboEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [theNum] = getNumfromEdit(theEditHandle,editName)
if nargin<2
    editName = 'unkown';
end

strIn = get(theEditHandle,'String');

theNum = str2num(strIn);
if isempty(theNum) 
    msgbox(['Edit "' editName '" needs a number !' ],'Edit Error!','modal');
    theNum = -99.99;
end

function [rel] = autoDetect()


% --- Executes on button press in setToCurrButton.
function setToCurrButton_Callback(hObject, eventdata, handles)
% hObject    handle to setToCurrButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iModeOperation = getSelObject(handles);
if  iModeOperation == 0 || iModeOperation == 2% Left only mode or Both & Synchro mode
    set(handles.automationStartEdit,'String',num2str(handles.currFrame));
end
if  iModeOperation == 1% Right only mode
    set(handles.automationStartEdit,'String',num2str(handles.currFrame_Right));
end


% --- Executes on button press in loadTrackFileBtn.
function loadTrackFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to loadTrackFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iModeOperation = getSelObject(handles);
if  iModeOperation == 0 || iModeOperation == 2% Left only mode or Both & Synchro mode
    currTrackFileDir = getCurrTrackFileDir(handles);
    if(currTrackFileDir ~= 0)
        % Open File Browser in Current Track File Directory
        [fName,pName] = uigetfile('*.mat','Select Tracking File',currTrackFileDir);
    else
        % No Tracking File Avaialble.  Check Video File
        currVidDir = getCurrVidDir(handles);
        if(currVidDir ~= 0)
            [fName,pName] = uigetfile('*.mat','Select Tracking File',currVidDir);
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
end

if  iModeOperation == 1 || iModeOperation == 2% Right only mode
    currTrackFileDir_R = getCurrTrackFileDir(handles,0);
    if(currTrackFileDir_R ~= 0)
        % Open File Browser in Current Track File Directory
        [fName,pName] = uigetfile('*.mat','Select Tracking File',currTrackFileDir_R);
    else
        % No Tracking File Avaialble.  Check Video File
        currVidDir_R = getCurrVidDir(handles,0);
        if(currVidDir_R ~= 0)
            [fName,pName] = uigetfile('*.mat','Select Tracking File',currVidDir_R);
        else
            [fName,pName] = uigetfile('*.mat','Select Tracking File');
        end
    end

    if(fName == 0)
        return;
    end

    % handles_BackUp = handles;
    
    load(fullfile(pName,fName),'saveStruct');

    saveFields = {'mark','pc','neighborhoodSize','templateFrame','inst',...
    'algoInfo','vidName','backName','imSize','numFrames',...
    'backIm','state','currentLabel','label','currFrame','analysisFileName'};

    loadFields = {'mark_Right','pc_Right','neighborhoodSize_Right','templateFrame_Right','inst_Right',...
    'algoInfo_Right','vidName_Right','backName_Right','imSize_Right','numFrames_Right',...
    'backIm_Right','state_Right','currentLabel_Right','label_Right','currFrame_Right','analysisFileName_Right'};

    % Load all Save Fields
    for k = 1:numel(saveFields)
        handles.(loadFields{k}) = saveStruct.(saveFields{k});
    end

    % Setup Video Objects
    handles.vidObj_Right = VideoReader(handles.vidName_Right);
    handles.backObj_Right = VideoReader(handles.backName_Right);
    displayFun(handles,0);

    % Set Axis Dimensions to that of Image
    set(handles.dispAxesRight,'Units','pixels');
    pos = get(handles.dispAxesRight,'Position');
    pos(3:4) = handles.imSize_Right([2,1]);
    set(handles.dispAxesRight,'Position',pos);

    % Display Video File & Analysis File Names
    set(handles.rightVideoFileEdit,'String',handles.vidName_Right);
    set(handles.RightAnalysisFileEdit,'String',handles.analysisFileName_Right);

    % Set Label to waitForLabel
    handles.state_Right = 'waitForLabel';
    set(handles.statusEdit_Right,'String','Waiting for Label');

    guidata(hObject,handles);
    
end



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
% This function display the Threshold vector
function displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx,vecIn)
% handles    structure with handles and user data (see GUIDATA)
% bLeft :  0: display right Threshold vector; 
%          1: display left Threshold vector;
% bTop  :  0: display the Threshold vector to bottom edit line; 
%          1: display the Threshold vector to top edit line;
% nVecIdx: 0: display threshVec_ini or threshVec_ini_Right
%          1: display threshVec or threshVec_Right
%          2: display input vec
% vecIn :  the input vector to display
if nargin<4
    warning('Threshold edits display error. Input is not enough.');
    return;
end
if nargin==4
    if nVecIdx > 1
        warning('Threshold edits display error. Input is not enough.');
        return;
    end
end

if nVecIdx==0
    if bLeft == 1
        theVec = handles.threshVec_ini;
    else
        theVec = handles.threshVec_ini_Right;
    end
end

if nVecIdx==1 
    if bLeft == 1
        theVec = handles.threshVec;
    else
        theVec = handles.threshVec_Right;
    end
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
function [theVec] = updateThreshVecFromEdits(handles)

theVec(1) = str2num(get(handles.usedThreshVecEdit01,'String'));
theVec(2) = str2num(get(handles.usedThreshVecEdit02,'String'));
theVec(3) = str2num(get(handles.usedThreshVecEdit03,'String'));
theVec(4) = str2num(get(handles.usedThreshVecEdit04,'String'));
theVec(5) = str2num(get(handles.usedThreshVecEdit05,'String'));
theVec(6) = str2num(get(handles.usedThreshVecEdit06,'String'));
theVec(7) = str2num(get(handles.usedThreshVecEdit07,'String'));
% handles.threshVec = theVec;

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


% --- Executes on button press in applyThresholdToVedioBtn.
function applyThresholdToVedioBtn_Callback(hObject, eventdata, handles)
% hObject    handle to applyThresholdToVedioBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iModelOperation = getSelObject(handles);
theVec = updateThreshVecFromEdits(handles);

if get(handles.radiobuttonLeftThreshold,'Value')== 1
    bLeft = 1;
    handles.threshVec = theVec;
else
    bLeft = 0;
    handles.threshVec_Right = theVec;
end

% if Left Only
if iModelOperation == 0 
    handles.threshVec_ini = handles.threshVec;
end

% if Right Only
if iModelOperation == 1 
    handles.threshVec_ini_Right = handles.threshVec_Right;
end

% if in Both & Synchro mode
if iModelOperation == 2 
    handles.threshVec_ini = handles.threshVec;
    handles.threshVec_ini_Right = handles.threshVec_Right;
end
guidata(hObject,handles);


bTop = 1; nVecIdx = 0;
displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx) ;
bTop = 0; nVecIdx = 1;
displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx) ;
%--------------------------------------------------------------------------
% --- Executes on button press in saveThresholdToFileBtn.
function saveThresholdToFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to saveThresholdToFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iModelOperation = getSelObject(handles);
% if Left Only or if in Both & Synchro mode
if iModelOperation == 0 || iModelOperation == 2
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

        
end

% if Right Only or if in Both & Synchro mode
if iModelOperation == 1 || iModelOperation == 2
    if(currVidDir ~= 0)
        % Open File Browser in current video directory
        [fName,pName] = uiputfile({'*.mat','(*.mat) Mat Files'; '*.txt','(*.txt) Text Files'}...
            ,'Save Right Threshold Vector',currVidDir);
    else
        [fName,pName] = uiputfile({'*.mat','(*.mat) Mat Files'; '*.txt','(*.txt) Text Files'}...
            ,'Save Right Threshold Vector');
    end
    
    % Check for User Cancel
    if(fName == 0)
        return;
    end
    thresholdVectorFileName_Right = fullfile(pName,fName);
    saveVec = handles.threshVec_ini_Right;
    
    [pName,fName,extName] = fileparts(thresholdVectorFileName_Right);
    extName = lower(extName);
    if (strcmp(extName,'.mat')==1)
        % mat file save
        save(thresholdVectorFileName_Right, 'saveVec');
    else 
        %txt file save
        fID = fopen(thresholdVectorFileName_Right,'wt');
        fprintf(fID,'%f\t',saveVec);
        fclose(fID);        
    end    
end

% 


%--------------------------------------------------------------------------
% --- Executes on button press in loadThresholdFromFileBtn.
function loadThresholdFromFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to loadThresholdFromFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iModelOperation = getSelObject(handles);
% if Left Only or if in Both & Synchro mode
if iModelOperation == 0 || iModelOperation == 2
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
    displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx) ;
    bTop = 0;
    displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx) ;
end

% if Right Only if in Both & Synchro mode
if iModelOperation == 1 || iModelOperation == 2
    currVidDir = getCurrVidDir(handles);
    if(currVidDir ~= 0)
        % Open File Browser in current video directory
        [fName,pName] = uigetfile({'*.mat','(*.mat) Mat Files'; '*.txt','(*.txt) Text Files'}...
            ,'Load right Threshold Vector',currVidDir);
    else
        [fName,pName] = uigetfile({'*.mat','(*.mat) Mat Files'; '*.txt','(*.txt) Text Files'}...
            ,'Load right Threshold Vector');
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
    handles.threshVec_ini_Right = saveVec;
    bLeft = 0; bTop = 1; nVecIdx = 0;
    displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx) ;
    bTop = 0;
    displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx) ;
end

guidata(hObject,handles);

%--------------------------------------------------------------------------
% --- Executes on button press in autoThresholdCheckBox.
function autoThresholdCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to autoThresholdCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoThresholdCheckBox

%--------------------------------------------------------------------------
function [returnHandle] = setLowerThreshold(handles,fVec,bLeft)
typeVec = [0,0,0,0,0,1,1] ;
if bLeft == 1
    for n = 1 : length(typeVec)
        if typeVec(n) == 0;
            if handles.threshVec(n) < fVec(n)
                handles.threshVec(n) = fVec(n);
            else
                continue;
            end
        else
            if handles.threshVec(n) > fVec(n)
                handles.threshVec(n) = fVec(n);
            else
                continue;
            end
        end
    end
else
    for n = 1 : length(typeVec)
        if typeVec(n) == 0;
            if handles.threshVec_Right(n) < fVec(n)
                handles.threshVec_Right(n) = fVec(n);
            else
                continue;
            end
        else
            if handles.threshVec_Right(n) > fVec(n)
                handles.threshVec_Right(n) = fVec(n);
            else
                continue;
            end
        end
    end
    
end
returnHandle = handles;


% --- Executes on button press in radiobuttonLeftThreshold.
function radiobuttonLeftThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonLeftThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonLeftThreshold
set(hObject,'Value',1);
set(handles.radiobuttonRightThreshold,'Value',0);
if(~isfield(handles,'threshVec_ini')||isempty(handles.threshVec_ini))
    handles.threshVec_ini = handles.threshVec_org;
end
bLeft = 1; bTop = 1; nVecIdx = 0;
displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx) ;

if(~isfield(handles,'threshVec')||isempty(handles.threshVec))
    handles.threshVec = handles.threshVec_ini;
end
bTop = 0; nVecIdx = 1;
displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx) ;

% --- Executes on button press in radiobuttonRightThreshold.
function radiobuttonRightThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonRightThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonRightThreshold
set(hObject,'Value',1);
set(handles.radiobuttonLeftThreshold,'Value',0);
if(~isfield(handles,'threshVec_ini_Right')||isempty(handles.threshVec_ini_Right))
    handles.threshVec_ini = handles.threshVec_org;
end
bLeft = 0; bTop = 1; nVecIdx = 0;
displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx) ;

if(~isfield(handles,'threshVec_Right')||isempty(handles.threshVec))
    handles.threshVec_Right = handles.threshVec_ini_Right;
end
bTop = 0; nVecIdx = 1;
displayThreshholdVecToEdits(handles,bLeft,bTop,nVecIdx) ;

% this function reads in the data from the edit 


% --------------------------------------------------------------------
function saveTrackingRightMenu_Callback(hObject, eventdata, handles)
% hObject    handle to saveTrackingRightMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RTrackFileSelBtn_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------


% --------------------------------------------------------------------
function newTrackingRightMenu_Callback(hObject, eventdata, handles)
% hObject    handle to newTrackingRightMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RVideoSelBtn_Callback(hObject, eventdata, handles);
% --------------------------------------------------------------------
