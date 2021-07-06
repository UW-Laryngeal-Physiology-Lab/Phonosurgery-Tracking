function varargout = multiErrorDlg(varargin)
% MULTIERRORDLG MATLAB code for multiErrorDlg.fig
%      MULTIERRORDLG by itself, creates a new MULTIERRORDLG or raises the
%      existing singleton*.
%
%      H = MULTIERRORDLG returns the handle to a new MULTIERRORDLG or the handle to
%      the existing singleton*.
%
%      MULTIERRORDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTIERRORDLG.M with the given input arguments.
%
%      MULTIERRORDLG('Property','Value',...) creates a new MULTIERRORDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before multiErrorDlg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to multiErrorDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% How to use multiErrorDlg
% This is a modal dialog that allows the user to set a stop frame and
% select an error type.  Additionally, a start frame is displayed in the
% dialog.  To tell the GUI the start frame, pass the Name, Value pair
% ('StartFrame',DESIREDSTARTFRAME) as input arguments.  Otherwise the GUI
% assumes the startFrame = 1.  The GUI outputs a 2 element array of the
% form [ENDFRAME,ERRORTYPE].  ENDFRAME is the integer end frame selected by
% the user.  ERRORTYPE is an integer indicator of the error type selected
% 1->Blur, 2->Occlusion, 3->Out of Frame 4->Other.  If the user cancelled
% or closed [0,0] is returned.

% Edit the above text to modify the response to help multiErrorDlg

% Last Modified by GUIDE v2.5 14-Nov-2011 15:43:17
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @multiErrorDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @multiErrorDlg_OutputFcn, ...
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

% --- Executes just before multiErrorDlg is made visible.
function multiErrorDlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to multiErrorDlg (see VARARGIN)

% Choose default command line output for multiErrorDlg
% output = [stopFrame,errorType]
% Error Type : 1->Blur, 2->Occlusion, 3->Out of Frame 4->Other
handles.output = [0,0];

% Display StartFrame
startFrameIndex = find(strcmp('StartFrame',varargin),1);

% If no start frame given assume frame 1
if(isempty(startFrameIndex))
    handles.startFrame = 1;
else
    handles.startFrame = varargin{startFrameIndex+1};
end

set(handles.startFrameText,'String',handles.startFrame);

% Update handles structure
guidata(hObject, handles);

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

% Make the GUI modal
set(handles.figure1,'WindowStyle','modal')

% UIWAIT makes multiErrorDlg wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = multiErrorDlg_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in okButton.
function okButton_Callback(hObject, eventdata, handles)
% hObject    handle to okButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get Stop Frame & Error Type Selection
stopFrame = round(str2double(get(handles.stopFrameEdit,'String')));
errorType = get(handles.errorTypeMenu,'Value');

% Check For Valid Stop Frame
if(stopFrame > handles.startFrame)
    handles.output = [stopFrame,errorType];
    
    % Update handles structure
    guidata(hObject, handles);

    % Use UIRESUME instead of delete because the OutputFcn needs
    % to get the updated handles structure.
    uiresume(handles.figure1);
else
    errordlg('Invalid Stop Frame','Invalid Stop Frame');
end

% --- Executes on button press in cancelButton.
function cancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Use UIRESUME instead of delete because the OutputFcn needs
% to get the updated handles structure.
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

function stopFrameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to stopFrameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stopFrameEdit as text
%        str2double(get(hObject,'String')) returns contents of stopFrameEdit as a double


% --- Executes during object creation, after setting all properties.
function stopFrameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stopFrameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in errorTypeMenu.
function errorTypeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to errorTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns errorTypeMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from errorTypeMenu

% --- Executes during object creation, after setting all properties.
function errorTypeMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to errorTypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
