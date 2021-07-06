function varargout = lowAlgoConfDialog(varargin)
% LOWALGOCONFDIALOG MATLAB code for lowAlgoConfDialog.fig
%      LOWALGOCONFDIALOG, by itself, creates a new LOWALGOCONFDIALOG or raises the existing
%      singleton*.
%
%      H = LOWALGOCONFDIALOG returns the handle to a new LOWALGOCONFDIALOG or the handle to
%      the existing singleton*.
%
%      LOWALGOCONFDIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOWALGOCONFDIALOG.M with the given input arguments.
%
%      LOWALGOCONFDIALOG('Property','Value',...) creates a new LOWALGOCONFDIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lowAlgoConfDialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lowAlgoConfDialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lowAlgoConfDialog

% Last Modified by GUIDE v2.5 27-Oct-2015 11:38:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lowAlgoConfDialog_OpeningFcn, ...
                   'gui_OutputFcn',  @lowAlgoConfDialog_OutputFcn, ...
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


% --- Executes just before lowAlgoConfDialog is made visible.
function lowAlgoConfDialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lowAlgoConfDialog (see VARARGIN)

% Choose default command line output for lowAlgoConfDialog
handles.output = 'No';

% Update handles structure
guidata(hObject, handles);

% Get User Defined Location
if(nargin == 5)
    % Custom Dialog Position
    set(hObject,'Units','pixels');
    defPos = get(hObject,'Position');
    defPos(1) = varargin{1}; defPos(2) = varargin{2};
    set(hObject,'Position',defPos);
end

% UIWAIT makes lowAlgoConfDialog wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lowAlgoConfDialog_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in yesButton.
function yesButton_Callback(hObject, eventdata, handles)
% hObject    handle to yesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
registerYes(hObject,handles);

% REGISTERYES Tells GUI to register 'yes' press
%
% registerYes(hObject,handles) Registers a press of the 'Yes' button.  Sets
% the output parameter of handles structure to 'Yes' and sends uiresume
% command which initiates closing of the GUI.
function registerYes(hObject,handles)
handles.output = 'Yes';
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes on button press in noButton.
function noButton_Callback(hObject, eventdata, handles)
% hObject    handle to noButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
registerNo(hObject,handles);

% REGISTERNO Tells GUI to register 'No' press
%
% registerNo(hObject,handles) Registers a press of the 'No' button.  Sets
% the output parameter of handles structure to 'No' and sends uiresume
% command which initiats closing of the GUI.
function registerNo(hObject,handles)
handles.output = 'No';
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

% Register 'Yes' or 'No' if 'y' or 'n' was pressed.  Otherwise do nothing.
switch(eventdata.Key)
    case 'y'
        registerYes(hObject,handles);
    case 'n'
        registerNo(hObject,handles);
end
