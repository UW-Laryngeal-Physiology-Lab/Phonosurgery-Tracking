% This dialog is designed for setting task frame range information. The
% frames which are out of this range will be taken as 'error'.
% By : Jonathan 
% Date: 2016, Feb

function varargout = taskRangeSettingDlg(varargin)
% TASKRANGESETTINGDLG MATLAB code for taskRangeSettingDlg.fig
%      TASKRANGESETTINGDLG, by itself, creates a new TASKRANGESETTINGDLG or raises the existing
%      singleton*.
%
%      H = TASKRANGESETTINGDLG returns the handle to a new TASKRANGESETTINGDLG or the handle to
%      the existing singleton*.
%
%      TASKRANGESETTINGDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TASKRANGESETTINGDLG.M with the given input arguments.
%
%      TASKRANGESETTINGDLG('Property','Value',...) creates a new TASKRANGESETTINGDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before taskRangeSettingDlg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to taskRangeSettingDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help taskRangeSettingDlg

% Last Modified by GUIDE v2.5 25-Jan-2016 12:44:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @taskRangeSettingDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @taskRangeSettingDlg_OutputFcn, ...
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


% --- Executes just before taskRangeSettingDlg is made visible.
function taskRangeSettingDlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to taskRangeSettingDlg (see VARARGIN)

% Choose default command line output for taskRangeSettingDlg
handles.output = cell(1,2);
handles.output{1} = 1;
% load rangeInfo
handles.rangeInfo = varargin{2};
[rangeInfoRow,rangeInfoCol] = size(handles.rangeInfo);
if rangeInfoCol ~= 1
    if rangeInfoRow == 1
        displayString{1} = ['Current task has ' num2str(rangeInfoRow) 'range'];
    else
        displayString{1} = ['Current task has ' num2str(rangeInfoRow) 'ranges'];
    end
    
    displayString{2} = ['range Index     startFrame No     endFrame No' ];
    
    for n=3:rangeInfoRow+2
        displayString{n} = ['        ' num2str(handles.rangeInfo(n-2,1)) '               ' num2str(handles.rangeInfo(n-2,2)) '                       ' num2str(handles.rangeInfo(n-2,3))];
    end
    set(handles.listboxTaskRangeInfo,'String',displayString);
else
    displayString{1} = ['Current task has no range setting'];
    set(handles.listboxTaskRangeInfo,'String',displayString);
end

handles.currentRangeIdx = 0;
set(handles.currentRangeEdit,'String',num2str(handles.currentRangeIdx));
set(handles.editStartFrameNo,'String','0');
set(handles.editStopFrameNo,'String','0');

% load maxRange
handles.maxRange = varargin{4};
set(handles.editMaxRange,'String',num2str(handles.maxRange));

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

% UIWAIT makes taskRangeSettingDlg wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = taskRangeSettingDlg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes on selection change in listboxTaskRangeInfo.
function listboxTaskRangeInfo_Callback(hObject, eventdata, handles)
% hObject    handle to listboxTaskRangeInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxTaskRangeInfo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxTaskRangeInfo
idxClickSel = get(hObject,'value');
if idxClickSel <= 2
    return;
end

handles.currentRangeIdx = idxClickSel-2;
set(handles.currentRangeEdit,'String',num2str(idxClickSel-2));

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function listboxTaskRangeInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxTaskRangeInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in doneBtn.
function doneBtn_Callback(hObject, eventdata, handles)
% hObject    handle to doneBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

bCancel = 0;

handles.output{1} = bCancel;
handles.output{2} = handles.rangeInfo; 
% Update handles structure
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in cancelBtn.
function cancelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to cancelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% to get the updated handles structure.
bCancel = 1;
handles.output{1} = bCancel;    
% Update handles structure
guidata(hObject, handles);
uiresume(handles.figure1);


function editStartFrameNo_Callback(hObject, eventdata, handles)
% hObject    handle to editStartFrameNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editStartFrameNo as text
%        str2double(get(hObject,'String')) returns contents of editStartFrameNo as a double


% --- Executes during object creation, after setting all properties.
function editStartFrameNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editStartFrameNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editStopFrameNo_Callback(hObject, eventdata, handles)
% hObject    handle to editStopFrameNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editStopFrameNo as text
%        str2double(get(hObject,'String')) returns contents of editStopFrameNo as a double


% --- Executes during object creation, after setting all properties.
function editStopFrameNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editStopFrameNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveBtn.
function saveBtn_Callback(hObject, eventdata, handles)
% hObject    handle to saveBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.currentRangeIdx==0)
    return;
end

startFrameNo = str2num(get(handles.editStartFrameNo,'String'));
stopFrameNo = str2num(get(handles.editStopFrameNo,'String'));
if stopFrameNo<=startFrameNo
    msgbox('StartFrame No should less than StopFrame No. Fail to save');
    refreshRangInfo(handles);
    return;
end

if startFrameNo<=0 || startFrameNo> handles.maxRange
    msgbox('StartFrame is not valid. Fail to save');
    refreshRangInfo(handles);
    return;
end

if stopFrameNo<=0 || stopFrameNo> handles.maxRange 
    msgbox('StopFrame is not valid. Fail to save');
    refreshRangInfo(handles);
    return;
end

% Overlap check;
[rangeInfoRow,rangeInfoCol] = size(handles.rangeInfo);
for n = 1 : rangeInfoRow
    if n == handles.currentRangeIdx
        continue;
    end
    infoRow = handles.rangeInfo(n,:);
    if startFrameNo >= infoRow(2) && startFrameNo <= infoRow(3)
        msgbox('StartFrame is in other range. Overlap exists. Fail to save');
        refreshRangInfo(handles);
        return;
    end
    
    if stopFrameNo <= infoRow(3) && stopFrameNo >= infoRow(2)
        msgbox('StopFrame is in other range. Overlap exists. Fail to save');
        refreshRangInfo(handles);
        return;
    end
    
end

infoRow(1) = handles.currentRangeIdx;
infoRow(2) = startFrameNo;
infoRow(3) = stopFrameNo;

if rangeInfoCol == 1
    handles.rangeInfo = infoRow;
else
    handles.rangeInfo(handles.currentRangeIdx,:) = infoRow(:);
end

% Update handles structure
guidata(hObject, handles);

refreshRangInfo(handles);


% --- Executes on button press in appendRangeBtn.
function appendRangeBtn_Callback(hObject, eventdata, handles)
% hObject    handle to appendRangeBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[rangeInfoRow,rangeInfoCol] = size(handles.rangeInfo);
if rangeInfoCol ~= 1
    displayString{1} = ['Current task has ' num2str(rangeInfoRow+1) 'ranges'];    
    displayString{2} = ['range Index     startFrame No     endFrame No' ];
    
    handles.currentRangeIdx = rangeInfoRow + 1;
    % Update handles structure
    guidata(hObject, handles);
    set(handles.currentRangeEdit,'String',num2str(handles.currentRangeIdx));
    
    for n=3:rangeInfoRow+2
        displayString{n} = ['        ' num2str(handles.rangeInfo(n-2,1)) '               ' num2str(handles.rangeInfo(n-2,2)) '                       ' num2str(handles.rangeInfo(n-2,3))];
    end
    displayString{rangeInfoRow+3} = ['        ' num2str(handles.currentRangeIdx) '                 ' '0' '                         ' '0'];
    set(handles.listboxTaskRangeInfo,'String',displayString);
else
    displayString{1} = ['Current task has 1 range'];
    displayString{2} = ['range Index     startFrame No     endFrame No' ];    
    handles.currentRangeIdx = 1;
    % Update handles structure
    guidata(hObject, handles);
    set(handles.currentRangeEdit,'String',num2str(handles.currentRangeIdx));
    
    displayString{rangeInfoRow+2} = ['        ' num2str(handles.currentRangeIdx) '                 ' '0' '                         ' '0'];
    set(handles.listboxTaskRangeInfo,'String',displayString);
end

set(handles.editStartFrameNo,'String','0');
set(handles.editStopFrameNo,'String','0');


% --- Executes on button press in modifyRangBtn.
function modifyRangBtn_Callback(hObject, eventdata, handles)
% hObject    handle to modifyRangBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.currentRangeIdx==0)
    return;
end

infoRow = handles.rangeInfo(handles.currentRangeIdx,:);
set(handles.editStartFrameNo,'String',num2str(infoRow(2)));
set(handles.editStopFrameNo,'String',num2str(infoRow(3)));


% --- Executes on button press in deleteRangeBtn.
function deleteRangeBtn_Callback(hObject, eventdata, handles)
% hObject    handle to deleteRangeBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.currentRangeIdx==0)
    return;
end

handles.rangeInfo(handles.currentRangeIdx,:) = [];
handles.currentRangeIdx = 0;
% Update handles structure
guidata(hObject, handles);

set(handles.currentRangeEdit,'String',num2str(handles.currentRangeIdx));

set(handles.editStartFrameNo,'String','0');
set(handles.editStopFrameNo,'String','0');

refreshRangInfo(handles);



function currentRangeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to currentRangeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentRangeEdit as text
%        str2double(get(hObject,'String')) returns contents of currentRangeEdit as a double


% --- Executes during object creation, after setting all properties.
function currentRangeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentRangeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function refreshRangInfo(handles)
% This function refresh listbox to load latest range info
[rangeInfoRow,rangeInfoCol] = size(handles.rangeInfo);
if rangeInfoCol ~= 1
    if rangeInfoRow == 1
        displayString{1} = ['Current task has ' num2str(rangeInfoRow) 'range'];
    else
        displayString{1} = ['Current task has ' num2str(rangeInfoRow) 'ranges'];
    end
    
    displayString{2} = ['range Index     startFrame No     endFrame No' ];
    
    for n=3:rangeInfoRow+2
        displayString{n} = ['        ' num2str(handles.rangeInfo(n-2,1)) '               ' num2str(handles.rangeInfo(n-2,2)) '                       ' num2str(handles.rangeInfo(n-2,3))];
    end
    set(handles.listboxTaskRangeInfo,'String',displayString);
else
    displayString{1} = ['Current task has no range setting'];
    set(handles.listboxTaskRangeInfo,'String',displayString);
end



function editMaxRange_Callback(hObject, eventdata, handles)
% hObject    handle to editMaxRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMaxRange as text
%        str2double(get(hObject,'String')) returns contents of editMaxRange as a double


% --- Executes during object creation, after setting all properties.
function editMaxRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMaxRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


