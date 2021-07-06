function varargout = UGT_Parameter_Extract(varargin)
% UGT_PARAMETER_EXTRACT MATLAB code for UGT_Parameter_Extract.fig
%      UGT_PARAMETER_EXTRACT, by itself, creates a new UGT_PARAMETER_EXTRACT or raises the existing
%      singleton*.
%
%      H = UGT_PARAMETER_EXTRACT returns the handle to a new UGT_PARAMETER_EXTRACT or the handle to
%      the existing singleton*.
%
%      UGT_PARAMETER_EXTRACT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UGT_PARAMETER_EXTRACT.M with the given input arguments.
%
%      UGT_PARAMETER_EXTRACT('Property','Value',...) creates a new UGT_PARAMETER_EXTRACT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UGT_Parameter_Extract_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UGT_Parameter_Extract_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UGT_Parameter_Extract

% Last Modified by GUIDE v2.5 30-Nov-2015 14:10:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UGT_Parameter_Extract_OpeningFcn, ...
                   'gui_OutputFcn',  @UGT_Parameter_Extract_OutputFcn, ...
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


% --- Executes just before UGT_Parameter_Extract is made visible.
function UGT_Parameter_Extract_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UGT_Parameter_Extract (see VARARGIN)

% Choose default command line output for UGT_Parameter_Extract
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes UGT_Parameter_Extract wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = UGT_Parameter_Extract_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Executes on button press in setThreshVecBtn02.
function setThreshVecBtn02_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn02 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in setThreshVecBtn03.
function setThreshVecBtn03_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn03 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in setThreshVecBtn04.
function setThreshVecBtn04_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn04 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in setThreshVecBtn05.
function setThreshVecBtn05_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn05 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in setThreshVecBtn06.
function setThreshVecBtn06_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn06 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in setThreshVecBtn07.
function setThreshVecBtn07_Callback(hObject, eventdata, handles)
% hObject    handle to setThreshVecBtn07 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in applyThresholdToVedioBtn.
function applyThresholdToVedioBtn_Callback(hObject, eventdata, handles)
% hObject    handle to applyThresholdToVedioBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in saveThresholdToFileBtn.
function saveThresholdToFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to saveThresholdToFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in loadThresholdFromFileBtn.
function loadThresholdFromFileBtn_Callback(hObject, eventdata, handles)
% hObject    handle to loadThresholdFromFileBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function editVideoFile_Callback(hObject, eventdata, handles)
% hObject    handle to editVideoFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editVideoFile as text
%        str2double(get(hObject,'String')) returns contents of editVideoFile as a double


% --- Executes during object creation, after setting all properties.
function editVideoFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVideoFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTrackingFile_Callback(hObject, eventdata, handles)
% hObject    handle to editTrackingFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTrackingFile as text
%        str2double(get(hObject,'String')) returns contents of editTrackingFile as a double


% --- Executes during object creation, after setting all properties.
function editTrackingFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTrackingFile (see GCBO)
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


% --- Executes on button press in previousFrameButton.
function previousFrameButton_Callback(hObject, eventdata, handles)
% hObject    handle to previousFrameButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nextFrameButton.
function nextFrameButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextFrameButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function currFrameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to currFrameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currFrameEdit as text
%        str2double(get(hObject,'String')) returns contents of currFrameEdit as a double


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
