function fileList = writeFileList(fileExt)
%WRITEFILELIST Prompts User For File Names Than writes locations to .txt
%
%FILELIST = writeFileList(FILEEXT) 
%Inputs:
%FILEEXT - Extension of Type of file you are trying to make a list of.
%This is a string.  For example '.wav'.  If the file does not have an
%extension use '.*'.
%Outputs:
%FILELIST - Cell Array Containing all the user selected file names
%
%Description: 
%Continually prompts user for file names of type FILEEXT.  All
%the file names are accumulated into the cell array FILELIST.  Additionally
%after all the file names have been selected the user is prompted to save a
%.txt file.  This file contains each file name on an individual line.

%KS 10/18/10

fileList = {}; % Holds the string of all files

prevPath = '.';
userContinue = 1;
k = 0;
while userContinue
    % Get File Name
    [fileName pathName] = uigetfile(['*' fileExt],'Get File',prevPath,...
        'MultiSelect','on');
    if iscell(fileName)
        % Store Multiple File Names
        numElements = length(fileName);
        for fileNum = 1:numElements
            ffile = fullfile(pathName,fileName{fileNum});
            k = k + 1;
            fileList{k,1} = ffile;
            prevPath = ffile;
        end
    elseif ischar(fileName)
        % Store Single File
        ffile = fullfile(pathName,fileName);
        k = k + 1;
        fileList{k,1} = ffile;
        prevPath = ffile;
    end
    % Ask User if they have another file
    button = questdlg('Select another file?','Continue?');
    if ~strcmp(button,'Yes')
        % User Wants to Stop
        userContinue = 0;
    end
end

% Save the file list into a .txt
[fileName pathName] = uiputfile('*.txt','Save File List');
fName = fullfile(pathName,fileName);
if ischar(fileName)
    fid = fopen(fName,'w');
    % Write The Files Names One at a Time
    for k = 1:length(fileList)
        fprintf(fid,'%s\n',fileList{k,1});
    end
    fclose(fid);
end
