function renameRot(prefix,rotDeg,newPrefix,imType)
% RENAMEROT Renames & Rotates a set of pictures with given prefix
%
% renameRot(PREFIX,ROTDEG,NEWPREFIX,IMTYPE='.tiff') Finds all images in the
% current directory with file name prefix PREFIX, rotates them ROTDEG &
% saves a file with PREFIX replaced by NEWPREFIX & image file type IMTYPE.

%% Input Arguments
if(nargin == 3)
    imType = '.tif';
end

%% Read Directory
dirStruct = dir(cd);
fNames = {dirStruct.name};
matches = strfind(fNames,prefix);

len = numel(prefix);
for k = 1:numel(matches)
    if(~isempty(matches{k}))
        if(matches{k} == 1)
            % Load Image
            fName = fNames{k}; [~,name,~] = fileparts(fName);
            imData = imread(fName);
            
            % Create new filename & rotated image
            newName = [newPrefix name(len+1:end)];
            newIm = imrotate(imData,rotDeg);
            
            % Write File
            imwrite(newIm,[newName imType]);
        end
    end
end