function interpFrames = findInterpFrames(label)
% FINDINTERPFRAMES Finds frames that need interpolation based on UGT labels
%
% INTERPFRAMES = findInterpFrames(LABEL) Finds frame indices in LABEL whose
% 3D position and orientation should be interpolated.  Logical array
% INTERPFRAMES is returned.  A value of 1 indicates an error label frame
% that should be interpolated.
%
% Interpolation Frame Criteria
% 1) Single error label sandwiched between correct labels
% 2) Multiple blur labels sandwiched between correct labels

% Find Error Labels Proceeding Correct Labels (atartErrors)
errorLabel = ~isnan(label) & label > 0;
startError = 1 + find(diff(errorLabel) == 1);

% Identify interpFrames
% 1) Single error label sandwiched between correct labels
% 2) Multiple blur labels sandwiched between correct labels
interpFrames = false(numel(label),1);
for k = 1:numel(startError)
    idx = startError(k);
    
    % Do Nothing if it is the last frame
    if(idx == numel(label))
        continue;
    end
    
    % Check if error is sandwiched between correct labels
    if(label(idx + 1) == 0)
        interpFrames(idx) = 1;
        continue;
    end
    
    % Check proceeding frames.  If a non-blur error or unlabeld frame is
    % reached prior to a correct label do nothing.  If all proceeding
    % frames are blur labels followed by correct label, add sequence of
    % frames to interpFrames
    for m = (idx+1):numel(label)
        % Unlabelld or non-blur error label
        if(isnan(label(m)) || label(m) > 1)
            break;
        end
        
        % Correct Label
        if(label(m) == 0)
            interpFrames(idx:(m-1)) = 1;
            break;
        end     
    end
end