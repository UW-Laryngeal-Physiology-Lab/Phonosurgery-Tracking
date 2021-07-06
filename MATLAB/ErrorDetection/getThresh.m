function thresh = getThresh(nullHypData,recallVal,fromDir)
% GETTHRESH Gets threshold associated w/ recall value

idx = ceil(recallVal * numel(nullHypData));

switch fromDir
    case 'fromMin'
        sortData = sort(nullHypData);
        thresh = mean(sortData([0,1] + idx));
        
    case 'fromMax'
        sortData = sort(nullHypData,'descend');
        thresh = mean(sortData([0,1] + idx));
end
