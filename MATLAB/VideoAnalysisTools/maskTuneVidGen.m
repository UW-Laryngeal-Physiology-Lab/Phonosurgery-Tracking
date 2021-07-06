function maskVid = maskTuneVidGen(vidObj,frames,algoInterface,pStruct)
% MASKTUNEVIDGEN 
%
% MASKVID = maskTuneVidGen(VIDOBJ,FRAMES,ALGOINTERFACE,PSTRUCT)

% Video Info
fNums = frames(1):frames(2);
M = vidObj.Height; N = vidObj.Width;

% Preallocate overlayVid
maskVid = false([M,N,1,numel(fNums)]);

for k = 1:numel(fNums)
    % Grab Frame Data
    frameIm = vidObj.read(fNums(k));

    % Compute Algo Mask
    try
        maskVid(:,:,1,k) = algoInterface(frameIm(:,:,1),pStruct);
    catch ME
        fprintf('Error @ Frame %u\n %s \n',fNums(k),getReport(ME));
    end
end
end

       