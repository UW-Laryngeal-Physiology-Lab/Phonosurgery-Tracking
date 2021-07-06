function overlayVid = overlayTuneVidGen(vidObj,frames,algoInterface,...
    pStruct)
%OVERALYTUNEVIDGEN Summary of this function goes here
%   Detailed explanation goes here

% Video Info
fNums = frames(1):frames(2);
M = vidObj.Height; N = vidObj.Width;

% Preallocate overlayVid
overlayVid = zeros([M,N,3,numel(fNums)],'uint8');

for k = 1:numel(fNums)
    % Grab Frame Data
    frameIm = vidObj.read(fNums(k));

    % Compute Algo Mask
    try
        maskIm = algoInterface(frameIm(:,:,1),pStruct);

        % Generate OverlayFrame
        % Mask Pixels are Pure Red
        maskIdx = find(maskIm);
        frameIm(maskIdx) = 255;
        frameIm([maskIdx + M*N;maskIdx + 2*M*N]) = 0; 
    catch ME
        fprintf('Error @ Frame %u\n %s \n',fNums(k),getReport(ME));
    end
    overlayVid(:,:,:,k) = frameIm;
end
end

