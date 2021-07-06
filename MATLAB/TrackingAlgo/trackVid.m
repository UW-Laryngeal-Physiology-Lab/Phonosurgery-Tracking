function tVid = trackVid(vidObj,frames,tCorner,nCorner,wStats)
% TRACKVID Generates Video with Neighborhood & Template Window Overlay
%
% TVID = trackVid(VIDOBJ,FRAMES,TCORNER,NCORNER,WSTATS) Generates image
% sequence from VideoReader object VIDOBJ from FRAMES(1) to FRAMES(2) with
% square overlays of simple tracking template match window and neighborhood
% search window.  TCORNER and NCORNER are (numFrames x 2) matrices where
% each column corresponds to the upper left corner of the template or
% neighborhood window along FRAMES.  WSTATS is structure describing the
% size of the template and neighborhood window.  A RGB video sequence TVID
% is returned with red overlays.
%
% Inputs:
% VIDOBJ VideoReader object with video simple tracking was performed on
% FRAMES 2 element array [StartFrame EndFrame] defining the subset of
% frames analyzed
% TCORNER (numFrames x 2) matrix.  Each column corresponds to a tracked
% frame in the sequence and is the [x,y] location of the upper left hand
% corner of the best match template window.
% NCORNER (numFrames x 2) matrix.  Each column corresponds to a tracked
% frame in the sequence and is the [x,y] location of the upper left hand
% corner of the neighborhood search window.
% WSTATS structure with fields wSize & nx,ny describing the size of the
% template window and neighborhood window.
%
% Outputs:
% TVID Overlay video sequence (rows x cols x 3 x numFrames) with tracked
% windows.


numFrames = 1 + (frames(2) - frames(1));
tVid = vidObj.read(frames);
wSize = wStats.wSize;
nx = wStats.nx-1;
ny = wStats.ny-1;

[rows,cols,rest] = size(tVid);
imSize = [rows,cols];

for k = 1:numFrames
    % Generate Mask with Neighborhood & Template Window
    tMask = trackMask(tCorner(k,:),nCorner(k,:),wSize,nx,imSize);
    
    % Generate Overlay Image and Save
    tVid(:,:,:,k) = genOverlayIm(squeeze(tVid(:,:,1,k)),tMask);
     
    %{
    % Draw Template Window
    tVid(tCorner(k,2)+[0:wSize-1],tCorner(k,1)+[0,0],:,k) = 0;
    tVid(tCorner(k,2)+[0:wSize-1],tCorner(k,1)+[wSize-1,wSize-1],:,k) = 0;
    tVid(tCorner(k,2)+[0,0],tCorner(k,1)+[0:wSize-1],:,k) = 0;
    tVid(tCorner(k,2)+[wSize-1,wSize-1],tCorner(k,1)+[0:wSize-1],:,k) = 0;

    % Draw Neighborhood Window
    tVid(nCorner(k,2)+[0:ny],nCorner(k,1)+[0,0],:,k) = 0;
    tVid(nCorner(k,2)+[0:ny],nCorner(k,1)+[nx,nx],:,k) = 0;
    tVid(nCorner(k,2)+[0,0],nCorner(k,1)+[0:nx],:,k) = 0;
    tVid(nCorner(k,2)+[ny,ny],nCorner(k,1)+[0:nx],:,k) = 0;
    %}
end
