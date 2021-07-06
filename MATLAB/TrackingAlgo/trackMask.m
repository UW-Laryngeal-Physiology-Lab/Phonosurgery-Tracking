function tMask = trackMask(tCurr,nCurr,tSize,nSize,imSize)
% TRACKMASK Generates Binary mask with template and neighborhood windows
%
% TMASK = trackMask(TCURR,NCURR,TSIZE,NSIZE,IMSIZE)  Given the upper left
% hand corners [x,y] of a template window of size TSIZE and neighborhood
% search window size NSIZE a binary mask of IMSIZE is generated with
% squares representing the two windows in TMASK.

tMask = false(imSize);

% Compute Crosshair Image Indices
halfSize = (round((tSize-1)/2));
midPts = tCurr + halfSize;
%crossHairDex = sub2ind(imSize,...
%                [(0:tSize-1)' + tCurr(2);midPts(2)*ones(tSize,1)],...
%                [midPts(1)*ones(tSize,1);(0:tSize-1)' + tCurr(1)]);
vertMidlineDex = sub2ind(imSize,midPts(2)*ones(tSize,1),...
                                        (0:tSize-1)' + tCurr(1));

%tMask([getSquareDex(tCurr,tSize,imSize),...
%       getSquareDex(nCurr,nSize,imSize),...
%       crossHairDex']) = 1;
tMask([getSquareDex(tCurr,tSize,imSize),...
       getSquareDex(nCurr,nSize,imSize),...
       vertMidlineDex']) = 1;


end
function sDex = getSquareDex(wCorner,wSize,imSize)
uLeft = [max([1,wCorner(1)]),...
         max([1,wCorner(2)])];
lRight = [min([imSize(2),wCorner(1)+(wSize-1)]),...
          min([imSize(1),wCorner(2)+(wSize-1)])];
numHoriz = 1 + (lRight(1) - uLeft(1));
numVert = 1 + (lRight(2) - uLeft(2));
      
startDex = sub2ind(imSize,uLeft(2),uLeft(1));
horiz1 = imSize(1)*(0:numHoriz-1) + startDex;
horiz2 = horiz1 + (numVert-1);

vert1 = startDex + (0:numVert-1);
vert2 = imSize(1)*(numHoriz-1) + vert1;
sDex = [horiz1,horiz2,vert1,vert2];
%sDex = sDex(and(sDex >= 1,sDex <= imSize(1)*imSize(2)));
end




