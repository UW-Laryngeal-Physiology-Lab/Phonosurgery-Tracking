function cboardSave(PPI,xSquares,ySquares,maxLen_mm,fileName)

PPMM = PPI/25.4; pxsideLen_mm = 1/PPMM;
pxPerSquare = round(maxLen_mm /(max([xSquares,ySquares])*pxsideLen_mm));
cBoard = checkerboard(pxPerSquare,ySquares,xSquares);
cBoard = (cBoard(1:ySquares*pxPerSquare,1:xSquares*pxPerSquare) > 0.5);

sideLen_mm = pxPerSquare * pxsideLen_mm;
fprintf('Square Side Length %2.6f\n',sideLen_mm);
imwrite(cBoard,fileName,'Compression','None');