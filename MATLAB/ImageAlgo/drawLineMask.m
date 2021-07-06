function lineMask = drawLineMask(imSize,rho,theta)
% DRAWLINEMASK Draws Line Mask from Foot of Normal Line Parameterization
%
% LINEMASK = drawLineMask(IMSIZE,RHO,THETA) 

lineMask = false(imSize);

% Draw Lines
try
    for k = 1:numel(rho)
        cosVal = cos(theta(k));
        sinVal = sin(theta(k));
        if (abs(cosVal) < (1/sqrt(2)))
            xPts = [1,imSize(2)];
            %yPts = (repmat(rho(k),1,2) - xPts .* repmat(cosVal,1,2)) ./ ...
            %            repmat(sinVal,1,2);
            yPts = (repmat(rho(k),1,2) - (xPts-1) .* repmat(cosVal,1,2)) ./ ...
                        repmat(sinVal,1,2);
            yPts = yPts + 1;
        else
            yPts = [1,imSize(1)];
            %xPts = (repmat(rho(k),1,2) - yPts .* repmat(sinVal,1,2)) ./ ...
            %    repmat(cosVal,1,2);
            xPts = (repmat(rho(k),1,2) - (yPts-1) .* repmat(sinVal,1,2)) ./ ...
                repmat(cosVal,1,2);
            xPts = xPts + 1;
        end

        ils = intLineSeg([xPts(1),yPts(1)],[xPts(2),yPts(2)]);
        [xIdx,yIdx] = ils.getIntLocs(); 
        goodIdx = and(xIdx >= 1,xIdx <= imSize(2));
        lineMask(sub2ind(imSize,yIdx(goodIdx),xIdx(goodIdx))) = 1;
    end
catch ME
    disp('Line Drawing Error');
end
end