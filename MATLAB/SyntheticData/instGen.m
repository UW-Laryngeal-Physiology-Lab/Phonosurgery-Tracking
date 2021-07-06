function [instIm instMask] = instGen(imSize,instInfo,markerInfo,posInfo)
% INSTGEN Generates synthetic instrument image
%
% [INSTIM,INSTMASK] = instGen(IMSIZE,INSTINFO,MARKERINFO,POSINFO)
%
% INSTINFO Fields:
% mean : Mean Intensity of Instrument (Gaussian Model)
% std : Intensity Std of Instrument (Gaussian Model)
% size : [Rows,Cols] array describing the rectangular instrument size
% 
% MARKERINFO
% offset : Distance in pixels from bottom of instrument at which the marker
%          begins
% numStripes : Number of Stripes
% thickness : Thickness in pixels of each stripes
% darkMean : Mean intensity of dark stripes (gaussian model)
% darkStd : Intensity std of dark stripes (gaussian model)
% lightMean : Mean intensity of light stripes (gaussian model)
% lightStd : Intensity std of light stripes (gaussian model)
%
% POSINFO
% x : x position of instrument lower left hand corner
% y : y position of instrument lower left hand corner
% ang : Angle in radians of instrument axis where 0 is vertical and
% positive values result in counter-clockwise rotation

% Get Template of Instrument
iTemp = drawInstTemp(instInfo,markerInfo);

% Utilize Inverse Warp to draw instrument
[X,Y] = meshgrid(0:imSize(2)-1,0:imSize(1)-1);
[X_inst,Y_inst] = meshgrid(0:instInfo.size(2)-1,0:instInfo.size(1)-1);

R = [cos(posInfo.ang),sin(posInfo.ang);
     -sin(posInfo.ang),cos(posInfo.ang)];
sourcePts = R.' * [(X(:).' - posInfo.x);(Y(:).'- posInfo.y)] + ...
    repmat([0;instInfo.size(1)],1,numel(X));

instIm = interp2(X_inst,Y_inst,double(iTemp),sourcePts(1,:),...
    sourcePts(2,:));
instIm = reshape(instIm,imSize);
instMask = (~isnan(instIm)); instIm = uint8(instIm);

end

function instTemp = drawInstTemp(instInfo,mInfo)
% DRAWINSTTEMP Generates Instrument template in canonical position
%
% INSTTEMP = drawInstTemp(INSTINFO,MINFO) Draws instrument template based
% on parameter structures INSTINFO and MINFO.  See below for description of
% structures.  INSTTEMP is an image the size of the instrument with the
% marker pattern based on MINFO.  In this image the instrument is vertical.
%
% INSTINFO Fields:
% mean : Mean Intensity of Instrument (Gaussian Model)
% std : Intensity Std of Instrument (Gaussian Model)
% size : [Rows,Cols] array describing the rectangular instrument size
% 
% MINFO
% offset : Distance in pixels from bottom of instrument at which the marker
%          begins
% numStripes : Number of Stripes
% thickness : Thickness in pixels of each stripes
% darkMean : Mean intensity of dark stripes (gaussian model)
% darkStd : Intensity std of dark stripes (gaussian model)
% lightMean : Mean intensity of light stripes (gaussian model)
% lightStd : Intensity std of light stripes (gaussian model)

instTemp = uint8(repmat(instInfo.mean,instInfo.size) + ...
    instInfo.std*randn(instInfo.size));
instTemp = stripeMarkerOverlay(instTemp,mInfo);

end

function instIm = stripeMarkerOverlay(instTemp,mInfo)
% STRIPEMARKEROVERLAY Overlays stripe marker on instrument template
%
% INSTIM = stripMarkerOverlay(INSTTEMP,MINFO) Draws a stripe marker overlay
% on INSTTEMP using parameters in MINFO structure.  See below for
% description of structure.  INSTTEMP should be an instrument template
% image, meanings it is the instrument in its canonical position.  INSTIM
% is the template image with the strip marker overlay.
%
% MINFO Fields
% offset : Distance in pixels from bottom of instrument at which the marker
%          begins
% numStripes : Number of Stripes
% thickness : Thickness in pixels of each stripes
% darkMean : Mean intensity of dark stripes (gaussian model)
% darkStd : Intensity std of dark stripes (gaussian model)
% lightMean : Mean intensity of light stripes (gaussian model)
% lightStd : Intensity std of light stripes (gaussian model)

[~,numCols] = size(instTemp);

% Initialize Start Row as dark
startRow = size(instTemp,1) - mInfo.offset;
darkRow = 1;

for k = 1:1:mInfo.numStripes
    endRow = startRow - (mInfo.thickness-1);
    
    if(darkRow)
        % Draw Dark Row
        instTemp(endRow:startRow,:) = mInfo.darkMean + ...
            mInfo.darkStd*randn([mInfo.thickness,numCols]);
        darkRow = 0;
    else
        % Draw Light Row
        instTemp(endRow:startRow,:) = mInfo.lightMean + ...
            mInfo.lightStd*randn([mInfo.thickness,numCols]);
        darkRow = 1;
    end
    startRow = startRow - mInfo.thickness;
end

instIm = instTemp;

end