function singInstStatNoise(inst_l,inst_r,X_end,X_start,fName,sheetName)
% SINGINSTSTATNOISE Generates Figures & XLSheet for static noise eval
%
% singInstNoise(INST_L,INST_R,X_END,X_START,FNAME,SHEETNAME)

% Marker Tracking Static Noise Error
% eMarkL = diff(mark_l.subT,1,1);
% eMarkR = diff(mark_r.subT,1,1);
eInstL = inst_l.endPt - inst_l.startPt;
eInstR = inst_r.endPt - inst_r.startPt;
eWorld = X_end.' - X_start.';

% markL_eStats2 = errorStats2D(eMarkL);
% markR_eStats2 = errorStats2d(eMarkR);

% Instrument Tracking Static Noise Error
instL_eStats2 = errorStats2D(eInstL);
instR_eStats2 = errorStats2D(eInstR);

% 3D Static Noise Error
eStats3 = errorStats3D(eWorld);

% Plots
plot2DError(eInstL,'Left View');
plot2DError(eInstR,'Right View');
plot3DError(eWorld,'World');

% Save XLS Data
genSpreadSheet(fName,sheetName,instL_eStats2,instR_eStats2,eStats3)

function plot2DError(errSig,titlePrefix)
% PLOT2DERROR Generates Plots related to 2D Error Signal

figure();
subplot(2,1,1);
plot(errSig(:,1)); ylabel('x_{Err}');
title([titlePrefix ' x error signal']);

subplot(2,1,2);
plot(errSig(:,2)); ylabel('y_{Err}');
title([titlePrefix ' y error signal']);

figure();
hist(errSig(:,1)); title([titlePrefix ' 2D x Error Histogram']);
figure();
hist(errSig(:,2)); title([titlePrefix ' 2D y error histogram']);

function plot3DError(errSig,titlePrefix)
% PLOT3DERROR Generates Plots related to 3D Error Signal

figure();
subplot(3,1,1);
plot(errSig(:,1)); ylabel('X_{Err}');
title([titlePrefix ' X error signal']);

subplot(3,1,2);
plot(errSig(:,2)); ylabel('Y_{Err}');
title([titlePrefix ' Y error signal']);

subplot(3,1,3);
plot(errSig(:,3)); ylabel('Z_{Err}');
title([titlePrefix ' Z error signal']);

figure();
hist(errSig(:,1)); title([titlePrefix ' 3D X Error Histogram']);
figure();
hist(errSig(:,2)); title([titlePrefix ' 3D Y error histogram']);
figure();
hist(errSig(:,3)); title([titlePrefix ' 3D Z error histogram']);


function genSpreadSheet(fName,sheetName,eStats2_L,eStats2_R,eStats3)
xlData = cell(12,6);
xlData{1,1} = '2D Left View Error';
xlData(2:5,1:2)  = {'X Max',eStats2_L.x.Max;...
                    'X Mean',eStats2_L.x.Mean;...
                    'X STD',eStats2_L.x.Std;...
                    'X RMS',eStats2_L.x.RMS};
xlData(2:5,3:4) = {'Y Max',eStats2_L.y.Max;...
                    'Y Mean',eStats2_L.y.Mean;...
                    'Y STD',eStats2_L.y.Std;...
                    'Y RMS',eStats2_L.y.RMS};

xlData{6,1} = '2D Right View Error';
xlData(7:10,1:2)  = {'X Max',eStats2_R.x.Max;...
                    'X Mean',eStats2_R.x.Mean;...
                    'X STD',eStats2_R.x.Std;
                    'X RMS',eStats2_R.x.RMS};
xlData(7:10,3:4) = {'Y Max',eStats2_R.y.Max;...
                    'Y Mean',eStats2_R.y.Mean;...
                    'Y STD',eStats2_R.y.Std
                    'Y RMS',eStats2_R.y.RMS};

xlData{11,1} = '3D World Error';
xlData(12:15,1:2) = {'X Max',eStats3.x.Max;...
                     'X Mean',eStats3.x.Mean;...
                     'X Std',eStats3.x.Std
                     'X RMS',eStats3.x.RMS};
xlData(12:15,3:4) = {'Y Max',eStats3.y.Max;...
                     'Y Mean',eStats3.y.Mean;...
                     'Y Std',eStats3.y.Std
                     'Y RMS',eStats3.y.RMS};
xlData(12:15,5:6) = {'Z Max',eStats3.z.Max;...
                     'Z Mean',eStats3.z.Mean;...
                     'Z Std',eStats3.z.Std
                     'Z RMS',eStats3.z.RMS};
                 
xlswrite(fName,xlData,sheetName);
     
     
                 

                
                







