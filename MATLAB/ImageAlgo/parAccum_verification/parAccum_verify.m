% parAccum_verify.m
%% Generate Image
imSize = [600,400];
ang = 20;
l1 = [150,(pi/180)*ang];
l2 = [185,(pi/180)*ang];

mask_left = drawLineMask(imSize,l1(1),l1(2));
mask_right = drawLineMask(imSize,l2(1),l2(2));
mask = mask_left | mask_right;
%figure();imshow(mask);

[HPar_L,Hpar_R,Match_R,T_L,R_L] = parAccum(mask_left,mask_right,[30,40],...
                                                [-1.5,1.5],10,[-2,2] + ang);
%% Maximum Pair
[~,peakIdx] = max(HPar_L(:) + Hpar_R(:));
[peakRow,peakCol] = ind2sub(size(HPar_L),peakIdx);

rho_left = R_L(peakRow); theta_left = T_L(peakCol);
rho_right = Match_R(peakRow,peakCol,2); 
theta_right = Match_R(peakRow,peakCol,1);

disp([rho_left,theta_left]);
disp([rho_right,theta_right]);
                                            
