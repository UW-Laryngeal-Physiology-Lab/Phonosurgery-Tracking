function [tCurr,nCurr,maxScore,fitMSE] = track_temp(frameIm,pc,...
                                                tPrev,tPrev2,nPrev,nSize)
% TRACK_TEMP Tracks weighted template w/ constant velocity model
%
% [TCURR,NCURR,MAXSCORE,FITMSE] = track_temp(FRAMEIM,PC,TPREV,TPREV2,
%                                                              NPREV,NSIZE)
% Utilizes constant velocity model and NCC to track weighted template.
% FRAMEIM is the grayscale image in which the template is searched for.  A
% bounding search region is defined by NPREV the previous upper left hand
% corner of the neighoborhood search window of sidelength NSIZE and the
% previous two template locations TPREV & TPREV2.  PC is a precomputed
% struct with template window info.  The function returns TCURR the upper
% left hand corner [x,y] of the template in the current frame, NCURR the
% location of the neighborhood search window centered on the template,
% MAXSCORE the NCC score associated with the match and fitMSE the MSE of
% the subpixel template position quadradic fitting.  The marker and
% neighborhood location are at subpixel accuracy.

%% Determine nCorn w/ constant velocity model
if(any(isnan(tPrev2)))
    % Keep Neighborhood Fixed not info for velocity estimate
    nEst = nPrev;
else
    nEst = nPrev + (tPrev - tPrev2);
end

% Perform Weighted Template Matching
[tCurr,nCurr,maxScore,fitMSE] = findTempMatch_w(frameIm,pc,nEst,nSize);