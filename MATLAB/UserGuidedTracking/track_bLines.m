function [rho,theta,ais] = track_bLines(frameIm,backIm,tCorn,tPrev,...
                                            trackPtPrev,thetaPrev,paramStruct)
% DETECTION_BLINES Tracks single instrument boundary lines
%
% [RHO,THETA,AIS] = track_bLines(FRAMEIM,BACKIM,TCORN,TPREV,TRACKPTPREV,
%                                                               THETAPREV)
% Utilizes previous frame line estimates & current frame template to track
% instrument boundary lines.  FRAMEIM is the grayscale frame image from
% which boundary lines need to be detected. BACKIM is background image used
% for background subtraction.  TCORN is the upper left corner [x,y] of the
% template window in the frame, TPREV corresponds to the template in the
% previous frame, TRACKPTPREV is the track point in the previous frame
% [x,y], THETAPREV are the boundary line theta parameterizations (in
% radians) in the previous frame.
%
% [RHO,THETA,AIS] = track_bLines(FRAMEIM,BACKIM,TCORN,TPREV,TRACKPTPREV,
%                                                    THETAPREV,PARAMSTRUCT)
% PARAMSTRUCT is an optional structure used to set several algorithm
% parameters.  See below for decription, all fields are optional.
%
% PARAMSTRUCT Fields 
% Backgnd Subtraction : backThresh, w_back, w_dark, edgeThresh,
% Boundary Line Estimation : delRho, delTheta, thetaRes, inlierThresh

%% Input Arguments
if(nargin == 6)
    paramStruct = struct();
end

% Combine Input Parameters with Defaults
pStruct = checkParameters(paramStruct);
edgeStruct = pStruct;
%edgeStruct.backIm = backIm;

%edgeStructL = edgeStruct; edgeStructL.edgeName = 'left';
%edgeStructR = edgeStruct; edgeStructR.edgeName = 'right';

%% Prediction Mask based on previous estimate
tDisp = tCorn - tPrev;
pMask = predictMask(tCorn,tDisp,trackPtPrev,...
                [0,45],thetaPrev,size(frameIm));
            
% Determine subIm
horiz = any(pMask,1); vert = any(pMask,2);
subIm = [find(horiz,1,'first'),find(horiz,1,'last'),...
         find(vert,1,'first'),find(vert,1,'last')];

pStruct.thetaBounds = (180/pi * mean(thetaPrev)) + [-5,5];
            
%% Boundary Line Estimation
% Detect Edges
[leftEdge,rightEdge] = boundEdgeDetect(frameIm,backIm,edgeStruct,subIm);
%leftEdge = Itune_edgeDetect(frameIm,edgeStructL); 
leftEdge_1 = leftEdge .* pMask;

%rightEdge = Itune_edgeDetect(frameIm,edgeStructR); 
rightEdge_1 = rightEdge .* pMask;

% Compute Initial Estimates
[r_1,t_1,pld_ais] = parLineDetect(leftEdge_1,rightEdge_1,pStruct);

% Sort [left,right]
[r_1,sortIdx] = sort(r_1,'ascend'); t_1 = t_1(sortIdx);

% Refit Boundary Lines
[yL_1,xL_1] = find(leftEdge_1); [yR_1,xR_1] = find(rightEdge_1);
[rL_1,tL_1,refitMSEL,numInliersL] = refitBoundLine(r_1(1),t_1(1),xL_1,yL_1,pStruct.inlierThresh);
[rR_1,tR_1,refitMSER,numInliersR] = refitBoundLine(r_1(2),t_1(2),xR_1,yR_1,pStruct.inlierThresh);
rho = [rL_1,rR_1]; theta = [tL_1,tR_1];

if(nargout == 3)
    ais = struct('leftNumInliers',numInliersL,'leftRefitMSE',refitMSEL,...
        'rightNumInliers',numInliersR,'rightRefitMSE',refitMSER,...
        'votesLeft',pld_ais.votesLeft,'votesRight',pld_ais.votesRight);
end
end           

function [pStruct] = checkParameters(paramStruct)
% CHECKPARAMETERS Checks parameter structure
%
% PSTRUCT = checkParameters(PARAMSTRUCT) Checks PARAMSTRUCT for fields
% utilized by algorithms in this function.  If the field is missiing in
% PARAMSTRUCT it is set to the default values.

% Background Subtraction & Edgel Detection Related
paramStruct = checkField(paramStruct,'backThresh',30);
paramStruct = checkField(paramStruct,'w_back',1.5);
paramStruct = checkField(paramStruct,'w_dark',1);
paramStruct = checkField(paramStruct,'edgeThresh',10);

% Line Detection Related
%paramStruct = checkField(paramStruct,'delRho',[35,45]);
%paramStruct = checkField(paramStruct,'delTheta',[-1.5,1.5]);
paramStruct = checkField(paramStruct,'delRho',[30,45]);
paramStruct = checkField(paramStruct,'delTheta',[-3.5,3.5]);
paramStruct = checkField(paramStruct,'thetaRes',0.1);

% Refitting Related
paramStruct = checkField(paramStruct,'inlierThresh',10);

pStruct = paramStruct;
end


function paramStruct = checkField(paramStruct,fieldName,defaultVal)
% CHECKFIELD Checks for field in structure
%
% PARAMSTRUCT = checkField(PARAMSTRUCT,FIELDNAME,DEFAULTVAL) Checks if
% PARAMSTRUCT has field FIELDNAME.  If it does not it is added to
% PARAMSTRUCT with value DEFAULTVAL.

if(~isfield(paramStruct,fieldName))
    paramStruct.(fieldName) = defaultVal;
end
end
