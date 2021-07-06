function [lineMask,rho,theta,ais,tStruct] = blobToBound(instBlob,candIm,...
                                                                    params)
% BLOBTOBOUND Finds Boundary Lines for blob associated w/ single inst
%
% [LINEMASK,RHO,THETA] = blobToBound(INSTBLOB,CANDIM,PARAMS=default) Fits
% vertical boundary lines to binary image INSTBLOB representing pixels of a
% single instrument.  CANDIM is the full candidate image from the same
% source image INSTBLOB was generated from.  PARAMS is a struct of
% parameters controlling the fitting algorithm.  PARAMS itself and its
% fields are optional (default values are used when not set by the user).
% See below for a full description of the available parameters.  LINEMASK
% is a binary image with the boundary line and midline fitting.  RHO and
% THETA are 2 element arrays containing the boundary line parameterization
% in foot-of-normal form.
%
% PARAMS fields
% fattenR : (default=3) Radius of disk structuring element used for
% fattening boundary candidate points
% gradThresh : (default=20) Threshold applied to gradient magnitude of
% CANDIM to determine boundary line candidate points
% rhoResolution : (default=0.5) Resolution of Rho parameter used in Hough
% Space when determining parameterized boundary lines.
% thetaVals : (default=-45:0.1:45) Hough space theta values used to
% determine parameterized boundary lines.
% suppress : (default=[45,201]) Size of suppresion mask used to determine
% peaks in the Hough Space when determining parameterized boundary lines.
% It is in the form [Rho,Theta].

%% Algo Parameters
% Parameter Defaults
fattenR = 3;
gradThresh = 20;
rhoResolution = 0.5;
thetaVals = -45:0.1:45;
suppress = [45,201];
inlierThresh = 4; % inlierThresh = 4;

% Assign User Defined Parameters
if nargin == 3
    fNames = fieldnames(params);
    for k = 1:numel(fNames)
        switch(fNames{k})
            case 'fattenR'
                fattenR = params.fattenR;
            case 'gradThresh'
                gradThresh = params.gradThresh;
            case 'thetaVals'
                thetaVals = params.thetaVals;
            case 'rhoResolution'
                rhoResolution = params.rhoResolution;
            case 'suppress'
                suppress = params.suppress;
            case 'inlierThresh'
                inlierThresh = params.inlierThresh;
        end
    end
end

%% Run Algo
if(nargout >= 4)
    % Run Functions with ais & tStruct data (nargout == 5)  
    % Boundary Candidates
    [candMask,fbc_ais] = findBoundCand(find(instBlob),size(instBlob),fattenR);

    % Line Candidates
    [lineCandMask,lc_ais] = lineCands(candIm,candMask,gradThresh);

    % Boundary Lines
    if(nargout == 5)
        [rho_init,theta_init,bl_ais,houghStruct] = boundLines(...
        lineCandMask,thetaVals,rhoResolution,suppress);
    else
        [rho_init,theta_init,bl_ais] = boundLines(...
        lineCandMask,thetaVals,rhoResolution,suppress);
    end
    
    % Refit Boundary Lines
    [y,x] = find(lineCandMask);
    [r1,t1,refitMSE1] = refitBoundLine(rho_init(1),theta_init(1),x,y,inlierThresh);
    [r2,t2,refitMSE2] = refitBoundLine(rho_init(2),theta_init(2),x,y,inlierThresh);
    rho = [r1,r2]; theta = [t1,t2];
    
    % Build the ais
    ais = struct('numBoundCands',fbc_ais.numBoundCands,...
                 'numLineCands',lc_ais.numLineCands,...
                 'pks_accum',bl_ais.pks_accum,'refitMSE1',refitMSE1,...
                 'refitMSE2',refitMSE2);
else
    % Boundary Candidates
    candMask = findBoundCand(find(instBlob),size(instBlob),fattenR);

    % Line Candidates
    lineCandMask = lineCands(candIm,candMask,gradThresh);

    % Boundary Lines
    [rho_init,theta_init] = boundLines(lineCandMask,thetaVals,rhoResolution,...
                                 suppress);

    % Refit Boundary Lines
    [y,x] = find(lineCandMask);
    [r1,t1] = refitBoundLine(rho_init(1),theta_init(1),x,y,inlierThresh);
    [r2,t2] = refitBoundLine(rho_init(2),theta_init(2),x,y,inlierThresh);
    rho = [r1,r2]; theta = [t1,t2];
end

midRho = mean(rho); midTheta = mean(theta);
lineMask = drawLineMask(size(instBlob),[rho,midRho],[theta,midTheta]);

if(nargout == 5)
    % Running Function in Test Mode Build Test Structure
    tStruct = struct('candMask',candMask,'lineCandMask',lineCandMask,...
                     'houghStruct',houghStruct);
end
