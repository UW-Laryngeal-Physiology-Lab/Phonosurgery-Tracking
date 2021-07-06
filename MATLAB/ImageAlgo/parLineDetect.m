function [rho,theta,ais] = parLineDetect(leftEdge,rightEdge,pStruct)
% PARLINEDETECT Detects parallel instrument boundary line pair
%
% [RHO,THETA,AIS] = parLineDetect(LEFTEDGE,RIGHTEDGE,PSTRUCT) Detects
% dominant parallel line pair using left & right edges pixels in an image.
% LEFTEDGE and RIGHTEDGE are binary images indicating the location of edge
% pixels. PSTRUCT is an optional structure with parameters controlling the
% algorithm.  The Hough transform is utilized for line detection.
%
% PSTRUCT fields
% thetaBounds {[-45,45]} : Bounds on theta used for Hough Accumulator in 
% degrees.  Basing this on past boundary line orientation estimates can
% reduce processing time (smaller range) and utilize this function within a
% tracking framework.
%
% thetaRes {0.1} : Hough Accumulator Theta Resolution
%
% delRho {[30,40]}: Along with delTheta determines bounds on search space
% in Hough Acumulator for finding right line pair to left line.  Should be
% based on the range of expected "width" of the instrument.
%
% delTheta {[-3,3]} : See description of delRho, this is the other
% parameter that determines the search space it is based on the angle
% variation in boundary line pairs (deviation from parallel) expected.  
%
% voteThresh {10} : Pixel count used to vote for potential parallel line
% pair candidates

%% Input Arguments
if(nargin == 2)
    pStruct = struct();
end

if(~isfield(pStruct,'thetaRes'))
    pStruct.thetaRes = 0.1;
end

if(~isfield(pStruct,'delRho'))
    pStruct.delRho = [30,40];
end

if(~isfield(pStruct,'delTheta'))
    pStruct.delTheta = [-3,3];
end

if(~isfield(pStruct,'thetaBounds'))
    pStruct.thetaBounds = [-45,45];
end

if(~isfield(pStruct,'voteThresh'))
    pStruct.voteThresh = 10;
end

%% Run Algorithm
% Build Parallel Line Pair Accumulators
[H_L,H_R,match,T,R] = parAccum(leftEdge,rightEdge,pStruct.delRho,...
    pStruct.delTheta,pStruct.voteThresh,pStruct.thetaBounds);

% Find Peak in Accumulator
[~,maxIdx] = max(H_L(:) + H_R(:)); 
[maxRhoIdx,maxThetaIdx] = ind2sub(size(H_L),maxIdx);

% Determine Line Pair Hough Parameters
peakMatchRhoIdx = sub2ind(size(match),maxRhoIdx,maxThetaIdx,2);
peakMatchThetaIdx = sub2ind(size(match),maxRhoIdx,maxThetaIdx,1);
rho = [R(maxRhoIdx),match(peakMatchRhoIdx).'];
theta = (pi/180) * [T(maxThetaIdx),match(peakMatchThetaIdx).'];

if(nargout == 3)
    % Algorithm Info Struct
    ais = struct();
    ais.votesLeft = H_L(maxIdx);
    ais.votesRight = H_R(maxIdx);
end
