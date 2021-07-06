function [Hpar_L,Hpar_R,match_R,T_L,R_L] = ...
    parAccum(leftEdgeIm,rightEdgeIm,delRho,delTheta,voteThresh,thetaBounds)
% PARACCUM Builds parallel line Hough Accumulator
%
% [HPAR_L,HPAR_R,MATCH_R,T_L,R_L] = 
%   parAccum(LEFTEDGEIM,RIGHTEDGEIM,DELRHO,DELTHETA,VOTETHRESH,THETABOUNDS)
% Uses binary images of left and right edgels LEFTEDGEIM AND RIGHTEDGEIM to
% build a Hough accumulators HPAR_L & HPAR_R with parallel line pair votes.
% Values in HPAR_L and HPAR_R correspond to the number of votes for the
% left and right line in the pair respectively. Locations in HPAR_L
% correspond to parameterizations of the left line in the parallel pair.
% The rows of HPAR_L correspond to rho values in R_L and the columns theta
% values (in degrees) in T_L.  The parameterization for the corresponding
% right line in the pair is found in MATCH_R (numel(R_L) x numel(T_L) x 2).
% Where the first 2 dimensions provide correspondence to HPAR_L and the
% third dimension is [theta,rho].  The search for parallel line pairs is
% controlled by DELRHO, DELTHETA, and VOTETHRESH.  See below for
% description of these parameters.  THETABOUNDS [thetaMin,thetaMax]
% determines the range of possible line orientations in degrees and is
% optional.
%
% Algorithm Parameters
% DELRHO [rhoAddMin,rhoAddMax] The range over rho to add to detected left
% line parameterizaiton that defines the search space for the corresponding
% right pair.
% DELTHEATA [thetaAddMin,thetaAddMax] The range over theta to add to the
% detected left line parameterization that defines the search space for the
% corresponding right line pair.
% VOTETHRESH : Scalar parameter that determines line candidates.  Hough
% accumulator bins with more votes than VOTETHRESH are candidate
% parameterizations for line pairs.
%% Input Arguments
if(nargin == 5)
    thetaBounds = [-45,45];
end

%% Compute Hough Transforms
rhoRes = 1; thetaRes = 0.1;
theta = thetaBounds(1):thetaRes:thetaBounds(2);

[H_L,T_L,R_L] = hough(leftEdgeIm,'RhoResolution',rhoRes,...
                               'Theta',theta);
                                                    
[H_R,T_R,R_R] = hough(rightEdgeIm,'RhoResolution',rhoRes,...
                               'Theta',theta);

Hpar_L = zeros(size(H_L));
Hpar_R = zeros(size(H_L));
match_R = zeros([size(H_L) 2]);
                           
nonZeroIdx = find(H_L > voteThresh);
numNonZero = numel(nonZeroIdx);
%parLineList = zeros(numNonZero,6);

%% Assign Parallel Pairs
% Convert Search Delta to Indices
delRhoIdx = (round(delRho(1)/rhoRes):1:round(delRho(2)/rhoRes)).';
delThetaIdx = (round(delTheta(1)/thetaRes):1:round(delTheta(2)/thetaRes)).';

searchRow = zeros(numel(delRhoIdx),1);
searchCol = zeros(numel(delThetaIdx),1);
numElem = numel(searchRow)*numel(searchCol);

% Generate Padded Array for Search
padRows = numel(delRhoIdx); padCols = numel(delThetaIdx);
H_R = padAccum(H_R,padRows,padCols);

hSize = size(H_L);
[rowNonZero,colNonZero] = ind2sub(hSize,nonZeroIdx);
rowNonZero2 = rowNonZero + padRows;
colNonZero2 = colNonZero + padCols;


vPrimeVals = zeros(numNonZero,1);
maxIdxVals = zeros(numNonZero,1);
for k = 1:numNonZero
    row = rowNonZero2(k); col = colNonZero2(k);
    
    searchRow(:) = row + delRhoIdx;
    searchCol(:) = col + delThetaIdx;
    
    %searchRow = max([1,row+delRhoIdx(1)]):min([...
    %                                    hSize(1),row+delRhoIdx(2)]);
    %searchCol = max([1,col+delThetaIdx(1)]):min([...
    %                                  hSize(2),col+delThetaIdx(2)]);
                                  
    %searchSpace = H_R(searchRow,searchCol);
    %[Vprime,maxIdx] = max(searchSpace(:)); 
    %Vprime = Vprime(1); maxIdx = maxIdx(1);
    %[maxRelRow,maxRelCol] = ind2sub([numel(searchRow),numel(searchCol)],...
    %                                 maxIdx);
    [Vprime,maxIdx] = max(reshape(H_R(searchRow,searchCol),[numElem,1]));
    vPrimeVals(k,1) = Vprime;
    maxIdxVals(k,1) = maxIdx;
    
    %match_R(row,col,:) = [T_R(searchCol(maxRelCol)),R_R(searchRow(maxRelRow))];
    %Hpar_L(row,col) = H_L(row,col) + Vprime;
    
    %parLineList(k,:) = [R_L(row),T_L(col),R_R(searchRow(maxRelRow)),...
    %    T_R(searchCol(maxRelCol)),H_L(row,col),Vprime];
end

%% Build Output Matrices
Hpar_L(nonZeroIdx) = H_L(nonZeroIdx);
Hpar_R(nonZeroIdx) = vPrimeVals;

[maxRelRow,maxRelCol] = ...
    ind2sub([numel(delRhoIdx),numel(delThetaIdx)],maxIdxVals);
rhoMatchIdx = delRhoIdx(maxRelRow) + rowNonZero;
thetaMatchIdx = delThetaIdx(maxRelCol) + colNonZero;
valid = (rhoMatchIdx > 0) & (rhoMatchIdx <= hSize(1)) & ...
        (thetaMatchIdx > 0) & (thetaMatchIdx <= hSize(2));

match_R(nonZeroIdx(valid) + numel(H_L)) = R_R(rhoMatchIdx(valid));
match_R(nonZeroIdx(valid)) = T_R(thetaMatchIdx(valid));

function H_pad = padAccum(H,padRows,padCols)
% PADACCUM Pads accumulator with zeros
%
% H_PAD = padAccum(H,padRows,padCols) Pads accumulator matrix H (numRho x
% numTheta) by padRows and padCols along each appropriate edge with zeros
% such that H_pad has size (2*padRows + size(H,1) x 2*padCols + size(H,2)).
aSize = size(H);
H_pad = zeros(2*padRows + aSize(1),2*padCols + aSize(2));
H_pad(padRows + (1:aSize(1)),padCols + (1:aSize(2))) = H;







