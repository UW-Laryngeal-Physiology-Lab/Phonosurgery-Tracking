function parLineList = parLineFind(leftEdgeIm,rightEdgeIm,delRho,delTheta,voteThresh)
% PARLINEFIND Finds parallel line candidates
%
% PARLINELIST = parLinFind(LEFTEDGEIM,RIGHTEDGEIM,DELRHO,DELTHETA)

%% Compute Hough Transforms
rhoRes = 1; thetaRes = 0.1;
theta = -45:thetaRes:45;

[H_L,T_L,R_L] = hough(leftEdgeIm,'RhoResolution',rhoRes,...
                               'Theta',theta);
                                                    
[H_R,T_R,R_R] = hough(rightEdgeIm,'RhoResolution',rhoRes,...
                               'Theta',theta);

nonZeroIdx = find(H_L > voteThresh);
numNonZero = numel(nonZeroIdx);
parLineList = zeros(numNonZero,6);

%% Assign Parallel Pairs
hSize = size(H_L);

% Convert Search Delta to Indices
delRhoIdx = round(delRho/rhoRes);
delThetaIdx = round(delTheta/thetaRes);

for k = 1:numNonZero
    [row,col] = ind2sub(hSize,nonZeroIdx(k));
    
    searchRow = max([1,row+delRhoIdx(1)]):min([...
                                        hSize(1),row+delRhoIdx(2)]);
    searchCol = max([1,col+delThetaIdx(1)]):min([...
                                      hSize(2),col+delThetaIdx(2)]);
                                  
    searchSpace = H_R(searchRow,searchCol);
    [Vprime,maxIdx] = max(searchSpace(:)); 
    Vprime = Vprime(1); maxIdx = maxIdx(1);
    [maxRelRow,maxRelCol] = ind2sub([numel(searchRow),numel(searchCol)],...
                                     maxIdx);
     
    parLineList(k,:) = [R_L(row),T_L(col),R_R(searchRow(maxRelRow)),...
        T_R(searchCol(maxRelCol)),H_L(row,col),Vprime];
end


                           