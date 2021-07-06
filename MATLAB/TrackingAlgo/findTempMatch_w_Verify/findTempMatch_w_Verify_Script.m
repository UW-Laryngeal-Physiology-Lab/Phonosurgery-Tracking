% findTempMatch_w_Verify_Script.m
%% User Parameters
% Noise
noiseSigma = 3;

% Template Size
wSize = 79;

% Search Window
nMultip = 2;

%% Simulation Parameters
% Image
nRows = 500; nCols = 500;

% Gaussian
sigmaX = round(wSize/4); sigmaY = round(wSize/4);
X0 = round((nCols-1)/2); Y0 = round((nRows-1)/2);

%% Initialize Template
% Generate Initialization Image
gIm = uint8(noiseSigma*randn(nRows,nCols)) + ...
    gaussian2dIm(X0,Y0,sigmaX,sigmaY,nRows,nCols);

% Initialize Template
[tPos,tSize,pc] = initWeightedTemplate(gIm,79);

nSize = round(tSize*nMultip);
if(mod(nSize,2) == 0)
% Define Neighborhood Search Window
    nSize = nSize + 1;
end
nCorn = tPos - ((nSize-1)/2 - (tSize-1)/2);

%% Run Test
shiftDelta = 5.5;
shiftMax = floor(0.9 * ((nSize-1)/2 - (tSize-1)/2));
shift = -shiftMax:shiftDelta:shiftMax;
matchScore = zeros(numel(shift).^2,1);
fitMSE = zeros(numel(shift).^2,1);
errorSig = zeros(numel(shift).^2,2);

count = 0;
for xShift = shift
    disp(xShift);
    for yShift = shift;
        count = count + 1;
        
        % Generate Shifted Gaussian Image
        gIm = uint8(noiseSigma*randn(nRows,nCols)) + ...
            gaussian2dIm(X0 + xShift,Y0 + yShift,sigmaX,sigmaY,...
                          nRows,nCols);
                      
        % Perform Weighted Template Matching & Compute Estimated
        % Displacement
        [estPos,~,score,fMSE] = findTempMatch_w(gIm,pc,nCorn,nSize);
        estDisp = estPos - tPos;
        
        % Computet Error in Displacement Estimate
        errorSig(count,:) = estDisp - [xShift,yShift];
        matchScore(count) = score;
        fitMSE(count) = fMSE;
    end
end

%% Display Results
figure();
subplot(2,1,1); plot(errorSig(:,1)); title('X Error');
subplot(2,1,2); plot(errorSig(:,2)); title('Y Error');

figure();
subplot(2,1,1); plot(matchScore); title('Match Score');
subplot(2,1,2); plot(fitMSE); title('Fit MSE');


