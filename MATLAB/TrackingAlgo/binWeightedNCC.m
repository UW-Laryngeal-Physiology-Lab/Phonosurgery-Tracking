function nccScore = binWeightedNCC(I,varargin)
% WEIGHTEDNCC Performs binary weighted template matching
%
% NCCSCORE = weightedNCC(I,T,W) Computes weighted normalized cross
% correlation of 2D template T with weighting matrix of same size W
% over I.  NCCSCORE is a matrix with same dimensions as I.  W is assumed
% binary valued (0 and 1's).
%
% NCCSCORE = weightedNCC(I,PRECOMPUTESTRUCT) Utilizes precomputed
% elements in struct PC for weighted template matching.  See needed fields
% in PRECOMPUTESTRUCT below.  PC should be generated using function
% wMatchBuildPcStruct.
%
% PRECOMPUTESTRUCT fields:
% t : Template
% w : Weighting Mask
% w2 : Squared weighting mask
% wSum : Sum over weighting mask
% w2Sum : Sum over squared weighting mask
% tMean : Weighted mean over weighted template
% tw : w .* (t - tMean)
% tw2 : w2 .* (t - tMean)
% tw2Sum : Sum over tw2
% tw_sqSum : Sum over tw.^2

if(nargin == 3)
    pc = struct('w',[],'t',[],'w2',[],'wSum',0,'w2Sum',0,'tMean',0,...
        'tw',[],'tw2',[],'tw2Sum',0,'tw_sqSum',0);
    pc.t = varargin{1}; pc.w = varargin{2};
    pc.w2 = pc.w.^2;
    pc.wSum = sum(pc.w(:));
    pc.w2Sum = sum(pc.w2(:));
    pc.tMean = (1/pc.wSum)*sum(pc.w(:) .* pc.t(:));
    pc.tw = pc.w .* (pc.t - pc.tMean);
    pc.tw2 = pc.w2 .* (pc.t - pc.tMean);
    pc.tw2Sum = sum(pc.tw2(:));
    pc.tw_sqSum = sum(pc.tw(:).^2);
elseif(nargin == 2)
    pc = varargin{1};
end

%% Precompute
% XCorrs
X_I_tw_sq = imfilter(I,pc.tw2,'same','corr');
X_I_w = imfilter(I,pc.w,'same','corr');
X_I_sq_w_sq = imfilter(I.^2,pc.w2,'same','corr');
%X_I_w_sq = imfilter(I,pc.w2,'same','corr');
X_I_w_sq = X_I_w; % This can be done because pc.w is binary.

IMean = (1/pc.wSum)*X_I_w;

%% Compute NCC Score
numerator = (X_I_tw_sq - IMean.*pc.tw2Sum);
denominator = ((X_I_sq_w_sq - 2.*IMean.*X_I_w_sq + (IMean.^2).*pc.w2Sum) .* ...
    (pc.tw_sqSum)).^(1/2);

% Set Low RMS Patches to NCC = 0
numerator(denominator < 1) = 0;
denominator(denominator < 1) = 1;

nccScore = numerator ./ denominator;