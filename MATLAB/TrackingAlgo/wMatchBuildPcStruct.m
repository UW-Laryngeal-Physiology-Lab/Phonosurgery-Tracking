function pc = wMatchBuildPcStruct(t,w)
% WMATCHBUILDPCSTRUCT Builds precomputes variabled for windowed matching
%
% PC = wMatchBuildPcStruct(T,W) Utilizes (MxN) template T and {MxN)
% weighting window W to calculate all precomputable components for weighted
% Normalized Cross Correlation. W is a binary window the same size of T
% where 1's indicate pixels to be used for matching and 0's indicate pixels
% that should not be used for matching.  A structure PC is returned, see
% below for a description of this structure.
%
% PC fields: (all fields stored as type double)
% w : Weighting window 
% t : Template
% w2 : Weighting Window Squared 
% wSum : Sum of all elements in weighting window
% w2Sum : Sum of all elements in squared weighting window
% tMean : Weighted mean of template window
% tw : Weighted mean removed template
% tw2 : Mean removed template weighted by w2
% tw2Sum : Sum of all elements in tw2
% tw_sqSum : Sum of all elements in tw squared.

pc = struct('w',[],'t',[],'w2',[],'wSum',0,'w2Sum',0,'tMean',0,...
    'tw',[],'tw2',[],'tw2Sum',0,'tw_sqSum',0);
pc.w = double(w); pc.t = double(t);
pc.w2 = pc.w.^2;
pc.wSum = sum(pc.w(:));
pc.w2Sum = sum(pc.w2(:));
pc.tMean = (1/pc.wSum)*sum(pc.w(:) .* pc.t(:));
pc.tw = pc.w .* (pc.t - pc.tMean);
pc.tw2 = pc.w2 .* (pc.t - pc.tMean);
pc.tw2Sum = sum(pc.tw2(:));
pc.tw_sqSum = sum(pc.tw(:).^2);