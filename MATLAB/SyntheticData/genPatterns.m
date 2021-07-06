function [instPattern,backPattern] = genPatterns(instStats,backStats,imSize)

%% Sythesize Texture Images
% Background
backPattern = uint8(repmat(backStats.mean,imSize) + ...
    backStats.std*randn(imSize));

% Instrument
instPattern = uint8(repmat(instStats.mean,imSize) + ...
    instStats.std*randn(imSize));
