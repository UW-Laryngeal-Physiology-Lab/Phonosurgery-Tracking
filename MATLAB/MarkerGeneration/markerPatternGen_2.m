function markerPatternGen_2(fileName)
% MARKERPATTERNGEN_2 Generates image of marker pattern #2
%
% markerPatternGen_2(FILENAME) Generates image of Pattern 2

%% Parameters
cols = ceil(72 * (3.14 * pi * 0.1 * (1/2.54)));

pat = false(17,cols);
pat([6,8,10,12],:) = 1;

imwrite(pat,fileName);