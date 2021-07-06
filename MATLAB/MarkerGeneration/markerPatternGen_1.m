function markerPatternGen_1(fileName)
% MARKERPATTERNGEN_1 Generates image of marker pattern #1
%
% markerPatternGen_1(FILENAME) Generates image of Pattern 1

%% Parameters
cols = ceil(72 * (3.14 * pi * 0.1 * (1/2.54)));

pat = false(17,cols);
pat([6,9,12],:) = 1;

imwrite(pat,fileName);
