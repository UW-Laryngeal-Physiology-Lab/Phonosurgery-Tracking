function classifyVec = threshClassifyFun(featureVec,threshVec,typeVec)
% THRESHCLASSIFYFUN Classifies feature vector based on thresholds
%
% RESULT = threshClassifyFun(FEATUREVEC,THRESHVEC,TYPEVEC) Classifies
% individual elements of (N x 1) vector FEATUREVEC by applying thresholding
% given by (N x 1) vectors THRESHVEC and TYPEVEC.  THRESHVEC contains the
% actual threshold values.  TYPEVEC determines how each element is
% thresholded.  It is a binary vector.  A value of 0 corresponds to
% a null hypothesis of feature <= threshold, 1 corresponds to null
% hypothesis of feature >= threshold.  CLASSIFYVEC is a binary vectory
% whose value is 1 if the null hypothesis is violated.

% Check Type Vector
assert(all((typeVec == 0) | (typeVec == 1)),'Incorrect Type Vector Value(s)');

classifyVec = false(size(featureVec));
type0 = (typeVec == 0); type1 = (typeVec == 1);

classifyVec(type0) = featureVec(type0) > threshVec(type0);
classifyVec(type1) = featureVec(type1) < threshVec(type1);

% Register underconfidence for NaNs
classifyVec(isnan(featureVec)) = 1;


