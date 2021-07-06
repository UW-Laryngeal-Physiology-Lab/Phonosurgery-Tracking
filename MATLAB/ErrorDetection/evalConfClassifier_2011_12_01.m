% evalConfClassifier_2011_11_18.m
% Evaluates instrument tracking algorithm confidence classifier using left
% and right view tracking examples.  Classifier is based on thresholding of
% parameters calculated using tracked instrument features at the current
% and previous frame.  Script generates an overview spreadsheet of the
% evaluation.  

%% Classifier Thresholds
% Thresholds should have been selected using a training data set.  This
% script is written to be compatible with confClassifierDesign_2011_11_18.

t_delta_NCC = 0.033178;
t_delta_width = 1.7152;
t_delta_orient = 0.0098665;

%% Load Tracking Data and Get Confidence Parameters
% Load tracking files of a single instrument in the left and right view.
% The tracking should have been performed using the UGT (user-guided
% tracking) application.  Correct tracking and tracking error labelling
% should be manually verified. After loading the tracking data for the left
% and right view, the value of the confidence parameters at each correctly
% tracked frame and tracking error frame is extracted.  

% Left Tracking File
[fName_L,pName_L] = uigetfile('*.mat','Load Left Tracking File');
[corr_L,err_L] = getConfParamData2(fullfile(pName_L,fName_L));

% Right Tracking File
[fName_R,pName_R] = uigetfile('*.mat','Load Right Tracking File',pName_L);
[corr_R,err_R] = getConfParamData2(fullfile(pName_R,fName_R));

%% Build Left-Right View Joint Paramter Array
%  Arrays consisting of left and right view confidence parameter data are
%  built.  Each array consists of confidence parameter values for correctly
%  tracked frames (with suffix _corr) or tracking error frames (with suffix
%  _err).  

% Joint Correct Arrays ----------------------------------------------------
% Delta Correct Tracking Parameters
delta_NCC_corr = [corr_L.delta_NCC_corr;corr_R.delta_NCC_corr];
delta_width_corr = [corr_L.delta_width_corr;corr_R.delta_width_corr];
delta_orient_corr = [corr_L.delta_orient_corr;corr_R.delta_orient_corr];
delta_inliersDistL_corr = [corr_L.delta_inliersDistL_corr;corr_R.delta_inliersDistL_corr];
delta_inliersDistR_corr = [corr_L.delta_inliersDistR_corr;corr_R.delta_inliersDistR_corr];
%delta_inliers_corr = delta_inliersL_corr + delta_inliersR_corr;

% Raw Correct Parameters
inliersL_corr = [corr_L.inliersL_corr;corr_R.inliersL_corr];
inliersR_corr = [corr_L.inliersR_corr;corr_R.inliersR_corr];

% Joint Tracking Error Arrays ---------------------------------------------
% Delta Tracking Erro Parameters
delta_NCC_err = [err_L.delta_NCC_error;err_R.delta_NCC_error];
delta_width_err = [err_L.delta_width_error;err_R.delta_width_error];
delta_orient_err = [err_L.delta_orient_error;err_R.delta_orient_error];
delta_inliersDistL_err = [err_L.delta_inliersDistL_error;err_R.delta_inliersDistL_error];
delta_inliersDistR_err = [err_L.delta_inliersDistR_error;err_R.delta_inliersDistR_error];
%delta_inliers_err = delta_inliersL_err + delta_inliersR_err;

% Raw Error Parameters
inliersL_err = [err_L.inliersL_error;err_R.inliersL_error];
inliersR_err = [err_L.inliersR_error;err_R.inliersR_error];

%% Evaluate Classifier
% Evaluate the classifier using standard machine learning parameters :
% recall, precision, and TNR (true negative rate).  Validation of the
% classifier is given by a high TNR (should classify almost all, if not all
% of the tracking errors correctly). Recall describes how effective the
% classifier is.  Low precision means that the algorithm will require more
% user intervention. I usually hope to see around 95%.

% Determine the whether the classifier correctly or incorrectly classifies
% each correctly tracked frame
class_corr = delta_NCC_corr <= t_delta_NCC & ...
               delta_width_corr <= t_delta_width & ...
               delta_orient_corr <= t_delta_orient & ...
               inliersL_corr >= 10 & ...
               inliersR_corr >= 10;

% Determine the whether the classifier correctly or incorrectly classifies
% each tracking error frame
class_err = delta_NCC_err > t_delta_NCC |...
               delta_width_err > t_delta_width | ...
               delta_orient_err > t_delta_orient | ...
               inliersL_err < 10 | ...
               inliersR_err < 10;

% tp : True Positive Rate | fp : False Positive Rate
% tn : True Negative Rate | fn : False Negative Rate
tp = sum(class_corr); fp = sum(~class_err);
tn = sum(class_err); fn = sum(~class_corr);

recall = tp/(tp + fn);
precision = tp/(tp + fp);
TNR = tn/(tn + fp);

disp(['Recall : ' num2str(recall)]);
disp(['Precision : ' num2str(precision)]);
fprintf('TNR : %u/%u -> %1.4f\n', tn, tn + fp, TNR);

%% Generate Overview Spreadsheet 
% The user is given an option to save spreadsheet that
% saves classifier summary data along with the tracking error frame numbers
% that it misclassified as confident.
%
% It is good to examine the tracking error frames that it classifies as
% being confident about.  They may have been incorrectly labelled as
% tracking error, when in fact they are correctly tracked frames.  Remember
% the confidence classifier is designed based on smoothness.  

summaryPage = {date,[];
               'Left File',fullfile(pName_L,fName_L);...
               'Right File',fullfile(pName_R,fName_R);...
               'Correctly Tracked Frames Used',tp+fn;...
               'Tracking Error Frames Used',tn+fp};...

classifierOverview = {'Confidence Parameter','Threshold Value',nan;...
                      'delta_NCC',t_delta_NCC,nan;...
                      'delta_width',t_delta_width,nan;...
                      'delta_orient',t_delta_orient,nan;...
                      'inliersL',10,nan;...
                      'inliersR',10,nan;...
                      nan,nan,nan;...
                      'Recall',recall,nan;...
                      'Precision',precision,nan;...
                      'TNR (correct/total)',tn,tn+fp};

% Misclassied Tracking Errors ---------------------------------------------
% Get Frame Numbers associated with tracking errors
frames_error_L = err_L.frames_error; num_error_L = numel(frames_error_L);
frames_error_R = err_R.frames_error;

% Break tracking error classification into left and right view
class_err_L = class_err(1:num_error_L);
class_err_R = class_err((1+num_error_L):end);

% Get Frame Numbers corresponding to tracking error frames that are
% classified as confident.
getFrameNums = @(classArray,frameArray,label) frameArray(classArray == label);
misclassify_err_L = getFrameNums(class_err_L,frames_error_L,0);
misclassify_err_R = getFrameNums(class_err_R,frames_error_R,0);

errorMisclassify = [{'Left View Tracking'};num2cell(misclassify_err_L(:));...
                    {'Right View Tracking'}; num2cell(misclassify_err_R(:))];
                  
                  
% Save Spreadsheet --------------------------------------------------------
[xlFile,xlPath] = uiputfile('*.xls','Save Classifier Info ');

if(xlFile ~= 0)
    xlFileName = fullfile(xlPath,xlFile);
    xlswrite(xlFileName,summaryPage,'Summary');
    xlswrite(xlFileName,classifierOverview,'Classifier');
    xlswrite(xlFileName,errorMisclassify,'Error_Misclassify');
end
