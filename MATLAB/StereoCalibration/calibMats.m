function [PL,PR,F] = calibMats(calibParameters)
% CALIBMATS Computes camera calib. matrices from calib. parameters
%
% [PL,PR,F] = calibMats(CALIBPARAMETERS) Computes camera projection
% matrices PL & PR (left & right) and corresponding fundamental matrix F
% from stereo calibration parameters generated using the Caltech Camera
% Calibration toolbox (http://www.vision.caltech.edu/bouguetj/calib_doc/).
% All parameters should be stored in struct CALIBPARAMETERS.  See below for
% form of CALIBPARAMETERS structure.  The camera calibration matrices are
% computed such that the coordinate frame of the left camera is the world
% frame.
%
% CALIBPARAMETERS Form
% Intrinsic Calibration Fields :
% fc_left : 2 elements array with left camera focal length parameters
% fc_right : 2 elements array with right camera focal length parameters
% cc_left : 2 element array with left camera principal point parameters
% cc_right : 2 elements array with right camera principal point parameters
%
% Extrinsic Calibration Fields
% om : rodrigues rotation vector (3 element array) that defines the
% orientation component of the rigid body transformation taking a vector in
% the left camera frame to the right camera frame.
% T : 3 element translation vector that defines the translation component
% of the rigid body transformation taking a point/vector in the left camera
% frame to the right camera frame.

% Deal Individual Parameters
[fc_left,fc_right,cc_left,cc_right,om,T] = ...
    deal(calibParameters.fc_left,calibParameters.fc_right,...
         calibParameters.cc_left,calibParameters.cc_right,...
         calibParameters.om,calibParameters.T);
 
%% Build Projection Matrices
K1 = [fc_left(1),0,cc_left(1);...
      0,fc_left(2),cc_left(2);...
      0,0,1];
PL = K1 * [eye(3,3),zeros(3,1)];
  
R = rodrigues(om);
K2 = [fc_right(1),0,cc_right(1);...
      0,fc_right(2),cc_right(2);...
      0,0,1];
PR = K2 * [R,T];
%% Compute Fundamental Matrix
e2 = PR * [0;0;0;1];
e2_x = [0,-e2(3),e2(2);...
        e2(3),0,-e2(1);...
        -e2(2),e2(1),0];
F = e2_x * (PR * pinv(PL));
