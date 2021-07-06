# Phonosurgery-Tracking
Experimental software for tracking and analyzing movement of phono-micro-surgical instruments in simulated settings.

This folder contains files for running the Video Based Phonomicrosurgery Instrument Tracking System (V-PITS).  A description of the individual folders is given below


1) Calibration : Contains StereoScreenGrab.exe.  This application is used to capture frames simultaneously from the left and right camera of the stereo rig.  This tool is generally used during camera calibration to take images of the calibration grid.  The CamCalToolbox folder in this directory contains the Camera Calibration Toolbox for Matlab which is used to compute calibration parameters.

2) Developer Documentation : Contains documentation related to the mechanical setup of the Simulated Surgery Station and the electrical and network setup of the camera rig.

3) Installers : Contains installer for Huffyuv codec.  Videos using the system are encoded using Huffyuv which is a lossless codec.  Additionally, contains an installer for Basler Pylon.  This is the SDK for the cameras used in the camera rig.  This needs to be installed on the PC communicating with the cameras.

4) Matlab : Contains all Matlab routines used by V-PITS to process video data.  Primarily this is related to instrument tracking and 3D trajectory estimation.

5) User_Documentation : Contains documentation on using different parts of the V-PITS system.  DO NOT IGNORE THIS FOLDER

6) VideoCapture : Contains the vidCaptureGUI application.   This application is used to capture synchronized video from the stereo calibration rig.

7) VirtualDub : Contains the VirtualDub video player.  This is a very good program for viewing videos.
