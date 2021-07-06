% viewAnalysis_Script.m
%% Load & Parameters
clear all;
frames = [3000,3100];
[fName,pName] = uigetfile('*.mat','Load Analysis File');
load(fullfile(pName,fName),'saveStruct');

%% Grab Needed Structures
vidObj = mmreader(saveStruct.vidName);
mark = saveStruct.mark;
inst = saveStruct.inst;
label = saveStruct.label;
algoInfo = saveStruct.algoInfo;
pc = saveStruct.pc;
backIm = saveStruct.backIm;

%% Code We Want to Profile
profileMe;
