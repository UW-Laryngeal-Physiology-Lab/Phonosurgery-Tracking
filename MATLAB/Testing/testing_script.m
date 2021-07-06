%% Part 1
try
    trackFile_1 = 'C:\KS_Thesis\VideoDump\08_24_2011\analysis\08_24-2011_leftView_inst2.mat';
    tic();
    newTrackData_1 = trackWithOldData(trackFile_1);
    toc();
    save('newTrackData_1.mat','newTrackData_1');
catch
    disp('Error With TrackFile_1');
end

%{
try
    trackFile_2 = 'C:\KS_Thesis\VideoDump\08_24_2011\analysis\08_24-2011_rightView_inst2.mat';
    tic();
    newTrackData_2 = trackWithOldData(trackFile_2);
    toc();
    save('newTrackData_2',newTrackData_2);
catch
    disp('Error With TrackFile_2');
end
%}

%% Part 2
compareTrackData(newTrackData_1,trackFile_1);