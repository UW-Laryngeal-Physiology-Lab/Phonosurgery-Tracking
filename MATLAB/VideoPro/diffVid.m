function dVid = diffVid(vid)
% DIFFVID Computes Video of differences between frames
% 
% DVID = DIFFVID(VID) Computes a video of the difference between frames of
% (Height x Width x channels x frames) video matrix VID.  The result is a
% double that is returned as DVID.  Note differences can be negative or
% positive.

dVid = diff(double(vid),1,4);

end

