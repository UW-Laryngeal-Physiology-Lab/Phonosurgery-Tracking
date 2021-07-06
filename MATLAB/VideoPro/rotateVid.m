function rotVid = rotateVid(sourceVid,flipDeg)
%ROTATEVID Rotates all frames of a video file
%
% ROTVID = rotateVid(SOURCEVID,FLIPDEG) Treats N-D matrix FLIPDEG as a
% video and rotates its channels counter-clockwise by a multiple of 90.  If
% flipDeg is not +/- {90,180,270} SOURCEVID is returned with nothing done.

switch flipDeg
    case {90,-270}
        rotVid = permute(sourceVid,[2,1,3,4]);
        rotVid = flipdim(rotVid,1);
        
    case {180,-180}
        rotVid = flipdim(sourceVid,1);
        
    case {270,-90}
        rotVid = permute(sourceVid,[2,1,3,4]);
        rotVid = flipdim(rotVid,2);
        
    otherwise
        disp('Function Only Rotates by a 90 Degree Multiple :(');
end
end

