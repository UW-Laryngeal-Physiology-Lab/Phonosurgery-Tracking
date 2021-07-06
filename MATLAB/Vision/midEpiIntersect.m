function im2Pts = midEpiIntersect(im1Pts,F,im2Midlines)
% MIDEPIINTERSECT 
%
% IM2PTS = midEpiIntersect(im1Pts,F,im2Midlines)

epiLines = (F * [im1Pts';ones(1,size(im1Pts,1))])';

Cstart = cross(epiLines,[cos(im2Midlines(:,2)),sin(im2Midlines(:,2)),...
                         -im2Midlines(:,1)]);
im2Pts = Cstart(:,1:2) ./ repmat(Cstart(:,3),1,2);