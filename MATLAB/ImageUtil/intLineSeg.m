classdef intLineSeg
    %Immutable Properties
    properties (SetAccess = immutable)
        %EndPoint Locations [x,y]
        pt1; pt2;

        %Gradient Related
        dx; dy; domDir;   
    end

    %Constructors
    %intLineSeg(XY_1,XY_2)
    methods
        %Class Constructor
        function obj = intLineSeg(xy_1,xy_2)
            obj.pt1 = xy_1;
            obj.pt2 = xy_2;

            [obj.dx obj.dy,obj.domDir] = obj.analyzeEndPoints(xy_1,xy_2);
        end
    end

    %Static Methods Used By the Class
    %[DX DY DOMDIR] = analyzeEndPoints(XY_1,XY_2)
    methods(Static=true,Access=protected)
        function [dx dy domDir] = analyzeEndPoints(xy_1,xy_2)
        %[DX DY DOMDIR] = analyzeEndpoints(xy_1,xy_2)
            dx = xy_2(1) - xy_1(1);
            dy = xy_2(2) - xy_1(2);

            if abs(dx) >= abs(dy)
                domDir = 'x';
            else
                domDir = 'y';
            end
        end
    end

    %Methods Used by the Outside
    %[XIDX,YIDX] = getIntLocs()
    methods
        function [xIdx yIdx] = getIntLocs(obj)
        %[XIDX,YIDX] = GETINTLOCS Finds Integer (x,y) locations associated
        %with line.  The points are choosen such that the error at each
        %point is guaranteed w/in +/- 0.5 along the non-dominant direction.
        %KS 12/17/2010
            if strcmp(obj.domDir,'x')
                %Walk Along X Direction
                m = obj.dy/obj.dx;
                xIdx = (obj.pt1(1):sign(obj.dx):obj.pt2(1)).';
                yIdx = round(m*(xIdx-obj.pt1(1)) + obj.pt1(2));
            else
                %Walk Along Y Direction
                m = obj.dx/obj.dy;
                yIdx = (obj.pt1(2):sign(obj.dy):obj.pt2(2)).';
                xIdx = round(m*(yIdx-obj.pt1(2)) + obj.pt1(1));
            end
        end
        function [midParam xIdx yIdx] = getParameterization(obj)
        %[MIDPARAM XIDX YIDX] = GETPARAMETERIZATION Returns an array
        %MIDPARAM with an element corresponding to each pixel along the
        %midline returned by [XIDX YIDX].  Each element contains the
        %parameterization of the midline pixel on [0,1] corresponding to
        %midline length.  Therefore the first row index of [XIDX YIDX]
        %has parameterization 0 and the row index parameterization 1.  In
        %between these endpoints the parameterization is the distance from
        %the 0 point divided by the total length.
        
            % Get Locations
            [xIdx yIdx] = obj.getIntLocs();

            %Compute Adjacent Distances
            cumDistVec = cumsum(sqrt(sum(diff([xIdx yIdx],1,1).^2,2)));

            %Compute Parameterization
            midParam = [0;(cumDistVec./cumDistVec(end))];
        end
        function ortho = getOrthoDir(obj)
        %ORTHO = GETOORTHODIR Returns a vector [dx dy] representing the
        %direction orthogonal to the midline.
            ortho = [0 1;-1 0] * [obj.dx;obj.dy];
        end 
    end
end
