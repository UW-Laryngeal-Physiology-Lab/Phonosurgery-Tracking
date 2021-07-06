%% Compute
theta = -45:0.1:45;
theta = (pi/180)*theta;
imSize = [492,658];

maxDelta = zeros(imSize);
maxTheta = zeros(imSize);

for x = 0:imSize(2)-1
    for y = 0:imSize(1)-1
        dp_dt = -x*sin(theta) + y*cos(theta);
        [d1,max_idx] = max(abs(dp_dt));
        maxTheta(y+1,x+1) = theta(max_idx);
        maxDelta(y+1,x+1) = dp_dt(max_idx);
    end
end

%
figure();
[X,Y] = meshgrid(0:imSize(2)-1,0:imSize(1)-1);
mesh(X,Y,(0.1 * (pi/180))*maxDelta); view(0,90);
xlabel('X Pixel');ylabel('Y Pixel');
title('Max Delta @ \theta_{res} = 0.1');

figure();
[X,Y] = meshgrid(0:imSize(2)-1,0:imSize(1)-1);
mesh(X,Y,maxTheta); view(0,90);
xlabel('X Pixel');ylabel('Y Pixel');
title('Max Delta Theta Val');
%}

disp(['drho for dtheta = 0.1 degrees is ' num2str(0.1*(pi/180)*max(maxDelta(:))) ' pixels']);