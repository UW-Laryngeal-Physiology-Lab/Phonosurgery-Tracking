% fitLineTLS_verify.m
%% Simulation Parameters
N = 500;
rhoVals = zeros(N,1); thetaVals = zeros(N,1);
error_rho = zeros(N,1); error_theta = zeros(N,1);

%% Run Simulation
for n = 1:N
    % Generate Line Coordinates
    rho = (-1)^(randi([0,1])) * randi([0,300]);
    theta = (-1)^(randi([0,1])) * (pi/2 * rand(1));

    if(abs(sin(theta)) < abs(cos(theta)))
        y = -1000:1000;
        x = (rho - y*sin(theta))/cos(theta);
    else
        x = -1000:1000;
        y = (rho - x*cos(theta))/sin(theta);
    end

    [rho_est,theta_est] = fitLineTLS(x,y);
    rhoVals(n,1) = rho; thetaVals(n,1) = theta;
    error_rho(n,1) = rho_est - rho;
    error_theta(n,1) = theta_est - theta;
end

%% Display Results
figure();
plot(1:N,error_rho); xlabel('Trial'); ylabel('Error');
title('Rho Error');

figure();
plot(1:N,error_theta); xlabel('Trial'); ylabel('Error');
title('Theta Error');



