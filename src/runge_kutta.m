clear; close all; clc;
g = 9.80665; % m / s^2;

%%                               PHYSICAL PARAMETERS PENDULUMS
%%                           (- YOU CAN MODIFY THESE PARAMETERS -)

% 1ST PENDULUM
L1_1 = 2.2; L2_1 = 1.3; % m;
m1_1 = 1.2; m2_1 = 1.7; % kg;
% 2ND PENDULUM
L1_2 = 1.8; L2_2 = 1.0; % m;
m1_2 = 1.0; m2_2 = 1.5; % kg;
% 3RD PENDULUM
L1_3 = 2.5; L2_3 = 1.8; % m;
m1_3 = 1.5; m2_3 = 2.0; % kg;

%%                                   INITIAL CONDITIONS

% FIRST
theta10_1 = 50; theta20_1 = 70; % degrees; initial angles
theta10_1 = theta10_1 * pi / 180; theta20_1 = theta20_1 * pi / 180; % conversion to radians
OM10_1 =  - 30; OM20_1 = 20; % degrees / s; initial angular velocities
OM10_1 = OM10_1 * pi / 180; OM20_1 = OM20_1 * pi / 180; % conversion to rad / s

% SECOND
theta10_2 = 65; theta20_2 = 80; % degrees; initial angles
theta10_2 = theta10_2 * pi / 180; theta20_2 = theta20_2 * pi / 180; % conversion to radians
OM10_2 =  - 20; OM20_2 = 25; % degrees / s; initial angular velocities
OM10_2 = OM10_2 * pi / 180; OM20_2 = OM20_2 * pi / 180; % conversion to rad / s

% THIRD
theta10_3 = 70; theta20_3 = 90; % degrees; initial angles
theta10_3 = theta10_3 * pi / 180; theta20_3 = theta20_3 * pi / 180; % conversion to radians
OM10_3 = 10; OM20_3 =  - 15; % degrees / s; initial angular velocities
OM10_3 = OM10_3 * pi / 180; OM20_3 = OM20_3 * pi / 180; % conversion to rad / s

%%                           DURATION: frequency and period components 

omega1_1 = sqrt(g / L1_1); omega2_1 = sqrt(g / L2_1);
T1_1 = 2 * pi / omega1_1; T2_1 = 2 * pi / omega2_1;
omega1_2 = sqrt(g / L1_2); omega2_2 = sqrt(g / L2_2);
T1_2 = 2 * pi / omega1_2; T2_2 = 2 * pi / omega2_2;
omega1_3 = sqrt(g / L1_3); omega2_3 = sqrt(g / L2_3);
T1_3 = 2 * pi / omega1_3; T2_3 = 2 * pi / omega2_3;

T = max([T1_1,T2_1,T1_2,T2_2, T1_3,T2_3]); % characteristic time of the double pendulum motion
ti = 0; tf = 5 * T; N = 200000; t = linspace(ti,tf,N); dt = t(2) - t(1); % discrete time

% Pre - allocation and initial values for the first pendulum:
theta1_1 = zeros(1,N); theta2_1 = theta1_1;
OM1_1 = zeros(1,N); OM2_1 = OM1_1;

theta1_1(1) = theta10_1; theta2_1(1) = theta20_1;
theta1_1(2) = theta10_1 + OM10_1 * dt; theta2_1(2) = theta20_1 + OM20_1 * dt;
OM1_1(1) = OM10_1; OM2_1(1) = OM20_1; 

% Pre - allocation and initial values for the second pendulum:
theta1_2 = zeros(1,N); theta2_2 = theta1_2;
OM1_2 = zeros(1,N); OM2_2 = OM1_2;

theta1_2(1) = theta10_2; theta2_2(1) = theta20_2;
theta1_2(2) = theta10_2 + OM10_2 * dt; theta2_2(2) = theta20_2 + OM20_2 * dt;
OM1_2(1) = OM10_2; OM2_2(1) = OM20_2; 

% Pre - allocation and initial values for the third pendulum:
theta1_3 = zeros(1,N); theta2_3 = theta1_3;
OM1_3 = zeros(1,N); OM2_3 = OM1_3;

theta1_3(1) = theta10_3; theta2_3(1) = theta20_3;
theta1_3(2) = theta10_3 + OM10_3 * dt; theta2_3(2) = theta20_3 + OM20_3 * dt;
OM1_3(1) = OM10_3; OM2_3(1) = OM20_3;

% Helper notations: miu, r - coefficients, a - motion coefficients 
miu_1 = 1 + m1_1 / m2_1; % adimensional coefficient
r_1 = L2_1 / L1_1; % adimensional coefficient
a11_1 = miu_1; a22_1 = r_1; % main diagonal coefficients (constants)
miu_2 = 1 + m1_2 / m2_2; % adimensional coefficient
r_2 = L2_2 / L1_2; % adimensional coefficient
a11_2 = miu_2; a22_2 = r_2; % main diagonal coefficients (constants)
miu_3 = 1 + m1_3 / m2_3; % adimensional coefficient
r_3 = L2_3 / L1_3; % adimensional coefficient
a11_3 = miu_3; a22_3 = r_3; % main diagonal coefficients (constants)

tic; % Start recurrency time counting.

% Define functions for system of differential equations
f1_1 = @(t, theta1, theta2, OM1, OM2) OM1;
f2_1 = @(t, theta1, theta2, OM1, OM2) OM2;
f3_1 = @(t, theta1, theta2, OM1, OM2) (-g * (sin(theta1) + r_1 * sin(theta2)) - ...
    miu_1 * (OM1^2 * sin(theta1 - theta2) + r_1 * OM2^2 * sin(theta1 - theta2)) ...
    * cos(theta1 - theta2)) / (L1_1 * (miu_1 - cos(theta1 - theta2)^2));
f4_1 = @(t, theta1, theta2, OM1, OM2) (g * miu_1 * (sin(theta1) + r_1 * sin(theta2)) + ...
    miu_1 * (OM1^2 * sin(theta1 - theta2) + r_1 * OM2^2 * sin(theta1 - theta2)) ...
    * cos(theta1 - theta2)) / (L2_1 * (miu_1 - cos(theta1 - theta2)^2));

f1_2 = @(t, theta1, theta2, OM1, OM2) OM1;
f2_2 = @(t, theta1, theta2, OM1, OM2) OM2;
f3_2 = @(t, theta1, theta2, OM1, OM2) (-g * (sin(theta1) + r_2 * sin(theta2)) - ...
    miu_2 * (OM1^2 * sin(theta1 - theta2) + r_2 * OM2^2 * sin(theta1 - theta2)) ...
    * cos(theta1 - theta2)) / (L1_2 * (miu_2 - cos(theta1 - theta2)^2));
f4_2 = @(t, theta1, theta2, OM1, OM2) (g * miu_2 * (sin(theta1) + r_2 * sin(theta2)) + ...
    miu_2 * (OM1^2 * sin(theta1 - theta2) + r_2 * OM2^2 * sin(theta1 - theta2)) ...
    * cos(theta1 - theta2)) / (L2_2 * (miu_2 - cos(theta1 - theta2)^2));

f1_3 = @(t, theta1, theta2, OM1, OM2) OM1;
f2_3 = @(t, theta1, theta2, OM1, OM2) OM2;
f3_3 = @(t, theta1, theta2, OM1, OM2) (-g * (sin(theta1) + r_3 * sin(theta2)) - ...
    miu_3 * (OM1^2 * sin(theta1 - theta2) + r_3 * OM2^2 * sin(theta1 - theta2)) * ...
    cos(theta1 - theta2)) / (L1_3 * (miu_3 - cos(theta1 - theta2)^2));
f4_3 = @(t, theta1, theta2, OM1, OM2) (g * miu_3 * (sin(theta1) + r_3 * sin(theta2)) + ...
    miu_3 * (OM1^2 * sin(theta1 - theta2) + r_3 * OM2^2 * sin(theta1 - theta2)) * ...
    cos(theta1 - theta2)) / (L2_3 * (miu_3 - cos(theta1 - theta2)^2));

% Solve the system of ODE 4 with Runge-kutta 
for n = 2:N - 1
    % First pendulum
    % Computes the slope at (t(n), theta1(n), theta2(n), OM1(n), OM2(n)), scaled by the step size dt.
    K1_1 = dt * f1_1(t(n), theta1_1(n), theta2_1(n), OM1_1(n), OM2_1(n));
    L1_1 = dt * f2_1(t(n), theta1_1(n), theta2_1(n), OM1_1(n), OM2_1(n));
    M1_1 = dt * f3_1(t(n), theta1_1(n), theta2_1(n), OM1_1(n), OM2_1(n));
    N1_1 = dt * f4_1(t(n), theta1_1(n), theta2_1(n), OM1_1(n), OM2_1(n));
    
    % Computes slope at midpoint, using k1 adjusted state.
    K2_1 = dt * f1_1(t(n) + dt / 2, theta1_1(n) + K1_1 / 2, theta2_1(n) + L1_1 / 2, OM1_1(n) + M1_1 / 2, OM2_1(n) + N1_1 / 2);
    L2_1 = dt * f2_1(t(n) + dt / 2, theta1_1(n) + K1_1 / 2, theta2_1(n) + L1_1 / 2, OM1_1(n) + M1_1 / 2, OM2_1(n) + N1_1 / 2);
    M2_1 = dt * f3_1(t(n) + dt / 2, theta1_1(n) + K1_1 / 2, theta2_1(n) + L1_1 / 2, OM1_1(n) + M1_1 / 2, OM2_1(n) + N1_1 / 2);
    N2_1 = dt * f4_1(t(n) + dt / 2, theta1_1(n) + K1_1 / 2, theta2_1(n) + L1_1 / 2, OM1_1(n) + M1_1 / 2, OM2_1(n) + N1_1 / 2);
    
    % Computes slope at midpoint, using k2 adjusted state.
    K3_1 = dt * f1_1(t(n) + dt / 2, theta1_1(n) + K2_1 / 2, theta2_1(n) + L2_1 / 2, OM1_1(n) + M2_1 / 2, OM2_1(n) + N2_1 / 2);
    L3_1 = dt * f2_1(t(n) + dt / 2, theta1_1(n) + K2_1 / 2, theta2_1(n) + L2_1 / 2, OM1_1(n) + M2_1 / 2, OM2_1(n) + N2_1 / 2);
    M3_1 = dt * f3_1(t(n) + dt / 2, theta1_1(n) + K2_1 / 2, theta2_1(n) + L2_1 / 2, OM1_1(n) + M2_1 / 2, OM2_1(n) + N2_1 / 2);
    N3_1 = dt * f4_1(t(n) + dt / 2, theta1_1(n) + K2_1 / 2, theta2_1(n) + L2_1 / 2, OM1_1(n) + M2_1 / 2, OM2_1(n) + N2_1 / 2);
    
    %  Computes slope at the next step, using k3 adjusted state.
    K4_1 = dt * f1_1(t(n) + dt, theta1_1(n) + K3_1, theta2_1(n) + L3_1, OM1_1(n) + M3_1, OM2_1(n) + N3_1);
    L4_1 = dt * f2_1(t(n) + dt, theta1_1(n) + K3_1, theta2_1(n) + L3_1, OM1_1(n) + M3_1, OM2_1(n) + N3_1);
    M4_1 = dt * f3_1(t(n) + dt, theta1_1(n) + K3_1, theta2_1(n) + L3_1, OM1_1(n) + M3_1, OM2_1(n) + N3_1);
    N4_1 = dt * f4_1(t(n) + dt, theta1_1(n) + K3_1, theta2_1(n) + L3_1, OM1_1(n) + M3_1, OM2_1(n) + N3_1);
    
    % Updates state using weighted average of slopes.
    theta1_1(n + 1) = theta1_1(n) + (K1_1 + 2 * K2_1 + 2 * K3_1 + K4_1) / 6;
    theta2_1(n + 1) = theta2_1(n) + (L1_1 + 2 * L2_1 + 2 * L3_1 + L4_1) / 6;
    % Updates time to the next step.
    OM1_1(n + 1) = OM1_1(n) + (M1_1 + 2 * M2_1 + 2 * M3_1 + M4_1) / 6;
    OM2_1(n + 1) = OM2_1(n) + (N1_1 + 2 * N2_1 + 2 * N3_1 + N4_1) / 6;

    % Second pendulum - ompute the increments using the Runge-Kutta method
    K1_2 = dt * f1_2(t(n), theta1_2(n), theta2_2(n), OM1_2(n), OM2_2(n));
    L1_2 = dt * f2_2(t(n), theta1_2(n), theta2_2(n), OM1_2(n), OM2_2(n));
    M1_2 = dt * f3_2(t(n), theta1_2(n), theta2_2(n), OM1_2(n), OM2_2(n));
    N1_2 = dt * f4_2(t(n), theta1_2(n), theta2_2(n), OM1_2(n), OM2_2(n));
    
    K2_2 = dt * f1_2(t(n) + dt / 2, theta1_2(n) + K1_2 / 2, theta2_2(n) + L1_2 / 2, OM1_2(n) + M1_2 / 2, OM2_2(n) + N1_2 / 2);
    L2_2 = dt * f2_2(t(n) + dt / 2, theta1_2(n) + K1_2 / 2, theta2_2(n) + L1_2 / 2, OM1_2(n) + M1_2 / 2, OM2_2(n) + N1_2 / 2);
    M2_2 = dt * f3_2(t(n) + dt / 2, theta1_2(n) + K1_2 / 2, theta2_2(n) + L1_2 / 2, OM1_2(n) + M1_2 / 2, OM2_2(n) + N1_2 / 2);
    N2_2 = dt * f4_2(t(n) + dt / 2, theta1_2(n) + K1_2 / 2, theta2_2(n) + L1_2 / 2, OM1_2(n) + M1_2 / 2, OM2_2(n) + N1_2 / 2);
    
    K3_2 = dt * f1_2(t(n) + dt / 2, theta1_2(n) + K2_2 / 2, theta2_2(n) + L2_2 / 2, OM1_2(n) + M2_2 / 2, OM2_2(n) + N2_2 / 2);
    L3_2 = dt * f2_2(t(n) + dt / 2, theta1_2(n) + K2_2 / 2, theta2_2(n) + L2_2 / 2, OM1_2(n) + M2_2 / 2, OM2_2(n) + N2_2 / 2);
    M3_2 = dt * f3_2(t(n) + dt / 2, theta1_2(n) + K2_2 / 2, theta2_2(n) + L2_2 / 2, OM1_2(n) + M2_2 / 2, OM2_2(n) + N2_2 / 2);
    N3_2 = dt * f4_2(t(n) + dt / 2, theta1_2(n) + K2_2 / 2, theta2_2(n) + L2_2 / 2, OM1_2(n) + M2_2 / 2, OM2_2(n) + N2_2 / 2);
    
    K4_2 = dt * f1_2(t(n) + dt, theta1_2(n) + K3_2, theta2_2(n) + L3_2, OM1_2(n) + M3_2, OM2_2(n) + N3_2);
    L4_2 = dt * f2_2(t(n) + dt, theta1_2(n) + K3_2, theta2_2(n) + L3_2, OM1_2(n) + M3_2, OM2_2(n) + N3_2);
    M4_2 = dt * f3_2(t(n) + dt, theta1_2(n) + K3_2, theta2_2(n) + L3_2, OM1_2(n) + M3_2, OM2_2(n) + N3_2);
    N4_2 = dt * f4_2(t(n) + dt, theta1_2(n) + K3_2, theta2_2(n) + L3_2, OM1_2(n) + M3_2, OM2_2(n) + N3_2);
    
    theta1_2(n + 1) = theta1_2(n) + (K1_2 + 2 * K2_2 + 2 * K3_2 + K4_2) / 6;
    theta2_2(n + 1) = theta2_2(n) + (L1_2 + 2 * L2_2 + 2 * L3_2 + L4_2) / 6;
    OM1_2(n + 1) = OM1_2(n) + (M1_2 + 2 * M2_2 + 2 * M3_2 + M4_2) / 6;
    OM2_2(n + 1) = OM2_2(n) + (N1_2 + 2 * N2_2 + 2 * N3_2 + N4_2) / 6;

    % Third pendulum- Compute the increments using the Runge-Kutta method
    K1_3 = dt * f1_3(t(n), theta1_3(n), theta2_3(n), OM1_3(n), OM2_3(n));
    L1_3 = dt * f2_3(t(n), theta1_3(n), theta2_3(n), OM1_3(n), OM2_3(n));
    M1_3 = dt * f3_3(t(n), theta1_3(n), theta2_3(n), OM1_3(n), OM2_3(n));
    N1_3 = dt * f4_3(t(n), theta1_3(n), theta2_3(n), OM1_3(n), OM2_3(n));
    
    K2_3 = dt * f1_3(t(n) + dt / 2, theta1_3(n) + K1_3 / 2, theta2_3(n) + L1_3 / 2, OM1_3(n) + M1_3 / 2, OM2_3(n) + N1_3 / 2);
    L2_3 = dt * f2_3(t(n) + dt / 2, theta1_3(n) + K1_3 / 2, theta2_3(n) + L1_3 / 2, OM1_3(n) + M1_3 / 2, OM2_3(n) + N1_3 / 2);
    M2_3 = dt * f3_3(t(n) + dt / 2, theta1_3(n) + K1_3 / 2, theta2_3(n) + L1_3 / 2, OM1_3(n) + M1_3 / 2, OM2_3(n) + N1_3 / 2);
    N2_3 = dt * f4_3(t(n) + dt / 2, theta1_3(n) + K1_3 / 2, theta2_3(n) + L1_3 / 2, OM1_3(n) + M1_3 / 2, OM2_3(n) + N1_3 / 2);
    
    K3_3 = dt * f1_3(t(n) + dt / 2, theta1_3(n) + K2_3 / 2, theta2_3(n) + L2_3 / 2, OM1_3(n) + M2_3 / 2, OM2_3(n) + N2_3 / 2);
    L3_3 = dt * f2_3(t(n) + dt / 2, theta1_3(n) + K2_3 / 2, theta2_3(n) + L2_3 / 2, OM1_3(n) + M2_3 / 2, OM2_3(n) + N2_3 / 2);
    M3_3 = dt * f3_3(t(n) + dt / 2, theta1_3(n) + K2_3 / 2, theta2_3(n) + L2_3 / 2, OM1_3(n) + M2_3 / 2, OM2_3(n) + N2_3 / 2);
    N3_3 = dt * f4_3(t(n) + dt / 2, theta1_3(n) + K2_3 / 2, theta2_3(n) + L2_3 / 2, OM1_3(n) + M2_3 / 2, OM2_3(n) + N2_3 / 2);
    
    K4_3 = dt * f1_3(t(n) + dt, theta1_3(n) + K3_3, theta2_3(n) + L3_3, OM1_3(n) + M3_3, OM2_3(n) + N3_3);
    L4_3 = dt * f2_3(t(n) + dt, theta1_3(n) + K3_3, theta2_3(n) + L3_3, OM1_3(n) + M3_3, OM2_3(n) + N3_3);
    M4_3 = dt * f3_3(t(n) + dt, theta1_3(n) + K3_3, theta2_3(n) + L3_3, OM1_3(n) + M3_3, OM2_3(n) + N3_3);
    N4_3 = dt * f4_3(t(n) + dt, theta1_3(n) + K3_3, theta2_3(n) + L3_3, OM1_3(n) + M3_3, OM2_3(n) + N3_3);
    
    theta1_3(n + 1) = theta1_3(n) + (K1_3 + 2 * K2_3 + 2 * K3_3 + K4_3) / 6;
    theta2_3(n + 1) = theta2_3(n) + (L1_3 + 2 * L2_3 + 2 * L3_3 + L4_3) / 6;
    OM1_3(n + 1) = OM1_3(n) + (M1_3 + 2 * M2_3 + 2 * M3_3 + M4_3) / 6;
    OM2_3(n + 1) = OM2_3(n) + (N1_3 + 2 * N2_3 + 2 * N3_3 + N4_3) / 6;
end

%%                              Plot the results

figure;

subplot(2,1,1);
plot(t, theta1_1 * 180 / pi);
hold on;
plot(t, theta1_2 * 180 / pi);
plot(t, theta1_3 * 180 / pi);
hold off;
title('Theta_1 vs. Time');
xlabel('Time (s)');
legend('L1/L2 = 1', 'L1/L2 = 2', 'L1/L2 = 3');

subplot(2,1,2);
plot(t, theta2_1 * 180 / pi);
hold on;
plot(t, theta2_2 * 180 / pi);
plot(t, theta2_3 * 180 / pi);
hold off;
title('Theta_2 vs. Time');
xlabel('Time (s)');
legend('L1/L2 = 1', 'L1/L2 = 2', 'L1/L2 = 3');
