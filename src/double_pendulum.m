clear; close all; clc;
g = 9.80665; % m / s^2;

% 1ST PENDULUM
L1_1 = 2.2; L2_1 = 1.3; % m;
m1_1 = 1.2; m2_1 = 1.7; % kg;
% 2ND PENDULUM
L1_2 = 1.8; L2_2 = 1.0; % m;
m1_2 = 1.0; m2_2 = 1.5; % kg;
% 3RD PENDULUM
L1_3 = 2.5; L2_3 = 1.8; % m;
m1_3 = 1.5; m2_3 = 2.0; % kg;

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

%%                   DURATION: frequency and period components 

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

for i = 2:N - 1 % recurrencies cycle
    aux_1 = theta2_1(i) - theta1_1(i);
    a21_1 = cos(aux_1); a12_1 = a21_1 * r_1; % secondary diagonal

    OM1_1(i) = (theta1_1(i) - theta1_1(i - 1)) / dt; % angular velocity 1st body
    OM2_1(i) = (theta2_1(i) - theta2_1(i - 1)) / dt; % angular velocity 2nd body

    % Free terms
    b1_1 = r_1 * OM2_1(i)^2 * sin(aux_1) - g / L1_1 * miu_1 * sin(theta1_1(i));
    b2_1 =  - OM1_1(i)^2 * sin(aux_1) - g / L1_1 * sin(theta2_1(i));

    % System matrix A and free terms column B
    A_1 = [a11_1,a12_1;a21_1,a22_1];
    B_1 = [b1_1;b2_1];
    E_1 = A_1\B_1;
    eps1_1 = E_1(1); eps2_1 = E_1(2); % angular accelerations

    % Recurrence of the II order:
    theta1_1(i + 1) = 2 * theta1_1(i) - theta1_1(i - 1) + dt^2 * eps1_1;
    theta2_1(i + 1) = 2 * theta2_1(i) - theta2_1(i - 1) + dt^2 * eps2_1;
end

for i = 2:N - 1 % recurrencies cycle
    aux_2 = theta2_2(i) - theta1_2(i);
    a21_2 = cos(aux_2); a12_2 = a21_2 * r_2; % secondary diagonal

    OM1_2(i) = (theta1_2(i) - theta1_2(i - 1)) / dt; % angular velocity 1st body
    OM2_2(i) = (theta2_2(i) - theta2_2(i - 1)) / dt; % angular velocity 2nd body
    
    % Free terms
    b1_2 = r_2 * OM2_2(i)^2 * sin(aux_2) - g / L1_2 * miu_2 * sin(theta1_2(i));
    b2_2 =  - OM1_2(i)^2 * sin(aux_2) - g / L1_2 * sin(theta2_2(i));

    % System matrix A and free terms column B
    A_2 = [a11_2,a12_2;a21_2,a22_2];
    B_2 = [b1_2;b2_2];
    E_2 = A_2\B_2;
    eps1_2 = E_2(1); eps2_2 = E_2(2); % angular accelerations

    % Recurrence of the II order:
    theta1_2(i + 1) = 2 * theta1_2(i) - theta1_2(i - 1) + dt^2 * eps1_2;
    theta2_2(i + 1) = 2 * theta2_2(i) - theta2_2(i - 1) + dt^2 * eps2_2;
end

for i = 2:N - 1 % recurrencies cycle
    aux_3 = theta2_3(i) - theta1_3(i);
    a21_3 = cos(aux_3); a12_3 = a21_3 * r_3;  % secondary diagonal

    OM1_3(i) = (theta1_3(i) - theta1_3(i - 1)) / dt; % angular velocity 1st body
    OM2_3(i) = (theta2_3(i) - theta2_3(i - 1)) / dt; % angular velocity 2nd body

    % Free terms
    b1_3 = r_3 * OM2_3(i)^2 * sin(aux_3) - g / L1_3 * miu_3 * sin(theta1_3(i));
    b2_3 =  - OM1_3(i)^2 * sin(aux_3) - g / L1_3 * sin(theta2_3(i));

    % System matrix A and free terms column B
    A_3 = [a11_3,a12_3;a21_3,a22_3];
    B_3 = [b1_3;b2_3];
    E_3 = A_3\B_3;
    eps1_3 = E_3(1); eps2_3 = E_3(2); % angular accelerations

    % Recurrence of the II order:
    theta1_3(i + 1) = 2 * theta1_3(i) - theta1_3(i - 1) + dt^2 * eps1_3;
    theta2_3(i + 1) = 2 * theta2_3(i) - theta2_3(i - 1) + dt^2 * eps2_3;
end

% ANGULAR VELOCITY:
OM1_1(N) = (theta1_1(N) - theta1_1(N - 1)) / dt;
OM2_1(N) = (theta2_1(N) - theta2_1(N - 1)) / dt;
OM1_2(N) = (theta1_2(N) - theta1_2(N - 1)) / dt;
OM2_2(N) = (theta2_2(N) - theta2_2(N - 1)) / dt;
OM1_3(N) = (theta1_3(N) - theta1_3(N - 1)) / dt;
OM2_3(N) = (theta2_3(N) - theta2_3(N - 1)) / dt;

toc; % Prints elapsed time, for the numerical solution.

% COORDONATES BODIES PENDULUMS:
x1_1 = L1_1 * sin(theta1_1); x2_1 = x1_1 + L2_1 * sin(theta2_1);
y1_1 =  - L1_1 * cos(theta1_1); y2_1 = y1_1 - L2_1 * cos(theta2_1);
x1_2 = L1_2 * sin(theta1_2); x2_2 = x1_2 + L2_2 * sin(theta2_2);
y1_2 =  - L1_2 * cos(theta1_2); y2_2 = y1_2 - L2_2 * cos(theta2_2);
x1_3 = L1_3 * sin(theta1_3); x2_3 = x1_3 + L2_3 * sin(theta2_3);
y1_3 =  - L1_3 * cos(theta1_3); y2_3 = y1_3 - L2_3 * cos(theta2_3);

% (T) Kinetic + (U) Potential = (H) Total Energy
% Hamiltonian of the systemtotal energy

T_1 = 1 / 2 * (m1_1 * L1_1^2 * OM1_1.^2 + ...
    m2_1 * (L1_1^2 * OM1_1.^2 + L2_1^2 * OM2_1.^2 + ...
    2 * L1_1 * L2_1 * OM1_1.* OM2_1.* cos(theta2_1 - theta1_1)));

U_1 =  - g * ((m1_1 + m2_1) * L1_1 * cos(theta1_1) + m2_1 * L2_1 * cos(theta2_1));
H_1 = T_1 + U_1;

T_2 = 1 / 2 * (m1_2 * L1_2^2 * OM1_2.^2 + ...
    m2_2 * (L1_2^2 * OM1_2.^2 + L2_2^2 * OM2_2.^2 + ...
    2 * L1_2 * L2_2 * OM1_2.* OM2_2.* cos(theta2_2 - theta1_2)));

U_2 =  - g * ((m1_2 + m2_2) * L1_2 * cos(theta1_2) + m2_2 * L2_2 * cos(theta2_2));
H_2 = T_2 + U_2;

T_3 = 1 / 2 * (m1_3 * L1_3^2 * OM1_3.^2 + ...
    m2_3 * (L1_3^2 * OM1_3.^2 + L2_3^2 * OM2_3.^2 + ...
    2 * L1_3 * L2_3 * OM1_3.* OM2_3.* cos(theta2_3 - theta1_3)));

U_3 =  - g * ((m1_3 + m2_3) * L1_3 * cos(theta1_3) + m2_3 * L2_3 * cos(theta2_3));
H_3 = T_3 + U_3;

% Length double - pendulum
Lmax_1 = L1_1 + L2_1;
Lmax_2 = L1_2 + L2_2; 
Lmax_3 = L1_3 + L2_3;

coef = 30; % Size coefficient for the objects

% Radius for the first pendulum
rg1_1 = coef * m1_1^(1 / 3); rg2_1 = coef * m2_1^(1 / 3);
rg1_2 = coef * m1_2^(1 / 3); rg2_2 = coef * m2_2^(1 / 3);
rg1_3 = coef * m1_3^(1 / 3); rg2_3 = coef * m2_3^(1 / 3);

% FOR EACH SYSTEM, ADJUST THE GRAPHICAL RADIUS OF EACH BODY
rg11 = coef * m1_1^(1 / 3);
rg21 = coef * m2_1^(1 / 3);
rg12 = coef * m1_2^(1 / 3);
rg22 = coef * m2_2^(1 / 3);
rg13 = coef * m1_3^(1 / 3);
rg23 = coef * m2_3^(1 / 3);

% Define the maximum number of iterations
max_iterations = 1000;

% Initialize arrays to store pendulum positions during the simulation
x1_1_trace = zeros(1, max_iterations);
y1_1_trace = zeros(1, max_iterations);
x1_2_trace = zeros(1, max_iterations);
y1_2_trace = zeros(1, max_iterations);
x1_3_trace = zeros(1, max_iterations);
y1_3_trace = zeros(1, max_iterations);
% Initialize arrays to store bottom objects' positions during the simulation
xb_1_trace = zeros(1, max_iterations);
yb_1_trace = zeros(1, max_iterations);
xb_2_trace = zeros(1, max_iterations);
yb_2_trace = zeros(1, max_iterations);
xb_3_trace = zeros(1, max_iterations);
yb_3_trace = zeros(1, max_iterations);

% Graphic simulation 3 double pendulums

trace_count = 1;
trace_count_bottom = 1;

figure('Name', '3 Systems Double Pendulum');
set(gcf, 'Position', [100, 100, 800, 800]);

tic; simt = 0; % Start the timer

while simt <= tf % Graphic simulation  
    j = abs(t - simt) == min(abs(t - simt)); % Find the index of the current time
    
    plot([0, x1_1(j), x2_1(j)], [0, y1_1(j), y2_1(j)], '-r', 'LineWidth', 3); hold on; % links for S1
    plot([0, x1_2(j), x2_2(j)], [0, y1_2(j), y2_2(j)], '-g', 'LineWidth', 3); hold on; % links for S2
    plot([0, x1_3(j), x2_3(j)], [0, y1_3(j), y2_3(j)], '-b', 'LineWidth', 3); hold on; % links for S3

    % Store positions of the pendulums for tracebacks
    x1_1_trace(trace_count) = x1_1(j);
    y1_1_trace(trace_count) = y1_1(j);
    x1_2_trace(trace_count) = x1_2(j);
    y1_2_trace(trace_count) = y1_2(j);
    x1_3_trace(trace_count) = x1_3(j);
    y1_3_trace(trace_count) = y1_3(j);
    % Store positions of the bottom objects for tracebacks
    xb_1_trace(trace_count_bottom) = x2_1(j);
    yb_1_trace(trace_count_bottom) = y2_1(j);
    xb_2_trace(trace_count_bottom) = x2_2(j);
    yb_2_trace(trace_count_bottom) = y2_2(j);
    xb_3_trace(trace_count_bottom) = x2_3(j);
    yb_3_trace(trace_count_bottom) = y2_3(j);

    % Plot tracebacks for each pendulum with lifted center
    plot(x1_1_trace, y1_1_trace, '-r', 'LineWidth', 1, 'Color', [1, 0, 0, 0.3]);
    plot(x1_2_trace, y1_2_trace, '-g', 'LineWidth', 1, 'Color', [0, 1, 0, 0.3]);
    plot(x1_3_trace, y1_3_trace, '-b', 'LineWidth', 1, 'Color', [0, 0, 1, 0.3]);
    % Plot tracebacks for each bottom object with lifted center
    plot(xb_1_trace, yb_1_trace, '--r', 'LineWidth', 1, 'Color', [1, 0, 0, 0.3]);
    plot(xb_2_trace, yb_2_trace, '--g', 'LineWidth', 1, 'Color', [0, 1, 0, 0.3]);
    plot(xb_3_trace, yb_3_trace, '--b', 'LineWidth', 1, 'Color', [0, 0, 1, 0.3]);
    
    % Draws the pendulum with lifted center
    plot([0 x1_1(j) x2_1(j)], [0 y1_1(j) y2_1(j)], ' -r', 'LineWidth', 1); hold on; % Links for S1
    plot(x1_1(j), y1_1(j), '.k', 'MarkerSize', rg1_1); % Body 1
    plot(x1_1(j), y1_1(j), '.r', 'MarkerSize', rg1_1 / 3); % Body 1
    plot(x2_1(j), y2_1(j), '.r', 'MarkerSize', rg2_1); % Body 2
    plot(x2_1(j), y2_1(j), '.w', 'MarkerSize', rg2_1  / 3); % Body 2
    % Draws the pendulum with lifted center
    plot([0 x1_2(j) x2_2(j)], [0 y1_2(j) y2_2(j)], ' -g', 'LineWidth', 1); % Links for S2
    plot(x1_2(j), y1_2(j), '.k', 'MarkerSize', rg1_2); % Body 1
    plot(x1_2(j), y1_2(j), '.g', 'MarkerSize', rg1_2  / 3); % Body 1
    plot(x2_2(j), y2_2(j), '.g', 'MarkerSize', rg2_2); % Body 2
    plot(x2_2(j), y2_2(j), '.w', 'MarkerSize', rg2_2  / 3); % Body 2
    % Draws the pendulum with lifted center
    plot([0 x1_3(j) x2_3(j)], [0 y1_3(j) y2_3(j)], ' -b', 'LineWidth', 1); % Links for S3
    plot(x1_3(j), y1_3(j), '.k', 'MarkerSize', rg1_3); % Body 1
    plot(x1_3(j), y1_3(j), '.b', 'MarkerSize', rg1_3 / 3); % Body 1
    plot(x2_3(j), y2_3(j), '.b', 'MarkerSize', rg2_3); % Body 2
    plot(x2_3(j), y2_3(j), '.w', 'MarkerSize', rg2_3 / 3); % Body 2

    xlabel('x / m'); ylabel('y / m');
    axis equal; % aspect ratio 1:1 for the axes
    axis([ - max([Lmax_1, Lmax_2, Lmax_3]), max([Lmax_1, Lmax_2, Lmax_3]), ...
           - max([Lmax_1, Lmax_2, Lmax_3]), max([Lmax_1, Lmax_2, Lmax_3])]); % axes limits

    annotation('textbox', [0.45, 0.5, 0.25, 0.15], 'String',...
        ['E1 = ', num2str(( - 1) * H_1(j)), ' J'], 'BackgroundColor', 'white', 'FitBoxToText', 'on');
    annotation('textbox', [0.45, 0.55, 0.25, 0.15], 'String',...
        ['E2 = ', num2str(( - 1) * H_2(j)), ' J'], 'BackgroundColor', 'white', 'FitBoxToText', 'on');
    annotation('textbox', [0.45, 0.6, 0.25, 0.15], 'String',...
        ['E3 = ', num2str(( - 1) * H_3(j)), ' J'], 'BackgroundColor', 'white', 'FitBoxToText', 'on');
    annotation('textbox', [0.45, 0.65, 0.25, 0.15], 'String',...
        ['t = ', num2str(round(t(j) * 100)), ' cs'], 'BackgroundColor', 'white', 'FitBoxToText', 'on');

    trace_count_bottom = trace_count_bottom + 1;
    trace_count = trace_count + 1;
    simt = simt + 1e-1; % cs; step for the simulation
    pause(1e-3);
    hold off;
end
