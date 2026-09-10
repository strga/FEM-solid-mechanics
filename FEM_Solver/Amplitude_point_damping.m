clear all
close all
clc

eps1 = 0.00005;
eps2 = 0.0006;
omega = 55.044002776095430;

d = -( eps1 + eps2 * omega * omega ) / 2;

omega_damped = sqrt( omega * omega - (eps1 + eps2 * omega * omega)^2 / 4 );

disp('Koeficient utlumu d:');
disp(d);

disp('Vlastni frekvence tlumeni:');
disp(omega_damped);

% time, x_displacement, y_displacement
% data = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Vibration_data/Bending_damping/New_Linear_point_displacement_damp.dat');
% data = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Vibration_data/New_Linear_point_displacement_damp_testing68.dat');

data = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Beam_diplom4_displacement_damping_100K_small.dat');

% Extract each column into separate variables
time = data(:,1);
u_x  = data(:,2);
u_y  = data(:,3);
u_z  = data(:,4);

A1 = max(abs(u_x));
A2 = max(abs(u_y));
A3 = max(abs(u_z));

y = exp(d * time);

plot(time, y)

num_plots = 3;

% Create a figure
figure;

% Plot the x-displacement over time
subplot(num_plots,1,1); % top subplot
plot(time, u_x, 'r-', 'LineWidth', 1.0);
hold on;
plot(time, A1 * y, 'k-', 'LineWidth', 2.0);
xlabel('Time (s)');
ylabel('X-Displacement (m)');
title('X-Displacement of Node Over Time');
hold off;
grid on;
xlim([0, 4.8]);

% Plot the y-displacement over time
subplot(num_plots,1,2); % bottom subplot
plot(time, u_y, 'b-', 'LineWidth', 1.0);
hold on;
plot(time, A2 * y, 'k-', 'LineWidth', 2.0);
xlabel('Time (s)');
ylabel('Y-Displacement (m)');
title('Y-Displacement of Node Over Time');
hold off;
grid on;
xlim([0, 4.8]);

% Plot the z-displacement over time
subplot(num_plots,1,3); % bottom subplot
plot(time, u_z, 'g-', 'LineWidth', 1.0);
hold on;
plot(time, A3 * y, 'k-', 'LineWidth', 2.0);
xlabel('Time (s)');
ylabel('Z-Displacement (m)');
title('Z-Displacement of Node Over Time');
hold off;
grid on;
xlim([0, 4.8]);

% === Separate figure for Y-displacement ===
figure;
plot(time, u_y, 'b-', 'LineWidth', 1.0); hold on;
plot(time, A2 * y, 'k-', 'LineWidth', 2.0);
xlabel('t [ s ]');
ylabel(' u [ m ]');
%title('Y-Displacement of Node Over Time');
grid on;
% xlim([0, 4.8]);
ylim([-0.05, 0.05]);
legend('u(t)', 'y(t)');
