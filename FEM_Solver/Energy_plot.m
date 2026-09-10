clear all;
clc;
close all;

% Load energy data
% data = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Vibration_data/Bending/New_Linear_energy_data.txt'); %...1A-3
data = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Vibration_data/Linear_energy_3D_cl021_11021e.txt');
time = data(:, 1);
kinetic_energy = data(:, 2);
potential_energy = data(:, 3);
total_energy = data(:, 4);

% Plot kinetic, potential, and total energy
figure(1);
plot(time, kinetic_energy, 'r', 'LineWidth', 0.1, 'DisplayName', 'Kinetic Energy');
hold on;
plot(time, potential_energy, 'b', 'LineWidth', 0.1, 'DisplayName', 'Potential Energy');
plot(time, total_energy, 'k--', 'LineWidth', 2, 'DisplayName', 'Total Energy');
xlabel('Time');
ylabel('Energy');
legend show;
title('Energy Conservation Over Time');
hold off;
xlim([0, 2.5]);
grid on;

%--------------------------------------------------------------------------
%{
% Load energy data1
data1 = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Vibration_data/Bending_1/New_Linear_energy_data_1.txt'); %...1A-3
time1 = data1(:, 1);
kinetic_energy1 = data1(:, 2);
potential_energy1 = data1(:, 3);
total_energy1 = data1(:, 4);

% Plot kinetic1, potential1, and total energy1
figure(2);
plot(time1, kinetic_energy1, 'r', 'LineWidth', 0.1, 'DisplayName', 'Kinetic Energy');
hold on;
plot(time1, potential_energy1, 'b', 'LineWidth', 0.1, 'DisplayName', 'Potential Energy');
plot(time1, total_energy1, 'k--', 'LineWidth', 2, 'DisplayName', 'Total Energy');
xlabel('Time');
ylabel('Energy');
legend show;
title('Energy Conservation Over Time');
hold off;
xlim([0, 2.5]);
grid on;

%--------------------------------------------------------------------------

% Load energy data2
data2 = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Vibration_data/Bending_2/New_Linear_energy_data_2.txt'); %...1A-3
time2 = data2(:, 1);
kinetic_energy2 = data2(:, 2);
potential_energy2 = data2(:, 3);
total_energy2 = data2(:, 4);

% Plot kinetic1, potential1, and total energy2
figure(3);
plot(time2, kinetic_energy2, 'r', 'LineWidth', 0.1, 'DisplayName', 'Kinetic Energy');
hold on;
plot(time2, potential_energy2, 'b', 'LineWidth', 0.1, 'DisplayName', 'Potential Energy');
plot(time2, total_energy2, 'k--', 'LineWidth', 2, 'DisplayName', 'Total Energy');
xlabel('Time');
ylabel('Energy');
legend show;
title('Energy Conservation Over Time');
hold off;
xlim([0, 2.5]);
grid on;

%--------------------------------------------------------------------------

% Load energy data3
data3 = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Vibration_data/Bending_3/New_Linear_energy_data_3.txt'); %...1A-3
time3 = data3(:, 1);
kinetic_energy3 = data3(:, 2);
potential_energy3 = data3(:, 3);
total_energy3 = data3(:, 4);

% Plot kinetic1, potential1, and total energy2
figure(4);
plot(time3, kinetic_energy3, 'r', 'LineWidth', 0.1, 'DisplayName', 'Kinetic Energy');
hold on;
plot(time3, potential_energy3, 'b', 'LineWidth', 0.1, 'DisplayName', 'Potential Energy');
plot(time3, total_energy3, 'k--', 'LineWidth', 2, 'DisplayName', 'Total Energy');
xlabel('Time');
ylabel('Energy');
legend show;
title('Energy Conservation Over Time');
hold off;
xlim([0, 2.5]);
grid on;

%--------------------------------------------------------------------------

% Load energy data4
data4 = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Vibration_data/Bending_4/New_Linear_energy_data_4.txt'); %...1A-3
time4 = data4(:, 1);
kinetic_energy4 = data4(:, 2);
potential_energy4 = data4(:, 3);
total_energy4 = data4(:, 4);

% Plot kinetic1, potential1, and total energy2
figure(5);
plot(time4, kinetic_energy4, 'r', 'LineWidth', 0.1, 'DisplayName', 'Kinetic Energy');
hold on;
plot(time4, potential_energy4, 'b', 'LineWidth', 0.1, 'DisplayName', 'Potential Energy');
plot(time4, total_energy4, 'k--', 'LineWidth', 2, 'DisplayName', 'Total Energy');
xlabel('Time');
ylabel('Energy');
legend show;
title('Energy Conservation Over Time');
hold off;
xlim([0, 2.5]);
grid on;

%--------------------------------------------------------------------------

% Load energy data5
data5 = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Vibration_data/Bending_5/New_Linear_energy_data_5.txt'); %...1A-3
time5 = data5(:, 1);
kinetic_energy5 = data5(:, 2);
potential_energy5 = data5(:, 3);
total_energy5 = data5(:, 4);

% Plot kinetic1, potential1, and total energy2
figure(6);
plot(time5, kinetic_energy5, 'r', 'LineWidth', 0.1, 'DisplayName', 'Kinetic Energy');
hold on;
plot(time5, potential_energy5, 'b', 'LineWidth', 0.1, 'DisplayName', 'Potential Energy');
plot(time5, total_energy5, 'k--', 'LineWidth', 2, 'DisplayName', 'Total Energy');
xlabel('Time');
ylabel('Energy');
legend show;
title('Energy Conservation Over Time');
hold off;
xlim([0, 2.5]);
grid on;
%}
