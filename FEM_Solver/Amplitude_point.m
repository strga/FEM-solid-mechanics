clear all;
clc;
close all;

% Load the displacement data from file
% data = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Vibration_data/Bending_5/New_Linear_point_displacement_5.dat'); %...1-3
data = load('/Users/strihavka/Documents/FEM_Diplom/FEM/Vibration_data/Linear_point_displacement_3D_6.dat'); %...1-3
time = data(:, 1);
amplitude = data(:, 2);

% Plot the amplitude over time
plot(time, amplitude);
xlabel('Time');
ylabel('Displacement Amplitude');
title('Displacement Amplitude at Point Over Time');
xlim([0, 1.0]);
grid on;