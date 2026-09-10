clear all;
clc;
close all;

% Load the displacement data from file
data = load('Linear_displacement_3D.dat');
time = data(:, 1);
amplitude = data(:, 2);

% Plot the amplitude over time
plot(time, amplitude);
xlabel('Time');
ylabel('Displacement Amplitude');
title('Displacement Amplitude at Point Over Time');
xlim([0, 1.0]);
grid on;
