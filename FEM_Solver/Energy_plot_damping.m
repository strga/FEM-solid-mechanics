clear all;
clc;
close all;

% Load dumping energy data
data_dump = load('Beam_energy_conserved_damping_100K_small.txt');
time_dump = data_dump(:, 1);
kinetic_energy = data_dump(:,2);
potential_energy = data_dump(:,3);
initial_energy = data_dump(:,4);
work_ref = data_dump(:,5);
total_energy_damp = data_dump(:, 6);

computed_total_energy = kinetic_energy + potential_energy;

% Plot total energy
figure(1);
% plot(time_dump, total_energy_dump, 'b-', 'DisplayName', 'Total Energy','LineWidth', 1.5);
plot(time_dump, total_energy_damp, 'b-', 'DisplayName', 'Celková energie','LineWidth', 2.0);
hold on;
plot(time_dump, kinetic_energy, 'r-', 'DisplayName', 'Kinetická energie', 'LineWidth', 1.3 );
plot(time_dump, potential_energy, 'k-', 'DisplayName', 'Potenciální energie', 'LineWidth', 1.3 );
xlabel('t [ s ]');
ylabel('E [ J ]');
legend show;
grid on;
%title('Energy Conservation Over Time');
hold off;
% xlim([0, 4.8]);
% Plot total energy - logarithmic
%{
figure(2);
loglog(time_dump, total_energy_dump, 'b-', 'DisplayName', 'Total Energy','LineWidth', 1.5);
xlabel('Time');
ylabel('Energy');
legend show;
title('Energy Conservation Over Time');
hold off;
xlim([0, 4.8]);
%}
