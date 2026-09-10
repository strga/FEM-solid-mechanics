clear;
clc;
format longG;

% Perform modal analysis on unscaled matrices
disp('Performing modal analysis on unscaled matrices...');
% K = spconvert(load('/Users/strihavka/Documents/FEM_Diplom/FEM/Matrices_modal/NUnscaled_Dyn_Stiff_damp.dat'));
% Mass = spconvert(load('/Users/strihavka/Documents/FEM_Diplom/FEM/Matrices_modal/NUnscaled_Dyn_Mass_damp.dat'));
% Dumping = spconvert(load('/Users/strihavka/Documents/FEM_Diplom/FEM/Matrices_modal/NUnscaled_Dyn_Dump_damp.dat'));
%-------------------------------------------------------------------------------------------------------------------%
K = spconvert(load('/Users/strihavka/Documents/FEM_Diplom/FEM/Matrices_modal/Stiff_matrix_damping_100K.dat'));
Mass = spconvert(load('/Users/strihavka/Documents/FEM_Diplom/FEM/Matrices_modal/Mass_matrix_damping_100K.dat'));
% Dumping = spconvert(load('/Users/strihavka/Documents/FEM_Diplom/FEM/Matrices_modal/Damping_matrix.dat'));

%------------------NEW_WING------------
% K = spconvert(load('/Volumes/MacOS/Orthotropic_Wing/Ortotropic_wing_-45_stiff_scaled.dat'));
% Mass = spconvert(load('/Volumes/MacOS/Orthotropic_Wing/Ortotropic_wing_-45_mass_scaled.dat'));

% Number of smallest eigenvalues to compute
% n = 2 * round(size(Mass, 1) / 8);

n = 90;

% Compute eigenvalues and eigenvectors
[V, D, flag] = eigs(K, Mass, n, 'smallestabs');

if flag ~= 0
   warning('Some eigenvalues did not converge.');
end

% Extract eigenvalues and angular frequencies
lambdas = diag(D); % Eigenvalues
omega = sqrt(lambdas); % Natural angular frequencies
frequencies = omega / (2 * pi); % Frequencies in Hz

% Filter out Dirichlet-associated eigenvalues
physical_indices = lambdas > 1 + 1e-6;
lambdas_filtered = lambdas(physical_indices);
omega_filtered = sqrt(lambdas_filtered);
smallest_omega = min(omega_filtered);
disp(smallest_omega);
% frequencies_filtered = omega_filtered / (2 * pi);
% V_filtered = V(:, physical_indices);

% Display filtered results
% disp('Filtered natural angular frequencies (rad/s):');
% disp(omega_filtered);

% frequecnies_filtered = omega_filtered / ( 2 * pi );
% disp('Filtered natural frequencies (Hz):');
% disp(frequencies_filtered);
disp('THE END');

%-------------------------------------------------------------------------%
%{

% Perform modal analysis on scaled matrices
disp('Performing modal analysis on scaled matrices...');
ScK = spconvert(load('/Users/strihavka/Documents/FEM_Diplom/FEM/Matrices_modal/NEW_Scaled_Dyn_Stiff.dat'));
ScMass = spconvert(load('/Users/strihavka/Documents/FEM_Diplom/FEM/Matrices_modal/NEW_Scaled_Dyn_Mass.dat'));
ScDumping = spconvert(load('/Users/strihavka/Documents/FEM_Diplom/FEM/Matrices_modal/NEW_Scaled_Dyn_Dump.dat'));

% Compute eigenvalues and eigenvectors
[V_scaled, D_scaled, flag] = eigs(ScK, ScMass, n, 'smallestabs');

if flag ~= 0
   warning('Some eigenvalues did not converge.');
end

% Extract eigenvalues and angular frequencies
lambdas_scaled = diag(D_scaled); % Eigenvalues
omega_scaled = sqrt(lambdas_scaled); % Natural angular frequencies
frequencies_scaled = omega_scaled / (2 * pi); % Frequencies in Hz

% Filter out Dirichlet-associated eigenvalues
physical_indices_scaled = lambdas_scaled > 1 + 1e-6;
lambdas_filtered_scaled = lambdas_scaled(physical_indices_scaled);
omega_filtered_scaled = sqrt(lambdas_filtered_scaled);
frequencies_filtered_scaled = omega_filtered_scaled / (2 * pi);
frequencies_filtered_unique = unique(frequencies_filtered_scaled);
V_filtered_scaled = V_scaled(:, physical_indices_scaled);

% Display filtered results
% disp('Filtered natural angular frequencies (rad/s):');
% disp(omega_filtered);

% frequecnies_filtered = omega_filtered / ( 2 * pi );
% disp('Filtered natural frequencies (Hz):');
% disp(frequencies_filtered);


%-------------------------------------------------------------------------%
% Export mode shapes to VTK files
% disp('Exporting mode shapes...');
% output_dir = 'Modal_vtk';
% if ~exist(output_dir, 'dir')
%     mkdir(output_dir);
% end
% 
% for i = 1:5%size(V_filtered_scaled, 2)
%     % Extract x and y displacements from eigenvectors
%     u_x = V_filtered_scaled(1:2:end, i); % X displacements
%     u_y = V_filtered_scaled(2:2:end, i); % Y displacements
% 
%     % Combine displacements into a single vector for VTK export
%     displacement = zeros(M.nbNod, 3); % Initialize [u_x, u_y, u_z=0]
%     displacement(:, 1) = u_x; % Assign X displacements
%     displacement(:, 2) = u_y; % Assign Y displacements
% 
%     % Save mode shape to VTK file
%     filename = fullfile(output_dir, sprintf('Test_Linear_Mode_Scaled_%d.vtk', i));
%     saveVTK(M, displacement, filename);
% end
% 
% disp('Mode shapes exported successfully.');

%}
