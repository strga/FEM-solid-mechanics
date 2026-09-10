clear all;
close all;
clc;

M = load_gmsh('/Users/strihavka/Documents/FEM_Diplom/FEM/Elasticity_Geo/AGARD_v2.msh');

mass = spconvert(load('/Volumes/MacOS/Orthotropic_Wing/Ortotropic_wing_0_mass_scaled.dat')); 
stiff = spconvert(load('/Volumes/MacOS/Orthotropic_Wing/Ortotropic_wing_0_stiff_scaled.dat'));

% M = load_gmsh('/Users/strihavka/Documents/FEM_Diplom/FEM/Elasticity_Geo/MESH_3D/BEAM_3D_CANTILEVER.msh');
% 
% mass = spconvert(load('/Users/strihavka/Documents/FEM_Diplom/FEM/Matrices_modal/Mass_matrix_TEST.dat')); 
% stiff = spconvert(load('/Users/strihavka/Documents/FEM_Diplom/FEM/Matrices_modal/Stiff_matrix_TEST.dat'));

%{
a = 10;
n = 15;

[V_scaled, D_scaled, flag] = eigs(stiff, mass, 2*n, 'smallestabs');

if flag ~= 0
   warning('Some eigenvalues did not converge.');
end

% --- Identify and zero out constrained DOFs ---
% fixedDOF = find(abs(diag(mass) - 1.0) < 1e-12);
% V_scaled(fixedDOF, :) = 0;  % enforce zero displacement at fixed DOFs

% --- Compute frequencies in Hz ---
freqs = sqrt(diag(D_scaled)) / (2*pi);
disp('Eigenfrequencies [Hz]:');
disp(freqs);

for i = a:n
    phi = reshape(V_scaled(:, i), 3, M.nbNod).';     % interleaved!
    % phi = phi ./ max(vecnorm(phi,2,2));             % normalize
    saveVTK3D(M, phi, sprintf('eigenmode_%d.vtk', i));
end
%}

n_all = 30;  
[V, D, flag] = eigs(stiff, mass, n_all, 'smallestabs');

if flag ~= 0
   warning('Some eigenvalues did not converge.');
end

% --- Identify and zero out constrained DOFs ---
fixedDOF = find(abs(diag(mass) - 1.0) < 1e-12);
V_scaled(fixedDOF, :) = 0;  % enforce zero displacement at fixed DOFs

% Convert eigenvalues to frequencies
freqs = sqrt(diag(D)) / (2*pi);

% Sort frequencies and corresponding eigenvectors
[freqs, sortIdx] = sort(freqs);
V = V(:, sortIdx);

% Find unique frequencies within a tolerance
tol = 1e-3; % tolerance for uniqueness (Hz)
uniqueIdx = [1; find(diff(freqs) > tol) + 1];
uniqueFreqs = freqs(uniqueIdx);

% Pick first n unique modes
n_unique = min(10, numel(uniqueIdx));  % e.g., first 10 unique frequencies
fprintf('Unique eigenfrequencies [Hz]:\n');
disp(uniqueFreqs(1:n_unique));

% Save the unique eigenmodes
for k = 1:n_unique
    idx = uniqueIdx(k);
    phi = reshape(V(:, idx), 3, M.nbNod).';  % reshape to n x 3
    % phi = phi / max(vecnorm(phi, 2, 2));     % normalize displacement
    saveVTK3D(M, phi, sprintf('eigenmode_unique_%d.vtk', k));
end
