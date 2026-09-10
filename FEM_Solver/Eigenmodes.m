clear all;
close all;
clc;

M = load_gmsh('AGARD_v2.msh');

mass = spconvert(load('Ortotropic_wing_0_mass_scaled.dat')); 
stiff = spconvert(load('Ortotropic_wing_0_stiff_scaled.dat'));

n_all = 30;  
[V, D, flag] = eigs(stiff, mass, n_all, 'smallestabs');

if flag ~= 0
   warning('Some eigenvalues did not converge');
end

% --- Identify and zero out constrained DOFs ---
fixedDOF = find(abs(diag(mass) - 1.0) < 1e-12);
V_scaled(fixedDOF, :) = 0;  % enforce zero displacement at fixed DOFs

% Convert eigenvalues to frequencies
freqs = sqrt(diag(D)) / (2*pi);

% Sort frequencies and corresponding eigenvectors
[freqs, sortIdx] = sort(freqs);
V = V(:, sortIdx);

% Find unique frequencies within tolerance
tol = 1e-3; % tolerance for uniqueness (Hz)
uniqueIdx = [1; find(diff(freqs) > tol) + 1];
uniqueFreqs = freqs(uniqueIdx);

% Pick first n unique modes
n_unique = min(10, numel(uniqueIdx));
fprintf('Unique eigenfrequencies [Hz]:\n');
disp(uniqueFreqs(1:n_unique));

% Save the unique eigenmodes
for k = 1:n_unique
    idx = uniqueIdx(k);
    phi = reshape(V(:, idx), 3, M.nbNod).';  % reshape to n x 3
    % phi = phi / max(vecnorm(phi, 2, 2));     % normalize displacement
    saveVTK3D(M, phi, sprintf('eigenmode_unique_%d.vtk', k));
end
