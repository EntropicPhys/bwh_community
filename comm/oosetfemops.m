function p = oosetfemops(p)
% oosetfemops — Sets problem-specific FEM matrices for the bwh model
%
% This function constructs the finite element mass and stiffness matrices
% needed for the spatial discretization of the system. It is placed in the 
% problem directory because the structure of the matrices is tightly linked 
% to the problem setup (e.g., number of components, coupling structure).

% === FEM Matrix Assembly for 1D Domain ===
[K, M, ~] = p.pdeo.fem.assema(p.pdeo.grid, 1, 1, 1);
% K = stiffness matrix (from ∇·(a∇u), here a=1)
% M = mass matrix (from cu, here c=1)
% The third output (for f) is unused

% === Store Stiffness Matrix ===
p.mat.K = K;

% === Block Diagonal Mass Matrix ===
% Construct mass matrix for the full system:
% - N+2 components (N species, 2 resources: w and h)
% - Each component has its own spatial mass matrix (M)
% - Combine into one big block-diagonal matrix
p.mat.M = kron(diag(ones(1, p.N + 2)), M);
