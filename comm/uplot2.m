function uplot2(p, wnr)
% uplot2 — Plot ⟨B(x, χ)⟩ over trait space χ
%
% This function computes the spatially averaged biomass ⟨B(x, χ)⟩
% for each trait χ and plots it as a function of χ.
%
% INPUT:
%   p   : PDE2PATH problem structure
%   wnr : figure window number

% === Unpack dimensions and parameters ===
N = p.N;                      % Number of traits
n = p.np;                     % Number of spatial grid points
par = p.u(p.nu+1:end);        % Parameter vector

chimin = par(19);             % Minimum trait value
chimax = par(20);             % Maximum trait value
chi = linspace(chimin, chimax, N);  % Trait grid

% === Retrieve biomass field B(x,χ) ===
b = p.u(1:N*n);               % Full biomass vector
B = reshape(b, n, N);         % Reshape into (space × trait)

% === Mass matrix and volume ===
M = p.mat.M(1:n, 1:n);        % Mass matrix (FEM spatial integration)
vol = p.vol;                  % Spatial domain volume

% === Compute ⟨B(x,χ)⟩ for each χ by spatial integration ===
ba = zeros(1, N);             % Preallocate
for i = 1:N
    ba(i) = sum(M * B(:, i)) / vol;  % Integrate over space
end

% === Plotting ===
figure(wnr); clf;
plot(chi, ba, 'LineWidth', 2);
xlabel('\chi'); ylabel('\langle B(x,\chi) \rangle');
title([p.file.pname, ' pt', num2str(p.file.count - 1)]);
axis tight; set(gca, 'FontSize', 12);
