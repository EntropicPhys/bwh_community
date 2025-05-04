function uplot3(p, wnr)
% uplot3 — Plots trait-integrated biomass ⟨B(.,chi)⟩ and trait-dependent coefficients
% Computes and visualizes:
%   - Trait biomass density ⟨B(.,chi)⟩
%   - Growth coefficient alpha(chi)
%   - Root uptake capacity K(chi)
%   - Mortality M(chi)
%   - Resource-saturated growth rate Lam(chi)

% === Problem dimensions and parameters ===
N = p.N;                  % Number of trait bins
n = p.np;                 % Number of spatial points
par = p.u(p.nu+1:end);    % Extract parameters from p.u

% Trait domain
chimin = par(19);
chimax = par(20);
delchi = (chimax - chimin) / N;
y = linspace(chimin, chimax, N);   % Trait mesh

% Trait shaping exponents
ap = par(21);  % exponent for K(chi)
bp = par(22);  % exponent for M(chi)

% Extract water and biomass
w = p.u(N*n+1);              % Scalar water value (assumed spatially uniform)
b = reshape(p.u(1:N*n), n, N);  % Biomass field reshaped to (x × chi)
M = p.mat.M(1:n,1:n);        % Mass matrix for integration

% Trait-aggregated biomass density
ba = zeros(1, N);
bt = 0;                     % Total biomass integrated over trait
for i = 1:N
    ba(i) = b(1, i);        % (assumption: spatially homogeneous, could use ba(i)=sum(M*b(:,i)))
    bt = bt + ba(i) * delchi;
end

% === Trait-dependent parameters ===
Ki    = zeros(1, N);
Mi    = zeros(1, N);
Lami  = zeros(1, N);
alpha = zeros(1, N);
Lam0  = par(2);
Kmin  = par(9);  Kmax = par(10);
Mmin  = par(11); Mmax = par(12);

for i = 1:N
    chii     = chimin + i * delchi;
    Ki(i)    = Kmax + chii^ap * (Kmin - Kmax);       % Root uptake capacity
    Mi(i)    = Mmax + chii^bp * (Mmin - Mmax);       % Mortality
    Lami(i)  = Lam0 * Ki(i) / (bt + Ki(i));          % Max growth rate
    alpha(i) = Lami(i) * w - Mi(i);                  % Effective trait growth rate
end

% === Plot 1: Biomass density vs chi ===
figure(wnr); clf;
plot(y, ba, '-*', 'LineWidth', 2);
xlabel('\chi'); ylabel('<B(.,\chi)>');
title([p.file.pname, mat2str(p.file.count - 1)]);
set(gca, 'FontSize', 12); axis tight;

% === Plot 2: Alpha(chi) ===
figure(11); clf;
plot(y, alpha, 'LineWidth', 2);
xlabel('\chi'); ylabel('\alpha(\chi)');
title('\alpha(\chi) = \Lambda(\chi) w - M(\chi)');
set(gca, 'FontSize', 12); axis tight; grid on;

% === Plot 3: Root uptake capacity K(chi) ===
figure(12); clf;
plot(y, Ki, 'LineWidth', 2);
xlabel('\chi'); ylabel('K(\chi)');
title('Trait-dependent uptake capacity');
set(gca, 'FontSize', 12); axis tight; grid on;

% === Plot 4: Mortality M(chi) ===
figure(13); clf;
plot(y, Mi, 'LineWidth', 2);
xlabel('\chi'); ylabel('M(\chi)');
title('Trait-dependent mortality');
set(gca, 'FontSize', 12); axis tight; grid on;

% === Plot 5: Growth rate Lambda(chi)*w ===
figure(14); clf;
plot(y, Lami * w, 'LineWidth', 2);
xlabel('\chi'); ylabel('\Lambda(\chi) w');
title('Saturated growth rate');
set(gca, 'FontSize', 12); axis tight; grid on;

% === Output total biomass for record-keeping ===
bt
