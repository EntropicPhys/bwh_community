function uplot1(p, wnr)
% uplot1 — Specialized visualization of trait-structured vegetation model
% Plots:
% 1. B(x,chi): biomass over space and trait space (as a surface plot)
% 2. w(x): soil water over space
% 3. h(x): surface water over space
% 4. ⟨B⟩(chi): trait-marginalized biomass (integrated over space)

% === Grid and Dimensions ===
N  = p.N;           % Number of traits/species
n  = p.np;          % Number of spatial grid points
par = p.u(p.nu+1:end);  % Extract parameter vector
chimin = par(19); chimax = par(20);  % Trait bounds

% === Trait-space and spatial mesh ===
y = linspace(chimin, chimax, N);  % Trait mesh (chi)
x = getpte(p);                    % Spatial mesh
[X, Y] = meshgrid(y, x);          % For 2D surface plot

% === Extract solution components ===
b = p.u(1:N*n);          % Biomass vector (concatenated)
B = reshape(b, n, N);    % Reshape to [space × trait]
bmax = 3;                % Max biomass value for plot scaling

w = p.u(N*n+1:(N+1)*n);              % Soil water
h = p.u((N+1)*n+1:(N+2)*n);          % Surface water

% === Figure layout (9x6 grid) and subpanel indices ===
mclf(wnr);                       % Clear and set figure window
p1 = [1:4 7:10 13:16 19:22 25:28];  % Surface plot
p2 = 37:40;                         % Water w(x)
p3 = 49:52;                         % Surface water h(x)
p4 = 6:6:28;                        % Trait-marginal biomass ⟨B⟩(chi)

x1 = min(x); x2 = max(x);          % Spatial domain bounds

% === 1. B(x, chi): Biomass over space and trait space ===
subplot(9, 6, p1);
surf(X, Y, B);
caxis([0 bmax]);
view(90, 270); shading interp; axis tight;
colormap(flipud(summer));
title([p.file.pname, mat2str(p.file.count - 1)]);
colorbar('southoutside');
set(gca, 'fontsize', 12);

% === 2. w(x): Soil water profile ===
subplot(9, 6, p2);
plot(x, w, 'linewidth', 2);
axis tight; title('w');
set(gca, 'fontsize', 12);

% Auto y-ticks for w
try
    z1 = round(min(w) * 105) / 100;
    z2 = round(max(w) * 95) / 100;
    yticks([z1 z2]);
catch
    try
        z1 = 0.9 * min(w); z2 = 1.1 * max(w);
        axis([x1 x2 z1 z2]);
        yticks(round(50 * (z1 + z2)) / 100);
    catch
    end
end

% === 3. h(x): Surface water profile ===
subplot(9, 6, p3);
plot(x, h, 'linewidth', 2);
axis tight; title('h');
set(gca, 'fontsize', 12);

% Auto y-ticks for h
try
    z1 = round(min(h) * 105) / 100;
    z2 = round(max(h) * 95) / 100;
    yticks([z1 z2]);
catch
    try
        z1 = 0.9 * min(h); z2 = 1.1 * max(h);
        axis([x1 x2 z1 z2]);
        yticks(round(50 * (z1 + z2)) / 100);
    catch
    end
end

% === 4. ⟨B⟩(chi): Biomass integrated over space ===
M = p.mat.M(1:n, 1:n);   % Mass matrix (single component)
vol = p.vol;             % Domain volume (length)
ba = zeros(size(y));     % Preallocate spatial integrals
for i = 1:N
    ba(i) = sum(M * B(:, i)) / vol;
end

subplot(9, 6, p4);
plot(ba, y, 'linewidth', 2);
xlim([0 bmax]);
xlabel('<B(·,χ)>'); ylabel('\chi');
set(gca, 'fontsize', 12);
grid on;
