function [p, ds] = bbdns(p, ds, nst)
% BBDNS — Time integration with continuously varying precipitation
% 
% Inputs:
%   p   — PDE2PATH problem structure (includes parameters and initial state)
%   ds  — Data storage array [P; ||B||; chi_max; k_fft; k_patch; k_patch_corrected]
%   nst — Number of outer steps (precipitation change iterations)
%
% Outputs:
%   p   — Updated PDE2PATH structure
%   ds  — Updated data with diagnostic values per step

% Extract initial precipitation
P = p.u(p.nu+1); 

% Spatial setup
N = p.N; n = p.np; nx = n;
x = getpte(p)'; lx = max(x);
xn = 2 * pi * x / lx;   % Rescaled space for periodic functions

% Time integration setup
dt = 0.005; 
nt = floor(p.T / dt);  % Number of timesteps per unit T
pmod = nt; 
smod = 1e9;            % No intermediate saves

% === Main time-stepping loop with precipitation variation ===
for i = 1:nst
    % Update precipitation value
    P = P + p.incr;  
    p.u(p.nu+1) = P;

    % === Perturb surface water field H before integration ===
    if 0
        % Optionally add random noise (unused here)
        p.u((N+1)*n+1 : (N+2)*n) = p.u((N+1)*n+1 : (N+2)*n) + rand([n,1]); 
    else
        % Add smooth perturbation using cosine modes (trait-like pattern)
        pert = zeros(size(xn));
        pwn = [0.5, 1:2:10];       % Frequencies (cos modes)
        pav = 0.1 * (1:6);         % Amplitudes
        for l = 1:6
            pert = pert + pav(l) * cos(pwn(l) * xn);
        end
        H = p.u((N+1)*n+1 : (N+2)*n);
        p.u((N+1)*n+1 : (N+2)*n) = H .* (1 + p.pa * pert);  % Apply perturbation
    end

    % === Time integration ===
    nc = 0; 
    t1 = 0; ts = [];
    [p, t1, ts, nc] = tintxs(p, t1, ts, dt, nt, nc, pmod, smod, @sGdns);

    % === Post-integration diagnostics ===
    out = sgbra(p, p.u);  % [||B||, ..., chi, ...]
    jjj = p.u(1:N*n); 
    br = reshape(jjj, n, N); 
    b = sum(br, 2)';       % Total biomass (spatial)

    % Extract surface water
    h = p.u(n*(N+1)+1 : n*(N+2));
    bdiff = max(b) - min(b);  % For detecting biomass collapse

    % === PATCH COUNTING: Based on surface water (k_patch ~ km2) ===
    threshold = min(h) + p.threshp * max(h);
    above_threshold = abs(h) > threshold;
    start_patches = find(diff([false; above_threshold(:)]) == 1);
    km2 = 2 * pi * (length(start_patches) / lx);  % Simple patch wavenumber

    % === FREQUENCY ANALYSIS (FFT) on surface water ===
    h = h - mean(h); hh = abs(fft(h)); 
    hh = hh(1:floor(nx/2)+1);
    hm = movmean(hh, 2); 
    [~, idx] = max(hm); 
    fs = nx / p.vol; 
    fr = (0:floor(nx/2)) * fs / nx;
    km = 2 * pi * fr(idx);     % Peak frequency mode from FFT

    % === PATCH COUNTING: Based on biomass (k_patch_corrected ~ km3) ===
    threshold = min(b) + p.threshp * max(b);
    above_threshold = abs(b) > threshold;
    start_patches = find(diff([false; above_threshold(:)]) == 1);

    if above_threshold(1) && above_threshold(end)
        km3 = 2 * pi * ((length(start_patches) - 1) / lx);  % Two boundary overlaps
    elseif above_threshold(1) || above_threshold(end)
        km3 = 2 * pi * ((length(start_patches) - 0.5) / lx);  % One boundary overlap
    else
        km3 = 2 * pi * (length(start_patches) / lx);         % Full patches only
    end

    % === Append results: [P; ||B||; chi; k_fft; k_patch; k_patch_corrected] ===
    bt = out(1); chi = out(3); 
    ds = [ds, [P; bt; chi; km; km2; km3]];

    % Save current state and check for extinction
    p.file.count = p.file.count + 1;
    p.fuha.savefu(p); 

    if bdiff < 1e-2
        disp('Biomass collapse detected. Ending run.');
        break;
    end
end
