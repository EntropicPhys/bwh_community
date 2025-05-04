% bpcontT1.m — Bifurcation point continuation from subcritical Turing bifurcation
% This script tracks the location of a subcritical bifurcation point (bpt6)
% on the T1 branch in a 2D parameter space, forming a Busse balloon domain.

close all; keep p2phome;

%% === Initialization: Bifurcation Point Continuation ===

% Initialize continuation from bifurcation point 'bpt6' on branch s1/T1
% We continue with respect to parameter indices 1 and 19
%  - 1 = precipitation (P)
%  - 19 = l (domain length rescaling or trait scaling)
% Output folder: 'bps1a'

p = bpcontini('s1/T1', 'bpt6', 19, 'bps1a', 3e-3);  % Start from bif point with ds = 3e-3

% Plotting setup
p.plot.bpcmp = 5;      % Plot 5th component during continuation
p.sw.spjac = 1;        % Use sparse Jacobian for efficiency
p.nc.dsmax = 1;        % Max step size
p.nc.dsmin = 1e-5;     % Min step size
p.nc.del   = 4e-3;     % Arclength increment constraint (step size control)
p.nc.lammin = 0;       % Lower bound for continuation parameter

% Jacobian for bifurcation point continuation
p.fuha.spjac = @bpjac;

% Disable spectrum calculation and bifurcation detection
p.sw.spcalc = 0;
p.sw.bifcheck = 0;

% Additional numerical tuning
p.nc.almine = 0.4;      % Minimum arclength weight for continuation
p.nc.tol = 8e-7;        % Newton tolerance
p.file.smod = 10;       % Save every 10th solution point

% Clean up old files if present and start continuation
huclean(p);             % Remove any previous output files in the folder
p = cont(p, 600);       % Continue for 600 steps

%% === Reverse Continuation ===

% Load initial point from the just-computed branch
p = loadp('bps1a', 'pt0');
p.sol.ds = -1.3 * p.sol.ds;  % Reverse direction with slightly larger step
p = setfn(p, 'bps1b');       % Set new output folder for reverse run
p = cont(p, 200);            % Continue for 200 steps in reverse

%% === Load Data and Extract Bifurcation Path ===

% Load branch point continuation results
p0 = loadpp('bps1a'); 
pp0 = p0.branch(11,:);   % Parameter 1 values (e.g. precipitation P)
l0  = p0.branch(4,:);    % Parameter 2 values (e.g. domain length l)

p1 = loadpp('bps1b'); 
pp1 = p1.branch(11,:);
l1  = p1.branch(4,:);

%% === Plotting the Busse Balloon ===

% Compute wavenumber from domain length:
% Assume 1D domain of length lx = 80 with 6.5 wavelengths → wl = 80 / 6.5
% Then k = 2π / wl = 2π / (80 / 6.5) ≈ 0.51051

kfactor = 0.51051;  % Converts domain rescaling l into wavenumber k = kfactor * l

figure;
hold on;box on 
plot(pp0, kfactor * l0, 'LineWidth', 2, 'Color', 'b');  % Forward BP curve
plot(pp1, kfactor * l1, 'LineWidth', 2, 'Color', 'b');  % Reverse BP curve

% Axes and labels
xlabel('Precipitation');
ylabel('Wavenumber');
title('Busse Balloon/Single Species');
set(gca, 'FontSize', 14);
