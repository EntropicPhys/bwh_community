% Bifurcation point continuation from subcritical Turing bifurcation
% tracks the location of bifurcation point (bpt5) on the T1 branch in a 
% 2D parameter space, yielding (part of) the bdry of the Busse balloon 
close all; keep p2phome;
%% === Init: continue 'bpt5' on branch s1/T1 wrt to parameters 1 and 19
%  - 1=precipitation (P)
%  - 19=l (domain length rescaling or trait scaling)
% Output folder: 's1/bpa'
p=bpcontini('s1/T1','bpt6',19,'s1/bpa',1e-3); % Start from BP5
p.plot.bpcmp=5;      % Plot 5th component during continuation
p.nc.dsmax=0.1; p.nc.dsmin=1e-5;  % Max and min step size 
p.nc.del  =0.2;     % finite difference discretization of extra terms
p.nc.lammin=0;       % Lower bound for continuation parameter
p.sw.spcalc=0; p.sw.bifcheck=0; % Disable spectral comp and bif detection
p.nc.tol=1e-8;        % Newton tolerance
p.sw.spjac=0;         % use numjac for 2nd directional derivatives 
hucl; p=cont(p, 60);  % clear windows and run
%% === Reverse Continuation ===
p=loadp('s1/bpa','pt0','s1/bpb'); % Load initial point from bpa
p.sol.ds=-1.3*p.sol.ds; p=cont(p,20); % Reverse direction and go 
%% === Load Data and Extract Bifurcation Path ===
% Load branch point continuation results
p0=loadpp('s1/bpa'); 
pp0=p0.branch(11,:);   % Parameter 1 values (e.g. precipitation P)
l0 =p0.branch(4,:);    % Parameter 2 values (e.g. domain length l)
p1=loadpp('s1/bpb'); pp1=p1.branch(11,:); l1 =p1.branch(4,:);
%% === Plotting the Busse Balloon ===
% Compute wavenumber from domain length:
% Assume 1D domain of length lx=80 with 6.5 wavelengths → wl=80 / 6.5
% Then k=2π / wl=2π / (80 / 6.5) ≈ 0.51051
kfactor=0.51051;  % Converts domain rescaling l into wavenumber k=kfactor * l
mclf(1); hold on; box on 
plot(pp0, kfactor * l0, 'LineWidth', 2, 'Color', 'b');  % Forward BP curve
plot(pp1, kfactor * l1, 'LineWidth', 2, 'Color', 'b');  % Reverse BP curve
% Axes and labels
xlabel('Precipitation'); ylabel('Wavenumber'); title('Upper Busse Balloon chi=1');
set(gca, 'FontSize', 14);
