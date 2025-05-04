function [kv, pv] = bfbb(kv, pv, dir, k0)
% bfbb — Brute Force Busse Balloon
%
% Appends stable solutions along a continuation branch to build
% a Busse balloon curve (i.e., (P, k) values where patterns are stable).
%
% Inputs:
%   kv   - existing vector of wavenumbers (can be empty [])
%   pv   - existing vector of parameter values (e.g., precipitation P)
%   dir  - folder name of the solution branch (e.g., 'comm/T1')
%   k0   - wavenumber associated with this branch (user-supplied or estimated)
%
% Outputs:
%   kv   - extended wavenumber vector (only stable points included)
%   pv   - extended parameter vector for those stable solutions

% === Load branch data from directory ===
p = loadpp(dir);           % Load full solution structure
pv0 = p.branch(4, :);      % Parameter values (e.g., precipitation)
inv = p.branch(3, :);      % Stability indicator: # of unstable eigenvalues

% === Loop through all solution points on the branch ===
lb = length(pv0);
for i = 1:lb
    if inv(i) < 1        % stable solution: no unstable modes
        kv = [kv, k0];   % append wavenumber
        pv = [pv, pv0(i)]; % append parameter
    end
end
