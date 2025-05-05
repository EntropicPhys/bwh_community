close all; clc; clear all;
keep p2phome;
global p2pglob;
p2pglob.gu = [];
p2pglob.nzi = 0;
p2pglob.ps = 1;

%% === Brute-force Busse Balloon Construction (chi = 0.88) ===
% This script scans multiple Turing branches and collects stable points
% for plotting the stability region (Busse Balloon)

kvcm0 = [];  % Wavenumber values
pvcm0 = [];  % Precipitation (or control parameter) values

% List of folders and associated wavenumbers
branchList = {
    '11/T1',  0.58643;
    '11/T2',  0.56549;
    '11/T3',  0.54454;
    '11/T4',  0.62832;
    '11/T5',  0.60737;
    '11/T6',  0.5236;
    '11/T7',  0.64926;
    '11/T8',  0.67021;
    '11/T9',  0.73304;
    '11/T10', 0.77493;
    '11/T11', 0.81681;
    '11/T12', 0.87965;
    '11/T13', 0.90059;
    '11/T14', 0.46077;
    '11/T15', 0.43982;
    '11/T16', 0.41888;
    '11/T17', 0.39794;
    '11/T18', 0.79587;
    '11/T19', 0.83776;
    '11/T20', 0.8587;
    '11/T21', 0.75398;
    '11/T22', 0.71209;
    '11/T23', 0.48171;
    '11/T24', 0.69115;
    '11/T25', 0.50265;
    '11/T26', 0.35605;
    '11/T27', 0.33510;
    '11/T28', 0.31416;
    '11/T29', 0.27227;
    '11/T30', 0.23038;
    '11/T31', 0.18850;
    '11/T32', 0.14661;
    '11/T33', 0.10472;
    '11/T34', 0.083776;
    '11/T35', 0.20944;
    '11/T36', 0.16755;
    '11/T37', 0.12566;
    '11/T38', 0.25133;
    '11/T39', 0.29322;
    '11/T40', 0.37699;
};

% === Loop through all branches and collect stable solutions
for i = 1:size(branchList, 1)
    dir = branchList{i, 1};
    k0  = branchList{i, 2};
    [kvcm0, pvcm0] = bfbb(kvcm0, pvcm0, dir, k0);
end

%% === Plot Busse Balloon
mclf(17);
hold on;
plot(pvcm0, kvcm0, 'x', 'Color', 'm');  % Stable points from brute force

% Optional: overlay bifurcation curves (must be defined elsewhere)
% Example: overlay continuation of bifurcation points (Turing threshold)
% plot(pp0, 0.07854 * l0, 'LineWidth', 2, 'Color', 'b');
% plot(pp1, 0.07854 * l1, 'LineWidth', 2, 'Color', 'r');

% Axes settings
ylim([0, 1.2]);
xlabel('Precipitation $P$', 'Interpreter', 'latex');
ylabel('Wavenumber $k$', 'Interpreter', 'latex');
title('Busse Balloon (Stable Pattern Region)', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

