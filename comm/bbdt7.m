%% === Transient Dynamics under Precipitation Changes ===
% This script simulates vegetation responses to time-varying precipitation,
% capturing transitions, stability loss, and pattern evolution using PDE2PATH.

%% === INITIALIZATION: Short Time Integration with Decreasing Precipitation ===
p = loadp('11/T9', 'pt80', '14/dns20_rate1_1yr_k9');
p.pa = 0.01;     % Precipitation step (mm/yr)
p.T = 1;         % Duration (years) per step
p.incr = -1;     % Direction of change: decrease
p.threshp = 0.03;
p.file.count = 0;
nst = 400; ds1 = [];

[p, ds1] = bbdns(p, ds1, nst);  % Time integration

%% === Short Increasing Precipitation (Rewetting Experiment) ===
p = loadp('14/dns15_rate1_1yr_kc', 'pt284', '14/dnsR15_rate1_1yr_kc');
p.pa = 0.1; p.T = 1; p.incr = 1; p.threshp = 0.03; nst = 3000; ds2 = [];
[p, ds2] = bbdns(p, ds2, nst);

%% === Medium Duration Increase ===
p = loadp('14/dnsR5_rate1_k2', 'pt402', '14/dnsR11_r1');
p.pa = 0.1; p.T = 10; p.incr = 0.5; p.threshp = 0.03; nst = 300; ds3 = [];
[p, ds3] = bbdns(p, ds3, nst);

%% === Long-Term Precipitation Increase (Mixed Modes Region) ===
p = loadp('9/T2', 'pt40', '9/dns1');
p.pa = 0.5; p.T = 10; p.incr = 1; p.threshp = 0.05; nst = 10; ds3 = [];
[p, ds3] = bbdns(p, ds3, nst);

%% === Reload Simulation and Continue with New Step ===
p = loadp('9/dns1', 'pt40', '9/dns1');
p.pa = 0.5; p.T = 10; p.incr = 1;
p.threshp = 0.05; p.file.count = 0; nst = 10; ds4 = [];
[p, ds4] = bbdns(p, ds4, nst);

%% === SAVE DATA ===
save('bb1', 'ds1', 'ds2', 'ds3', 'ds4');

%% === LOAD DATA ===
load('bb1', 'ds3');  % Reload just one dataset (optional)

%% === QUICK PLOT: Chi-max over Precipitation ===
mclf(20); 
plot(ds1(1,:), ds1(3,:), 'b'); hold on;
plot(ds2(1,:), ds2(3,:), 'r');
xlabel('P (mm/yr)'); ylabel('\chi_m'); ylim([0.4 0.65]);

%% === BUSSE BALLOON PLOTS ===
% kcmp: how to compute k (4=FFT, 5=patch count, 6=hybrid)
% cc: quantity to plot (2=||B||, 3=chi_m)
wnr = 14; kc = 2*pi / (p.vol / 14); cc = 3; kcmp = 6;
mclf(wnr); hold on;
plotds(ds1, kcmp, kc, cc, 0, 0); 
plotds(ds2, kcmp, kc, cc, 0, 0); 
plotds(ds3, kcmp, kc, cc, 0, 2); 

xlabel('P (mm/yr)'); ylabel('k/k_c');
colormap('turbo'); clim([0.4 0.65]);
title(['Precipitation change: \pm' num2str(abs(p.incr)) ' mm/yr every ' num2str(p.T) ' yrs']);

%% === ADD BRUTE-FORCE BUSSE BALLOON POINTS ===
fig(10); hold on;
plot(pv, kv / kc, '*');

kv = []; pv = [];
[kv, pv] = bfbb2(kv, pv, '9/T1');
[kv, pv] = bfbb2(kv, pv, '7/T4');
[kv, pv] = bfbb2(kv, pv, '7/T5');
[kv, pv] = bfbb2(kv, pv, '7/T6');
[kv, pv] = bfbb2(kv, pv, '7/T8');
[kv, pv] = bfbb2(kv, pv, '7/T9');

%% === PATCH COUNTING AND FREQUENCY ANALYSIS (Loop over Time Series) ===
p.threshp = 0.1; ds5 = []; kc = 2*pi / (p.vol / 14);
frames = struct();  % to store video frames

for ii = 285:569
    ff = sprintf('pt%1d', ii);
    p = loadp('14/dnsR10_01_rate1_kc', ff);
    P = p.u(p.nu+1);

    % Extract biomass B(x, χ) and compute spatial profile
    N = p.N; n = p.np;
    x = getpte(p)';
    br = reshape(p.u(1:N*n), n, N);
    b = sum(br, 2)';

    % Extract and FFT surface water
    h = p.u(n*(N+1)+1 : n*(N+2)); 
    h = h - mean(h); 
    hh = abs(fft(h)); 
    hh = hh(1:floor(n/2)+1);
    hm = movmean(hh, 2);
    [~, idx] = max(hm);
    fs = n / p.vol;
    fr = (0:floor(n/2)) * fs / n;
    km = 2*pi * fr(idx);

    % Patch counting
    threshold = min(b) + p.threshp * max(b);
    above_thresh = abs(b) > threshold;
    start_patches = find(diff([false; above_thresh(:)]) == 1);

    if above_thresh(1) && above_thresh(end)
        km2 = 2*pi * ((length(start_patches) - 1) / p.vol);
    else
        km2 = 2*pi * (length(start_patches) / p.vol);
    end

    ds5 = [ds5, [P; km; km2]];

    % === Diagnostic Plotting ===
    mclf(20);
    subplot(4,2,1); imagesc(br'); colormap(flipud(summer));
    xticks([1 n-1]); yticks([chimin chimax]); xticklabels([0 p.vol]); 
    yticklabels([chimin chimax]); xlabel('x'); ylabel('\chi');

    subplot(4,2,2); plot(2*pi*fr(:), hh(:)); xlabel('f'); ylabel('A(f)');
    subplot(4,2,[3,4]); plot(h); ylabel('H');
    subplot(4,2,[5,6]); plot(b); xlabel('x'); ylabel('B');
    subplot(4,2,[7,8]); plot(ds4(1,:), ds4(3,:)); xlim([50 450]); ylim([0 1]); xlabel('P'); ylabel('k'); hold on;

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12);
    frames(ii) = getframe(gcf);  % Save for video
    pause(0.01);

    % Early stop if vegetation disappears
    if max(b) - min(b) < 1e-2
        disp('Vegetation collapsed'); break;
    end
end

%% === EXPORT VIDEO FROM FRAMES ===
v = VideoWriter('myMovie.avi');
v.FrameRate = 10;
open(v);
for k = 285:569
    writeVideo(v, frames(k));
end
close(v);
