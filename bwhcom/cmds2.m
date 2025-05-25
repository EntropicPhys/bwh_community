%% === Brute-force Busse Balloon Construction 
% scan multiple Turing branches and collects stable points
% for plotting the stability region (Busse Balloon)
kvcm0=[];  % Wave numbers 
pvcm0=[];  % Precipitation P values
% List of folders and associated wavenumbers
branchList={'1/T1',0.59; '1/T2',0.52; '1/T1b',0.55; '1/T1bb',0.55;
    '1/T2b',0.42; '1/T2bb',0.42; '1/T2C',0.664; '1/T3D',0.838};
for i=1:size(branchList, 1) % Loop through list and collect stable solutions
    dir=branchList{i, 1}; k0 =branchList{i, 2};
    [kvcm0, pvcm0]=bfbb(kvcm0, pvcm0, dir, k0);
end
%% === Plot Busse Balloon
mclf(10); hold on;
plot(pvcm0, kvcm0, 'x', 'Color', 'm');  % Stable points from sampling 
ylim([0, 0.85]); xlabel('Precipitation $P$', 'Interpreter', 'latex');
ylabel('Wavenumber $k$', 'Interpreter', 'latex');
title('Busse balloon sampling', 'Interpreter', 'latex');
set(gca, 'FontSize', 14); box on; 
%% === DNS through BB with change of P 
p=loadp('1/T3D','pt60','1/1yra'); % load start point and set up output folder
p.u(p.nu+1)=230; % set P to some convenient starting value
p.incr=-0.5; p.T=1;  % increment and time-interval in precipitaion.  
p.threshp=0.1;  % threshold to define deviation from baresoil, i.e., patches
nst=400; p.file.count=0;  % #time steps, and counter for saved output files
p.pa=0.05;         % amplitude of perturbations in H 
ds1=[];            % for output data
[p,ds1]=bbdns(p,ds1,nst); % go! (repeat if useful) 
%% backward: load point before dropping to bare soil, reverse rate
dir='1/1yra'; lastpt=['pt' mat2str(max(getlabs(dir))-10)];
p=loadp(dir,lastpt,'1/1yrb'); p.T=1; p.incr=-p.incr; nst=100; ds2=[]; 
%% go
[p,ds2]=bbdns(p,ds2,nst);           % (repeat if useful) 
%save('bb1','ds1','ds2'); % save results (useful for large-scale long DNS) 
%load('bb1','ds1');                  % Load specific results
%% plotting the data; 
wnr=10; mclf(wnr); kc=1;  % new figure, and wave-nr normalization;
%kc=2*pi/(p.vol/14);    % uncomment for normalizing k with kc
cc=3; ctit='\chi_c'; % component for colormap 2=||B||, 3=chicent, 4=chicent
%cc=2; ctit='||B||'; % uncomment to choose B 
kcmp=5;   % choose the "wavenumber k--variable" to plot 
% kcmp: 1=P; 2=||B||; 3=chi_center; 4=k from fft; 5=k from patch counting; 
plotds(ds1,kcmp,kc,cc,0,0);   % plot
h=plotds(ds2(:,1:end-1),kcmp,kc,cc,0,0);  % Run 2
xlabel('P (mm/yr)'); ylabel('k'); incr=1; T=1; 
%clim([0.4 0.65]);    % fixing colorscale, uncomment if desired
title(['P changes by \pm' mat2str(abs(incr)) 'mm/yr every ' mat2str(T) ' year']);
title(h,ctit,'fontsize',14); set(gca,'fontsize',14); 
plot(pvcm0, kvcm0, 'x', 'Color', 'm'); % overlay BB