close all; keep p2phome; % script for 1D 
%% continuation of the subcritical BP for T1 up
p=bpcontini('s1/T1','bpt4',19,'bps1a',3e-3); % bif. point. cont
p.plot.bpcmp=5; p.sw.spjac=0; p.nc.dsmax=0.1;p.nc.lammin=0;p.nc.del=4e-3;
p.nc.dsmin=1e-5; p.fuha.spjac=@bpjac; p.sw.spcalc=0; p.sw.bifcheck=0;
p.nc.almine=0.4;p.nc.tol=9e-9; p.file.smod=10; 
huclean(p); p=cont(p,20);
%% continuation in the reverse direction 
p=loadp('bps1a','pt0');p.sol.ds=-3*p.sol.ds;p=setfn(p,'bps1b');p=cont(p,20);
%% branch point continuation data
p0=loadpp('bps1a'); pp0=p0.branch(11,:);l0=p0.branch(4,:);
p1=loadpp('bps1b'); pp1=p1.branch(11,:);l1=p1.branch(4,:);

%% plot of branch point cont. results
% wl=lx/n; k=2*pi/wl; our case wl=90/7;k=2*pi/wl=0.48869

hold on
plot(pp0,0.48869*l0,'LineWidth',2,'Color','b');
plot(pp1,0.48869*l1,'LineWidth',2,'Color','b');
ylim([0,0.85])
xlabel('P'); ylabel('k'); title('Busee-Balloon Single Species','Interpreter','latex'); 
set(gca,'fontsize',14); 