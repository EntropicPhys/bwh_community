close all; clc;clear all; 
keep p2phome; global p2pglob; p2pglob.gu=[]; p2pglob.nzi=0;p2pglob.ps=1;
%% ------- brute force Busse ballon, extracting stability info from computed branches, chi=0.88
% here passing wave-nr to bfbb; to improve: compute more T-branches, 
kvcm0=[]; pvcm0=[];

[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T1',0.58643) ;[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T2',0.56549);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T3',0.54454) ;[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T4',0.62832); 
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T5',0.60737) ;[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T6',0.5236);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T7',0.64926) ;[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T8',0.67021);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T9',0.73304) ;[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T10',0.77493);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T11',0.81681);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T12',0.87965);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T13',0.90059);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T14',0.46077);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T15',0.43982);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T16',0.41888);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T17',0.39794);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T18',0.79587);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T19',0.83776);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T20',0.8587);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T21',0.75398);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T22',0.71209);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T23',0.48171);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T24',0.69115);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T25',0.50265);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T26',0.35605);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T27',0.3351); [kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T28',0.31416);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T29',0.27227);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T30',0.23038);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T31',0.1885); [kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T32',0.14661);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T33',0.10472);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T34',0.083776);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T35',0.20944);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T36',0.16755);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T37',0.12566);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T38',0.25133);
[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T39',0.29322);[kvcm0 pvcm0]=bfbb(kvcm0,pvcm0,'11/T40',0.37699);

%% plot the BB 
mclf(17); 
hold on
plot(pvcm0,kvcm0,'x','Color','m'); 
plot(pp0,0.07854*l0,'LineWidth',2,'Color','b');
plot(pp1,0.07854*l1,'LineWidth',2,'Color','r');

ylim([0,1.2])
xlabel('P'); ylabel('k'); title('BB','Interpreter','latex'); 
set(gca,'fontsize',14); 
