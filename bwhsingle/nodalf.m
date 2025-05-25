function r=nodalf(p,u) % for bwh, N=1, i.e., only space, here used as a dummy
n=p.np; b=u(1:n); w=u(n+1:2*n); h=u(2*n+1:3*n); 
par=u(p.nu+1:end); pp=par(1); Lam0=par(2); Ga=par(3); A=par(4); R=par(5); L0=par(6); 
f=par(7); Q=par(8); Kmin=par(9); Kmax=par(10); Mmin=par(11); Mmax=par(12); 
Ymin=par(13); Ymax=par(14); Db=par(15); Dw=par(16); Dh=par(17); chi=par(18); l=par(19);
K=p.mat.K; M=p.mat.M(1:n,1:n); ov=ones(n,1);

Yi=Ymax+chi*(Ymin-Ymax); 
Mi=Mmax+chi*(Mmin-Mmax); 
Ki=Kmax+chi*(Kmin-Kmax);
Lam=Lam0*ov*Ki./(b+Ki); 
I=A*(Yi*b+f*Q)./(Yi*b+Q); 
L=L0*ov./(1+R*b);  
r1=M*(Lam.*w.*b-Mi*b); 
r2=M*(I.*h-L.*w-Ga*w.*b); 
r3=M*(pp-I.*h);   
r=[r1;r2;r3]; 