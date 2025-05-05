function duGuph=bpjac(p,u) % second derivative for BP continuation 
n=p.np; b=u(1:n);w=u(n+1:2*n);h=u(2*n+1:3*n);par=u(p.nu+3*p.np+1:end);
pp=par(1); Lam0=par(2); Ga=par(3); A=par(4); R=par(5); L0=par(6);
f=par(7); Q=par(8); Kmin=par(9); Kmax=par(10); Mmin=par(11); Mmax=par(12);
Ymin=par(13); Ymax=par(14); Db=par(15); Dw=par(16); Dh=par(17); chi=par(18);
ov=ones(n,1);Yi=Ymax+chi*(Ymin-Ymax);Ki=Kmax+chi*(Kmin-Kmax);Mi=Mmax+chi*(Mmin-Mmax);

f1bb=-2*(Lam0*w*Ki^2)./(b+Ki).^3; 
f1wb=(Lam0*Ki^2)*ov./(b+Ki).^2;
f2bb=-2*Lam0*R^2*w./(1+R*b).^3 + 2*A*Q*Yi^2*(-1+f)*h./(Q+Yi*b).^3;
f2bw=-ov*Ga + L0*R*ov./(1+R*b).^2; 
f2bh=-A*Q*Yi*(-1+f)*ov./(Q+Yi*b).^2;
f3bb=-2*A*(-1+f)*Q*Yi^2*h./(Q+Yi*b).^3;
f3bh=A*(-1+f)*Q*Yi*ov./(Q+Yi*b).^2;

ph1=u(p.nu+1:p.nu+p.np); ph2=u(p.nu+p.np+1:p.nu+2*p.np);ph3=u(p.nu+2*p.np+1:p.nu+3*p.np);

M1=spdiags(f1bb.*ph1+f2bb.*ph2+f3bb.*ph3,0,n,n);
M2=spdiags(f1wb.*ph1+f2bw.*ph2,0,n,n);
M3=spdiags(f2bh.*ph2+f3bh.*ph3,0,n,n);

duGuph = - [[M1 M2 M3];[M2 0*M2 0*M2];[M3 0*M3 0*M3]]*p.mat.M;





