function Gu=sGjac(p,u) % Jac of single species sG
n=p.np; par=u(p.nu+1:end); Db=par(15); Dw=par(16); Dh=par(17);l=par(19);
[f1b,f1w,f1h,f2b,f2w,f2h,f3b,f3w,f3h]=njac(p,u); % derivative of nonlinearity 
Fu=[[spdiags(f1b,0,n,n),spdiags(f1w,0,n,n),spdiags(f1h,0,n,n)];
    [spdiags(f2b,0,n,n),spdiags(f2w,0,n,n),spdiags(f2h,0,n,n)];
    [spdiags(f3b,0,n,n),spdiags(f3w,0,n,n),spdiags(f3h,0,n,n)]];
Gu=kron([l^2*Db,0,0; 0,l^2*Dw,0; 0,0,l^2*Dh],p.mat.K)... % Gu=diffusion terms 
    - p.mat.M*Fu;   %                                          -local terms 
end

function [f1b,f1w,f1h,f2b,f2w,f2h,f3b,f3w,f3h]=njac(p,u)
%Jacobian, nodal version
n=p.np; b=u(1:n); w=u(n+1:2*n); h=u(2*n+1:3*n);
par=u(p.nu+1:end); pp=par(1); Lam0=par(2); Ga=par(3); A=par(4); R=par(5); L0=par(6);
f=par(7); Q=par(8); Kmin=par(9); Kmax=par(10); Mmin=par(11); Mmax=par(12);
Ymin=par(13); Ymax=par(14); Db=par(15); Dw=par(16); Dh=par(17); chi=par(18);
ov=ones(n,1);
Yi=Ymax + chi*(Ymin-Ymax);
Ki=Kmax + chi*(Kmin-Kmax);
Mi=Mmax + chi*(Mmin-Mmax);
I=A*(Yi*b+f*Q)./(Yi*b+Q);
L=L0*ov./(1+R*b);
Lam=Lam0*Ki*ov./(b+Ki);
dI=A*Yi*ov./(Yi*b+Q) - A*Yi*(Yi*b+f*Q)./(Yi*b+Q).^2;
dLam=-Lam0*Ki*ov./(b+Ki).^2;
dL =-L0*R*ov./(1+R*b).^2;
f1b=dLam.*w.*b + Lam.*w - Mi; f1w=Lam.*b; f1h =0*b;
f2b=dI.*h-Ga*w-dL.*w; f2w=-L - Ga*b; f2h=I;
f3b=-dI.*h; f3w=0*b ; f3h=-I;
end
