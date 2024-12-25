function J=sGjac(p,u)
N=p.N;n=p.np;ng=(N+2)*n;   % field assign. 
w=u(N*n+1:(N+1)*n);h=u((N+1)*n+1:(N+2)*n); 
par=u(p.nu+1:end);pp=par(1);Lam0=par(2);Ga=par(3); % syst. params
A=par(4);R=par(5);L0=par(6);f=par(7);Q=par(8);
Kmin=par(9);Kmax=par(10);Mmin=par(11);Mmax=par(12); 
Ymin=par(13);Ymax=par(14);Db=par(15);Dw=par(16);
Dh=par(17);Dchi=par(18); chimin=par(19); chimax=par(20);  
K=p.mat.K; M=p.mat.M(1:n,1:n); % Stiffness and Mass matrices
bt=zeros(n,1); btt=bt; Yi=zeros(N,1); % b-tilde and b-tilde-tilde init.
Ki=Yi;Mi=Yi; chii=linspace(chimin,chimax,N);dchi=chii(2)-chii(1);
for i=1:N % b-tilde and b-tilde-tilde assig.
    Bi=u((i-1)*n+1:i*n);
    bt=bt+Bi; 
    Yi(i)=Ymax+chii(i)*(Ymin-Ymax); 
    Ki(i)=Kmax+chii(i)*(Kmin-Kmax); 
    Mi(i)=Mmax+chii(i)*(Mmin-Mmax); 
    btt=btt+Yi(i)*Bi; 
end
ov=ones(n,1);I=A*(btt+f*Q)./(btt+Q); % infil. and evapo. funcs and dev.  
L=ov*L0./(1+R*bt);dL=-L0*R*ov./(1+R*bt).^2; 
if p2pglob.nzi==0; J=zeros(ng); else J=p2pglob.gu; end % for allocation! 
for i=1:N  % finite-difference scheme for the Jacobian.
    Bi    = u((i-1)*n+1:i*n); 
    Lami  = Lam0*Ki(i)*ov./(bt+Ki(i)); 
    dLami = -Lam0*Ki(i)*ov./(bt+Ki(i)).^2;     
    dI    = (A*Yi(i)*ov./(btt+Q) - A*Yi(i)*(btt+f*Q)./(btt+Q).^2);
    for j=1:N % D_j f_i    
     djfi=dLami.*w.*Bi+(Lami.*w-Mi(i))*(i==j);   
     J((i-1)*n+1:i*n,(j-1)*n+1:j*n)=-M*spdiags(djfi,0,n,n)+(i==j)*Db*K;   
    end
    dwfi=Lami.*Bi; % w-derivative of fi  
    J((i-1)*n+1:i*n,ng-2*n+1:ng-n)=-M*spdiags(dwfi,0,n,n);          
    dbfw=dI.*h-dL.*w-Ga*w; 
    J(N*n+1:(N+1)*n,(i-1)*n+1:i*n)=-M*spdiags(dbfw,0,n,n); % b-der. of fw   
    dbfh=-dI.*h; 
    J((N+1)*n+1:(N+2)*n,(i-1)*n+1:i*n)=-M*spdiags(dbfh,0,n,n); % b-der. of fh
    switch i  % add chi-diffusion, distinguish bulk from boundaries 
  case 1; J(1:n,1:n)=J(1:n,1:n)-(Dchi/dchi^2)*M*spdiags(ov,0,n,n); 
      J(1:n,n+1:2*n)=J(1:n,n+1:2*n)+(Dchi/dchi^2)*M*spdiags(ov,0,n,n); 
  case N;  J((N-1)*n+1:N*n,(N-2)*n+1:(N-1)*n)=...
          J((N-1)*n+1:N*n,(N-2)*n+1:(N-1)*n)-(Dchi/dchi^2)*M*spdiags(ov,0,n,n); 
      J((N-1)*n+1:N*n,(N-1)*n+1:N*n)=...
          J((N-1)*n+1:N*n,(N-1)*n+1:N*n)+(Dchi/dchi^2)*M*spdiags(ov,0,n,n); 
  otherwise; 
     J((i-1)*n+1:i*n,(i-2)*n+1:(i-1)*n)=...
      J((i-1)*n+1:i*n,(i-2)*n+1:(i-1)*n)-(Dchi/dchi^2)*M*spdiags(ov,0,n,n); 
     J((i-1)*n+1:i*n,(i-1)*n+1:i*n)=...
      J((i-1)*n+1:i*n,(i-1)*n+1:i*n)+2*(Dchi/dchi^2)*M*spdiags(ov,0,n,n); 
     J((i-1)*n+1:i*n,i*n+1:(i+1)*n)=...
      J((i-1)*n+1:i*n,i*n+1:(i+1)*n)-(Dchi/dchi^2)*M*spdiags(ov,0,n,n); 
     end   
end     
dwfw=-(L+Ga*bt); 
J(N*n+1:(N+1)*n,N*n+1:(N+1)*n)=Dw*K-M*spdiags(dwfw,0,n,n); % w-der of fw
dhfw=I; 
J(N*n+1:(N+1)*n,(N+1)*n+1:(N+2)*n)=-M*spdiags(dhfw,0,n,n); % w-der of fw
dhfh=-I; 
J((N+1)*n+1:(N+2)*n,(N+1)*n+1:(N+2)*n)=Dh*K-M*spdiags(dhfh,0,n,n); % h-der of f
if p2pglob.nzi==0; J=sparse(J); p2pglob.gu=J; p2pglob.nzi=1; end 
end 

