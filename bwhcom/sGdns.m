function r=sGdns(p,u) % for DNS of community model; essentially rhs without 
% spatial diffusion, i.e., chi-diffusion included here 
N=p.N; n=p.np; b=reshape(u(1:N*n),n,N); w=u(N*n+1:(N+1)*n); h=u((N+1)*n+1:(N+2)*n); 
par=u(p.nu+1:end); pp=par(1); Lam0=par(2); Ga=par(3); A=par(4); R=par(5); L0=par(6); 
f=par(7); Q=par(8); Kmin=par(9); Kmax=par(10); Mmin=par(11); Mmax=par(12); 
Ymin=par(13); Ymax=par(14); Dchi=par(18); chimin=par(19); chimax=par(20); 
M=p.mat.M(1:n,1:n); dchi=(chimax-chimin)/(N-1); r=zeros(N*n,1); bt=zeros(n,1); btt=bt; 
for i=1:N 
    chii=chimin+(i-1)*dchi; bt=bt+b(:,i); Yi=Ymax+chii*(Ymin-Ymax); 
    btt=btt+b(:,i)*Yi;
end

for i=1:N
    chii=chimin+(i-1)*dchi; 
    Mi=Mmax+chii*(Mmin-Mmax); Ki=Kmax+chii*(Kmin-Kmax); % pause
    Lami=Lam0*Ki./(bt+Ki); I=A*(btt+f*Q)./(btt+Q); L=L0./(1+R*bt);   
    bi=b(:,i);    
    switch i
        case 1; bcc=(b(:,2)-2*b(:,1))/dchi^2; 
        case N;  bcc=(-2*b(:,N)+b(:,N-1))/dchi^2; 
        otherwise; bcc=(b(:,i+1)-2*b(:,i)+b(:,i-1))/dchi^2; 
    end    
    r1=-M*(Lami.*w.*bi-Mi*bi+Dchi*bcc);%+Db*K*bi; 
    r((i-1)*n+1:i*n)=r1; 
end
r2=-M*(I.*h-L.*w-Ga*w.*bt); 
r3=-M*(pp-I.*h); 
r=-[r;r2;r3]; 