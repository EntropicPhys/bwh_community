function r=sG(p,u) % RHS of bwh-community
N=p.N; n=p.np; b=reshape(u(1:N*n),n,N);     % field assign. 
w=u(N*n+1:(N+1)*n); h=u((N+1)*n+1:(N+2)*n);  
par=u(p.nu+1:end);pp=par(1);Lam0=par(2);Ga=par(3); % syst. params
A=par(4);R=par(5);L0=par(6);f=par(7);Q=par(8);
Kmin=par(9);Kmax=par(10); Mmin=par(11); Mmax=par(12); 
Ymin=par(13);Ymax=par(14);Db=par(15);Dw=par(16); 
Dh=par(17);Dchi=par(18);chimin=par(19);chimax=par(20); 
K=p.mat.K; M=p.mat.M(1:n,1:n); % Stiffness and Mass matrices 
ov=ones(n,1);bt=zeros(n,1); btt=bt;% b-tilde and b-tilde-tilde
chii=linspace(chimin,chimax,N);dchi=chii(2)-chii(1);  
for i=1:N % fill bt and btt 
    bt=bt+b(:,i); 
    Yi=Ymax+chii(i)*(Ymin-Ymax); 
    btt=btt+b(:,i)*Yi; 
end 
I=A*(btt+f*Q)./(btt+Q); L=ov*L0./(1+R*bt); % infil. and evapo. funcs 
for i=1:N    %loop that assemble the finite-difference scheme
    Ki=Kmax+chii(i)*(Kmin-Kmax);  
    Mi=Mmax+chii(i)*(Mmin-Mmax); 
    Lami=Lam0*Ki./(bt+Ki); 
    bi=b(:,i);    
    switch i  % chi-diffusion terms, i=1 and i=N with Neumann BCs 
        case 1; bcc=(b(:,2)-2*b(:,1))/dchi^2;
        case N;  bcc=(-2*b(:,N)+b(:,N-1))/dchi^2; 
        otherwise; bcc=(b(:,i+1)-2*b(:,i)+b(:,i-1))/dchi^2; 
    end    
    r1=-M*(Lami.*w.*bi-Mi*bi+Dchi*bcc)+Db*K*bi; 
    r((i-1)*n+1:i*n)=r1; 
end
r2=-M*(I.*h-L.*w-Ga.*bt.*w)+Dw*K*w;
r3=-M*(pp-I.*h)+Dh*K*h;   
r=[r;r2;r3]; 