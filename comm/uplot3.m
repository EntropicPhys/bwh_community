function uplot3(p,wnr) % also plot al(chi) for B_t=al(chi)*B+D_chi*B''
N=p.N;  par=p.u(p.nu+1:end); chimin=par(19); chimax=par(20); delchi=(chimax-chimin)/N; 
ap=par(21); bp=par(22); 
n=p.np; y=linspace(chimin,chimax,N); ba=0*y; 
w=p.u(N*n+1); % was wrong! 
b=p.u(1:n*N); b=reshape(b,n,N); M=p.mat.M(1:n,1:n); 
bt=0; al=zeros(1,N); par=p.u(p.nu+1:end); Lam0=par(2); 
Kmin=par(9); Kmax=par(10); Mmin=par(11); Mmax=par(12); 
for i=1:N
    %ba(i)=sum(M*b(:,i)); 
    ba(i)=b(1,i); 
    bt=bt+ba(i)*delchi;   
end
%ppp=[chimin,chimax,delchi,w,bt,Kmin,Kmax,Mmin,Mmax]; ppp
for i=1:N; 
     chii=chimin+i*delchi;      
     Ki(i)=Kmax+chii^ap*(Kmin-Kmax);  
     Mi(i)=Mmax+chii^bp*(Mmin-Mmax); 
     Lami(i)=Lam0*Ki(i)/(bt+Ki(i)); 
     al(i)=Lami(i)*w-Mi(i);     
end 
figure(wnr); plot(y,ba,'-*','linewidth',2); set(gca,'fontsize',12);  
xlabel('\chi'); ylabel('<B(.,\chi)>'); axis tight; 
tits=[p.file.pname mat2str(p.file.count-1)]; title(tits); 

figure(11); clf; plot(y,al,'linewidth',2); set(gca,'fontsize',12);  axis tight; 
xlabel('\chi'); ylabel('\alpha'); title(''); grid on

figure(12); clf; plot(y,Ki,'linewidth',2); set(gca,'fontsize',12);  axis tight; 
xlabel('\chi'); ylabel('K'); title(''); grid on

figure(13); clf; plot(y,Mi,'linewidth',2); set(gca,'fontsize',12);  axis tight; 
xlabel('\chi'); ylabel('M'); title(''); grid on
figure(14); clf; plot(y,Lami*w,'linewidth',2); set(gca,'fontsize',12);  axis tight; 
xlabel('\chi'); ylabel('Lam'); title(''); grid on
bt

