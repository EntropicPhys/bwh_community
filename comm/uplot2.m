function uplot2(p,wnr) % only plot <B(.,chi)>  over chi 
N=p.N; n=p.np; par=p.u(p.nu+1:end); chimin=par(19); chimax=par(20); 
y=linspace(chimin,chimax,N); vol=p.vol; ba=0*y; 
b=p.u(1:N*n); B=reshape(b,n,N); M=p.mat.M(1:n,1:n); 
for i=1:N
    ba(i)=sum(M*B(:,i))/vol; 
end
figure(wnr); plot(y,ba,'linewidth',2); set(gca,'fontsize',12);  
xlabel('\chi'); ylabel('<B(.,\chi)>'); axis tight; 
tits=[p.file.pname mat2str(p.file.count-1)]; title(tits); 

