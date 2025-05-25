function uplot2(p,wnr) % plot B and W
N=p.N; n=p.np; x=getpte(p); mclf(wnr); x1=min(x); x2=max(x); 
par=p.u(p.nu+1:end); chimin=par(19); chimax=par(20); 
y=linspace(chimin,chimax,N); ba=0*y; [X,Y]=meshgrid(y,x); 
p1=[1:4 7:10 13:16 19:22 25:28]; 
p1=[1:4 7:10 13:16 19:22 25:28 31:34]; 
p2=37:40; p4=6:6:28; 
b=p.u(1:N*n); B=reshape(b,n,N); 
w=p.u(N*n+1:(N+1)*n); h=p.u((N+1)*n+1:(N+2)*n); bmax=1.7;
subplot(7,6,p2); plot(x,w,'linewidth',2); axis tight; title('w'); 
set(gca,'fontsize',12); 
try; z1=round(min(w)*105)/100; z2=round(max(w)*95)/100; yticks([z1 z2]); 
catch; try; z1=0.9*min(w); z2=1.1*max(w); axis([x1 x2 z1 z2]); zt=round(50*(z1+z2))/100; 
    yticks(zt); catch; end 
end 
M=p.mat.M(1:n,1:n); vol=p.vol; for i=1:N; ba(i)=sum(M*B(:,i))/vol; end
subplot(7,6,p4); plot(ba,y,'linewidth',2);xlim([0 bmax]); 
set(gca,'fontsize',12);  ylabel('\chi'); xlabel('<B(.,\chi)>'); grid on; 
subplot(7,6, p1); surf(X,Y,B);caxis([0 bmax]); view(90,270);  shading interp; 
axis tight; colormap(flipud(summer)); tits=[p.file.pname mat2str(p.file.count-1)]; 
colorbar('southoutside'); title(tits); set(gca,'fontsize',12); 