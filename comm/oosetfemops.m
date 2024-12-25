 function p=oosetfemops(p) % in problem-dir, since highly problem dependent
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); 
p.mat.K=K; p.mat.M=kron(diag(ones(1,p.N+2)),M); % scalar Lapl., full M