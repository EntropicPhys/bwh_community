function p=oosetfemops(p) % set FEM operators, homog. Neuman BC by default
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); % FEM matrices
p.mat.K=K; p.mat.M=kron(diag([1,1,1]),M); % scalar Lapl., full M
