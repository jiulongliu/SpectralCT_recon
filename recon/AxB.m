function Y=AxB(X,var)
B=var.B;
nx=var.nx;
ny=var.ny;
nb=var.nb;
% size(X)
% nb
X=reshape(X,[nx*ny,nb])*B;
X=X(:);
Y=Amx3d(X,var);

% AtmX=Atx_cone(X,var);