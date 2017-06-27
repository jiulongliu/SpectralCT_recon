function Y=AtxB(X,var)
B=var.B;
nx=var.nx;
ny=var.ny;
nbin=var.nbin;
Y=reshape(Atmx3d(X,var),[nx*ny,nbin])*B';

Y=Y(:);
% AtmX=Atx_cone(X,var);