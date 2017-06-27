function x=Atmx3d(y,para)
N=para.nx*para.ny;
M=para.nd*para.nv;%para.nd=256;
nbin=para.nbin;
x=zeros(nbin*N,1,'single');
for i=1:nbin
x((i-1)*N+1:i*N)=Atx_fan_mf(y((i-1)*M+1:i*M),para);
end