function y=Amx3d(x,para)
N=para.nx*para.ny;
M=para.nd*para.nv;%para.nd=256;
nbin=para.nbin;
y=zeros(nbin*M,1,'single');
for i=1:nbin
%     i
%     length(x)
%     length(x((i-1)*N+1:i*N))
y((i-1)*M+1:i*M)=Ax_fan_mf(x((i-1)*N+1:i*N),para);
end