function Z=ReconZNLTV(y0,x0,B,mu0,lambda,GPU)
[nd,nv0,nbin]=size(y0);
[nb,nbin2]=size(B);

nx=512;ny=512;n=nx*ny;
dx=0.02;
US=1;
SO=965;OD=1055-965;
dy_det=0.02;
%dx=dx*SO/(SO+OD);
%dx=0.1;

sd_phi=single((0:(nv0-1))*(2*pi/nv0));
V=single(US:nv0);
sd_phi=sd_phi(V);nv=numel(V);

y_os=single(0*dy_det);

y_det=single(((-nd/2:nd/2-1)+0.5)*dy_det)+y_os;
scale=dx;
SO=SO/scale;OD=OD/scale;y_det=y_det/scale;dy_det=dy_det/scale;y_os=y_os/scale;

if 1
    % single frame %
    nt=1;
    Nv(1)=nv;
    id_X=0*ones(1,nv);
    tmp_size=nv;
    id_Y=zeros(tmp_size,nt);
    id_Y(1:Nv(1),1)=uint32(0:nv-1);
else
    % multiple frames %
    nt=300;
    dt=1/(nt/2);
    I=sin(pi/2*(1:-dt:0));
    I=unique(sort([I -I]));
%     dI=2/(nt-1);
%     I=[-1 -1+0.5*dI:dI:1-0.5*dI 1];
    
    id_X=-1*ones(1,nv);
    Nv=zeros(1,nt);
    id_Y=zeros(nv,nt);
    for i=1:nt
        tmp=intersect(find(s>=I(i)),find(s<=I(i+1)));
        id_X(tmp)=uint32(i-1);
        Nv(i)=numel(tmp);
        id_Y(1:Nv(i),i)=uint32(tmp-1);
    end
    tmp_size=max(Nv);
    id_Y=id_Y(1:tmp_size,:);
end
%%%%%%%%%%%%%%%
% nbin=5;

para=struct('version',uint32(1),'GPU',uint32(GPU),'SO',single(SO),'OD',single(OD),'scale',single(scale),'nx',uint32(nx),'ny',uint32(ny),'nv',uint32(nv),...
    'sd_phi',sd_phi,'y_det',y_det,'id_X',uint32(id_X),'nt',uint32(nt),'dy_det',dy_det,'y_os',y_os/dy_det,...
    'id_Y',uint32(id_Y),'Nv',uint32(Nv),'tmp_size',uint32(tmp_size),...
    'cos_phi',[],'sin_phi',[],'cos_det',[],'sin_det',[],'nbin',uint32(nbin),'nb',uint32(nb),'B',single(B),'nd',uint32(nd));
para.cos_phi=single(cos(sd_phi));para.sin_phi=single(sin(sd_phi));
angle_det=atan2(y_det,SO+OD);para.cos_det=single(cos(angle_det));para.sin_det=single(sin(angle_det));
Wght=exp(-y0(:));
Wght=Wght/mean(Wght);
Wght(Wght<0.2)=0.2;
para.Wght=single(Wght);

% tic;maxAtA=maxAtA_fan_mf_v2(para);toc;
maxAtA=single(0.1402);

N_iter=18;
nplot=1;
var_plot=struct('nx',nx,'ny',ny,'nt',nt);
plotx=@plotx_2d_sf;
% AmX=@Ax_fan_mf_v2;
% AtmX=@Atx_fan_mf_v2;
% AmX=@Amx3d;
% AtmX=@Atmx3d;
AmX=@AxB;
AtmX=@AtxB;
% lambda=lambda*maxAtA;
isC=0;
para.isC=uint32(isC);


varW1=computeweight(single(reshape(x0(:,:,1),[nx ny])));
varW2=computeweight(single(reshape(x0(:,:,2),[nx ny])));
varW3=computeweight(single(reshape(x0(:,:,3),[nx ny])));
wps=struct('isC',uint32(0),'nx',uint32(nx),'ny',uint32(ny),'nt',uint32(nb*10),'nb',uint32(nb),'B',single(B));

wps.W=zeros([size(varW1.W),nb],'single');
wps.Y=zeros([size(varW1.Y),nb],'single');
wps.SY=zeros([size(varW1.SY),nb],'single');
wps.W(:,:,:,1)=varW1.W;
wps.Y(:,:,:,1)=varW1.Y;
wps.SY(:,:,1)=varW1.SY;
wps.W(:,:,:,2)=varW2.W;
wps.Y(:,:,:,2)=varW2.Y;
wps.SY(:,:,2)=varW2.SY;
wps.W(:,:,:,3)=varW3.W;
wps.Y(:,:,:,3)=varW3.Y;
wps.SY(:,:,3)=varW3.SY;
wps.VecParametersNLTV=varW1.VecParametersNLTV;

Wxs=@Wznltvm3d;Wtxs=@Wtznltvm3d;


% lambda
% lambda=0.06;
mu=mu0*maxAtA;
N_iter=20;
cg_tol=1e-4;
cg_iter=10;
nplot=5;

% wps=struct('isC',uint32(0),'nx',uint32(nx),'ny',uint32(ny),'nt',uint32(nbin),'nbin',uint32(nbin));
% % wps.W=varW.W;
% % wps.Y=varW.Y;
% % wps.SY=varW.SY;
% % wps.VecParametersNLTV=varW.VecParametersNLTV;

% Wxs=@Wxnltvm;Wtxs=@Wtxnltvm;
% Wxs=@wx_tv3d;Wtxs=@wtx_tv3d;
Sxs=@shrink_tv3d;
JtJ='JtJ_wtw';
plotx=@plotx_3d_dynamic;
var_plot=struct('nx',nx,'ny',ny,'nt',nt);

ip=struct('lambda_xs',lambda,'mu_xs',mu,'isC',uint32(0),...
'Min_iter',20,'threshold',1e-3,'N_iter',N_iter,'N',nx*ny*nb,'nplot',nplot,'plotx',plotx,'var_plot',var_plot,...
'cg_iter',20,'cg_tol',1e-4,'JtJ',JtJ,'Wxs',Wxs,'Wtxs',Wtxs,'wps',wps,'Sxs',Sxs);
var_CG=struct('id_x12',uint32(0),'id_xs',uint32(1),'id_xr',uint32(0),...
'AmX',AmX,'AtmX',AtmX,'var_AtA',para,'Wxs',Wxs,'Wtxs',Wtxs,'wps',wps,'isC',uint32(0));
tic;[Z n_iter d Ind Ind2]=bregman_tv([],y0,ip,var_CG);toc;  
Z=reshape(Z,[nx,ny,nb]);
% save('Xrec0.mat','X');
% Xr=[Xr X];
% Xr=Xr(:,1:3*256);
% save('recX.mat','Xr');
% imshow(X,[]);
% % Xbin4=X;
% Xbin4=X;
% save('Xbin4n.mat','Xbin4');
% Xr=Xr(:,1:3*256);