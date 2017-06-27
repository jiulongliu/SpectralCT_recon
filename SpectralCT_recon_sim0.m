% *----------------------------------------------
% 	Author Contact Information:
% Jiulong Liu and Hao Gao
% jiulong.liu@gmail.com
% Institute of Natural Sciences, Shanghai Jiao Tong University, China
% If you find this code useful, you may cite the following reference:
% Jiulong Liu, Huanjun Ding, Sabee Molloi, Xiaoqun Zhang, Hao Gao, Ticmr: Total image constrained material reconstruction via nonlocal total variation regularization for spectral ct. IEEE transactions on medical imaging 35.12 (2016): 2578-2586.

addpath([pwd '\recon']);
addpath([pwd '\CT-fp']); 
addpath([pwd '\NLTV']); 
addpath([pwd '\data']); 

load SpectralCTData.mat mx0 nbin pc_sino_bins_ut Borig

y0=zeros(size(pc_sino_bins_ut),'single');
for i=1:nbin
    y0(:,:,i)=single(pc_sino_bins_ut(:,:,i)/mx0(i));
    B(:,i)=Borig(:,i)/mx0(i);
end


% sctmea=sum(exp(-pc_sino_bins_ut(:,:,1:5)),3);% figure;imshow(sctmea(1:10,1:20),[]);
% sino_3dtat=-log(sctmea/mean(mean(sctmea(1:10,1:20))));
sino_3dtat=sum(pc_sino_bins_ut(:,:,1:5),3);






GPU=1; % 1 for CUDA GPU, 0 for CPU
%% TV Reconstruction for total image Xref
mu=6.0;
lambda=0.15;
Xref=ReconTV(single(sino_3dtat),mu,lambda,GPU);
figure;imshow([Xref],[]);



%% TV Reconstruction for Z


mu=1;
lambda=0.08;
Ztv=ReconZTV(y0,B,mu,lambda,GPU);
figure;imshow([Ztv(:,:,1) Ztv(:,:,2) Ztv(:,:,3)],[0 1]);



%% Nonlocal TV Reconstruction for Z

mu=0.18;
lambda=0.1;
[nb,~]=size(B);
Zref=repmat(Xref,[1,1,nb]);
Znltv=ReconZNLTV(y0,Zref,B,mu,lambda,GPU);
figure;imshow([Znltv(:,:,1) Znltv(:,:,2) Znltv(:,:,3)],[0 1]);