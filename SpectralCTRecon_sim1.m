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


load SpectralCTData.mat
y0=zeros(size(pc_sino_bins_ut),'single');
for i=1:nbin
    y0(:,:,i)=single(pc_sino_bins_ut(:,:,i)/mx0(i));
    B(:,i)=Borig(:,i)/mx0(i);
end

[nd,nv,~]=size(pc_sino_bins_ut);
[nb,~]=size(B);
% sctmea=sum(exp(-pc_sino_bins_ut(:,:,1:5)),3);% figure;imshow(sctmea(1:10,1:20),[]);
% sino_3dtat=-log(sctmea/mean(mean(sctmea(1:10,1:20))));
sino_3dtat=sum(pc_sino_bins_ut(:,:,1:5),3);
y0b=reshape(y0,[nd*nv,nbin])*B'/(B*B');
y0b=reshape(y0b,[nd,nv,nb]);


GPU=1; % 1 for CUDA GPU, 0 for CPU
%% TV Reconstruction for total image Xref
mu=6.5;
lambda=0.15;
Xref=ReconTV(single(sino_3dtat),mu,lambda,GPU);
figure;imshow([Xref],[]);





%% TV Reconstruction Z

 for mu=[10 30 50]
    for lambda=[0.05 0.15 0.25]
        Ztv1=ReconTV(y0b(:,:,1),mu,lambda,GPU);                
        method=['base1tv_'  num2str(mu) '_'  num2str(lambda)];
        imwrite(Ztv1,['output/', method '.png']);
    end
end
 for mu=[10 30 50]
    for lambda=[0.05 0.15 0.25]
        Ztv2=ReconTV(y0b(:,:,2),mu,lambda,GPU);                
        method=['base2tv_'  num2str(mu) '_'  num2str(lambda)];
        imwrite(Ztv2,['output/', method '.png']);
    end
end
 for mu=[10 30 50]
    for lambda=[0.05 0.15 0.25]
        Ztv3=ReconTV(y0b(:,:,3),mu,lambda,GPU);                
        method=['base3tv_'  num2str(mu) '_'  num2str(lambda)];
        imwrite(Ztv3,['output/', method '.png']);
    end
end
%% Nonlocal TV Reconstruction Z

 for mu=[4 8 12 ]
    for lambda=[0.1 0.2 0.3]
        Znltv1=ReconNLTV(y0b(:,:,1),Xref,mu,lambda,GPU);
        method=['base1nltv_'  num2str(mu) '_'  num2str(lambda)];
        imwrite(Znltv1,['output/', method '.png']);
    end
 end
 for mu=[4 8 12 ]
    for lambda=[0.1 0.2 0.3]
        Znltv2=ReconNLTV(y0b(:,:,2),Xref,mu,lambda,GPU);
        method=['base2nltv_'  num2str(mu) '_'  num2str(lambda)];
        imwrite(Znltv2,['output/', method '.png']);
    end
  end
 for mu=[4 8 12 ]
    for lambda=[0.1 0.2 0.3]
        Znltv3=ReconNLTV(y0b(:,:,3),Xref,mu,lambda,GPU);
        method=['base3nltv_'  num2str(mu) '_'  num2str(lambda)];
        imwrite(Znltv3,['output/', method '.png']);
    end
 end

