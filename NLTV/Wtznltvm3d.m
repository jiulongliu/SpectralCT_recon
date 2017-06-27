function utt=Wtznltvm3d(gradu,para)
% %  [Im0,map] =imread('images/barbara512.tif'); 
% % Im0=double(Im0)/255;
% % u=Im0;
% % gradu=Wxnltvm(u,para);
% %  utt=Wtxnltvm(gradu,para);
% 
% % 
% % function test_nonlocalTV
% % 
% % close all;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Global constants
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% YES = 0;
% NO = 1;
% 
% Im0=para.Iref;
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Comment and un-comment experiments
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%
% % Barbara266
% %%%%%%%%%%%%%%%
% % [Im0,map] = imread('images/barbara512.tif'); 
% f=1.5; % scale for function truesize
% [Ny,Nx] = size(Im0); % image size
% % Im0 = double(Im0); Im0 = Im0/ max(Im0(:)); % normalize image
% % k = 0.05; % noise level
% % ImOrig = Im0; Im0 = Im0 + k*randn(size(Im0)); % add noise
% Min1 = min(Im0(:)); Max1 = max(Im0(:));
% Im0 = ( Im0-Min1 )/ ( Max1-Min1 ); % noisy image
% % ImOrig = ( ImOrig-Min1 )/ ( Max1-Min1 ); % original image
% 
% % parameters for Weight function
% m = 5; % patch size
% w = 11; % window search size
% h = (0.25)^2; % scale parameter for weight function
% NbNeigh = 10; % number of neighbors
% includeCloseNeigh = YES; % include 4 closest neighbors
% includeCloseNeigh = NO;
% VecParametersW = [ Ny; Nx; m; w; h; NbNeigh; includeCloseNeigh; ];
% 
% % parameters for Split-Bregman Nonlocal TV
% mu = 70; % 
% lambdaInv = 0.05; %  
% lambda = 1/lambdaInv;
% nbIters = 4; % nomber of outer iterations
% NbInnerIter = 2; % number of inner iterations
% VecParametersNLTV = [ Ny; Nx; m; w; NbNeigh; lambda; mu;...
%     nbIters; NbInnerIter; includeCloseNeigh;];
% 
% % other parameters
% cpt_fig = 0; % figure number
% 
% 
% 
% 
% 
% % %%%%%%%%%%%%%%%
% % % Barbara512
% % %%%%%%%%%%%%%%%
% % [Im0,map] = imread('images/barbara512.tif'); 
% % f=1.5; % scale for function truesize
% % [Ny,Nx] = size(Im0); % image size
% % Im0 = double(Im0); Im0 = Im0/ max(Im0(:)); % normalize image
% % k = 0.05; % noise level
% % ImOrig = Im0; Im0 = Im0 + k*randn(size(Im0)); % add noise
% % Min1 = min(Im0(:)); Max1 = max(Im0(:));
% % Im0 = ( Im0-Min1 )/ ( Max1-Min1 ); % noisy image
% % ImOrig = ( ImOrig-Min1 )/ ( Max1-Min1 ); % original image
% % 
% % % parameters for Weight function
% % m = 5; % patch size
% % w = 11; % window search size
% % h = (0.25)^2; % scale parameter for weight function
% % NbNeigh = 10; % number of neighbors
% % includeCloseNeigh = YES; % include 4 closest neighbors
% % includeCloseNeigh = NO;
% % VecParametersW = [ Ny; Nx; m; w; h; NbNeigh; includeCloseNeigh; ];
% % 
% % % parameters for Split-Bregman Nonlocal TV
% % mu = 70; % 
% % lambdaInv = 0.05; %  
% % lambda = 1/lambdaInv;
% % nbIters = 4; % nomber of outer iterations
% % NbInnerIter = 2; % number of inner iterations
% % VecParametersNLTV = [ Ny; Nx; m; w; NbNeigh; lambda; mu;...
% %     nbIters; NbInnerIter; includeCloseNeigh;];
% % 
% % % other parameters
% % cpt_fig = 0; % figure number
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % compute NL-Weights
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %file = 'weights.mat';
% [W,Y,SY] = compute_fastNLWeights_mex(single(Im0),...
%     single(VecParametersW));
% %save(file,'W','Y','SY');
% %load(file,'W','Y','SY');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Perform NL-TV
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Initialization
% NbNeigh2 = VecParametersNLTV(5);
% includeCloseNeigh2 = VecParametersNLTV(10);
% if includeCloseNeigh2==YES; i=4; else; i=0; end
% % d = zeros(Nx,Ny,2*(NbNeigh2+i));
% % b = zeros(Nx,Ny,2*(NbNeigh2+i));
% % u = Im0;
% 
% % gradu=Wxnltv(single(u),single(W),int32(Y),int32(SY),single(VecParametersNLTV));

VecParametersNLTV=para.VecParametersNLTV;
nb=para.nb;
utt=zeros(para.ny,para.nx,para.nb,'single');
[ws1,ws2,ws3]=size(para.W(:,:,:,1));
ws3=ws3/2;
for i=1:nb
     W=para.W(:,:,:,i);
Y=para.Y(:,:,:,i);
SY=para.SY(:,:,i);   
utt(:,:,i)=Wtxnltv(single(gradu((i-1)*para.ny*para.nx*ws3+1:(i)*para.ny*para.nx*ws3)),single(W),int32(Y),int32(SY),single(VecParametersNLTV));
end
utt=utt(:);
% utt=utt(:);