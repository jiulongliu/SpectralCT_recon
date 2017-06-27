%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: 	Xavier Bresson (xbresson@math.ucla.edu)
% Description: Nonlocal TV Minimization for Image Denoising
% For more details, see the report:
% X. Bresson, "A Short Note for Nonlocal TV Minimization"
% See also these reports: 
% X. Zhang, M. Burger, X. Bresson, and S. Osher
% "Bregmanized Nonlocal Regularization for Deconvolution and Sparse Reconstruction"
% CAM Report 09-03, 2009
% Last version: April 28, 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% 
% function test_nonlocalTV

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

YES = 0;
NO = 1;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comment and un-comment experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%
% Barbara266
%%%%%%%%%%%%%%%
[Im0,map] = imread('images/barbara266.tif'); 
f=1.5; % scale for function truesize
[Ny,Nx] = size(Im0); % image size
Im0 = double(Im0); Im0 = Im0/ max(Im0(:)); % normalize image
k = 0.05; % noise level
ImOrig = Im0; Im0 = Im0 + k*randn(size(Im0)); % add noise
Min1 = min(Im0(:)); Max1 = max(Im0(:));
Im0 = ( Im0-Min1 )/ ( Max1-Min1 ); % noisy image
ImOrig = ( ImOrig-Min1 )/ ( Max1-Min1 ); % original image

% parameters for Weight function
m = 5; % patch size
w = 11; % window search size
h = (0.25)^2; % scale parameter for weight function
NbNeigh = 10; % number of neighbors
includeCloseNeigh = YES; % include 4 closest neighbors
includeCloseNeigh = NO;
VecParametersW = [ Ny; Nx; m; w; h; NbNeigh; includeCloseNeigh; ];

% parameters for Split-Bregman Nonlocal TV
mu = 70; % 
lambdaInv = 0.05; %  
lambda = 1/lambdaInv;
nbIters = 4; % nomber of outer iterations
NbInnerIter = 2; % number of inner iterations
VecParametersNLTV = [ Ny; Nx; m; w; NbNeigh; lambda; mu;...
    nbIters; NbInnerIter; includeCloseNeigh;];

% other parameters
cpt_fig = 0; % figure number





% %%%%%%%%%%%%%%%
% % Barbara512
% %%%%%%%%%%%%%%%
% [Im0,map] = imread('images/barbara512.tif'); 
% f=1.5; % scale for function truesize
% [Ny,Nx] = size(Im0); % image size
% Im0 = double(Im0); Im0 = Im0/ max(Im0(:)); % normalize image
% k = 0.05; % noise level
% ImOrig = Im0; Im0 = Im0 + k*randn(size(Im0)); % add noise
% Min1 = min(Im0(:)); Max1 = max(Im0(:));
% Im0 = ( Im0-Min1 )/ ( Max1-Min1 ); % noisy image
% ImOrig = ( ImOrig-Min1 )/ ( Max1-Min1 ); % original image
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








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute NL-Weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%file = 'weights.mat';
[W,Y,SY] = compute_fastNLWeights_mex(single(Im0),...
    single(VecParametersW));
%save(file,'W','Y','SY');
%load(file,'W','Y','SY');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform NL-TV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
NbNeigh2 = VecParametersNLTV(5);
includeCloseNeigh2 = VecParametersNLTV(10);
if includeCloseNeigh2==YES; i=4; else; i=0; end
d = zeros(Nx,Ny,2*(NbNeigh2+i));
b = zeros(Nx,Ny,2*(NbNeigh2+i));
u = Im0;

gradu=Wxnltv(single(u),single(W),int32(Y),int32(SY),single(VecParametersNLTV));

utt=Wtxnltv(single(gradu),single(W),int32(Y),int32(SY),single(VecParametersNLTV));

% gradu=reshape(gradu,[512 512 2*(NbNeigh2+i)]);
% imshow(gradu(:,:,5),[]);

% utt=reshape(utt,[512 512]);
% imshow(utt,[]);
[u_new,d_new,b_new,Temp] = SBNLTV_mex(single(u),single(d),single(b),...
    single(Im0),single(W),int32(Y),int32(SY),single(VecParametersNLTV));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var_u = var(u_new(:));
var_u_uRef = var(u_new(:)-ImOrig(:));
SNR = 10*log10(var_u/var_u_uRef)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cpt_fig = cpt_fig + 1;
figure(cpt_fig); clf;
subplot(1,2,1); imagesc(Im0,[0 1]); colormap(gray); title('Im0'); hold on; %colorbar;
subplot(1,2,2); imagesc(u_new,[0 1]); colormap(gray); title(['u, SNR= ',num2str(SNR)]); hold on; %colorbar;
truesize(cpt_fig,[round(f*Ny) round(f*Nx)]);


















