
function [X n_iter d Ind Ind2]=bregman_tv(x0,f0,ip,var)
N_iter=ip.N_iter;cg_tol=ip.cg_tol;cg_iter=uint32(ip.cg_iter);
Wxs=ip.Wxs;Wtxs=ip.Wtxs;Sxs=ip.Sxs;wps=ip.wps;
AmX=var.AmX;AtmX=var.AtmX;JtJ=ip.JtJ;
N=ip.N;
ns=ip.wps.nx*ip.wps.ny;
nb=ip.wps.nt/10;
Min_iter=ip.Min_iter;

Wght=var.var_AtA.Wght;
X=zeros(N,1,'single');
if ~isempty(x0)
    X=x0;
end
if ip.isC==1
    X=complex(X,X);
end

lambda_xs=ip.lambda_xs;
mu_xs=ip.mu_xs;
var.mu_xs=mu_xs;
d_xs=zeros(size(Wxs(X(:,1),wps)),'single');v_xs=d_xs;

% Bregman iteration
f=f0(:);
Xt=X;

n_Ind=2;ind=0;Ind=zeros(n_Ind,1);
n_Ind2=1;ind2=0;Ind2=zeros(n_Ind2,1);
d=zeros(N_iter,1);

for n_iter=1:N_iter
    tic;
    Xi=Xt;
%     muWt=Wtxs(d_xs+v_xs,wps);
%     for i=1:nb
%         muWt(1+(i-1)*ns:ns+(i-1)*ns)=mu_xs(i)*muWt(1+(i-1)*ns:ns+(i-1)*ns);
%     end
%     g=AtmX(Wght.*f,var.var_AtA)+muWt;

    g=AtmX(Wght.*f,var.var_AtA)+mu_xs*Wtxs(d_xs+v_xs,wps);
    
    [X,cg_err,cg_n]=CG033114(X,g,JtJ,var,cg_tol,cg_iter);
    
    disp([num2str(n_iter) ' -- CG: N=' num2str(cg_n) ' error=' num2str(cg_err)]);

    temp_xs=Wxs(X,wps)-v_xs;

    d_xs=Sxs(temp_xs,lambda_xs,wps);

    v_xs=d_xs-temp_xs;

    Xt=X;

 
    if (cg_err>5)
        break;
    end
    d(n_iter)=sum(abs(Xi-Xt))/sum(abs(Xt));
    if(n_iter>Min_iter&&(d(n_iter)>d(n_iter-1)))
        ind=ind+1;Ind(ind)=n_iter;
        disp(['Ind ' num2str(ind) '=' num2str(n_iter)]);
    end
    if(n_iter>Min_iter&&(d(n_iter)>d(n_iter-1)&&d(n_iter-1)>d(n_iter-2)))
        ind2=ind2+1;Ind2(ind2)=n_iter;
        disp(['Ind2 ' num2str(ind2) '=' num2str(n_iter)]);
    end

    if(n_iter>Min_iter&&(d(n_iter)<ip.threshold||ind>=n_Ind||ind2>=n_Ind2))
        break
    end


    toc
end
