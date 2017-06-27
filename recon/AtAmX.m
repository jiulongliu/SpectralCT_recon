function Y=AtAmX(X,AmX,AtmX,var)
% var
L=2*length(X);
if L==length(AmX(X,var))
    Wght=ones(L,1,'single');
elseif 10*L==length(AmX(X,var))
    Wght=ones(10*L,1,'single');
else
    Wght=var.Wght;
end
% size(Wght)
% figure;plot(AmX(X,var));
% size(AmX(X,var))
Y=AtmX(Wght.*AmX(X,var),var);
