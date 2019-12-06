function res = gradf(x,y)
res = exp(-x.^2-y.^4).*[ 3.*x.^2-2.*x.^4 ; -4.*x.^3.*y.^3 ];
end
