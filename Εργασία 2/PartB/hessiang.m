function res = hessiang(x,y)
res = [ 12.*x.^2+0.8*pi^2.*sin(2*pi.*x) , 0 ;
           0    , 2+1.2*pi^2.*cos(2*pi.*y) ];
end