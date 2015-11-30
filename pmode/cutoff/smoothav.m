function [ spres, srho ] = smoothav( consts, pres, rho, nsmoothsteps )
%smoothav smooth an array of pressure and density values
%   Detailed explanation goes here

[nr,nc]=size(pres);

spres=zeros(nr,nc);
srho=zeros(nr,nc);

for j=1:nsmoothsteps
   pres(nr+j)=pres(nr);
   rho(nr+j)=rho(nr);
end

for i=1:nr
    ptot=pres(i);
    rhotot=rho(i);
    for j=1:nsmoothsteps-1
        ptot=ptot+pres(i+j);
        rhotot=rhotot+rho(i+j);      
    end
    spres(i)=ptot/nsmoothsteps;
    srho(i)=rhotot/nsmoothsteps;
end
    
    


end

