

//lagrange interpolation

//xval is the value of x at which we wish to interpolate
//f is the series of values
function y=lagrange_interp(xval,f,x,i)  
    t1=(xval-x(i))*(xval-x(i+1))/((x(i-1)-x(i))*(x(i-1)-x(i+1)));
    t2=(xval-x(i-1))*(xval-x(i+1))/((x(i)-x(i-1))*(x(i)-x(i+1)));
    t3=(xval-x(i-1))*(xval-x(i))/((x(i+1)-x(i-1))*(x(i+1)-x(i)));
    y=t1*f(i-1)+t2*f(i)+t3*f(i+1);
endfunction

//3 and 5 point functions to compute derivatives
//need to use extrapolation to determine derivative at boundary


function diff5p=diff5p(y,i,h)
    diff5p=(y(i-2)-8*y(i-1)+8*y(i+1)-y(i+2))/(12*h);
endfunction

function diff3p=diff3p(y,i,h)
    diff3p=(y(i+1)-y(i-1))/(2*h);
endfunction


function diffdd5p=diffdd5p(y,i,h)
    diffdd5p=(-y(i-2)+16*y(i-1)+16*y(i+1)-y(i+2)-30*y(i))/(12*h^2);
endfunction

function diffdd3p=diffdd3p(y,i,h)
    diffdd3p=(y(i+1)+y(i-1)-2*y(i))/(h^2);
endfunction


function vcs=cs (vars,p,rho)
    vcs=sqrt(vars.fgamma*p/rho);
endfunction

//pressure scale height
function vlambda=lambda (vars,p,rho)
    vlambda=-p/(vars.g*rho);
endfunction

//z, p and rho are arrays of values


function [vomega,lam0,vcs]=omegaiso(z,vars,p,rho)
    
    sz=size(p);
    nrows=sz(1);
    h=z(2)-z(1);
    
    for i=1:nrows
      vomega(i)=0;
      vlam0dash(i)=0;
      vcs(i)=cs(vars,p(i),rho(i));
      lam0(i)=lambda(vars,p(i),rho(i));
    end
    for i=3:nrows-2
     ;// vlam0dash(i)=diff5p(lam0,i,h);
    end
 
    for i=1:nrows

      vomega(i)=sqrt(1+2*vlam0dash(i))*(vcs(i)/(2*lam0(i)));
    end    
endfunction

function [vomega,lam0,vcs]=omega(z,vars,p,rho)
    
    sz=size(p);
    nrows=sz(1);
    h=z(2)-z(1);
    
    for i=1:nrows
      vomega(i)=0;
      vlam0dash(i)=0;
      vcs(i)=cs(vars,p(i),rho(i));
      lam0(i)=lambda(vars,p(i),rho(i));
    end
    for i=3:nrows-2
      vlam0dash(i)=diff5p(lam0,i,h);
    end
 
    for i=1:nrows

     // vomega(i)=sqrt(1+2*vlam0dash(i))*(vcs(i)/(2*lam0(i)));
     vomega(i)=(vars.fgamma*vars.g/(4*%pi*vcs(i)))*sqrt(1+2*vlam0dash(i));
    end    
endfunction


function [vomega,lam0,vcs]=omega1(z,vars,p,rho)
    
    sz=size(p);
    nrows=sz(1);
    h=z(2)-z(1); 
    
    for i=1:nrows
      vomega(i)=0;
      vlam0dash(i)=0;
      vcs(i)=cs(vars,p(i),rho(i));
       vcssq(i)=vcs(i).*vcs(i);
      lam0(i)=lambda(vars,p(i),rho(i));
    end
    
    for i=3:nrows-2
      vlam0dash(i)=diff5p(vcssq,i,h);
    end
 
    for i=1:nrows
      om0(i)=-(vars.g*vars.fgamma)./(2.*vcs(i))    
      vomega(i)=sqrt(om0(i).*om0(i)+(om0(i).*vlam0dash(i)./vcs(i)));
    end    
endfunction

function [vomega,lam0,vcs]=omega2(z,vars,p,rho)
    
    sz=size(p);
    nrows=sz(1);
    h=z(2)-z(1);
    
    for i=1:nrows
      vomega(i)=0;
      vlam0dash(i)=0;
      vcs(i)=cs(vars,p(i),rho(i));
      lam0(i)=lambda(vars,p(i),rho(i));
    end
    for i=3:nrows-2
      vlam0dash(i)=diff5p(lam0,i,h);
    end
 
    for i=1:nrows

      vomega(i)=sqrt(1+2*vlam0dash(i))*(vcs(i)/(2*lam0(i)));
    end    
endfunction

function ombv=brunt(z,vars,p,rho)

    for i=1:nrows
      logp(i)=0;
      logrho(i)=0;
      ombv(i)=0;
    end
    
    for i=1:nrows
      logp(i)=log(p(i));
      logrho(i)=log(rho(i));
    end
    for i=3:nrows-2      
      dlogp=diff5p(logp,i,h);
      dlogrho=diff5p(logrho,i,h);
      ombv(i)=sqrt(((vars.g)*((dlogp/vars.fgamma))-dlogrho));
    end
  
endfunction

function vbfreq=bfreq(z,vars,p,rho,bv,vcs)

    for i=3:nrows-2
      vbfreq(i)=0;
    end

    for i=3:nrows-2
        k0=vars.vomega/vcs(i);
        vbfreq(i)=1.0/sqrt((k0^2/(k0^2+(1)^2))*(bv(i).*bv(i)));
    end

  
endfunction

function vfshift=fshift(vars,vomega,p,rho)    
    vfshift=1-vars.fgamma*vars.g*vars.g*rho/(8*vomega*vomega*p);
endfunction
