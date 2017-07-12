function [res]=par3(x,x0,A)

    A=A/(x0.^2);
    
  

if (x <= -2*x0)
    res=0;
elseif ((x >= -2*x0) & (x <= -x0)) 
    res=A*(x+3*x0)*(x+x0)+A*x0.^2;
elseif ((x >= -x0) & (x <= x0)) 
    res=-A*(x+x0)*(x-x0)+A*x0.^2;
elseif ((x >= x0) & (x <= 2*x0)) 
    res=A*(x-3*x0)*(x-x0)+A*x0.^2;
elseif (x >= 2*x0) 
    res=0;
end

%return,res
end


