function [res]=par5(x,x0,A)
    res=A.*exp(-1.0*((x+abs(min(x)))/x0)^0.95);
end



