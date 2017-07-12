function [res]=inte(f,dx)

nsiz=size(f);

if nsiz(1)>nsiz(2)
    nel=nsiz(1);
else
    nel=nsiz(2);
end

res=0.0;
if (nel > 1) 
    if (nel == 2) 
        res=dx*0.5*(f(1)+f(0));
    end

    if (nel > 2)
      for k=1:nel-1 
          res=res+0.5*(f(k-1)+f(k))*dx;
      end
    end
end %outer

%return,res
end


