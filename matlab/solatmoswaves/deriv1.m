function [res]=deriv1(f,x)
nel=size(f);
nel1=size(x);

if nel(1)>nel(2)
    num=nel(1);
else
    num=nel(2);
end

% if (nel ne nel1) then begin
%  print,'Inconsistant input, stop.'
%  stop
% endif

res=zeros(num);
for i=2:num-3 
    res(i)=(1./12.0./(x(i+1)-x(i)))*(8.0*f(i+1)-8.0*f(i-1)-f(i+2)+f(i-2));
end
%for i=1,nel-2 do res(i)=(1.d0/2.d0/(x(i+1)-x(i)))*(f(i+1)-f(i-1))
res(0)=res(2);
res(1)=res(2);
res(num-1)=res(num-3);
res(num-2)=res(num-3);

end



