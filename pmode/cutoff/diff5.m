function [ diff ] = diff5( y,i,h )
%diff5 5 point derivative

diff=(y(i-2)-8*y(i-1)+8*y(i+1)-y(i+2))/(12*h);

end

