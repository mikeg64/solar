xy = -2.5 + 5*gallery('uniformdata',[200 2],0);
x = xy(:,1);
y = xy(:,2);
v = x.*exp(-x.^2-y.^2);
% [xq,yq] = meshgrid(-2:.2:2, -2:.2:2);
[xq,yq] = meshgrid(-2:.2:2, -2);
vq = griddata(x,y,v,xq,yq);


% xy = height;
% 
% x = xy(:,1);
% y = 1:0.1:1;
% v = temp;
% % [xq,yq] = meshgrid(-2:.2:2, -2:.2:2);
% [xq,yq] = meshgrid(-2:.2:2, -2);
% vq = griddata(x,y,v,xq,yq);
