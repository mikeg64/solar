m=.05;  %m=pi.*a.^2.*I/c ring current I radius a
len=1.0;
b0=1.0;

x=0:0.05:1.0;
y=0.1:0.1:2.0;

[X,Y]=meshgrid(x,y);






bx=b0*cos(pi*X./len).*exp(-pi*Y./len);
by=-b0*sin(pi*X./len).*exp(-pi*Y./len);

bmag=sqrt(bx.^2+by.^2);
contour(X,Y,bmag,[0.005 0.006 0.007 0.008 0.009 0.01 0.05 0.1 0.15 0.2 0.25 0.3 .35 .4 .45 2]);
hold on
quiver(X(8:20,:),Y(8:20,:),bx(8:20,:),by(8:20,:),3);
hold off
