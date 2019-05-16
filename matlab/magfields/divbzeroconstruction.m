m=.05;  %m=pi.*a.^2.*I/c ring current I radius a
len=0.1;
b0=1.0;
wid=1.0;
shift=0.5;

x=-1.0:0.05:1.0;
y=0.1:0.1:2.0;

[X1,Y1]=meshgrid(x-shift,y);
[X2,Y2]=meshgrid(x+shift,y);


%bx1=-0.1*(X1./sqrt(X1.^2))-2.0.*b0*X1.*exp(-wid.*(X1.^2)./len);
bx1=-2.0.*b0*X1.*exp(-wid.*(X1.^2)./len);
by1=+b0.*exp(-wid.*(X1.^2)./len);

%bx2=0.1*(X2./sqrt(X2.^2))+2.0.*b0*X2.*exp(-wid.*(X2.^2)./len);
bx2=+2.0.*b0*X2.*exp(-wid.*(X2.^2)./len);
by2=-b0.*exp(-wid.*(X2.^2)./len);


bx=bx1+bx2;
by=by1+by2;


bmag=sqrt(bx.^2+by.^2);
contour(X,Y,bmag,[0.005 0.006 0.007 0.008 0.009 0.01 0.05 0.1 0.15 0.2 0.25 0.3 .35 .4 .45 2]);
hold on
quiver(X(8:20,:),Y(8:20,:),bx(8:20,:),by(8:20,:),3);
hold off
