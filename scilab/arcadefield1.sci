


b0=1;  //m=pi.*a.^2.*I/c ring current I radius a
l=1;

x=-1.0:0.1:1.0;
y=0.1:0.1:2.0;

[X,Y]=meshgrid(x,y);



bx=b0.*cos(%pi*X./l).*exp(-%pi.*Y./l);
by=-b0.*sin(%pi*X./l).*exp(-%pi.*Y./l);


bmag=sqrt(bx.^2+by.^2);
contour(x,y,bmag',[0.005 0.006 0.007 0.008 0.009 0.01 0.05 0.1 0.15 0.2 0.25 0.3 .35 .4 .45 2]);
//contour(x,y,bmag',5);

//hold on
//quiver(X(8:20,:),Y(8:20,:),bx(8:20,:),by(8:20,:),3);
//hold off

