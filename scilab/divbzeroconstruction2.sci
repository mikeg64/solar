

exec('../vtk/vacscalar2vtk3d.sce',-1)

m=.05;  //m=%pi.*a.^2.*I/c ring current I radius a
len=0.1;
b0=1.0;  //1kG
bz0=15/1000; //G
wid=1.1;
shift=0.0;
Z0=3.1

x=0.1:0.1:2.0;
y=0.1:0.1:3.0;
z=0.1:0.1:2;


[X1,Y1]=meshgrid(x-shift,z);
[X2,Y2]=meshgrid(x+shift,z);
X2=0.5*X2;
Y2=Z0-1.*Y2;
g=(-exp(-Y2)-Y2);
//R=sqrt(X2.^2+Y2.^2);
R=X2;
f=R.^2.*(-wid.*exp(-Y2)-Y2);
br=0.5*b0.*R.*exp(f);
bz=b0.*exp(f)+bz0

//br=-2.0.*b0*X1.*exp(-wid.*(X1.^2)./len);
//bx1=-0.1*(X1./sqrt(X1.^2))-2.0.*b0*X1.*exp(-wid.*(X1.^2)./len);
//bx1=-2.0.*b0*X1.*exp(-wid.*(X1.^2)./len);
//by1=+b0.*exp(-wid.*(X1.^2)./len);

//bx2=0.1*(X2./sqrt(X2.^2))+2.0.*b0*X2.*exp(-wid.*(X2.^2)./len);
//bx2=+2.0.*b0*X2.*exp(-wid.*(X2.^2)./len);
//by2=-b0.*exp(-wid.*(X2.^2)./len);


//bx=bx1+bx2;
//by=by1+by2;


bmag=sqrt(br.^2+bz.^2);
//contour(x,y,bmag',[0.005 0.006 0.007 0.008 0.009 0.01 0.05 0.1 0.15 0.2 0.25 0.3 .35 .4 .45 2]);
contour(x,z,bmag',7);
//hold on
//quiver(X(8:20,:),Y(8:20,:),bx(8:20,:),by(8:20,:),3);
//hold off

//vacscalar2vtk3d( pict,wd,xx,yy,zz,field,vecsize,filename,directory )
field=0
vecsize=0
nbmag=zeros(20,20,20)
nbmag(:,:,1)=bmag;
vacscalar2vtk3d( 1,nbmag,x',y',z',field,vecsize,'tests','.' )
