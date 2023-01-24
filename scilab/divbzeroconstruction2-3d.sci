

exec('../vtk/vacscalar2vtk3d.sce',-1)

m=.05;  //m=%pi.*a.^2.*I/c ring current I radius a
len=0.1;
b0=1.0;  //1kG
bz0=15/1000; //G
wid=1.1;
shift=0.0;
Z0=3.1
bmag3=zeros(size(x'),size(y'),size(z'));
bx3=zeros(size(x'),size(y'),size(z'));
by3=zeros(size(x'),size(y'),size(z'));
bz3=zeros(size(x'),size(y'),size(z'));
x=0.1:0.1:2.0;
y=0.1:0.1:3.0;
z=0.1:0.1:2;

x=x-shift;
y=y-shift;

for i=1:size(x')
    for j=1:size(y')
        for k=1:size(z')
            r=sqrt(x(i)^2+y(j)^2);
            g=(-exp(-z(k))-z(k));
            f=r^2.*(-wid.*exp(-z(k))-z(k));
            br=0.5*b0.*r.*exp(f);
            bz=b0.*exp(f)+bz0;
            bthe=0;
            bmag=sqrt(br*br+bz*bz);
            bx=br*x(i)/r;
            by=br*y(j)/r;
            bmag3(i,j,k)=bmag;
            bx3(i,j,k)=bx;
            by3(i,j,k)=by;
            bz3(i,j,k)=bz;                        
        end
    end
end


//vacscalar2vtk3d( pict,wd,xx,yy,zz,field,vecsize,filename,directory )
field=0
vecsize=0


vacscalar2vtk3d( 1,bmag3,x',y',z',field,vecsize,'tests','.' )
