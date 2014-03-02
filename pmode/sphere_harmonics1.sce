//Spherical harmonics
//solution for perturbed pressure for a spherical hydrodynamical system

//Using scilab
//http://www.scilab.org/
//Scilab is free and open source software for numerical computation providing a powerful computing environment 
//for engineering and scientific applications. 








//The initial surface definition 
//----------------------
nsize=100;
phi=linspace(-%pi,%pi,nsize);
the=linspace(0,2*%pi,nsize);

x=zeros(1,nsize);
y=x;



nmx=8;
nmy=8;
modeamps=100*rand(nmx,nmy);
modewav=10*rand(nmx,nmy);

x=modeamps(1,1)*cos(the).*sin(phi);
y=modeamps(1,1)*sin(the).*sin(phi);

modeamps=100;
modewav=10;

modeamps=zeros(nmx,nmy);
modewavx=zeros(nmx,nmy);
modewavy=zeros(nmx,nmy);
modeamps(1,1)=20;
modewavx(1,1)=1;
modewavy(1,1)=1;
modeamps(2,2)=8;
modewavx(2,2)=4;
modewavy(2,2)=4;

modeamps(1,2)=5;
modewavx(1,2)=16;
modewavy(1,2)=16;

[X,Y] = meshgrid(x,y);
myones=ones(nsize,nsize);
//Z=sqrt(()^2-x^2-y^2)
//[mmx,mmy]=meshgrid(x,y);

for i=1:size(X,1)
  for j=1:size(X,2)
    Z(i,j) = sqrt((modeamps(1,1))^2-(X(i,j)^2+Y(i,j)^2));
  end
end



//Z=2*max(modeamps)*sin(x)'*cos(y);
//Z=sqrt(()^2-x^2-y^2)

myselec=1;
//global freq;
//global fchange;
global imxsel;
global imysel;

mfreq=0*ones(nmx,nmy);
fchange=500;
imxsel=1;
imysel=1;

// set a new colormap
//-------------------
//cmap= curFig.color_map; //preserve old setting


// plot of a sphere using facets computed by eval3dp 
deff("[x,y,z]=sph(alp,tet)",["x=r*cos(alp).*cos(tet)+orig(1)*ones(tet)";.. 
     "y=r*cos(alp).*sin(tet)+orig(2)*ones(tet)";.. 
     "z=r*sin(alp)+orig(3)*ones(tet)"]); 
     
     
//deff("[x,y,amp]=harmonic(alp,tet)",["x=r*cos(alp).*cos(tet)+orig(1)*ones(tet)";.. 
//     "y=r*cos(alp).*sin(tet)+orig(2)*ones(tet)";..
//     "amp=legendre(l,m,alp').*cos(m.*tet')"]); 


deff("[x,y,amp]=harmonic(alp,tet)",["x=r*cos(alp).*cos(tet)+orig(1)*ones(tet)";.. 
     "y=r*cos(alp).*sin(tet)+orig(2)*ones(tet)";..
     "amp=(-1)^m.*cos(m.*tet).*legendre(l,m,cos(alp))"]); 



//mfreq(1,1)=200; 
//l=10; 
//m=4; 
//shift=10;  //adjust zero of the colour scale
//scale=1.5*10^(-2); //adjust the range of the scale

//mfreq(1,1)=200; 
//l=1; 
//m=0; 
//shift=45;  //adjust zero of the colour scale
//scale=20; //adjust the range of the scale

//mfreq(1,1)=200; 
//l=2; 
//m=0; 
//shift=45;  //adjust zero of the colour scale
//scale=20; //adjust the range of the scale

//mfreq(1,1)=200; 
//l=2; 
//m=2; 
//shift=45;  //adjust zero of the colour scale
//scale=20; //adjust the range of the scale

//mfreq(1,1)=200; 
//l=4; 
//m=2; 
//shift=45;  //adjust zero of the colour scale
//scale=20; //adjust the range of the scale

//mfreq(1,1)=200; 
//l=4; 
//m=4; 
//shift=45;  //adjust zero of the colour scale
//scale=1; //adjust the range of the scale

//mfreq(1,1)=200; 
//l=10; 
//m=4; 
//shift=45;  //adjust zero of the colour scale
//scale=1*10^(-2); //adjust the range of the scale

//mfreq(1,1)=200; 
//l=20; 
//m=0; 
//shift=40;  //adjust zero of the colour scale
//scale=1*10^(2); //adjust the range of the scale

mfreq(1,1)=200; 
l=20; 
m=2; 
shift=5;  //adjust zero of the colour scale
scale=0.015*10^(2); //adjust the range of the scale

 ampones=ones(4,6241);  
r=1; orig=[0 0 0]; 


curFig             = scf(100001);
clf(curFig,"reset");
//curFig.color_map = hotcolormap(64);
//curFig.color_map = autumncolormap(64);
curFig.color_map = jetcolormap(64);

//colorbar();
[xx,yy,zz]=eval3dp(sph,linspace(-%pi/2,%pi/2,80),linspace(0,%pi*2,80));
[xn,yn,amp]=eval3dp(harmonic,linspace(-%pi/2,%pi/2,80),linspace(0,%pi*2,80));
//s=plot3d1(xx,yy,zz,35,45,"X@Y@Z",[-2,2,3]);
//s.color_flag=1 ; //assign facet color according to Z value
plot3d1(xx,yy,list(zz,shift*ampones+scale*amp),35,45,"X@Y@Z",[-2,2,3]);
//colorbar(42*min(amp),42*max(amp));
//plot3d1(xx,yy,zz,35,45,"X@Y@Z",[-2,2,3]);
s=gce();
 for i=1:2000000   



//s=plot3d1(xx,yy,zz,35,45,"X@Y@Z",[-2,2,3]);

//s=plot3d1(xx,yy,list(zz,amp),35,45,"X@Y@Z",[-2,2,3]);
//curFig.color_map = jetcolormap(64);

//s=gce(); //the handle on the surface
//curFig.color_flag=1 ; //assign facet color according to Z value
//curFig=plot3d1(xx,yy,zz,35,45,"X@Y@Z",[-2,2,3]);
//curFig=plot3d1(xx,yy,zz,35,45,"X@Y@Z",[-2,2,3]);
//curFig             = scf(100001);
//clf(curFig,"reset");
//demo_viewCode("membrane.sce");

//drawlater();
//s=surf(xx,yy,zz);
//xselect(); //raise the graphic window
//s.color_flag=2 ; //assign facet color according to Z value
//s.color_map = jetcolormap(64);
//s.color_flag=1 ; //assign facet color according to Z value
pamp=shift*ampones+scale*amp.*sin(mfreq(1,1)*%pi*i/50000);
//plot3d1(xx,yy,list(zz,pamp));
//disp(sin(mfreq(1,1)*%pi*i/10000),i);


//s.data.z=list(zz,pamp);
//s.data.z=pamp;
 
//clf();//plot3d(xx,yy,zz) 


// for i=1:200

//r=sin(i*2*%pi/10);
//[xx,yy,zz]=eval3dp(sph,linspace(-%pi/2,%pi/2,40),linspace(0,%pi*2,20)); 
//s.data.x=xx;
//s.data.y=yy;
//s.data.z=zz;

s.data.color=pamp;


end






