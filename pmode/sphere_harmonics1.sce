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





myones=ones(nsize,nsize);




//Z=2*max(modeamps)*sin(x)'*cos(y);
//Z=sqrt(()^2-x^2-y^2)

myselec=1;
//global freq;
//global fchange;
global imxsel;
global imysel;
nmx=1;
nmy=1;
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
     

deff("[x,y,amp]=harmonic(alp,tet)",["x=r*cos(alp).*cos(tet)+orig(1)*ones(tet)";.. 
     "y=r*cos(alp).*sin(tet)+orig(2)*ones(tet)";..
     "amp=(-1)^m.*cos(m.*tet).*legendre(l,m,cos(alp))"]); 

//different options for frequecy and spherical hormonic
//l increase number of modes from pole to pole
//m increases number of mdes observed  arounfd the axis

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
l=2; 
m=0; 
shift=45;  //adjust zero of the colour scale
scale=20; //adjust the range of the scale

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

mfreq(1,1)=400; 
//l=20; 
//m=2; 
//shift=5;  //adjust zero of the colour scale
//scale=0.015*10^(2); //adjust the range of the scale

 ampones=ones(4,6241);  
r=1; orig=[0 0 0]; 


curFig             = scf(100001);
clf(curFig,"reset");
//curFig.color_map = hotcolormap(64);
//curFig.color_map = autumncolormap(64);
curFig.color_map = jetcolormap(64);


[xx,yy,zz]=eval3dp(sph,linspace(-%pi/2,%pi/2,80),linspace(0,%pi*2,80));
[xn,yn,amp]=eval3dp(harmonic,linspace(-%pi/2,%pi/2,80),linspace(0,%pi*2,80));

plot3d1(xx,yy,list(zz,shift*ampones+scale*amp),35,45,"X@Y@Z",[-2,2,3]);
//colorbar(42*min(amp),42*max(amp));

s=gce();
sf=gcf();
sf.figure_size=[493,576];
//for i=1:2000000
for i=1:200   
    pamp=shift*ampones+scale*amp.*sin(mfreq(1,1)*%pi*i/50000);
    s.data.color=pamp;
    
    //uncomment lines below to save images to jpg
    //imfile=msprintf('images/im_%d.jpg',i);
    //xs2jpg(sf,imfile,1);
end






