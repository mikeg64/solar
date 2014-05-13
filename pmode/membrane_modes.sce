
// =============================================================================
// Animation of a vibrating membrane
// ============================================================================


//Using scilab
//http://www.scilab.org/
//Scilab is free and open source software for numerical computation providing a powerful computing environment 
//for engineering and scientific applications.

curFig             = scf(100001);
clf(curFig,"reset");


drawlater();

// set a new colormap
//-------------------
cmap= curFig.color_map; //preserve old setting
curFig.color_map = jetcolormap(64);


cs=1343.5;
a00=5;
lx=2;
ly=2;

n1=3;
n2=3;
n1o=n1;
n2o=n2;
anm=2*a00/(n1*n1+n2*n2+2*(n1+n2)+2);

n1_1=0;
n2_1=1;
n1=n1_1;
n2=n2_1;
anm1=2*a00/(n1*n1+n2*n2+2*(n1+n2)+2);

n1_2=1;
n2_2=2;
n1=n1_2;
n2=n2_2;
anm2=2*a00/(n1*n1+n2*n2+2*(n1+n2)+2);

n1_3=2;
n2_3=3;
n1=n1_3;
n2=n2_3;
anm3=2*a00/(n1*n1+n2*n2+2*(n1+n2)+2);


n1_4=3;
n2_4=4;
n1=n1_4;
n2=n2_4;
anm4=2*a00/(n1*n1+n2*n2+2*(n1+n2)+2);

n1=n1o;
n2=n2o;

anm2=0;
anm3=0;
anm4=0;

ommn=%pi*cs*sqrt(((n1+1)/lx).^2+((n2+1)/ly).^2);
ommn1=%pi*cs*sqrt(((n1_1+1)/lx).^2+((n2_1+1)/ly).^2);
ommn2=%pi*cs*sqrt(((n1_2+1)/lx).^2+((n2_2+1)/ly).^2);
ommn3=%pi*cs*sqrt(((n1_3+1)/lx).^2+((n2_3+1)/ly).^2);
ommn4=%pi*cs*sqrt(((n1_4+1)/lx).^2+((n2_4+1)/ly).^2);

//ommn=950*%pi;

//The initial surface definition 
//----------------------
nsize=100;
x=linspace(-%pi,%pi,nsize);
y=x;
Z=15*sin(x)'*cos(y);

myones=ones(nsize,nsize);
[mmx,mmy]=meshgrid(x,y);
//Creates and set graphical entities which represent the surface
//--------------------------------------------------------------
plot3d1(x,y,Z,35,45,"X@Y@Z",[-2,2,3]);
s=gce(); //the handle on the surface
s.color_flag=1 ; //assign facet color according to Z value
title("evolution of a 3d surface","fontsize",3)

I=4000:-0.1:1;
realtimeinit(0.1);;//set time step (0.1 seconds)  and date reference


drawnow();

for i=1:max(size(I))
  realtime(i); //wait till date 0.1*i seconds
 
  s.data.z = anm1*sin((ommn1*i*myones/max(size(I)))).*(sin(((n1_1+1)*mmx/lx)).*sin(((n2_1+1)*mmy/ly)))+..
             anm2*sin((ommn2*i*myones/max(size(I)))).*(sin(((n1_2+1)*mmx/lx)).*sin(((n2_2+1)*mmy/ly)))+..
             anm3*sin((ommn3*i*myones/max(size(I)))).*(sin(((n1_3+1)*mmx/lx)).*sin(((n2_3+1)*mmy/ly)))+..
             anm4*sin((ommn4*i*myones/max(size(I)))).*(sin(((n1_4+1)*mmx/lx)).*sin(((n2_4+1)*mmy/ly)));

    //uncomment lines below to save images to jpg
    //imfile=msprintf('images/im_%d.jpg',i);
    //xs2jpg(sf,imfile,1);
  
end
