%script to compute coronal loop intensities

imnames1=importdata('img2_2014_0312to_0313/img_193/img.dat');
imnames2=importdata('img2_2014_0312to_0313/img_171/img.dat');

i2=45;
i=i2;
i1=(2*i)-1;

%read in the loop properties
loopproperties_img2;

istart1=83;
istart2=50;

%event at 54, 89
%list of corresponding times for each sequence of images
%the index below corresponds to the index of the image in the above arrays

iia1=[50 51 52 53 54 55 56 57 58 59 60 61 62];
iia2=[83 85 87 88 89 90 91 92 93 95 96 97 99];

% i2= 62;   %50,51,52,53,54,55,56,57,58,59,60,61,62
% i1= 99;   %83,85,87,88,89,90,91,92,93,95,96,97,99
% ii1= 13;

for i=1:13

ii1=i;
ii2=i;
i1=iia1(ii1);
i2=iia2(ii2);

impath1=['img2_2014_0312to_0313/crop2img_193/',imnames1{i1}];
impath2=['img2_2014_0312to_0313/crop2img_171/',imnames2{i2}];

li1=imread(impath1);
li2=imread(impath2);

ig1=rgb2gray(li1);
ig2=rgb2gray(li2);



%iprofile;
%impixel;
%imhist; used to detemine intensity along the line



bgprof1=improfile(ig1,lxbg193(ii1,:),lybg193(ii1,:));
bgprof2=improfile(ig2,lxbg171(ii1,:),lybg171(ii1,:));
prof1=improfile(ig1,lx193(ii1,:),ly193(ii1,:));
prof2=improfile(ig2,lx171(ii1,:),ly171(ii1,:));

sumbg1=sum(bgprof1);
sumbg2=sum(bgprof2);

meanbg1=mean(bgprof1);
stdbg1=std(bgprof1);
meanbg2=mean(bgprof2);
stdbg2=std(bgprof2);

mean1=mean(prof1);
std1=std(prof1);
mean2=mean(prof2);
std2=std(prof2);

sum1=sum(prof1);
sum2=sum(prof2);


%compute line length
np=2;
for i=1:np-1
   bgll=(lxbg193(ii1,i+1)-lxbg193(ii1,i)).^2+(lybg193(ii1,i+1)-lybg193(ii1,i)).^2;
   bgll=sqrt(bgll); 
end
bgll193=bgll;

np=npos;
for i=1:np-1
   ll=(lx193(ii1,i+1)-lx193(ii1,i)).^2+(ly193(ii1,i+1)-ly193(ii1,i)).^2;
   ll=sqrt(ll); 
end
ll193=ll;


np=2;
for i=1:np-1
   bgll=(lxbg171(ii1,i+1)-lxbg171(ii1,i)).^2+(lybg171(ii1,i+1)-lybg171(ii1,i)).^2;
   bgll=sqrt(bgll); 
end
bgll171=bgll;

np=npos;
for i=1:np-1
   ll=(lx171(ii1,i+1)-lx171(ii1,i)).^2+(ly171(ii1,i+1)-ly171(ii1,i)).^2;
   ll=sqrt(ll); 
end
ll171=ll;





%from Isothermal and Multithermal Analysis of Coronal Loops Observed with AIA,
%Schmelz, 2011ApJ...731...49S

% where the fluxes are given in units of data numbers per second (DN s?1) 
%after normalizing the counts (DN) by the exposure time, which is different
%in each filter.
% 
% Here is the technique for computing the DN's from ref 12
% After co-aligning the data, we chose 10 pixels along the spine of each 
%loop and 10 pixels from a clean background area. We averaged the loop- and
%background-pixel values and calculated standard deviations. 
%Then we subtracted the background and propagated the errors. 
%The background-subtracted averages and uncertainties were then normalized
%by the exposure time for the appropriate filter, resulting in units of 
%Data Numbers s?1. A background region is adjacent to the loop under consideration.

%compute the DNs using the above method

%EXPTIME171=0.069883000; %secEXPTIME: 1.9995570
%EXPTIME193= 2.0000320;  %2.0000630

et171=mean(exptime171);
et193=mean(exptime193);


lintensity193(ii1)=((sum1/ll193  )-(sumbg1/bgll193))/et193;
lstd193(ii1)=abs(sqrt(((std1/ll193  ).^2+(stdbg1/bgll193).^2))/et193);


lintensity171(ii2)=((sum2/ll171  )-(sumbg2/bgll171))/et171;
lstd171(ii2)=abs(sqrt(((std2/ll171  ).^2+(stdbg2/bgll171).^2))/et171);


end

