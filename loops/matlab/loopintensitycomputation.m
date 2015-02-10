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
ii1=1;
ii2=1;
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



bgprof1=improfile(I,lxbg193(ii1,:),lybg193(ii1,:));
bgprof2=improfile(I,lxbg171(ii1,:),lybg171(ii1,:));
prof1=improfile(I,lx193(ii1,:),ly193(ii1,:));
prof2=improfile(I,lx171(ii1,:),ly171(ii1,:));


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


