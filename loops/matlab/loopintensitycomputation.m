%script to compute coronal loop intensities

imnames1=importdata('img2_2014_0312to_0313/img_193/img.dat');
imnames2=importdata('img2_2014_0312to_0313/img_171/img.dat');

i2=45;
i=i2;
i1=(2*i)-1;

var1=imnames1{1};
stim1=var1(10:15);
sdat1=var1(1:8);
htim1=str2double(stim1(1:2));
mtim1=str2double(stim1(3:4));
sectim1=str2double(stim1(5:6));
daydat1=str2double(var1(7:8));
startsectim1=(daydat1-12)*24*3600+(htim1*3600+mtim1*60+sectim1);

var2=imnames2{1};
stim2=var2(10:15);
sdat2=var2(1:8);
htim2=str2double(stim2(1:2));
mtim2=str2double(stim2(3:4));
sectim2=str2double(stim2(5:6));
daydat2=str2double(var2(7:8));
startsectim2=(daydat2-12)*24*3600+(htim2*3600+mtim2*60+sectim2);

%read in the loop properties
loopproperties_partialloops_img2;

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

for j=1:13
%for i=1:1

ii1=j;
ii2=j;
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

lx=lx193{ii1};
ly=ly193{ii1};
prof1=improfile(ig1,lx(:),ly(:));

lx=lx171{ii2};
ly=ly171{ii2};
prof2=improfile(ig2,lx(:),ly(:));

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

sz=size(lx193{ii1});
lx=lx193{ii1};
ly=ly193{ii1};
%np=npos;
np=sz(2);
for i=1:np-1
   ll=(lx(i+1)-lx(i)).^2+(ly(i+1)-ly(i)).^2;
   ll=sqrt(ll); 
end
ll193=ll;


np=2;
for i=1:np-1
   bgll=(lxbg171(ii1,i+1)-lxbg171(ii1,i)).^2+(lybg171(ii1,i+1)-lybg171(ii1,i)).^2;
   bgll=sqrt(bgll); 
end
bgll171=bgll;

sz=size(lx171{ii2});
lx=lx171{ii2};
ly=ly171{ii2};
%np=npos;
np=sz(2);

for i=1:np-1
   ll=(lx(i+1)-lx(i)).^2+(ly(i+1)-ly(i)).^2;
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


f171_193=lintensity171./lintensity193;
figure;
plot(startsectim1+tsec193 ,lintensity193, '-', startsectim2+tsec171  ,lintensity171, '--');
% hold on
% plot(startsectim1+tsec193 ,lintensity193)


%flux_171/flux_193=resp_171/resp_193
%use tabulated response functions to compute
resp_171_193=zeros(30,1);
temp_171_193=zeros(30,1);
 for ii=1:30
   logtemp= 5.55+2.25*(ii-1)/29; 
   temp_171_193(ii)=logtemp;
   
   iiresp=1;
   iresploop=1;
   sz=size(respxlogt193);
   while iresploop==1  && (iiresp+1)<=sz(2)
       iiresp=iiresp+1;
       if logtemp>=reslogt193(iiresp-1) && logtemp<reslogt193(iiresp)
           iresploop=0;
           xval=logtemp;
           x=reslogt193;
           %lagrange interpolate logDN
%            t1=(xval-x(1,iiresp)).*(xval-x(1,iiresp+1))/((x(1,iiresp-1)-x(1,iiresp)).*(x(1,iiresp-1)-x(1,iiresp+1)));
%            t2=(xval-x(1,iiresp-1)).*(xval-x(1,iiresp+1))/((x(1,iiresp)-x(1,iiresp-1)).*(x(1,iiresp)-x(1,iiresp+1)));
%            t3=(xval-x(1,iiresp-1)).*(xval-x(1,iiresp))/((x(1,iiresp+1)-x(1,iiresp-1)).*(x(1,iiresp+1)-x(1,iiresp)));
         
           %logdn193=t1*reslogdn193(iiresp-1)+t2*reslogdn193(iiresp)+t3*reslogdn193(iiresp+1); 
           
           t1=(xval-x(1,iiresp))/(x(1,iiresp-1)-x(1,iiresp));
           t2=(xval-x(1,iiresp-1))/(x(1,iiresp)-x(1,iiresp-1));
           logdn193=t1*reslogdn193(iiresp-1)+t2*reslogdn193(iiresp);   
       end    
   end
   
   iresploop=1;
   iiresp=1;
   sz=size(respxlogt171);
   logtemp= 5.55+2.25*(ii-1)/29;
  
    while iresploop==1  && (iiresp+1)<=sz(2)
       iiresp=iiresp+1;
       if logtemp>=reslogt171(iiresp-1) && logtemp<reslogt171(iiresp)
           iresploop=0;
           xval=logtemp;
           x=reslogt171;
           %lagrange interpolate logDN
%            t1=(xval-x(1,iiresp)).*(xval-x(1,iiresp+1))/((x(1,iiresp-1)-x(1,iiresp)).*(x(1,iiresp-1)-x(1,iiresp+1)));
%            t2=(xval-x(1,iiresp-1)).*(xval-x(1,iiresp+1))/((x(1,iiresp)-x(1,iiresp-1)).*(x(1,iiresp)-x(1,iiresp+1)));
%            t3=(xval-x(1,iiresp-1)).*(xval-x(1,iiresp))/((x(1,iiresp+1)-x(1,iiresp-1)).*(x(1,iiresp+1)-x(1,iiresp)));
         
           %logdn171=t1*reslogdn171(iiresp-1)+t2*reslogdn171(iiresp)+t3*reslogdn171(iiresp+1); 
           
           t1=(xval-x(1,iiresp))/(x(1,iiresp-1)-x(1,iiresp));
           t2=(xval-x(1,iiresp-1))/(x(1,iiresp)-x(1,iiresp-1));
           logdn171=t1*reslogdn171(iiresp-1)+t2*reslogdn171(iiresp);   
       end    
   end  
   
   
   
   
   
   
   resp_171_193(ii)=(10.^logdn171)/(10.^logdn193);
   
     
 end
 
 figure;
 plot(temp_171_193 ,resp_171_193);
 
 %read off where the 
f171_193=lintensity171./lintensity193;

finaltemp=[ 5.535 5.524 5.518 5.57 5.53 5.519 5.51 5.515 5.517 5.528 5.522 5.51 5.51  ];


