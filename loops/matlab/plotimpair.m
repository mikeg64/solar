%This script is used to generate and image plots which are used to measure
%coronal llop prperties
%these data values are then cellected in the file lop properties_img2.m

%view of full solar disk
%http://helioviewer.org/?date=2014-03-12T22:36:37.000Z&imageScale=2.4204409&centerX=30.25551125&centerY=-50.8292589&imageLayers=%5BSDO,AIA,AIA,171,1,49%5D,%5BSDO,AIA,AIA,193,1,52%5D&eventLayers=&eventLabels=true



imnames1=importdata('img2_2014_0312to_0313/img_193/img.dat');
imnames2=importdata('img2_2014_0312to_0313/img_171/img.dat');

i2=45;
i=i2;
i1=(2*i)-1;

istart1=83;
istart2=50;

%event at 54, 89
%list of corresponding times for each sequence of images
%the index below corresponds to the index of the image in the above arrays

i2= 62;   %50,51,52,53,54,55,56,57,58,59,60,61,62
i1= 99;   %83,85,87,88,89,90,91,92,93,95,96,97,99
ii1= 13;
impath1=['img2_2014_0312to_0313/crop2img_193/',imnames1{i1}];
impath2=['img2_2014_0312to_0313/crop2img_171/',imnames2{i2}];

var1=imnames1{1};
stim1=var1(10:15);
sdat1=var1(1:8);
htim1=str2double(stim1(1:2));
mtim1=str2double(stim1(3:4));
sectim1=str2double(stim1(5:6));
daydat1=str2double(var1(7:8));
startsectim1=(daydat1-12)*24*3600+(htim1*3600+mtim1*60+sectim1);

var1=imnames1{i1};
stim1=var1(10:15);
sdat1=var1(1:8);
htim1=str2double(stim1(1:2));
mtim1=str2double(stim1(3:4));
sectim1=str2double(stim1(5:6));
daydat1=str2double(var1(7:8));
totsectime1=(daydat1-12)*24*3600+(htim1*3600+mtim1*60+sectim1)-startsectim1;
strtim1=var1(7:15);

var2=imnames2{1};
stim2=var2(10:15);
sdat2=var2(1:8);
htim2=str2double(stim2(1:2));
mtim2=str2double(stim2(3:4));
sectim2=str2double(stim2(5:6));
daydat2=str2double(var2(7:8));
startsectim2=(daydat2-12)*24*3600+(htim2*3600+mtim2*60+sectim2);

var2=imnames2{i2};
stim2=var2(10:15);
sdat2=var2(1:8);
htim2=str2double(stim2(1:2));
mtim2=str2double(stim2(3:4));
sectim2=str2double(stim2(5:6));
daydat2=str2double(var2(7:8));
totsectime2=(daydat2-12)*24*3600+(htim2*3600+mtim2*60+sectim2)-startsectim2;
strtim2=var2(7:15);

loopproperties_img2;
    
li1=imread(impath1);
li2=imread(impath2);

ig1=rgb2gray(li1);
ig2=rgb2gray(li2);

%balance the contrast to make measuremetns easier
ige1=histeq(ig1);
ige2=histeq(ig2);

figure
imshow(ige1);
x=lx193(ii1,:);
y=ly193(ii1,:);
line(x,y);
%plot(x,y,'+');


figure
imshow(ige2);
x=lx171(1+i2-istart2,:);
y=ly171(1+i2-istart2,:);
line(x,y);