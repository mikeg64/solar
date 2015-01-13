Processing imagery data from SDO.


Step 1.

Need to correct for suns rotation
translatelopsequence.m

Redirect a list of image files to a file of filenames img.dat

Translate image relative to first image in the set. Extract date and time from file name

    transtime=(daydat-5)*24*3600+(htim*3600+mtim*60+sectim)-startsectim;

    %solar rotation period 24.47 days
    % 360deg rotation takes 2114208sec
    % 1deg rotation takes 5872.8 sec



    %solar radius=6.955x10^5km
    %lh edge=405pixels
    %rh edge=3709 pixels
    %each pixel = 421004.84262m
    %rotation corresponds to translation of 2066.9467 m/s
    % equivalent to 0.00491 pixels/sec

    translation=transtime*0.00491;

	trans_image=imtranslate(li2, [-translation, 0],'FillValues',255);



Step 2.
Select regions of interest from stabilised images
use script
croploopsequence.m

	    xmin=550;
	    ymin=2000;
	    width=1300;
	    height=1200;
	    crop_image=imcrop(li2, [ xmin ymin width height]);



Step 3.
Feature tracking and extraction
using the detect features
points1 = detectMinEigenFeatures(intensi1);
points2 = detectBRISKFeatures(intensi1);
points3 = detectFASTFeatures(intensi1);
points4 = detectHarrisFeatures(intensi1);
points5 = detectSURFFeatures(intensi1);
points6 = extractHOGFeatures(intensi1);

Step 4
Most effective method is to click and select points from features, for example
clicked points along a flux loop

%line feature from from1 crop1img
x=[848 865 888 913 947 984 1007 1005 1005 1005];
y=[543 518 493 470 453 447 472 502 532 575];

we plot these points on subsequent plots to detect changes
e.g. see comparepairs script

f1=['crop1img/',imnames{1}];
f2=['crop1img/',imnames{5}];
f3=['crop1img/',imnames{7}];
f4=['crop1img/',imnames{15}];
f5=['crop1img/',imnames{5}];

li1=imread(f1);
li2=imread(f2);
li3=imread(f3);
li4=imread(f4);
li5=imread(f5);

figure
subplot(3,2,1);
%imshowpair(li1,li2);
imshow(li1);

subplot(3,2,2);
%imshowpair(li1,li3);
imshow(li1);
hold on
line(x,y);

subplot(3,2,3);
%imshowpair(li1,li4);
imshow(li2);





