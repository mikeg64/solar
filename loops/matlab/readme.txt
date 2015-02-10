Processing imagery data from SDO.

img2_2014_0312to_0313
http://sdo.gsfc.nasa.gov/gallery/main/item/499

The Sun unleashed a M-9.3 flare, just short of an X class (the largest) from an active region right at the Sun's edge (Mar. 12-13, 2014). The bright flash is the tell tale sign of a flare. The brightness of the flare causes very bright saturation and blooming above and below the flare region on the CCD detector and caused extended diffraction patterns to spread out across the SDO imager. The video clip shows a smaller flare preceded this one as well. The video covers about 15 hours. The still shows the peak of the flare at 22:38 UT on Mar. 12. Images were taken in extreme ultraviolet light, showing ionized iron at 10 million degrees. Credit: Solar Dynamics Observatory/NASA.

Search Tag(s): 131, flare, active region


%event at 54, 89
i2=60;  %52,53,54,55,56,57,58,59,60
i1=96;  %87,88,89,90,91,92,93,95,96





Step 1.

Need to correct for suns rotation
translatelopsequence.m

Redirect a list of image files to a file of filenames img.dat

Translate image relative to first image in the set. Extract date and time from file name


    transtime=(daydat-12)*24*3600+(htim*3600+mtim*60+sectim)-startsectim;

    %solar rotation period for different latitudes
    %http://en.wikipedia.org/wiki/Solar_rotation
    % omega=A+Bsin^2(phi)+Csin^4(phi)
    %omega angular velocity in degrees per day
    % phi = solar latitude
     solA=14.713;%deg/day (p/m 0.0491)
     solB=-2.396;% deg/day (p/m 0.188)
     solC=-1.787;% deg/day (p/m 0.253)
    
    %approximate expression for longitude
    %see https://delicious.com/mikeg64/search/solar-coordinates
    %view image on helioviewer and see information tab on a selected image
    x1=4096/2;
    y1=4096/2;
    r=(3865-421)/2;
    x2=2863;
    y2=1800;
    long=sqrt((x2-x1)^2)/r;
    longdeg=360*long/(2*pi);
    
    solomega=solA+solB*(sin(longdeg))^2+solC*(sin(longdeg))^4;
    
    
    %solar rotation period 24.47 days
    % 360deg rotation takes 2114208sec
    % 1deg rotation takes 5872.8 sec



    %solar radius=6.955x10^5km
    %lh edge=405pixels
    %rh edge=3709 pixels
    %each pixel = 421004.84262m
    %rotation corresponds to translation of 2066.9467 m/s
    % equivalent to 0.00491 pixels/sec
    
    
    solomega=solomega/(24*3600); %convert from deg per day to deg per sec
    solrad=solomega*2*pi/360;
    transrate=solrad*6.955*(10^8)/421004.84262; %pixels per sec

    translation=transtime*transrate;


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





Step 4
An initial examination of the loops was made we prepared the comparpairs script to compare 
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

step 5

we used the histeq to balance the histogram to make it easier to contrast the positions for the loops

Measure the loop properties and co-ordinates along selected loops at the selected range of time steps


step 6

Using the grayscale images measure the intensities of the loops using te 

iprofile;
impixel;
imhist; used to detemine intensity along the line


Step 6.
Feature tracking and extraction
using the detect features
points1 = detectMinEigenFeatures(intensi1);
points2 = detectBRISKFeatures(intensi1);
points3 = detectFASTFeatures(intensi1);
points4 = detectHarrisFeatures(intensi1);
points5 = detectSURFFeatures(intensi1);
points6 = extractHOGFeatures(intensi1);





