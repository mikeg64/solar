%Translate sequence of so that region of interest is fixed using solar
%rotation and feature position on the solar disc

%see also for solar rotation corrections
%Solar Corona Loop Studies with AIA: I. Cross-Sectional Temperature Structure
%Markus J. Aschwanden, Paul Boerner, arXiv.org > astro-ph > arXiv:1103.0228, 
%(ADS 2011ApJ...732...81A)

%imnames=importdata('img2_2014_0312to_0313/img_193/img.dat');
imnames=importdata('img2_2014_0312to_0313/img_171/img.dat');
%loopimfile='20130105_000012_4096_0171.jpg';
loopimfile='20140312_002947_4096_HMI171';
%loopimfile='20140312_004007_4096_0193';
sres='4096';
%sinst='0171';
%sinst='0193';
sinst='HMI171';
%sdat='20130105';
sdat='20140312';

%startsectim=12;
%startsectim=4007;
startsectim=2947;
tim=startsectim;
%stim=['0000',int2str(tim)];
stim=['00',int2str(tim)];

%loopimfile=['img2_2014_0312to_0313/img_193/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
loopimfile=['img2_2014_0312to_0313/img_171/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
li1=imread(loopimfile);

% tim=tim+1512;
%for i=1:1
for i=1:121 
    var1=imnames{i};
    stim=var1(10:15);
    sdat=var1(1:8);


    htim=str2double(stim(1:2));
    mtim=str2double(stim(3:4));
    sectim=str2double(stim(5:6));
    tim=31436;
    %stim=[int2str(htim),int2str(mtim),int2str(sectim)];
    %stim=[int2str(htim),'00',int2str(sectim)];
    %stim=[int2str(htim),int2str(mtim),'00'];
    %stim=[int2str(htim),'00','00'];
    %stim=['0',int2str(htim),int2str(mtim),int2str(sectim)];
    %stim=['0',int2str(htim),'0',int2str(mtim),int2str(sectim)];
    %stim=['0',int2str(htim),'00',int2str(sectim)];
    %stim=['0',int2str(htim),int2str(mtim),'00'];
    %stim=['0',int2str(htim),'00','00'];
    %tim=154424;
    daydat=str2double(var1(7:8));
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

    %stim=int2str(tim);

    %looprotimfile=['img2_2014_0312to_0313/transimg_193/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
    looprotimfile=['img2_2014_0312to_0313/transimg_171/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
    %loopimfile=['img2_2014_0312to_0313/img_193/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
    loopimfile=['img2_2014_0312to_0313/img_171/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
    li2=imread(loopimfile);

    %rotate image so that rotational axis is in the same direction as y-axis of the image 
    %rot_image=imrotate(li2,7.25);
    trans_image=imtranslate(li2, [-translation, 0],'FillValues',255);
    %rotate the image back
    %trans_image=imrotate(temptrans_image,-7.25);
    imshowpair(li2,trans_image,'diff');
    %imshowpair(li2,rot_image,'blend','Scaling','joint');
    %imshow(rot_image);
    imwrite(trans_image, looprotimfile);

end