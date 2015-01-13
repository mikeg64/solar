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





