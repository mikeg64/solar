imnames=importdata('img.dat');
loopimfile='20130105_000012_4096_0171.jpg';
sres='4096';
sinst='0171';
sdat='20130105';

startsectim=12;
tim=startsectim;
stim=['0000',int2str(tim)];

loopimfile=['img/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];

li1=imread(loopimfile);

% tim=tim+1512;
%for i=1:1
for i=1:232   
    var1=imnames{i};
    stim=var1(10:15);
    sdat=var1(1:8);


    htim=str2double(stim(1:2));
    mtim=str2double(stim(3:4));
    sectim=str2double(stim(5:6));
    

    %stim=int2str(tim);

    loopcropimfile=['crop2img/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
    loopimfile=['transimg/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
    li2=imread(loopcropimfile);
    
    intensi1=rgb2ind(li2,0.061);
    bw1=edge(crop_intensi,'sobel');
    figure
    subplot(2,2,1);
    imshow(loopcropimfile);
    subplot(2,2,2);
    imshow(bw1);
    subplot(2,2,3);
    bw1=edge(crop_intensi,'canny');
    imshow(bw1);
    subplot(2,2,4);
    imshowpair(crop_image,bw1);
    
    
    xmin=550;
    ymin=2000;
    width=1300;
    height=1200;
    crop_image=imcrop(li2, [ xmin ymin width height]);

    %imshowpair(li2,trans_image,'diff');
    %imshowpair(li2,rot_image,'blend','Scaling','joint');
    imshow(crop_image);
    imwrite(crop_image, loopcropimfile);

end