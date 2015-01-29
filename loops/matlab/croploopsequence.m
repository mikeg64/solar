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
startsectim=2947;
%startsectim=4007;

tim=startsectim;
%stim=['0000',int2str(tim)];
stim=['00',int2str(tim)];

%loopimfile=['img/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
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
    

    %stim=int2str(tim);

    loopcropimfile=['img2_2014_0312to_0313/crop2img_171/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
    loopimfile=['img2_2014_0312to_0313/transimg_171/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
    li2=imread(loopimfile);

%     xmin=550;
%     ymin=2000;
%     width=1300;
%     height=1200;
    
%    xmin=1769;
%    ymin=1653;
%    width=900;
%    height=860;

    xmin=2213;
    ymin=1400;
    width=1300;
    height=800;

    crop_image=imcrop(li2, [ xmin ymin width height]);

    %imshowpair(li2,trans_image,'diff');
    %imshowpair(li2,rot_image,'blend','Scaling','joint');
    imshow(crop_image);
    imwrite(crop_image, loopcropimfile);

end