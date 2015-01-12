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
for i=98:232
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

    transtime=(htim*3600+mtim*60+sectim)-startsectim;

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

    %stim=int2str(tim);

    looprotimfile=['transimg/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
    loopimfile=['img/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
    li2=imread(loopimfile);

    trans_image=imtranslate(li2, [-translation, 0],'FillValues',255);

    imshowpair(li2,trans_image,'diff');
    %imshowpair(li2,rot_image,'blend','Scaling','joint');
    %imshow(rot_image);
    imwrite(trans_image, looprotimfile);

end