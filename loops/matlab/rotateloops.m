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

htim=1;
mtim=15;
sectim=1;
tim=31436;
%stim=['0',int2str(htim),int2str(mtim),int2str(sectim)];
stim=['0',int2str(htim),int2str(mtim),'0',int2str(sectim)];
%stim=['0',int2str(htim),'00','00'];
%tim=154424;

rottime=(htim*3600+mtim*60+sectim)-startsectim;

%solar rotation period 24.47 days
% 360deg rotation takes 2114208sec
% 1deg rotation takes 5872.8 sec
rotation=rottime/5872.8;

%stim=int2str(tim);

looprotimfile=['rotimg/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
loopimfile=['img/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
li2=imread(loopimfile);

rot_image=imrotate(li2, rotation,'nearest', 'crop');

imshowpair(li2,rot_image,'diff');
%imshowpair(li2,rot_image,'blend','Scaling','joint');
%imshow(rot_image);
imwrite(rot_image, looprotimfile);