loopimfile='20130105_000012_4096_0171.jpg';
sres='4096';
sinst='0171';
sdat='20130105';

tim=12;
stim=['0000',int2str(tim)];

loopimfile=['img/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];

li1=imread(loopimfile);

% tim=tim+1512;
tim=31436;
stim=['0',int2str(tim)];
%tim=154424;

%stim=int2str(tim);

loopregimfile=['regimg/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
loopimfile=['img/',sdat, '_', stim, '_', sres,'_',sinst,'.jpg'];
li2=imread(loopimfile);

[optimizer,metric] = imregconfig('monomodal');
ttype='translation';
fixed  = rgb2gray(li1);
moving = rgb2gray(li2);
moving_reg = imregister(moving,fixed,ttype,optimizer,metric);

imshow(moving_reg);
imwrite(moving_reg, loopregimfile);