
imnames1=importdata('img2_2014_0312to_0313/img_193/img.dat');
imnames2=importdata('img2_2014_0312to_0313/img_171/img.dat');




ntot=82;
%ntot=10;

for i=1:ntot
    
    i2=(i+1)-1;
    i1=(2*i)-1;
    impath1=['img2_2014_0312to_0313/crop2img_193/',imnames1{i1}];
    impath2=['img2_2014_0312to_0313/crop2img_171/',imnames2{i2}];
    
    st1=imnames1{i1};
    st2=imnames2{i2};

    sdat=st1(1:8);
    stim1=st1(10:15);
    stim2=st2(10:15);
    
    
    
    outpath=['img2_2014_0312to_0313/crop2impairs/imp',sdat,'_',stim1,'_193_',stim2,'_HMI171_',int2str(i),'.jpg'];
    li1=imread(impath1);
    li2=imread(impath2);
    
    newfig=figure;
    subplot(1,2,1);
    %imshowpair(li1,li2);
    imshow(li1);

    subplot(1,2,2);
    %imshowpair(li1,li3);
    imshow(li2);
    
    %imwrite(newfig,outpath);
    saveas(newfig,outpath);
    close(gcf);
    %clear newfig;

    
end