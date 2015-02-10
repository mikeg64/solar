%imnames=importdata('img.dat');

%view of full solar disk
%http://helioviewer.org/?date=2014-03-12T22:36:37.000Z&imageScale=2.4204409&centerX=30.25551125&centerY=-50.8292589&imageLayers=%5BSDO,AIA,AIA,171,1,49%5D,%5BSDO,AIA,AIA,193,1,52%5D&eventLayers=&eventLabels=true


%line feature from from1 crop1img
%x=[848 865 888 913 947 984 1007 1005 1005 1005];
%y=[543 518 493 470 453 447 472 502 532 575];

x=[405 409 453 519 547 553 543];
y=[259 227 173 147 179 223 255 ];
f1=['img2_2014_0312to_0313/crop2img/',imnames{1}];
f2=['img2_2014_0312to_0313/crop2img/',imnames{5}];
f3=['img2_2014_0312to_0313/crop2img/',imnames{7}];
f4=['img2_2014_0312to_0313/crop2img/',imnames{10}];
f5=['img2_2014_0312to_0313/crop2img/',imnames{40}];

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
plot(x,y,'+');

subplot(3,2,3);
%imshowpair(li1,li4);
imshow(li4);


subplot(3,2,4);
%imshowpair(li1,li5);
imshow(li4);
hold on
%line(x,y);
plot(x,y,'+');


subplot(3,2,5);
%imshowpair(li1,li4);
imshow(li5);


subplot(3,2,6);
%imshowpair(li1,li5);
imshow(li5);
hold on
%line(x,y);
plot(x,y,'+');