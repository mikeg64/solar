%imnames=importdata('img.dat');

%line feature from from1 crop1img
x=[848 865 888 913 947 984 1007 1005 1005 1005];
y=[543 518 493 470 453 447 472 502 532 575];
f1=['crop1img/',imnames{1}];
f2=['crop1img/',imnames{5}];
f3=['crop1img/',imnames{7}];
f4=['crop1img/',imnames{15}];
f5=['crop1img/',imnames{40}];

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
imshow(li4);


subplot(3,2,4);
%imshowpair(li1,li5);
imshow(li4);
hold on
line(x,y);


subplot(3,2,5);
%imshowpair(li1,li4);
imshow(li5);


subplot(3,2,6);
%imshowpair(li1,li5);
imshow(li5);
hold on
line(x,y);