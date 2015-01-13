figure
subplot(2,2,1);
imshow(crop_image);
subplot(2,2,2);
imshow(crop_image);
hold on
plot(points1.selectStrongest(60));
subplot(2,2,3);
imshow(crop_image);
hold on
plot(points2.selectStrongest(60));
subplot(2,2,4);
imshow(crop_image);
hold on
plot(points3.selectStrongest(60));