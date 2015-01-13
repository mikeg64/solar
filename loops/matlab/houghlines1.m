%intensi1=rgb2ind(crop_image,0.061);
[H,T,R] = hough(intensi1,'RhoResolution',0.1,'Theta',-90:0.1:89.5);
%[H,T,R] = hough(intensi1);
%P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
%lines = houghlines(intensi1,T,R,P,'FillGap',5,'MinLength',7);

figure, imshow(crop_image), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end