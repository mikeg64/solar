
val3f=fopen('val3.dat','w');
for i=128:-1:1
   fprintf(val3f,'%g %g %g %g \n', xmin+dx*(i-1),densg(i),tempg(i),presg(i) ); 
    
    
end
fclose(val3f);