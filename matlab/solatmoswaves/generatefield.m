function [simparams, simgridinfo, simdata]=generatefield(simparams, simgridinfo, simdata, mode)

%Generate magnetic field configuration
%simple fluxtube using self similarity and hydrostatic pressure correction



%Generate field
% calculate hydrostatic pressure update


%cases
%1.  Vertical field
%2. Horizontal field
%3. Inclined (off vertical field)
%4. Flux tube array

if strcmp(mode,'fluxtube')
    
    
    %tube width =100km
    %footpoint intensity 1kG ala thin photospheric flux tubes
    nb=4;
    mb=4;
    
    nx1=simgridinfo.grid_dimensions(1);
    nx2=simgridinfo.grid_dimensions(2);
    nx3=simgridinfo.grid_dimensions(3);
    
    dx=(simparams.domain_left_edge(1)-simparams.domain_right_edge(1))/simgridinfo.grid_dimensions(1);
    dy=(simparams.domain_left_edge(2)-simparams.domain_right_edge(2))/simgridinfo.grid_dimensions(2);
    dz=(simparams.domain_left_edge(3)-simparams.domain_right_edge(3))/simgridinfo.grid_dimensions(3);

    bx=zeros(nx1,nx2,nx3);
    by=zeros(nx1,nx2,nx3);
    bz=zeros(nx1,nx2,nx3);
    b0z=zeros(nx1);
    
    xf=zeros(nx1,nx2,nx3);
    
    d_z=1.5; % width of Gaussian in Mm
    z_shift= 0.0; % shift in Mm
    A=0.45; % amplitude
    scale=1.0e6;
    
    x=zeros(nx1);
    y=zeros(nx2);
    z=zeros(nx3);
    
    for j=1:nx2
        y(j)=simparams.domain_left_edge(2)+dy*(j-1);
    end
    for k=1:nx3
        z(k)=simparams.domain_left_edge(3)+dz*(k-1);
    end
    
    for i=1:nx1
        x(i)=simparams.domain_left_edge(1)+dx*(i-1);
        b0z(i)=(par4((x(i)/scale-z_shift),d_z,A)).^2;
    end
    
    
    b0z_top=0.08;
    b0z=b0z./max(b0z);
    f0=2.0e6; %tube opening factor

    Ab0z=20.d0; % bz - amplitude
    b0z=Ab0z.*b0z+b0z_top;
    
    dbz=deriv1(b0z,x);

%     x=x-max(x)/2.d0
%     y=y-max(y)/2.d0

    
    
    
    for ib=0:nb-1
        for jb=0:mb-1
            
            if nb>1
                ybp=simparams.domain_left_edge(2)+(ib+1)*nx2*dy*(simgridinfo.grid_dimensions(2)/((nb+1)))-nx2*dy*(simgridinfo.grid_dimensions(2)/(2*(nb+1)));
            else
                ybp=simparams.domain_left_edge(2)+(ib+1)*nx2*dy*(simgridinfo.grid_dimensions(2)/((nb+1)));                
            end
            
            if mb>1
                zbp=simparams.domain_left_edge(3)+(jb+1)*nx3*dz*(simgridinfo.grid_dimensions(3)/((mb+1)))-nx3*dz*(simgridinfo.grid_dimensions(3)/(2*(mb+1)));                
            else
                zbp=simparams.domain_left_edge(3)+(jb+1)*nx3*dz*(simgridinfo.grid_dimensions(3)/((mb+1)));
            end
            
  

for k=1:nx3
for j=1:nx2
for i=1:nx1

f=b0z(i)*sqrt((y(j)-ybp).^2+(z(k)-zbp).^2);

xf(i,j,k)=(par4(f,f0,0.5)).^2;



end
end
end
xf=xf/max(xf);

for k=1:n3 
for j=1:n2
for i=1:n1

bz(i,j,k)=bz(i,j,k)+(b0z(i)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k));
bx(i,j,k)=bx(i,j,k)-(dbz(i)*(x(j)-ybp)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k));
by(i,j,k)=by(i,j,k)-dbz(i)*(y(k)-zbp)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k);

end
end
end           
                       
        end  %end of loop
    end
    
end   %building fluxtube


function [res]=par4(x,x0,A)
    res=A.*exp(-1.0*x.^2./(x0.^2));

%end


%compute magnetostatics pressure correction




%update the fields and save to output









function [res]=deriv1(f,x)
nel=size(f);
nel1=size(x);

if nel(1)>nel(2)
    num=nel(1);
else
    num=nel(2);
end

% if (nel ne nel1) then begin
%  print,'Inconsistant input, stop.'
%  stop
% endif

res=zeros(num);
for i=2:num-3 
    res(i)=(1./12.0./(x(i+1)-x(i)))*(8.d0*f(i+1)-8.d0*f(i-1)-f(i+2)+f(i-2));
end
%for i=1,nel-2 do res(i)=(1.d0/2.d0/(x(i+1)-x(i)))*(f(i+1)-f(i-1))
res(0)=res(2);
res(1)=res(2);
res(num-1)=res(num-3);
res(num-2)=res(num-3);

%end



