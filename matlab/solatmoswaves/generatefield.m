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
    nb=1;
    mb=1;
    
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
    b0z_top=0.08;
   
    f0=2.0e6; %tube opening factor

    Ab0z=20.d0; % bz - amplitude
    
    
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
        %b0z(i)=(par3((x(i)/scale-z_shift),d_z,A));
    end
    
    dyb=(simparams.domain_right_edge(2)-simparams.domain_left_edge(2))/nb;
    dzb=(simparams.domain_right_edge(3)-simparams.domain_left_edge(3))/mb;
    
    
    

    
    
    
    
    b0z=b0z./max(b0z);
   
    b0z=Ab0z.*b0z+b0z_top;
    
    dbz=deriv1(b0z,x);


  xold=x;
  yold=y;
  zold=z;
    
    
 ib=0;
 jb=0;
    for ib=0:nb-1
        for jb=0:mb-1





            
          %  if nb>1
          %      ybp=simparams.domain_left_edge(2)+(ib+1)*nx2*dy*(simgridinfo.grid_dimensions(2)/((nb+1)))-nx2*dy*(simgridinfo.grid_dimensions(2)/(2*(nb+1)));
          %  else
          %      ybp=simparams.domain_left_edge(2)+(ib+1)*nx2*dy*(1/((nb+1)));                
          %  end
            
          %  if mb>1
          %      zbp=simparams.domain_left_edge(3)+(jb+1)*nx3*dz*(simgridinfo.grid_dimensions(3)/((mb+1)))-nx3*dz*(simgridinfo.grid_dimensions(3)/(2*(mb+1)));                
          %  else
          %      zbp=simparams.domain_left_edge(3)+(jb+1)*nx3*dz*(1/((mb+1)));
          %  end
            
  
          %  z=zold-max(z)/2.d0;
          %  y=yold-max(y)/2.d0;

            
            z=zold-(jb-1)*dzb-(dzb/2);
            y=yold-(ib-1)*dyb-(dyb/2);           
            
            
            
            
for k=1:nx3
for j=1:nx2
for i=1:nx1

f=b0z(i)*sqrt((y(j)).^2+(z(k)).^2);

xf(i,j,k)=(par4(f,f0,0.5)).^2;



end
end
end
xf=xf./max(max(max(xf)));

for k=1:n3 
for j=1:n2
for i=1:n1

% bz(i,j,k)=bz(i,j,k)+(b0z(i)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k));
% bx(i,j,k)=bx(i,j,k)-(dbz(i)*(x(j)-ybp)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k));
% by(i,j,k)=by(i,j,k)-dbz(i)*(y(k)-zbp)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k);
bz(i,j,k)=bz(i,j,k)+(b0z(i)/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k));
bx(i,j,k)=bx(i,j,k)-(dbz(i)*(x(j))/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k));
by(i,j,k)=by(i,j,k)-dbz(i)*(y(k))/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k);



end
end
end           
                       
         end  %end of loop
     end
%     
 end   %building fluxtube
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%checked to here


%compute magnetostatics pressure correction
dbzdz=zeros(nx1,nx2,nx3);
dbxdz=zeros(nx1,nx2,nx3);
dbydz=zeros(nx1,nx2,nx3);

dbzdx=zeros(nx1,nx2,nx3);
dbxdx=zeros(nx1,nx2,nx3);
dbydx=zeros(nx1,nx2,nx3);

dbzdy=zeros(nx1,nx2,nx3);
dbxdy=zeros(nx1,nx2,nx3);
dbydy=zeros(nx1,nx2,nx3);

br=zeros(nx1,nx2,nx3);
Bvarix=zeros(nx1,nx2,nx3);
Bvariy=zeros(nx1,nx2,nx3);
Bvar=zeros(nx1,nx2,nx3);
dbvardz=zeros(nx1,nx2,nx3);

for k=1:nx3
for j=1:nx2
 dbzdz(:,j,k)=deriv1(bz(:,j,k),x);
 dbxdz(:,j,k)=deriv1(bx(:,j,k),x);
 dbydz(:,j,k)=deriv1(by(:,j,k),x);
end
end


for k=1:nx3
for i=1:nx1
 dbzdx(i,:,k)=deriv1(bz(i,:,k),y);
 dbxdx(i,:,k)=deriv1(bx(i,:,k),y);
 dbydx(i,:,k)=deriv1(by(i,:,k),y);
end
end


for j=1:nx2
for i=1:nx1
 dbzdy(i,j,:)=deriv1(bz(i,j,:),z);
 dbxdy(i,j,:)=deriv1(bx(i,j,:),z);
 dbydy(i,j,:)=deriv1(by(i,j,:),z); 
end
end

divb=dbzdz+dbxdx+dbydy;


for i=1:nx1
for j=1:nx2
for k=1:nx3
  br(i,j,k)=(bx(i,j,k)+bz(i,j,k)).*sqrt(x(j).^2+y(k).^2)/(x(j)+y(k));
end
end
end



for i=1:nx1
for j=1:nx2
for k=1:nx3
  dbzdr(i,j,k)=dbzdx(i,j,k).*(x(j)./sqrt(x(j).^2+y(k).^2))+dbzdy(i,j,k).*(y(k)./sqrt(x(j).^2+y(k).^2))
end
end
end


% ***** dbrdz
for k=1:nx3
for j=1:nx2
 dbrdz(:,j,k)=deriv1(br(:,j,k),z);
end
end



F=bz.*(dbrdz-dbzdr)
G=br.*(dbrdz-dbzdr)


bvar=bz*br

for k=1:nx3
for j=1:nx2
 dbvardz(:,j,k)=deriv1(bvar(:,j,k),z);
end
end



%***** Bvar+br^2/r
for i=1:nx1
for j=1:nx2
for k=1:nx3
  Bvar(i,j,k)=dbvardz(i,j,k)+br(i,j,k).^2/sqrt(x(j).^2+y(k).^2);
end
end
end

for i=1,nx1
for k=1,nx3
  for j=1,nx2
   sum=inte((Bvar(i,1:j,k).*x(j)./sqrt(x(j).^2+y(k).^2)),x(1)-x(0))
   Bvarix(i,j,k)=sum
  end
end
end

for i=1,nx1
for j=1,nx2
  for k=1,nx3
   sum=inte((Bvar(i,j,0:k).*y(k)./sqrt(x(j).^2+y(k).^2)),y(1)-y(0));
  Bvariy(i,j,k)=sum
  end
end

end


Bvari=Bvarix+Bvariy-bz^2./2+br.^2./2;


for i=1:n2 
    bvaridz(:,i)=deriv1(bvari(:,i),z);
end


%intFdz=-bvaridz/ggg


for j=1:n1
  for i=1:n2
  rho1(j,i)=(bvaridz(j,i)+G(j,i))./ggg;
  end
end
%update the fields and save to output






rho1=rho+rho1
p=Bvari+simdata.w(:,:,:,5)*(consts.fgamma-1.0);





%rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
%set background density
%set background energy


%update the background energy and magnetic fields
simdata.w(:,:,:,10)=rho1;
simdata.w(:,:,:,9)=p./((consts.fgamma-1.0))+0.5*(bx.*bx+bz.*bz+by.*by)
simdata.w(:,:,:,11)=bx;
simdata.w(:,:,:,12)=by;
simdata.w(:,:,:,13)=bz;









function [res]=inte(f,dx)

nsiz=size(f);

if nsiz(1)>nsiz(2)
    nel=nsiz(1);
else
    nel=nsiz(2);
end

res=0.0;
if (nel > 1) 
    if (nel == 2) 
        res=dx*0.5*(f(1)+f(0));
    end

    if (nel > 2)
      for k=1:nel-1 
          res=res+0.5*(f(k-1)+f(k))*dx;
      end
    end
end %outer

%return,res
end




function [res]=par4(x,x0,A)
    res=A.*exp(-1.0*x.^2./(x0.^2));

end


function [res]=par3(x,x0,A)

    A=A/(x0.^2);
    
  

if (x <= -2*x0)
    res=0;
elseif ((x >= -2*x0) & (x <= -x0)) 
    res=A*(x+3*x0)*(x+x0)+A*x0.^2;
elseif ((x >= -x0) & (x <= x0)) 
    res=-A*(x+x0)*(x-x0)+A*x0.^2;
elseif ((x >= x0) & (x <= 2*x0)) 
    res=A*(x-3*x0)*(x-x0)+A*x0.^2;
elseif (x >= 2*x0) 
    res=0;
end

%return,res
end



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
    res(i)=(1./12.0./(x(i+1)-x(i)))*(8.0*f(i+1)-8.0*f(i-1)-f(i+2)+f(i-2));
end
%for i=1,nel-2 do res(i)=(1.d0/2.d0/(x(i+1)-x(i)))*(f(i+1)-f(i-1))
res(0)=res(2);
res(1)=res(2);
res(num-1)=res(num-3);
res(num-2)=res(num-3);

end



