
exec('cutofffuns.sce');

[fd,sst,sheetnames,sheetpos]=xls_open('atmos.xls');
[value,textind] = xls_read(fd,sheetpos(1));
mclose(fd);

//[fd1,sst1,sheetnames1,sheetpos1]=xls_open('timeperiods_solareigen1.xls');
//[shifts,textind1] = xls_read(fd1,sheetpos1(2));
//mclose(fd1);


vomega0=0.020944;

vars.mu=0.6e0;
vars.R=8.31e3;
vars.fgamma=1.66666667e0;
vars.ggg=-274.0e0;
vars.g=vars.ggg;
vars.mu=4.d0*%pi/1.0e7;
vars.vomega=vomega0;


//plot(value(:,1),value(:,2));
sz=size(value);
nrows=sz(1);

for i=1:nrows
      shift(i) =fshift(vars,vomega0,value(i,4),value(i,3));
end

ombv=zeros(nrows,1);
vomega=zeros(nrows,1);
lam0=zeros(nrows,1);
vcs=zeros(nrows,1);
//omegas(i) =omega(value(:,1),vars,value(:,4),value(:,3));
    h=value(2,1)-value(1,1);
    
    
    [vomega,lam0,vcs]=omega(value(:,1),vars,value(:,3),value(:,4));
    
    ombv=brunt(value(:,1),vars,value(:,3),value(:,4));
    vbfreq=bfreq(value(:,1),vars,value(:,3),value(:,4),ombv,vcs)
//    for i=1:nrows
//      vomega(i)=0;
//      vomegadash(i)=0;
//    end
//    for i=3:nrows-2
//      vomegadash(i)=diff5p(vomega,i,h);
//    end
 
//    for i=1:nrows
//      vcs(i)=cs(vars,value(i,4),value(i,3));
//      lam0=lambda(vars,value(i,4),value(i,3));
//      vomega(i)=(1+2*vomegadash(i))*(vcs(i)/(2*lam0))^2;
//      tperiod(i)=2*%pi*((2*lam0)/vcs(i))*sqrt( 1/(1+2*vomegadash(i)));
//    end   
    
    
//    deltaom=sqrt((vomega0)^2-(vomega).^2);
//    deltat=2*%pi./deltaom;
    
    
    
      
