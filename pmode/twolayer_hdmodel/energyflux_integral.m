%compute enrgy flux integral

%run the script twolayer_hdmodel.m to initiailise the variables
%use runge-kutta to compute solution at top boundary?
nt=period*3/dt;
for it=1:nt;
    velocity=zeros(nx,ny,nz);
    pressure=zeros(nx,ny,nz);
    locflux=zeros(nx,ny,nz);
    intflux=zeros(nz);
    
    %eflux=(gamma-1).*(perturbedenergy-(mom_z./2)).*(mom_z./(density_bkg+density_perturb));
    %compute energy-flux out each cros-section position for different heights
    for iz=1:nz;
     
    z=consts.dz*iz;

         %use function and definition to compute p0        
         %use function and definition to compute rho0
         if z>consts.l1
            %corona 
             [temp, lrho0, lp0]=solar_profiles_corona(z,consts.rho0c, consts.pres0c, consts.t0c,consts.hc);
         else
            %chromosphere 
             [temp, lrho0, lp0]=solar_profiles_phot(z,consts.rho0p, consts.pres0p, consts.t0p,consts.hp);
         end
    
    
    for ix=1:nx;    
     for iy=1:ny;
         
         
         %because of different lower boundary conditiondifferent points
         %have different values for constants
         
         %compute epl1 and epl2
         %values of lagrange displacement at base of layer 1 and base of
         %layer 2
         
         tpres=0;
         [epl1,tpres]=lowerboundary(ix,iy,1,it, consts,consts.rho0p);
         [epl2,tpres]=lowerboundary(ix,iy,1,it, consts,consts.rho0c);
         [ d1c, d1p, d2c, d2p ] = dconsts( omega, consts, rho0c, rho0p, pres0c, pres0p, epl1, epl2, l1, l2 );

        consts.d1c=d1c;
        consts.d2c=d2c;
        consts.d2c=d2c;
        consts.d2p=d2p;

         
         
         if iz==1
             
            if iz>l1
               [tt, trho0, tpres]=solar_profiles_corona(z,consts.rho0c, consts.pres0c, consts.t0c,consts.hc);
            else
               [tt, trho0, tpres]=solar_profiles_phot(z,consts.rho0p, consts.pres0p, consts.t0p,consts.hp);
            end
             
            [lagdisp,pres]=lowerboundary(ix,iy,iz,it, consts,trho0); 
            lagdispm=lagdisp;
         elseif iz==nz
             pres=0;
             lagdisp=0;
             lagdispm=lagdisp;
         else
             %use expression for lagarngian pressure to determine pressure
             %perturbation
         %compute lagrange displacement (differential with respect to time
         %to get velpcity
         
         %compute pressure
         lagdisp=lagrange_disp(z,consts, epl1, epl2, l1, l2, d1c, d1p,d2c,d2p );
         
         ddispdz=(lagdisp-lagdispm)/(consts.dz);
         
         

         
         
         pres=consts.g*lrho0*lagdisp-consts.gamma*lp0*ddispdz;
         lagdispm=lagdisp;
         
         
             
         end
         
         %compute flux
         locflux(ix,iy,iz)=lagdisp*pres*consts.omega;   
         intflux(iz)=intflux(iz)+locflux(ix,iy,iz);
         
     end
    end
     
    totflux(it,iz)=intflux(iz);
   
        
    end %compute energy-flux out each cros-section position for different heights
    
    savname=['fluxout_',num2str(it),'.mat'];
    save(savname, 'it','velocity', 'pressure', 'locflux', 'intflux');
    
    
end