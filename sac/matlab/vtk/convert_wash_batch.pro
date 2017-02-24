tarr=dblarr(1)
maxa=fltarr(1)
mina=fltarr(1)
cuta=fltarr(2000,50)

;DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
;WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
;PRINT, 'Date:      ', systime(0)
;PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
;PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


ii=1

nnp=1

if (ii eq 1) then begin
;loadct,4
;mixct
endif else begin
loadct,0
tek_color
endelse


;timesteps = 78
timesteps=1999
;timesteps = 3

npic=timesteps

time_nn=dblarr(1+timesteps)


pic=1
nn=0;
for ipic=pic,npic do begin

mass=dblarr(1)
egas=dblarr(1)
tm=dblarr(1)
dtt=dblarr(1)

ia=1.0

headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)


;nn=0


close,1

;openr,1,'/data/ap1vf/3D/torsional_driver_puls_long/3D_tube_196_100_100.out',/f77_unf

;openr,1,'/data/ap1vf/3D/torsional_driver_cont_2min_LA/3D_tube_196_100_100.out',/f77_unf

;openr,1,'/data/ap1vf/3D_tube_196_100_100_multidriver_lower.out',/f77_unf

;directory='/data/cs1mkg/smaug_spicule1/spicule4b0_b4_3d/'
;directory='/fastdata/cs1mkg/smaug/spic4b0_b4_3d/'
;directory='/fastdata/cs1mkg/smaug/spic6p7a_0_0_3d/'
;directory='/fastdata/cs1mkg/smaug/spic2p3a_0_3_3d/'
;directory='/fastdata/cs1mkg/smaug/spic4p3a_0_1_3d/'
;directory='/fastdata/cs1mkg/smaug/spic2p3a_0_3_3d/'
;directory='/fastdata/cs1mkg/smaug/washing_mach/'
directory='/fastdata/cs1mkg/smaug/washmc_2p5_2p5_12p5_180_kg/'

;name='zeroOT_'
;name='3D_em_t1_bin_np010808_00'
;name='3D_em_f1_bin_'
name='washmc__'



;while not(eof(1)) do begin

;picid=ipic*5+4
picid=ipic*1000L
;picid=ipic;
outfile=directory+name+strtrim(string(picid),2)+'.out'
print,'ipic=',ipic
;openr,1,outfile,/f77_unf
openr,1,outfile








readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
tarr=[tarr,time]
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname



xout=dblarr(3)
yout=dblarr(3)


n1=nx(0)
n2=nx(1)
n3=nx(2)
x=dblarr(n1,n2,n3,ndim)

wi=dblarr(n1,n2,n3)

w=fltarr(n1,n2,n3,nw)

readu,1,x
for iw=0,nw-1 do begin
 readu,1,wi
  w(*,*,*,iw)=wi
endfor


if nnp eq 1 then begin

print,time, nn

ttime=dblarr(1)


zmin=0
zmax=110

zmax=180
zmax=127


xmin=19
xmax=127

xmin=0
xmax=127


ymin=19
ymax=79

ymin=0
ymax=127


ttime[0]=time

;///////////////////////////////////////
;-- VAPOR
;//////////////////////////////////////
dim = [zmax-zmin+1,xmax-xmin+1,ymax-ymin+1]

num_levels = 0



if (ipic eq pic) then  begin

mfd = vdf_create(dim,num_levels)


vdf_setnumtimesteps, mfd,timesteps

time_nn[nn]=time

;strTrim(string(time),1)+'  nn = '+strTrim(string(nn),1)

;VDF_SETTSTR,mfd,nn,'UserTimeStampString','0.1' ;strcompress(ss,/remove_all)
  
  FOR i=0,timesteps-1 DO BEGIN
    VDF_SETTSTR,mfd,i,'UserTimeStampString',strTrim(string(FORMAT='(6F10.2)', time_nn[i]),2)+'  n = '+strTrim(string(i),1)
;print, strTrim(string(time_nn[i]),1)
  ENDFOR ;i
  
  ;strTrim(string(FORMAT='(6F10.2)', time),2)+' it ='+strTrim(string(it),1)
  
varnames = ['bx','by','bz','vx','vy','vz','vv','rhot','rho']

vdf_setvarnames, mfd, varnames

extents = [x(zmin,0,0,0), x(0,xmin,0,1), x(0,0,ymin,2), x(zmax,0,0,0), x(0,xmax,0,1), x(0,0,ymax,2) ]
;print, x(zmin,0,0,0), x(0,xmin,0,1), x(0,0,ymin,2), x(zmax,0,0,0), x(0,xmax,0,1), x(0,0,ymax,2)
;stop
vdf_setextents, mfd, extents

;
; Set a global comment
;

vdf_setcomment, mfd, 'This is my SAC data'

attribute_name = 'MyMetadata'
f = findgen(100)
Vdf_setdbl,mfd,attribute_name, f

;vdffile = '/data/ap1vf/3D/torsional_driver_cont_2min_LA/3D_tube_196_100_100_cont_traj_2min.vdf'
;vdffile = '/data/ap1vf/3D/torsional_driver_puls_long/3D_tube_196_100_100_puls_long_1.vdf'
;vdffile = '/data/cs1mkg/smaug_spicule1/spicule4b0_b4_3d/vap/vapt.vdf'
;vdffile = '/fastdata/cs1mkg/smaug/spic5b0_3d/vap/vapt.vdf'

;vdffile = '/fastdata/cs1mkg/smaug/spic6p7a_0_0_3d/vap/vapt.vdf'
;vdffile = '/fastdata/cs1mkg/smaug/spic2p3a_0_3_3d/vap/vapt.vdf'
;vdffile = '/fastdata/cs1mkg/smaug/washing_mach/vap/vapt.vdf'
vdffile = '/fastdata/cs1mkg/smaug/washmc_2p5_2p5_12p5_180_kg/vap/vapt.vdf'

vdf_write, mfd, vdffile
vdf_destroy, mfd

endif


;vdffile = '/data/ap1vf/3D/torsional_driver_cont_2min_LA/3D_tube_196_100_100_cont_traj_2min.vdf'
;vdffile = '/data/ap1vf/3D/torsional_driver_puls_long/3D_tube_196_100_100_puls_long_1.vdf'
;vdffile = '/data/cs1mkg/smaug_spicule1/spicule4b0_b4_3d/vap/vapt.vdf'
;vdffile = '/fastdata/cs1mkg/smaug/spic4b0_b4_3d/vap/vapt.vdf'
;vdffile = '/fastdata/cs1mkg/smaug/spic6p7a_0_0_3d/vap/vapt.vdf'
;vdffile = '/fastdata/cs1mkg/smaug/spic2p3a_0_3_3d/vap/vapt.vdf'
;vdffile = '/fastdata/cs1mkg/smaug/washing_mach/vap/vapt.vdf'
vdffile = '/fastdata/cs1mkg/smaug/washmc_2p5_2p5_12p5_180_kg/vap/vapt.vdf'

dfd = vdc_bufwritecreate(vdffile)


; Get the data volume that we wish to store.


sac_rho=reform(w(zmin:zmax,xmin:xmax,ymin:ymax,0))


sac_rhot=reform(w(zmin:zmax,xmin:xmax,ymin:ymax,0)+ $ 
               w(zmin:zmax,xmin:xmax,ymin:ymax,9))

sac_vx=reform(w(zmin:zmax,xmin:xmax,ymin:ymax,2)/(w(zmin:zmax,xmin:xmax,ymin:ymax,0)+ $
              w(zmin:zmax,xmin:xmax,ymin:ymax,9)))
sac_vy=reform(w(zmin:zmax,xmin:xmax,ymin:ymax,3)/(w(zmin:zmax,xmin:xmax,ymin:ymax,0)+ $ 
              w(zmin:zmax,xmin:xmax,ymin:ymax,9)))
sac_vz=reform(w(zmin:zmax,xmin:xmax,ymin:ymax,1)/(w(zmin:zmax,xmin:xmax,ymin:ymax,0)+ $ 
              w(zmin:zmax,xmin:xmax,ymin:ymax,9)))

sac_bx=reform(w(zmin:zmax,xmin:xmax,ymin:ymax,6)+w(zmin:zmax,xmin:xmax,ymin:ymax,11))
sac_by=reform(w(zmin:zmax,xmin:xmax,ymin:ymax,7)+w(zmin:zmax,xmin:xmax,ymin:ymax,12))
sac_bz=reform(w(zmin:zmax,xmin:xmax,ymin:ymax,5)+w(zmin:zmax,xmin:xmax,ymin:ymax,10))
;sac_bx=reform(w(zmin:zmax,xmin:xmax,ymin:ymax,6))
;sac_by=reform(w(zmin:zmax,xmin:xmax,ymin:ymax,7))
;sac_bz=reform(w(zmin:zmax,xmin:xmax,ymin:ymax,5))


;sac_bt=sqrt((w(zmin:zmax,xmin:xmax,ymin:ymax,5)+w(zmin:zmax,xmin:xmax,ymin:ymax,10))^2.0+  $
;            (w(zmin:zmax,xmin:xmax,ymin:ymax,6)+w(zmin:zmax,xmin:xmax,ymin:ymax,11))^2.0+  $
;	    (w(zmin:zmax,xmin:xmax,ymin:ymax,7)+w(zmin:zmax,xmin:xmax,ymin:ymax,12))^2.0)


bx = sac_bx
by = sac_by
bz = sac_bz
vx = sac_vx
vy = sac_vy
vz = sac_vz

vv=sqrt(vx*vx+vy*vy)

;bt = sac_bt
rho = sac_rho
rhot = sac_rhot


; Prepare the data set for writing. We need to identify
; the time step and the name of the variable that
; we wish to store. In this case, the first time step,
; zero, and the variable named ÔvxÕ
;

dim= [zmax-zmin+1,xmax-xmin+1,xmax-xmin+1]

vdc_openvarwrite, dfd, nn, 'bx'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, bx[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'by'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, by[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'bz'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, bz[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'vx'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, vx[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'vy'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, vy[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'vz'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, vz[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'vv'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, vv[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'rhot'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, rhot[*,*,z]
endfor
vdc_closevar, dfd

vdc_openvarwrite, dfd, nn, 'rho'
for z = 0, dim[2]-1 do begin
vdc_bufwriteslice, dfd, rho[*,*,z]
endfor
vdc_closevar, dfd

;vdc_openvarwrite, dfd, nn, 'bt'
;for z = 0, dim[2]-1 do begin
;vdc_bufwriteslice, dfd, bt[*,*,z]
;endfor
;vdc_closevar, dfd
                          ;
; Write (transform) the volume to the data set one
; slice at a time


;////////////////////////////////////////
;////////////////////////////////////////

nn=nn+1

nnp=0

endif 

nnp=nnp+1
;vdc_bufwritedestroy, dfd

endfor
;endwhile


;-- CLOSE VAPOR

;An Overview of VAPOR Data Collections 12
; Close the currently opened variable/time-step. We're
; done writing to it
;
;
; Destroy the "buffered write" data transformation
; object. We're done with it.
;
vdc_bufwritedestroy, dfd

end
