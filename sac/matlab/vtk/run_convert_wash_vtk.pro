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
;timesteps=1999
timesteps = 3
timesteps = 480
npic=timesteps

time_nn=dblarr(1+timesteps)


pic=1
pic=80
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
directory='/fastdata/cs1mkg/smaug/washmc_2p5_2p5_12p5_mov4_kg/'

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

for iw=1,3 do begin
  w(*,*,*,iw)=w(*,*,*,iw)/(w(*,*,*,0)+w(*,*,*,9))
endfor

;w(*,*,*,0)=w(*,*,*,0)+w(*,*,*,9)
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


;vac2vtk3d,ipic,w,x,3,6,name
vac2vtk3d,ipic,w,x,5,3,'b',directory

; bfield
vac2vtk3d,ipic,w,x,10,3,'bb',directory


; vfield
vac2vtk3d,ipic,w,x,1,3,'v',directory

;density perturbation
vacscalar2vtk3d,ipic,w,x,0,1,'dens',directory

nnp=nnp+1

endfor
;endwhile



end
