pro line2_b_classic, bx,bz,x, z, jj, scale ; jj = 0  plot, jj=1 oplot

x(*)=x(*)/scale
z(*)=z(*)/scale

n1=n_elements(z)
n2=n_elements(x)

xl=fltarr(1)
yl=fltarr(1)
dx=fltarr(1)

xline=fltarr(1)
yline=fltarr(1)

xl[0]=x[0]
yl[0]=z[0]


coeff=dblarr(n2)

nn=30
xtec=x[0]
i=0
bmax=max(sqrt(bx(0,*)^2.d0+bz(0,*)^2.d0))

coeff(*)=sqrt(bx(0,*)^2.d0+bz(0,*)^2.d0)/bmax


while xtec le x[n2-1] do begin
 coefftec=INTERPOL(coeff,x, xtec)

 if i eq 0 then dx(i)=(x(n2-1)-x(0))/(2.0*coefftec+1.d0)/nn else $ 
                dx=[dx,(x(n2-1)-x(0))/(2.0*coefftec+1.d0)/nn]
 
 
 xtec=xtec+dx[i]
 xl=[xl,xtec]
 yl=[yl,z[0]]
 print, xl[i],yl[i],i
 i=i+1
endwhile
 
 dl=1.0d4/scale

bxbz=sqrt(bx^2.d0+bz^2.d0)

for i=3, n_elements(xl)-1 do begin
 xxline=xl[i]
 yyline=yl[i]

 xle=fltarr(1)
 yle=fltarr(1)
 
 xle=xxline
 yle=yyline

while (xxline le x[n2-1]) and (xxline ge x[0]) and $ 
      (yyline le z[n1-1]) and (yyline ge z[0]) do begin
  
     xc=interpol(dindgen(n2),x,xxline)
     yc=interpol(dindgen(n1),z,yyline)

  Bbx=interpolate(bx, yc, xc)
  Bby=interpolate(bz, yc, xc)    
  Bb=interpolate(bxbz, yc, xc)    

 ddx=dl*Bbx/Bb
 ddy=dl*Bby/Bb
 
 xxline=xxline+ddx
 yyline=yyline+ddy     
 
 
 xle=[xle,xxline]
 yle=[yle,yyline]
 
 print, '###', xxline, yyline     
endwhile
     if n_elements(xle) gt 2 then begin
   if jj eq 0 then  begin 
                        plot,yle, xle, xrange=[z[0],z[n1-1]], yrange=[x[0],x[n2-1]]
                        jj=1
                    endif  else begin
                	oplot,yle, xle
                    endelse
  endif		    	   
 
endfor 

end

function deriv1,f,x
nel=n_elements(f)
nel1=n_elements(x)
if (nel ne nel1) then begin
 print,'Inconsistant input, stop.'
 stop
endif
res=dblarr(nel)
for i=2,nel-3 do res(i)=(1.d0/12.D0/(x(i+1)-x(i)))*(8.d0*f(i+1)-8.d0*f(i-1)-f(i+2)+f(i-2))
;for i=1,nel-2 do res(i)=(1.d0/2.d0/(x(i+1)-x(i)))*(f(i+1)-f(i-1))
res(0)=res(2)
res(1)=res(2)
res(nel-1)=res(nel-3)
res(nel-2)=res(nel-3)
return,res
end

function par3,x,x0

if (x le -2.d0*x0) then res=0.d0
if ((x ge -2.d0*x0) and (x le -x0)) then res=(x+3.d0*x0)*(x+x0)+x0^2.d0
if ((x ge -x0) and (x le x0)) then res=-(x+x0)*(x-x0)+x0^2.d0
if ((x ge x0) and (x le 2.d0*x0)) then res=(x-3.d0*x0)*(x-x0)+x0^2.d0
if (x ge 2.d0*x0) then res=0.d0

return,res
end

function par4,x,x0

res=exp(-1.d0*(x/x0)^1.0d0)

return,res
end


function inte ,f,dx
nel=n_elements(f)

res=0.d0
if (nel gt 1) then begin
 
if (nel eq 2) then res=dx*0.5d0*(f(1)+f(0))

if (nel gt 2) then begin 

  nel1=1.0d0*nel
  f2=congrid(f,nel1)

; res=int_tabulated(dindgen(nel)*dx,f,/double) 

  for k=1,nel1-1 do res=res+0.5d0*(f2(k-1)+f2(k))*(dx/1.d0)
;  sum=sum+dx*(f(0)+f(nel-1))/3.d0

endif
endif

return,res
end

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)


window, 0,xsize=1025,ysize=300,XPOS = 700, YPOS = 900 
window, 1,xsize=600,ysize=300,XPOS = 100, YPOS = 900 
window, 2,xsize=600,ysize=600,XPOS = 50, YPOS = 500 
window, 4,xsize=900,ysize=600,XPOS = 300, YPOS = 500 
;***************
headline='                                                                               '
it=long(1)
ndim=long(1)
neqpar=long(1)
nw=long(1)
varname='                                                                               '
time=double(1)
dum=long(1)
dumd=long(1)
nn=0
close,1
openr,1,'/data/cs1mkg/smaug_pmode/configs/3D_128_spic_bin.ini',/f77_unf
;openr,1,'/fastdata/cs1mkg/smaug/configs/3D_256_4Mm_bin2.ini'


;*****************************************************

readu,1,headline
readu,1,it,time,ndim,neqpar,nw
gencoord=(ndim lt 0)
ndim=abs(ndim)
nx=lonarr(ndim)
readu,1,nx
print,'tuta', neqpar
eqpar=dblarr(neqpar)
readu,1,eqpar
readu,1,varname

n1=nx(0)
n2=nx(1)
n3=nx(2)

x_code=dblarr(n1,n2,n3,ndim)
w=dblarr(n1,n2,n3,nw)

wi=dblarr(n1,n2,n3)

readu,1,x_code
for iw=0,nw-1 do begin
 print, iw
 readu,1,wi
  w(*,*,*,iw)=wi
endfor
print, n1,n2,n3

;*************************************************
;stop
yy=49

gamma=1.66666667d0
ggg=-274.0d0
mu=4.d0*!PI/1.0d7


y=reform(x_code(0,0,*,2))
x=reform(x_code(0,*,0,1))
z=reform(x_code(*,0,0,0))

wset,0
!p.multi = [0,3,0,0,1]

plot,x, title='x', charsize=2.0
plot,y, title='y', charsize=2.0
plot,z, title='z', charsize=2.0

scale=1.0d6

;*************** start pressure ******************

p=dblarr(n1,n2,n3)

p[*,*,*]=w[*,*,*,4]+w[*,*,*,8]
p[*,*,*]=p[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/ $
         (w[*,*,*,0]+w[*,*,*,9])/2.0
p[*,*,*]=p[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0 $
          +(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.0
p[*,*,*]=(gamma-1.d0)*p[*,*,*]


wset,1
!p.multi = [0,2,0,0,1]

plot, z(*)/1.0d6,alog10(p(*,yy,yy)), title='alog(p_gas)', charsize=1.2, xtitle='[Mm]'

;*************** end pressure ******************


bx=dblarr(n1,n2,n3)
by=dblarr(n1,n2,n3)
bz=dblarr(n1,n2,n3)


rho=reform(w(*,*,*,0)+w(*,*,*,9))
rho1=dblarr(n1,n2,n3)

plot, z(*)/1.0d6,alog10(rho(*,yy,yy)), title='log(rho)', charsize=1.2, xtitle='[Mm]'


b0z=dblarr(n1)


x=x-(max(x)-min(x))/2.d0 ;+5000.d0
y=y-(max(y)-min(y))/2.d0 ;+5000.d0


d_z=0.2d0 ; width of Gaussian in Mm
z_shift=0.d0 ; shift in Mm
A=0.05d0 ; amplitude

for i=0,n1-1 do b0z[i]=par4((z[i]/scale-z_shift),d_z)

;for i=0,n1-1 do b0z[i]=par3((z[i]/scale-z_shift),d_z)

;zc=0.2d0
;for i=0,n1-1 do b0z[i]=1.d0-((atan((z[i]/max(z)-zc)*15.d0))+!Pi/2.d0)/!Pi


wset,2
!p.multi = [0,2,2,0,1]

xf=dblarr(n1,n2,n3)


;xr=0.20d5
;yr=0.20d5


;xr=1.0d6
;yr=1.0d6

xr=0.2d6
yr=0.2d6


R2=(xr^2.d0+yr^2.d0)

A=R2/2.d0

for k=0,n3-1 do begin 
for j=0,n2-1 do begin
for i=0,n1-1 do begin

f=(x[j]^2.d0+y[k]^2.d0)/R2

xf[i,j,k]=exp(-f)

endfor
endfor
endfor







plot, z/scale, b0z, title='b0z', charsize=1.2, xtitle='[Mm]'

;for i=n1-2,0,-1 do b0z[i]=b0z[i]+b0z[i+1]


;Bmax=0.0075d0  ; mag field Tesla 75G
;Bmax=0.0005d0  ; mag field Tesla 50G
;Bmax=0.0025d0  ; mag field Tesla 25G
;Bmax=0.009d0  ; mag field Tesla 90G

;Bmax=0.01d0  ; mag field Tesla 100G
;Bmax=0.005d0  ; mag field Tesla 50G
Bmax=0.02d0  ; mag field Tesla 200G

;Bmax=0.00010d0  ; mag field Tesla 10G
;Bmax=0.010d0  ; mag field Tesla
;Bmin=0.0006d0  ; mag field Tesla
;Bmin=0.0005d0  ; mag field Tesla 50G
Bmin=0.00001d0  ; mag field Tesla 1G

; ************** convert to VAC magnetic field

bx(*,*,*)=bx(*,*,*)/sqrt(mu)
by(*,*,*)=by(*,*,*)/sqrt(mu)
bz(*,*,*)=Bmax*bz(*,*,*)/sqrt(mu)

;*********************************************

;tvframe, divb(20,*,*), title='divb',/bar, charsize=1.0


;tvframe,p-(gamma-1.d0)*(bx^2.0+bz^2.0)/2.0, title='delta p', charsize=1.8, /bar
;print, 'min p=',min(p-(gamma-1.d0)*(bx^2.0+bz^2.0)/2.0)


rho1=dblarr(n1,n2,n3)


rho1=0
p=dblarr(n1,n2,n3)

rho1=rho+rho1

print, 'p'
for i=0,n1-1 do begin
for j=0,n2-1 do begin
for k=0,n3-1 do begin

p(i,j,k)=p(i,j,k)+w(i,j,k,8)*(gamma-1.d0)
bz(i,j,k)=Bmax*xf[i,j,k]/sqrt(mu)

endfor
endfor
endfor


beta=(bx^2.d0+by^2.d0+bz^2.d0)/2.d0/p

tvframe,beta(*,*,yy), charsize=1.2, title='1/beta', /bar;CT='VpVd'
contour, beta(*,*,yy), LEVELS = [0.001,0.01,0.1,1.0, 10.0], $
         C_Annotation = ['1000.0','100.0','10.0','1.0','0.1'], /overplot,/follow
	 
tvframe,rho1(*,*,yy), charsize=1.2, title='rho', /bar	 
tvframe,p(*,*,yy), charsize=1.2, title='p', /bar	 	 

mu_gas=1.2d0
R=8.3e+003

T=mu_gas*p/R/rho1

tvframe,T(*,*,yy), charsize=1.2, title='T', /bar	 	 


;stop

;lower boundary

for ix_1=3,2,-1 do begin
  for ix_2=0,n2-1 do begin
  for ix_3=0,n3-1 do begin  
         p_2=rho1(ix_1,ix_2,ix_3)*ggg
         p(ix_1-1,ix_2,ix_3) = (z(1)-z(0))*p_2+p(ix_1,ix_2,ix_3)
  endfor  
  endfor
 endfor


;upper boundary

for ix_1=n1-3,n1-2 do begin
   for ix_2=0,n2-1 do begin
   for ix_3=0,n3-1 do begin   
           p_2=rho1(ix_1,ix_2,ix_3)*ggg
           p(ix_1+1,ix_2,ix_3) = -(z(1)-z(0))*p_2+p(ix_1,ix_2,ix_3)
   endfor	   
   endfor
endfor


e=p/(gamma-1.d0)+0.5d0*(bx*bx+by*by+bz*bz)

;'h m1 m2 m3 e b1 b2 b3 eb rhob bg1 bg2 bg3'

w(*,*,*,9)=rho1
w(*,*,*,8)=e
;w(*,*,*,10)=rotate(bz,-1)
;w(*,*,*,11)=rotate(by,-1)
;w(*,*,*,12)=rotate(bx,-1)

w(*,*,*,10)=bz
w(*,*,*,11)=bx
w(*,*,*,12)=by

;w(*,*,*,5)=bz
;w(*,*,*,6)=bx
;w(*,*,*,7)=by



tvframe, w(*,*,yy,0)+w(*,*,yy,9), /bar, charsize=1.5, title='rho_t' 
;tvframe, w(*,*,yy,10), title='bz', /bar, charsize=1.5
;tvframe, w(*,*,yy,11), title='bx', /bar, charsize=1.5
;tvframe, w(*,yy,*,12), title='by', /bar, charsize=1.5
tvframe, w(*,*,yy,5), title='bz', /bar, charsize=1.5
tvframe, w(*,*,yy,6), title='bx', /bar, charsize=1.5
tvframe, w(*,yy,*,7), title='by', /bar, charsize=1.5



; PROVERIT !!!!!
;dpdz=dblarr(n1,n2)
;dpdx=dblarr(n1,n2)
;for i=0,n2-1 do dpdz(*,i)=deriv1(p(*,i),z)
;for i=0,n1-1 do dpdx(i,*)=deriv1(p(i,*),x)
;diff1=-dpdz-G+rho1*ggg
;diff2=-dpdx+F

;tvframe, diff1, title='diff1', /bar, charsize=2.0
;tvframe, diff2, title='diff2', /bar, charsize=2.0


;qq=dbzdz+dbxdx
;goto, wwww
close,1
;openw,1,'/data/cs1mkg/smaug_realpmode/configs/magvert/3D_128_spic_bvertbg500G_bin.ini',/f77_unf
;openw,1,'/fastdata/cs1mkg/smaug/configs/3D_256_4Mm_bvertbg200Gn_bin.ini',/f77_unf
openw,1,'/data/cs1mkg/smaug_realpmode/configs/magvert/3D_128_spic_bvertbg200Gn_bin.ini',/f77_unf


writeu,1,headline
writeu,1,it,time,ndim,neqpar,nw
writeu,1,nx
writeu,1,eqpar
writeu,1,varname
writeu,1,x_code
for iw=0,nw-1 do begin
wi=w(*,*,*,iw)
writeu,1,wi
endfor


 
close,1

;wwww :
print, 'complete'
end





