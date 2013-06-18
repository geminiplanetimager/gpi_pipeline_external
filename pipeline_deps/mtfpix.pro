; Pixel MTF function
; By Christian Marois
; used in FFTSCALE


;Genere la mtf des pixels (FLOAT) centree en (0,0) par defaut. Voir Astronomical optics p309.

function mtfpix,D,lambda,dim,res,xc=xc,yc=yc,scfact=scfact

if not keyword_set(xc) then xc = 0.
if not keyword_set(yc) then yc = 0.
if not keyword_set(scfact) then scfact=1.0

fc=dim*res*1./((lambda/D)*(10^(-6.))*180.*(1/(!pi))*3600.)

f=fc/dim
if keyword_set(fact) then f=f*fact

x=findgen(dim)#replicate(1.,dim)

y=transpose(x)

x=x-fix(dim/2)
y=y-fix(dim/2)

tmp=scfact*sqrt(x^2.+y^2.)/(fc)

mtf=tmp*0.0

i=where(tmp gt 0.)
mtf(i) = sin(!pi*f*tmp(i))/(!pi*f*tmp(i))
;if fc lt dim/2. then i=where(tmp gt 1.)
;if fc lt dim/2. then mtf(i) = 1.0

i=where(tmp eq 0.,count)

if count ne 0 then mtf(i) = 1.0

return,shift(mtf,(xc-fix(dim/2)),(yc-fix(dim/2)))

end

