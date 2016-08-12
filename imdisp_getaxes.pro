;+
; NAME: imdisp_getaxes 
; PURPOSE:
;	create arrays with axes values for use with IMDISP 
; NOTES:
;
;	By default, if the axes span more than 120 units, they are divided
;	by 60 (i.e. converted into arcmin instead of arcsec)
;
; INPUTS:
; 	image		an image
; INPUT KEYWORDS:
; 	/center		if set, forces the axes origin to be the exact center, OR ELSE
; 	center=		2d array containing the pixel coords for the axes zero point
; 				(if not present, default is brightest pixel in the image)
; 	platescale=		platescale for the image. Can be 1 element or a 2d array
; 					for differing x and y scales. Default assumption is arcsec
; 	pixelscale=	synonym for platscale.
; 			
;
; 	/reversex	flip x axis sign (e.g. for RA)
;	/arcsec		display output in arcsec
;	/arcmin		display output in arcmin (i.e. divide by 60)
; 	
; OUTPUTS:
; 	xr=			xrange array suitable for using with imdisp
; 	yr=			yrange array suitable for using with imdisp
; 	x=			array of x axis coords converted to arcsec
; 	y=			array of y axis coords converted to arcsec
;
; HISTORY:
; 	Began 2005-10-22 17:24:51 by Marshall Perrin 
; 	2006-06-07		Added documentation
; 	2013-02-19		More careful handling of subpixel center for even/odd sizes.
;					Default pixel scale is now 1.0 if not specified. -MP
;-

PRO imdisp_getaxes,image,center=center,xr,yr,platescale=platescale0,$
		pixelscale=pixelscale, $
		x=x,y=y,reversex=reversex,$
		arcsec=arcsec, arcmin=arcmin

		if keyword_set(pixelscale) then platescale0=pixelscale

	if keyword_set(platescale0) then platescale=platescale0
	if ~(keyword_set(platescale)) then platescale=1.0 ; 0.0754
		;print,platescale
	if n_elements(platescale) lt 2 then platescale = [platescale,platescale]

	sz = size(image)
	if ~(keyword_set(center)) then begin
		whereismax,image[*,*,0],mx,my,/single,/silent
	endif else begin
		if n_elements(center) eq 1 then begin
			mx=(sz[1]-1)/2.
			my=(sz[2]-1)/2.
		endif else begin
			mx = center[0]
			my = center[1]
		endelse
 	endelse

	if keyword_set(reversex) then coeff=-1 else coeff=1
	x = (findgen(sz[1]) - mx) *platescale[0]*coeff
	y = (findgen(sz[2]) - my) *platescale[1]

	xr=[x[0],x[n_elements(x)-1]]
	yr=[y[0],y[n_elements(y)-1]]


	if (abs( xr[0]-xr[1]) gt 120 or abs( yr[0]-yr[1]) gt 120) $
			and ~(keyword_set(arcsec)) or keyword_set(arcmin ) then begin
		xr /=60
		yr /= 60
		x /=60
		y /= 60
	endif
			



end
