;+
; NAME:  imdisp_with_contours
; PURPOSE:
; 		imdisp an image, and overplot contours using subtle colors a la Tufte
;
; INPUTS:
; 	image
; REQUIRED KEYWORDS: (unless you use the pixelscale= keyword)
;	xrange,yrange	X and Y range for whole array
;	xvals,yvals		X and Y for each pixel
;	contourlevels	array giving levels for contours 
;
; OPTIONAL KEYWORDS
; 	/axis			draw axis
; 	
;	/alog, /asinh	what scaling to use?
;	/noscale		display image exactly as provided
;	min=, max=		min and max for scaling
;	/negative		black-on-white display 
;
;	contourcolors	colors for contours
; 	contour_step	offset between regular plot colors and contour colors
; 					default = 40
; 	contour_image	plot contours of this image instead of main image
; 					(useful for overplots or smoothed versions)
; 	contoursmooth=  median-smooth the image by this much before contouring.
;	/nocontours		don't plot contours! (kind of silly give the routine name, 
;					but sometimes useful...)
;
;	force_contourcolors=	numerical color value override (can be array)
;	contourstring=			string to call fsc_color with 
;							(useful for pstscript due to IDL brokenness
;							requiring fsc_color call *after* image display)
;
;
; OUTPUTS:
; 	out_position	output position used to plot (see imdisp.pro)
;
; HISTORY:
; 	Began 2007-01-04 14:59:55 by Marshall Perrin 
; 	2007-06025	all the bugs out.
; 	2009-04-16  added ccols_out keyword
;-

PRO imdisp_with_contours, image, $
	; outputs from imdisp_getaxes:
	xrange=xr,yrange=yr,xvals=xvals, yvals=yvals,$
	contourlevels=contourlevels, contourcolors=contourcolors,$
	contoursmooth=contoursmooth, $
	out_position=out_pos,negative=negative,$
	_extra=_extra,axis=axis,$
	alog=alog, asinh=asinh, max=max,min=min, $ ;,beta=beta,$
	;autoasinh = autoasinh, $
	nocontours=nocontours, $
	contour_image = contour_image,$
	contour_step=constep, $
	pixelscale=pixelscale, center=center, $
	ccols_out =ccols_out, $
	force_contourcolors=force_contourcolors, contourstring=contourstring, $
	noscale=noscale

;---- some defaults

	sz = size(image)
	if ~(keyword_set(xvals)) then xvals = indgen(sz[1])
	if ~(keyword_set(yvals)) then yvals = indgen(sz[2])

	if keyword_set(pixelscale) and (~(keyword_set(xr)) or ~(keyword_set(yr))) then begin
		imdisp_getaxes, image, center=center, xr,yr,platescale=pixelscale,$
			x=xvals, y=yvals,/reversex
	endif

	; automatically determine good levels to use
	if ~(keyword_set(contourlevels)) then begin
		mx = max(image,/nan)
		contourlevels = reverse( mx / [10, 31.6, 100,316] )
		print, contourlevels
	endif
    if ~(keyword_set(constep)) then constep=40

;--- the imdisp part ----------

	if keyword_set(asinh) then begin
		scaled_image = asinhscl(image,min=min,max=max,beta=beta) ;,auto=autoasinh )
		contourcolors = asinhscl(contourlevels,  min=min,max=max,beta=beta) ;,beta=beta)
	endif else if keyword_set(alog) then begin
		scaled_image = alogscale(image,min,max)
		contourcolors =alogscale(contourlevels, min,max)

	endif else if keyword_set(noscale) then begin
		; really no scaling!
		scaled_image=image
		contourcolors = contourlevels
	
	endif else begin
		; bytscl suffers tremendously from roundoff error on integer input, so
		; make everything a float before calling bytscl!
		scaled_image = bytscl(image*1.0,max=max,min=min)
		contourcolors =bytscl(contourlevels*1.0,  min=min,max=max)
	endelse
	
	imdisp,scaled_image, out_pos=out_pos,_extra=_extra,xr=xr,yr=yr,/axis,negative=negative,/keep,/noscale

;--- the contours part ----------

if ~(keyword_set(nocontours)) then begin

    if keyword_set(negative) then begin
		contourcolors = 255L-contourcolors
    	contourcolors = (contourcolors>constep ) - constep ; in negative mode, make the contours DARKER
	endif else begin
    	contourcolors = (contourcolors< (255-constep) ) + constep ; in regular mode, make them LIGHTER
	endelse
    if !d.name eq 'PS' then ccols = contourcolors else ccols = contourcolors*256L+contourcolors*65536L+contourcolors
    loadct,0


    if keyword_set(contoursmooth) then begin
		contour_image = median(image,contoursmooth) 
		wnan = where(~finite(image),nanct)
		if nanct gt 0 then contour_image[wnan]=!values.f_nan
	endif

	; overrides for the contour colors
	if keyword_set(force_contourcolors) then ccols=force_contourcolors
	if keyword_set(contourstring) then ccols = fsc_color(contourstring)

	if keyword_set(contour_image) then begin
		contour, contour_image, xvals, yvals, /over,position=out_pos,$
			levels=contourlevels, c_colors=ccols, _extra=_extra; contourcolors
	endif else begin
		contour, image, xvals, yvals, /over,position=out_pos,$
			levels=contourlevels, c_colors=ccols, _extra=_extra; contourcolors
	endelse
	if arg_present(ccols_out) then ccols_out=ccols

endif

end
