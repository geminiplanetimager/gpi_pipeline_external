;+
; NAME: mrecenter 
;
; PURPOSE:
;	Find precise center for a star, like RECENTER, but using
;	MPFITPEAK to do the peak fitting.
;
; INPUTS:
;	xold, yold		Initial guess for the center of some peak in an image
;
; KEYWORDS:
; 	/moffat, /lorentz, /gauss		what type of function to fit?
; 									see MPFIT2DPEAK for details
; 	/silent			suppress printed output
; 	/nodisp			suppress image display on screen
; 	box=			box size to use in peak fitting. Default=20 pixels
; 	/all			don't use a box subimage, fit to the whole image
; 	/stop			have a 'stop' breakpoint after execution (For debugging/test)
;
; OUTPUTS:
;	xnew, ynew		Improved estimates for the center of that peak
;	width			Estimate for FWHM of that peak
;
; HISTORY:
; 	Began 2006-08-03 16:41:26 by Marshall Perrin 
;-

PRO mrecenter, image, xold, yold, xnew, ynew, box=box,stop=stop,$
	silent=silent,nodisp=nodisp,$ ; for drop-in replacement to recenter
	moffat=moffat, lorentz=lorentz, gauss=gauss, width=width,all=all

	compile_opt defint32, strictarr, logical_predicate

	if ~(keyword_set(moffat)) and ~(keyword_set(lorentz)) and ~(keyword_set(gauss)) then moffat=1

	if ~(keyword_set(box)) then box = 20.0

	if keyword_set(all) then begin
		subim=image 
		x0=0
		y0=0
	endif else  subim = imcut(image,box, xold, yold, x0=x0,y0=y0,/silent)

	sz = size(subim)
	x = indgen(sz[1])+x0
	y = indgen(sz[2])+y0

	fwhm = mfwhm(subim)
	estimates = [0.0, max(subim), fwhm,fwhm, xold, yold, 0.0, 0.0]

	;yfit = mpfit2dpeak(subim, params, X, Y,estimates=estimates,/lor)
	yfit = mpfit2dpeak(subim, params, X, Y,estimates=estimates,$
		moffat=moffat, lorentz=lorentz, gauss=gauss)
	; lorentzian is more peaky, works better
	; moffat is maybe a bit still better, but not sure if it'll make
	; an appreciable difference. Oh, what the heck. 

	;atv, [[[subim]],[[yfit]]],/bl
	xnew = params[4]
	ynew = params[5]
	if arg_present(width) then width = 2.35*(params[2]+params[3])/2

if ~(keyword_set(silent)) then begin
	print, "  OLD      NEW     DIFF"
	print, xold, xnew, xold-xnew
	print, yold, ynew, yold-ynew
endif

	if keyword_set(stop) then stop

end
