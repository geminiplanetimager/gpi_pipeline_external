;From tmurphy@mop.caltech.edu Tue Jun 30 09:48:54 1998
;From: Tom Murphy <tmurphy@mop.caltech.edu>
;Date: Tue, 30 Jun 1998 09:50:24 -0700
;To: larkin@astro.ucla.edu
;Subject: noise & fixpix
;X-Status: ;;
;;
;
;Hi James,
;
;The routines noise.pro and fixpix.pro follow in separate e-mails.
;The noise.pro routine is universal, meaning that it will work on any
;;array input.  But fixpix currently has some d80-dependent stuff.
;
;Fixpix is mostly independent of the array characteristics, evaluating
;the noise on the fly (with noise.pro).  But sometimes you run into a
;feature (a star, spectrum, skyline, fringe, etc.).  In these cases, the
;range of pixel values is dominated by the shape of the feature itself,;
;and I resort to a calculation of the expected noise based on the array
;characteristics (read plus shot).  Before doing this, I check to see
;if the data has been run through D80P (quad float and mparm2 adjust).
;If not, I divide by mparm2 (usu. 16).  So these are obviously d80-
;dependent things, but not difficult to replace with your own array
;characteristics.  Just look for the string mparm to show you where to
;change things.
;
;I also use two different rejection criteria for pixels in the Gaussian
;background and those thought to exist on a (positive or negetive) feature.
;Currently, I use 5-sigma on the background, and 3 sigma on the feature.
;This may sound backwards, but my sigma estimates on features are on the
;high side, so that rejection here is unlikely.
;
;Anyway, if you give me your FAX number, I'll also send you some supplemental
;info that will better enable you to modify fixpix to suit your needs.  As it
;is, it works quite well on a variety of data types: subtracted spectra, raw
;spectra, G-stars, images, etc.  So I think if you have a comparable array
;noise estimate, you will end up with an equally robust algorithm.
;
;Tom
;
;
;
;
;
;----- End Included Message -----
;
;
;----- Begin Included Message -----
;
;From tmurphy@mop.caltech.edu Tue Jun 30 10:32:34 1998
;From: Tom Murphy <tmurphy@mop.caltech.edu>
;Date: Tue, 30 Jun 1998 10:34:05 -0700
;To: larkin@astro.ucla.edu
;Subject: noise.pro
;X-Status: 
;;
;
;(I first tried larkin@astro.caltech.edu a half hour ago, and went away
;thinking it was done!  opps.)
;
function noise,im

; function to find the std dev of a block of pixels

; sort the array

imsort = im(sort(im))

; find the number of pixels

npix = (size(imsort))[1]

; pick the integer that is closest to 2/3 of the number of pixels
; this is an integer

num = 2*npix/3 + 1

; find the noise as a the minimum difference between
; the averages of two pairs of points approximately 2/3 of the 
; entire array length apart

noise = min(imsort(num-1:npix-1)+imsort(num-2:npix-2)-imsort(0:npix-num)-imsort(1:npix-num+1),/nan)/2.0

; now correct the noise for the fact that it's not exactly over 2/3 of the 
; elements

noise = noise / (1.9348+(float(num-2)/float(npix) - 0.66667)/.249891)

; now correct for the finite number of pixels

if npix gt 100 then $ 
	noise = noise/(1.0-0.6141*exp(-0.5471*alog(npix))) $
	else $
	noise = noise/(1.0-0.2223*exp(-0.3294*alog(npix)))

return, noise

end

;+
;function ns_fixpix, ims,h,maskarr,iter=iternum,silent=silent
;
; Clean up hot and cold pixels in images.
; Usage: imgs_fixed = ns_fixpix(imgs)
;
; KEYWORDS:
; 	iter		set number of cleaning iterations. default=3
; 	/NaN		Clean up NaNs in addition to hot and cold pixels.
; 	/silent		Suppress printed output
; 	
;
; HISTORY:
;
;	Originally by Tom Murphy, Caltech. June 1998
;	Modified by Marshall Perrin at UC Berkeley. June 2003
;
;	/NaN added 2003 July 1. M Perrin
;-
function ns_fixpix, im,h,maskarr,iter=iternum,silent=silent,NaN=NaN

; INITIALIZATIONS

; in case of an error return to the calling prog

;on_error, 2

; create a new image to be the cleaned up output
imnew = im

; set the maskarray to zero
maskarr = im*0.0
; determine the median value of the entire image
medimg = median(im, /even)
sz = size(imnew)
if sz[0] eq 3 then maximnum=sz[3]-1 else  maximnum=0 

; if the number of iterations is not explicitly declared set it to 3
if (not keyword_set(iternum)) then iternum =3


for imagenum=0L,maximnum do begin
	if maximnum gt 0 and not keyword_set(silent) then print,"Cleaning image "+strc(imagenum+1)+" of "+strc(maximnum+1)

; set mparm to 1 as a default
mparm = 1
iter=1

qf = (not keyword_set(silent))

; READ THE HEADER FILE

;if (not sxpar(h,"D80P")) then begin
;    if qf then print, "Warning: This image has not been run through D80P"
;    mparm = sxpar(h,"MPARM201")
;    if (mparm eq 0) then begin
;        if qf then print, "Warning: No MPARM found in the header file. Setting to default."
;        mparm = 1
;    endif
;    if qf then print, format='("mparm = ",i0)',mparm
;endif	


; DETERMINE THE AVERAGE NOISE IN THE ENTIRE IMAGE

; scan through 5x5 boxes to determine some average estimate of the
; image noise
sz = size(im)
xs = sz[1]
ys = sz[2]
;sigma = fltarr((xs-1)/5*(ys-1)/5)
sigma = fltarr((xs)/5*(ys)/5)

n=float(0)
;for i=2,xs-3,5 do begin
for i=2L,xs-4,5 do BEGIN  ; changed for a cropped image AMG

    for j=2,ys-3,5 do begin

        tmp=im[i-2:i+2,j-2:j+2,imagenum]

        srt=tmp(sort(tmp))

        sigma[n] = (min(srt(16:24)+srt(15:23)-srt(0:8)-srt(1:9))/2.0)/1.540

        n = n+1
    endfor

endfor

; define the median value of the sigmas, and the sigma of the sigmas

medsig = median(sigma,/even)
sigsig = noise(sigma) 

print, "Medsig, sigsig:", medsig, sigsig
; BEGIN SCANNING FOR HOT & COLD PIXELS

; start at (4,4) (so that any pixel in the box can have a 5x5 
; square centered on it.  find the hottest and coldest pixels

; loop through several iterations of replactments

;print,hoti,coldi
	; find NaN pixels
	tmp  = intarr(5,5) & tmp[1:3,1:3]=1 & tmp[2,2]=2
	ind8 = where(tmp eq 1)
	ind16 = where(tmp eq 0)
	
for iter=1L,iternum do begin

	if keyword_set(NaN) then begin
		wnan = where(finite(imnew[*,*,imagenum],/NaN),nanct)
		if (nanct gt 0) then begin
			sz = size(imnew)
			whereis,imnew[*,*,imagenum],wnan,xn,yn
			x0 = xn-2>0
			x1 = xn+2<(sz[1]-1)
			y0 = yn-2>0
			y1 = yn+2<(sz[2]-1)
			for i = 0L,nanct-1 do begin
				imnew[xn[i],yn[i],imagenum] = median(imnew[x0[i]:x1[i],y0[i]:y1[i],imagenum])
			endfor
	    	if qf then print, format='("     ",i5,"  NaN pixels in iter. ",i0)',nanct,iter
		endif
	endif



	; find hot and cold pixels
    hotcnt = 0                  ; variables to count the replacements of hot and cold pixels
    coldcnt = 0
    addon = ((iter+1) mod 2) *2

	sigcut = medsig + 2.0*sigsig

    for i=4L + addon,xs-5,5 do begin

        for j=4 + addon,ys-5, 5 do begin

            box = imnew(i-2:i+2,j-2:j+2,imagenum)
            
            hotval = max(box,hoti) ; hoti is returned from max = index of hot value in box
            hotx = hoti mod 5 + i-2 ; coords in original image
            hoty = hoti/5 + j-2

            coldval = min(box,coldi)
            coldx = coldi mod 5 + i-2
            coldy = coldi/5 + j-2

                                ; begin the decision process for the hottest pixel
            
            hot = imnew(hotx-2:hotx+2,hoty-2:hoty+2,imagenum)
			med8 = median(hot[ind8],/even)

            ;med8 = median([hot(1:3,1),hot(1:3,3),hot(1,2),hot(3,2)],/even) ; 8-pixel box centered on hot pixel
			med16 = median(hot[ind16],/even)
            ;med16 = median([hot(0:4,0),hot(0:4,4),transpose(hot(0,1:3)),transpose(hot(4,1:3))],/even)
   
            srt = hot(sort(hot))
            sig = (min(srt(16:24)+srt(15:23)-srt(0:8)-srt(1:9))/2.0)/1.540

                                ; decide from the noise in the box if we 
                                ; are on a feature or a gaussian background

            if sig gt sigcut then $
              sig = max([sigcut,sqrt(5.0+.210*abs(med8))])
              ;sig = max([medsig+2.0*sigsig,sqrt(5.0+.210*abs(med8)/mparm)*mparm])

                                ; decide whether to replace pixel

            if ((med8+2.0*med16)/3.0 -medimg) gt 2.0*sig then begin
                if (imnew(hotx,hoty,imagenum) gt (2.0*med8-medimg+3.0*sig)) then begin
					;maskb[hotx,hoty]=1
                    imnew(hotx,hoty,imagenum) = med8
                    if keyword_set(domask) then maskarr(hotx,hoty,imagenum)=1
                    hotcnt++
                endif
            endif else begin
                if ((imnew(hotx,hoty,imagenum)-(med8+2.0*med16)/3) gt 5.0*sig) then begin
					;maskb[hotx,hoty]=1
                    imnew(hotx,hoty,imagenum) = med8
                    if keyword_set(domask) then maskarr(hotx,hoty,imagenum)=1
                    hotcnt++
                endif
            endelse

                                ; begin the decision process for the coldest pixel
            
            cld = imnew(coldx-2:coldx+2,coldy-2:coldy+2,imagenum)
            ;med8 = median([cld(1:3,1),cld(1:3,3),cld(1,2),cld(3,2)],/even)
			med8 = median(cld[ind8],/even)
			med16 = median(cld[ind16],/even)
            ;med16 = median([cld(0:4,0),cld(0:4,4),transpose(cld(0,1:3)),transpose(cld(4,1:3))],/even)
            srt = cld(sort(cld))
            sig = (min(srt(16:24)+srt(15:23)-srt(0:8)-srt(1:9))/3.08)

                                ; decide from the noise in the box if we 
                                ; are on a feature or a gaussian background

            if sig gt sigcut then $
              sig = max([sigcut,sqrt(5.0+.210*abs(med8))])

                                ; decide whether to replace pixel

            if ((med8+2.0*med16)/3.0 -medimg) lt -2.0*sig then begin
                if (imnew(coldx,coldy,imagenum) lt (2.0*med8-medimg-3.0*sig)) then begin
                    imnew(coldx,coldy,imagenum) = med8
                    if keyword_set(domask) then maskarr(coldx,coldy,imagenum)=1
                    coldcnt++
                endif
            endif else begin
                if ((imnew(coldx,coldy,imagenum)-(med8+2.0*med16)/3) lt -5.0*sig) then begin
                    imnew(coldx,coldy,imagenum) = med8
                    if keyword_set(domask) then maskarr(coldx,coldy)=-1
                    coldcnt++
                endif
            endelse
            

        endfor		
        
    endfor


    if qf then print, format='(i5,i5,"  hot and cold pixels in iter. ",i0)',hotcnt,coldcnt,iter
endfor

; end of for loop for image number
endfor 
; add keyword to header file
;sxaddpar,h,'FIX_ITER',iternum

return, imnew

end




;----- End Included Message -----


