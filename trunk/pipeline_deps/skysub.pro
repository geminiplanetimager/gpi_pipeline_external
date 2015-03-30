pro skysub, imgs, sky, imgsub,  $
            badpix = badpix, $
            scale = scale,  $
            save = save,  $
            silent = silent,$
			exptimes=exptimes,$
			nosave=nosave, $
			err_imgs=err_imgs, err_sky=err_sky,err_imgsub=err_imgsub
			

;+
;pro Skysub, imgs, sky, imgsub,  $
;            badpix = badpix, $
;            scale = scale,  $
;            save = save,  $
;            silent = silent
;
; PURPOSE:
; subtract a 2d sky frame (or bias frame) from a 3d array of images
; by default, it will delete the input images unless /save
;
; KEYWORDS:
;
; /scale  match the median of sky image with the object before subtraction
; /save   preserve the input images (default is to delete them)
; exptimes	array of exposure times for the images. If set, the sky frame
; 			is assumed to be for exptime=1 and is scaled appropriately
; 			before being subtracted from each image.
;
; NOTES
; to save memory when subtracting large/numerous images, 
; should call this routine with the TEMPORARY() function, i.e.:
;       IDL> skysub, temporary(imgs), skyimg, imgsub
;
; HISTORY: Written by M. Liu 07/07/94
; 09/27/95 MCL: added /scale
; 10/25/96 MCL: added badpix option
; 07/11/98 MCL: copies BADVAL in input obj & sky images to output images
; 02/16/01 MCL: significant improvements in memory & some in speed
;                * default is to delete input images unless /save
;                - grouped scalar operations
;                - more careful use of memory
;               now ignores BADVAL pix when using /scale
;               added & organized comments
; 06/11/01 MCL: if only a single frame is passed, then /save
; 2004-04-13 MDP:  Added exptimes option. Don't overwrite by default!
; 2005-05-03 MDP:  Now is /save by default, unless /nosave.
; 2005-10-13 MDP:  Removed horrible evil 'retall on any error' code!
;-


if ~(keyword_set(nosave)) then save=1

BADVAL =  -1.e6
if n_params() lt 2 then begin
	print,'skysub,imgs,sky,[imgsub],[scale],[save],[badpix=],[silent]'
	return
endif

; sanity check
imsz = size(imgs)
if (imsz(0) eq 2) then begin
    nn = 1 
    save = 1
endif else  begin
  nn = imsz(3)
endelse
  
sksz = size(sky)
if (imsz(1) ne sksz(1)) and (imsz(2) ne sksz(2)) then begin
	message,'sky and image frames not same dimension!',/info
	stop
endif


; define output array 
if keyword_set(save) then  $
  imgsub = fltarr(imsz(1),imsz(2),nn)


; get list of BADVAL pixels which apply to all images
; (from 2-d badpix mask and sky frame to be subtracted)
if keyword_set(badpix) then $
  wbad = where(badpix eq 0 or sky eq BADVAL, nbad) $
else $
  wbad = where(sky eq BADVAL, nbad)


; also remember bad pixels in original 3-d image stack
wbad3 = where(imgs eq BADVAL, nbad3)


; get scaling factors
if keyword_set(scale) then begin
	message, '* scaling sky frame to match object *', /info
	msky = median(sky)
endif


; loop to do subtraction
if not(keyword_set(silent)) then print, format = '($,"image:  ")'
for i = 0, nn-1 do begin

    if not(keyword_set(silent)) then print, format = '($,A," ")', strc(i)

    ; do subtraction 
    if keyword_set(scale) then begin
        wg =  where(imgs[*, *, i] ne BADVAL, ng)
        if (ng ne 0) then begin
            mobj = median((imgs[*, *, i])[wg])
            tmp = imgs[*, *, i] - (mobj/msky)*sky
        endif 
    endif else begin
		if keyword_set(exptimes) then begin
      		tmp = imgs[*, *, i] - sky*exptimes[i]
		endif else $
      		tmp = imgs[*, *, i] - sky
    endelse    
    ; mask BADVAL pixels originating from the sky frame
    if (nbad gt 0) then $
      tmp[wbad] =  BADVAL

    ; save the result (using clever IDL memory access)
    if keyword_set(save) then  $
      imgsub[0, 0, i] = temporary(tmp) $
    else  $
      imgs[0, 0, i] = temporary(tmp)

endfor
if not(keyword_set(silent)) then print


; propagate any BADVAL pixels from original 3-d stack
if keyword_set(save) then if (nbad3 gt 0) then imgsub(wbad3) = BADVAL $
else if (nbad3 gt 0) then imgs(wbad3) = BADVAL


end

