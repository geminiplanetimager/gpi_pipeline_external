pro fourier_coreg,im1,im2,out,window=window,findshift=findshift
;+
; NAME:
;      FOURIER_COREG
;
; PURPOSE:
;      Return the parallactic angle of a source in degrees.
;
; CALLING SEQUENCE:
;      fourier_coreg,im1,im2,out,[/window,/findshift]
;
; INPUTS:
;      im1 - Image
;      im2 - Template
;
; KEYWORD PARAMETERS:
;      /window - Apply hanning window.
;      /findshift - Return offset between images
;
; OUTPUTS:
;      out - Either offset between images (if /findshift is set) or 
;      the convolution of the images.
;
; COMMON BLOCKS:
;      None.
;
; NOTES:
;      If input images have odd and even numbers of pixels, respectively, 
;      you will get half pixel shifts. 
;
; MODIFICATION HISTORY:
;      Written by Dmitry Savransky based partially on code by Lisa Poyneer.
;-

;;preserve original inputs
img1 = im1
img2 = im2

;;get sizes and types
typ = max([size(img1,/type),size(img2,/type)])
sz_img1 = size(img1,/dim)
sz_img2 = size(img2,/dim)
if total((sz_img1 mod 2.) - (sz_img2 mod 2.)) ne 0 then begin
   message,'WARNING: Input images have non-matching odd and even numbers of pixels in at least one dimension.',/continue
   message,'This will result in a half pixel shift.' /continue
endif

;;define final sizes and calculate centers
sz_fin = max([[sz_img1],[sz_img2]],dim=2)
cent1 = sz_img1/2. - 0.5
cent2 = sz_img2/2. - 0.5
cent = sz_fin/2. - 0.5

;;pad and cast images as necessary
tmp = make_array(sz_fin,type=typ)
tmp[cent[0]-cent1[0]:cent[0]+cent1[0],cent[1]-cent1[1]:cent[1]+cent1[1]] = img1
img1 = tmp
tmp = make_array(sz_fin,type=typ)
tmp[cent[0]-cent2[0]:cent[0]+cent2[0],cent[1]-cent2[1]:cent[1]+cent2[1]] = img2
img2 = tmp

if keyword_set(window) then begin
   wind = hanning(sz_fin[0],sz_fin[1])
   img1 = wind*img1
   img2 = wind*img2
endif

;;do transforms and convolve
ft1 = fft(img1)
ft2 = fft(img2)
cps = ft1 * conj(ft2)

if keyword_set(findshift) then begin
   invcps = real_part(fft(cps,1))
   
   ;;find largest value & surrounding indices
   val = max(invcps, ind)
   inds = array_indices(invcps,ind)  
   minds = inds - 1
   pinds = inds + 1

   ;;any +1 values get set to 0, any -1 values get set to n-1
   bad = where(pinds gt sz_fin-1,ct)
   if ct gt 0 then pinds[bad] = 0
   bad = where(minds lt 0,ct)
   if ct gt 0 then minds[bad] = sz_fin[1]-1

   out = make_array(2,type=typ)
   out[0] = -0.5*(invcps[minds[0],inds[1]] - invcps[pinds[0],inds[1]])/$
         (invcps[minds[0],inds[1]] + invcps[pinds[0],inds[1]] - 2*invcps[inds[0],inds[1]])
   out[1] = -0.5*(invcps[inds[0],minds[1]] - invcps[inds[0],pinds[1]])/$
         (invcps[inds[0],minds[1]] + invcps[inds[0],pinds[1]] - 2*invcps[inds[0],inds[1]])
   out += sz_fin-inds

   ;;wrap accordingly
   b = where(out gt sz_fin/2.,ct)
   if ct gt 0 then out[b] -= sz_fin[b]
   b = where(out lt -sz_fin/2.,ct)
   if ct gt 0 then out[b] += sz_fin[b]

endif else begin
   ;;generate pre-shifted grids (zero always comes first)
   ;;odd  # of pix: 0,1,...(n-1)/2,-(n-1)/2...-1
   ;;even # of pix: 0,1,...n/2-1,-n/2,-n/2+1...-1
   if (sz_fin[0] mod 2) eq 1 then $
      xs = [indgen((sz_fin[0]+1)/2,type=typ), indgen((sz_fin[0]-1)/2,type=typ) - (sz_fin[0]-1)/2] else $
         xs = [indgen((sz_fin[0]/2),type=typ), indgen((sz_fin[0]/2),type=typ) - (sz_fin[0]/2)]
   xs = xs  # (make_array(sz_fin[1],type=typ) + 1)
   if (sz_fin[1] mod 2) eq 1 then $
      ys = [indgen((sz_fin[1]+1)/2,type=typ), indgen((sz_fin[1]-1)/2,type=typ) - (sz_fin[1]-1)/2] else $
         ys = [indgen((sz_fin[1]/2),type=typ), indgen((sz_fin[1]/2),type=typ) - (sz_fin[1]/2)]
   ys = ys  ## (make_array(sz_fin[0],type=typ) + 1)

   ;;half image shift is a pi shift:
   shft = exp(complex(0,1)*!dpi*(xs+ys))
   out = real_part(fft(cps*shft,1))
endelse

end
