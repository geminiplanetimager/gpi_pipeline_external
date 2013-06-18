;+
;
; FUNCTION: fftshift
; PURPOSE:
;   shifts an image by (dx,dy) pixels using fourier transforms. 
; NOTES:
;   if
;   dx and dy are both nonzero, shifts in 2 dimensions. If dy is zero or not
;   present, the shift is performed in only the x dimension. 
;
;   Thus, calling fftshift,img,dx will perform a 1d shift of either 
;   a 1d series or a 2d image. Calling fftshift,img,dx,dy will perform
;   a 2d shift of a 2d image. 
;
;   based on the fourier transform shift theorem:
;   	f(x-dx,y-dy) <==> exp(-2pi i(u*dx+v*dy)) F(u,v)
;
; INPUTS:
; 	img		an image
; 	dx,dy	pixel shifts. Can be fractional. dy is optional.
; KEYWORDS:
; 	/edge_mask	Apply a windowing function around the edges to prevent
; 				abrupt ringing effects. 
; 	/silent		suppress informative statements
; 	/null		
; OUTPUTS:
; 	the shifted version of the image
;
; HISTORY:
;   2001-07-27 - Marshall Perrin
;   2012-12-18	 made edge mask optional, added /silent
;   			 docs cleaned up slightly
;
;-

FUNCTION fftshift,img,dx,dy,null=null,edge_mask=edge_mask, silent=silent
	sz=size(img)

	; calculate the  u indices. From IDL FFT help
	n=sz[1]
	u=findgen(n)
	n21=n/2+1
	if n mod 2 then c=1 else c=2
	u[n21] = n21-n+findgen(n21-c) 
								; FIXME should this coeff be 1 or 2?
								  ; I think it matters whether n is odd or not.
	u=u/n*dx; convert to frequency, mult by shift

	; replicate the indices to 2 dimensions
	uu=fltarr(sz[1],sz[2])
	for j=0,sz[2]-1 do uu[*,j]=u
			
	if keyword_set(dy) then begin
		; calcvlate the v indices. From IDL FFT help
		n=sz[2]
		v=findgen(n)
		n21=n/2+1
		if n mod 2 then c=1 else c=2
		v[n21] = n21-n+findgen(n21-c)
		v=v/n*dy ; convert to frequency, mult by shift
		
		; replicate the indices to 2 dimensions
		vv=fltarr(sz[1],sz[2])
		for j=0,sz[1]-1 do vv[j,*]=v

	endif else vv=0 ; no dy, so don't scale phases in v.
	
	if keyword_set(edge_mask) then begin
		if ~(keyword_set(silent)) then message,/info, "Now computing FFT shift with edge masking"
		if sz[0] eq 2 then begin
			edgemask=fltarr(sz[1],sz[2])+1
			edgemask[0,*]=0.
			edgemask[sz[1]-1,*]=.0
			edgemask[1,*]=.1
			edgemask[sz[1]-2,*]=.1
			edgemask[2,*]=.3
			edgemask[sz[1]-3,*]=.3
			edgemask[3,*]=.7
			edgemask[sz[1]-4,*]=.7
			edgemask[4,*]=.9
			edgemask[sz[1]-5,*]=.9

			edgemask[*,0]=0.
			edgemask[*,sz[2]-1]=.0
			edgemask[*,1]=edgemask[*,1]*.1
			edgemask[*,sz[2]-2]=edgemask[*,sz[2]-2]*.1
			edgemask[*,2]=edgemask[*,2]*.3
			edgemask[*,sz[2]-3]=edgemask[*,sz[2]-3]*.3
			edgemask[*,3]=edgemask[*,3]*.7
			edgemask[*,sz[2]-4]=edgemask[*,sz[2]-4]*.7
			edgemask[*,4]=edgemask[*,4]*.9
			edgemask[*,sz[2]-5]=edgemask[*,sz[2]-5]*.9
		endif else begin 
			edgemask=1
		endelse
		
		m=median(img) ; use the median not the mean because it's
					  ; better at grabbing the sky rather than the star
		f=fft((img-m)*edgemask) ; remove the median DC level before FFTing
		;stop
		i=complex(0,1)
		e=exp(-2*!pi*i*(uu+vv))

		if keyword_set(null) then begin
			en=e
			e[*,sz[2]/2]=0
			e[sz[1]/2,*]=0
		endif

		return,float(fft(e*f,/inverse))+m ; don't forget to add mean back in
	endif else begin
		; Simple shift, with no edge masking or DC offset removal.
		if ~(keyword_set(silent)) then message,/info, "Now computing FFT shift"
		f=fft((img)) 
		i=complex(0,1)
		e=exp(-2*!pi*i*(uu+vv))

		return,float(fft(e*f,/inverse)) 
	endelse 


end
