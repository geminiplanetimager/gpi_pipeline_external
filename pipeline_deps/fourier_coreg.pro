pro fourier_coreg,im1,im2,invcps

IMG1 = im1
IMG2 = im2

IMG1_ORIG = IMG1
IMG2_ORIG = IMG2
SZ_IMG1 = size(IMG1)
SZ_IMG2 = size(IMG2)
SZ_1 = MAX([SZ_IMG1[1],SZ_IMG2[1]])
SZ_2 = MAX([SZ_IMG1[2],SZ_IMG2[2]])
IF ( (SZ_IMG1[1] NE SZ_1) OR (SZ_IMG1[2] NE SZ_2)) THEN BEGIN
   IMG1_NEW = fltarr(SZ_1,SZ_2)
   xstart = (SZ_1 - SZ_IMG1[1])/2
   ystart = (SZ_2 - SZ_IMG1[2])/2
   IMG1_NEW[xstart:xstart+SZ_IMG1[1]-1, $
                ystart:ystart+SZ_IMG1[2]-1] = IMG1
   IMG1 = IMG1_NEW
   IMG1_ORIG = IMG1
ENDIF
IF ( (SZ_IMG2[1] NE SZ_1) OR (SZ_IMG2[2] NE SZ_2)) THEN BEGIN
   IMG2_NEW = fltarr(SZ_1,SZ_2)
   xstart = (SZ_1 - SZ_IMG2[1])/2
   ystart = (SZ_2 - SZ_IMG2[2])/2
   IMG2_NEW[xstart:xstart+SZ_IMG2[1]-1, $
            ystart:ystart+SZ_IMG2[2]-1] = IMG2
   IMG2 = IMG2_NEW
   IMG2_ORIG=IMG2
ENDIF


IMG1 = HANNING(SZ_1,SZ_2)*IMG1
IMG2 = HANNING(SZ_1,SZ_2)*IMG2
FFT1 = FFT(IMG1,-1)
FFT2 = FFT(IMG2,-1)
CPS = FFT1 * CONJ(FFT2) ;;/ ABS(( FFT1) * (FFT2))
INVCPS = FFT(CPS,1)

INVCPS = ROTATE(SHIFT(REAL_PART(INVCPS),SIZE(INVCPS,/dim)/2.),2)

end
