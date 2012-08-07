;+
; NAME: polcalplots
; PIPELINE PRIMITIVE DESCRIPTION: plots using the polarization spot calibration
;
;	
;	
;
; INPUTS: data-cube
;
;
; KEYWORDS:
;	/Save	Set to 1 to save the output image to a disk file.
; KEYWORDS:
; GEM/GPI KEYWORDS:
; DRP KEYWORDS: 
;
; OUTPUTS:  
;
; PIPELINE COMMENT: plots using the wavelength solution
; PIPELINE ARGUMENT: Name="ps_figure" Type="int" Range="[0,500]" Default="2" Desc="1-500: choose # of saved fig suffix name, 0: no ps figure "
; PIPELINE ARGUMENT: Name="suffix" Type="string"  Default="-spec" Desc="Enter output suffix (fits)"
; PIPELINE ORDER: 2.51
; PIPELINE TYPE: ALL-SPEC
; PIPELINE SEQUENCE: 
;
; HISTORY:
; 	
;   JM 2010-03 : created module.
;- 

function polcalplots, DataSet, Modules, Backbone
common PIP
COMMON APP_CONSTANTS

primitive_version= '$Id: polcalplots.pro 96 2010-10-20 13:47:13Z maire $' ; get version from subversion to store in header history
	;getmyname, functionname
	  @__start_primitive

   	; save starting time
   	T = systime(1)

  	polcal=*(dataset.currframe[0])

        band=strcompress(sxpar( *(dataset.headers[numfile]), 'FILTER',  COUNT=cc),/rem)
        if cc eq 1 then begin
          cwv=get_cwv(band)
          CommonWavVect=cwv.CommonWavVect
          lambda=cwv.lambda
          lambdamin=CommonWavVect[0]
          lambdamax=CommonWavVect[1]
          NLam=CommonWavVect[2]
        endif else begin
          NLam=0
          lambda=(indgen((size(main_image_stack))[3]))
        endelse
	
	thisModuleIndex = Backbone->GetCurrentModuleIndex()


;stop
 ;tilt of pol spots:
tilt3= (180./!dpi)*atan(reform(polcal[1,*,*,1]-polcal[1,*,*,0])/reform(polcal[0,*,*,1]-polcal[0,*,*,0]))
;spectral spacing:
w3=((reform(polcal[0,*,*,0])-shift(reform(polcal[0,*,*,0]),0,1)    ))
w3m=mean(w3)
;;dispersion distance
;dp=((sqrt( (reform(polcal[0,*,*,0])-shift(reform(polcal[0,*,*,0]),2,-1))^2 + (reform(reform(polcal[1,*,*,0]))-shift(reform(polcal[1,*,*,0]),2,-1))^2) ))
;dpm=mean(dp)
;;micro-pupil pattern:
;P3=sqrt(abs((dp-w3)/w3))
;;lenslet tilt:
;theta=(180./!dpi)*atan(1./P3)


tx=900
ty=900
plotc, tilt3, 3, tx,ty,'micro-lens','micro-lens','Tilt of spectra [degrees]'
plotc, w3, 4, tx,ty,'micro-lens','micro-lens','w [detector pixel]'
;plotc, p3, 5, tx,ty,'micro-lens','micro-lens','P micro-pupil pattern'
;plotc, ((180./!dpi)*p3 mod 180.), 6, tx,ty,'micro-lens','micro-lens','P micro-pupil pattern'
;plotc, dp, 7, tx,ty,'micro-lens','micro-lens','distance between spectra [detector pixel]'
;plotc, theta, 8, tx,ty,'micro-lens','micro-lens','distance between spectra [detector pixel]'
;;plotc, reform(wavcal[*,*,3]), 9, tx,ty,'micro-lens','micro-lens','coeef wavcal'


lamzem=readfits(getenv('GPI_IFS_DIR')+path_sep()+'dst'+path_sep()+'zemdispLam'+band+'.fits')

zem=readfits(getenv('GPI_IFS_DIR')+path_sep()+'dst'+path_sep()+'zemdisp_pol'+band+'.fits')
zemX=zem[*,*,*,0]+1024.
zemY=zem[*,*,*,1]+1024.

zemX2=zem[*,*,*,2]+1024.
zemY2=zem[*,*,*,3]+1024.
;shiftx=5 ;not needed
;shifty=-5  
;zemX0=zemX[*,*,0] ;zem locations at lambda_min ;shift(zemX[*,*,0], shiftx,shifty)
;zemY0=zemY[*,*,0] ;zem locations at lambda_min ;shift(zemY[*,*,0], shiftx,shifty)
;    zemtheoX=zemx
;      zemtheoY=zemy
      zemtilt=(180./!dpi)*atan((zemY[*,*,18]-zemY2[*,*,18])/(zemX[*,*,18]-zemX2[*,*,18]))
      zemw=(shift(zemy[*,*,0],0,1)-zemy[*,*,0])
;      dpzem=sqrt( (zemx[*,*,0]-shift(zemx[*,*,0],-1,2))^2 + ((zemy[*,*,0]-shift(zemy[*,*,0],-1,2))^2)) 
;      zemP=sqrt(abs((dpzem[*,*]-zemw[*,*])/zemw[*,*]))      
;      theta_zem=(180./!dpi)*atan(1./zemP)
;      pl=30
;      coefzem=(lambda[pl]-lambda[0])/(sqrt((zemX[*,*,pl]-zemX[*,*,0])^2+(zemY[*,*,pl]-zemY[*,*,0])^2))
;      
      indout=where(~finite(tilt3),cout)
      zemtilt2=rotate(transpose(zemtilt),1)
      zemtilt2[indout]= !VALUES.F_NAN
      
      zemw2=rotate(transpose(zemw),1)
      zemw2[indout]= !VALUES.F_NAN
      
;      dpzem2=rotate(transpose(dpzem),2)
;      dpzem2[indout]= !VALUES.F_NAN
;      
;      zemp2=rotate(transpose(zemp),2)
;      zemp2[indout]= !VALUES.F_NAN
;      
;      theta_zem2=rotate(transpose(theta_zem),2)
;      theta_zem2[indout]= !VALUES.F_NAN
;      
;      coefzem2=rotate(transpose(coefzem),2)
;      coefzem2[indout]= !VALUES.F_NAN
 ;plotc, zemtilt, 13, tx,ty,'Zemax micro-lens','Zemax  micro-lens','Tilt of spectra [degrees]'
 plotc, tilt3, 3, tx,ty,'micro-lens','micro-lens','Tilt of polarization spots [degrees]',valmin=-24.,valmax=-22. ;valmin=valmin1,valmax=valmax1
  plotc, zemtilt2, 23, tx,ty,'Zemax micro-lens','Zemax  micro-lens','Tilt of polarization spots [degrees]',valmin=-19.,valmax=-18. ;valmin=valmin1,valmax=valmax1

  plotc, zemw2, 24, tx,ty,'Zemax micro-lens','Zemax micro-lens','spectral spacing w [detector pixel]',valmin=4.7,valmax=4.9;,valmin=valmin2,valmax=valmax2
  plotc, w3, 4, tx,ty,'micro-lens','micro-lens','spectral spacing w [detector pixel]',valmin=4.4,valmax=4.6;,valmin=valmin2,valmax=valmax2

;  plotc, dpzem2, 27, tx,ty,'Zemax micro-lens','Zemax micro-lens','Dispersion distance [detector pixel]',valmin=valmin2,valmax=valmax2
;  plotc, dp, 7, tx,ty,'micro-lens','micro-lens','Dispersion distance [detector pixel]',valmin=valmin2,valmax=valmax2
;  
;  plotc, zemp2, 25, tx,ty,'Zemax micro-lens','Zemax micro-lens','P micro-pupil pattern';,valmin=valmin3,valmax=valmax3
;    plotc, p3, 5, tx,ty,'micro-lens','micro-lens','P micro-pupil pattern';,valmin=valmin3,valmax=valmax3
;  
;plotc, ((180./!dpi)*zemp2 mod 180.), 26, tx,ty,'Zemax micro-lens','Zemax micro-lens','P micro-pupil pattern',valmin=valmin4,valmax=valmax4
;plotc, ((180./!dpi)*p3 mod 180.), 6, tx,ty,'micro-lens','micro-lens','P micro-pupil pattern',valmin=valmin4,valmax=valmax4
;
;plotc, theta_zem2, 28, tx,ty,'Zemax micro-lens','Zemax micro-lens','lenslet tilt Theta [degree]',valmin=valmin6,valmax=valmax6
;plotc, theta, 8, tx,ty,'micro-lens','micro-lens','lenslet tilt Theta [degree]'  ,valmin=21.,valmax=22.
;
;
;plotc, coefzem2, 29, tx,ty,'Zemax micro-lens','Zemax micro-lens','slope of the linear dispersion [microms / detector pixel]',valmin=valmin9,valmax=valmax9
;plotc, reform(wavcal[*,*,3]), 9, tx,ty,'micro-lens','micro-lens','slope of the linear dispersion [microms / detector pixel]',valmin=valmin9,valmax=valmax9


hdr=*(dataset.headers[numfile])



	thisModuleIndex = Backbone->GetCurrentModuleIndex()
  

 
;drpPushCallStack, functionName

return, ok


end