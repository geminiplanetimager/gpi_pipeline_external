;+
;
; NAME:
;       ATV
; 
; PURPOSE: 
;       Interactive display of 2-D or 3-D images.
;
; CATEGORY: 
;       Image display.
;
; CALLING SEQUENCE:
;       atv [,array_name OR fits_file] [,min = min_value] [,max=max_value] 
;           [,/linear] [,/log] [,/histeq] [,/asinh] [,/block]
;           [,/align] [,/stretch] [,header = header]
;
; REQUIRED INPUTS:
;       None.  If atv is run with no inputs, the window widgets
;       are realized and images can subsequently be passed to atv
;       from the command line or from the pull-down file menu.
;
; OPTIONAL INPUTS:
;       array_name: a 2-D or 3-D data array to display
;          OR
;       fits_file:  a fits file name, enclosed in single quotes
;
; KEYWORDS:
;       min:        minimum data value to be mapped to the color table
;       max:        maximum data value to be mapped to the color table
;       linear:     use linear stretch
;       log:        use log stretch 
;       histeq:     use histogram equalization
;       asinh:      use asinh stretch
;       block:      block IDL command line until ATV terminates
;       align:      align image with previously displayed image
;       stretch:    keep same min and max as previous image
;       header:     FITS image header (string array) for use with data array
;       
; OUTPUTS:
;       None.  
; 
; COMMON BLOCKS:
;       atv_state:  contains variables describing the display state
;       atv_images: contains the internal copies of the display image
;       atv_color:  contains colormap vectors
;       atv_pdata:  contains plot and text annotation information
;
; RESTRICTIONS:
;       Requires IDL version 6.0 or greater.
;       Requires Craig Markwardt's cmps_form.pro routine.
;       Requires the GSFC IDL astronomy user's library routines.
;       Some features may not work under all operating systems.
;
; EXAMPLE:
;       To start atv running, just enter the command 'atv' at the idl
;       prompt, either with or without an array name or fits file name 
;       as an input.  Only one atv window will be created at a time,
;       so if one already exists and another image is passed to atv
;       from the idl command line, the new image will be displayed in 
;       the pre-existing atv window.
;
; MODIFICATION HISTORY:
;       Written by Aaron J. Barth, with contributions by 
;       Douglas Finkbeiner, Michael Liu, David Schlegel,
;       Wesley Colley, Jim Brauher, and Marshall Perrin.  
;       First released 17 December 1998.
;
;		 ***** DISCLAIMER ***** DISCLAIMER ***** DISCLAIMER *****
;
;			This is a highly modified version of atv! It's probably
;			buggy, and you use it at your own risk. Questions, comments
;			and complaints should go to Marshall Perrin,
;			mperrin@astro.berkeley.edu, and NOT to Aaron Barth.
;
;		 ***** DISCLAIMER ***** DISCLAIMER ***** DISCLAIMER *****
;
;
;       This version is 2.0pre1-MP, which contain's Marshall Perrin and Henry
;       Roe's modifications to ATV 1.3, merged into Aaron Barth's 
;       release 1.5 of 11 Aug 2004, and then with Jim Brauher's version 2.0.
;       2006-12-17: And then *again* partially merged with Aaron Barth's 2.0pre1
;       	(but not entirely...)
;       2008-10-06: Partial merge with Aaron Barth's 2.0pre4
;      
;      	The most notable additions by Marshall Perrin:
;      		- Preserves users colormap and !p.multi for external windows
;      		- Polarimetry vector field overplotting
;      		- 'Measure' mouse mode to measure distances
;      		- Better support for NaNs to indicate missing pixels
;      		- FITS reading code allows loading multiple files from disk to a
;      		  datacube (files must all be the same size)
;      		- Option for having a different title extra ("name") for each image
;      		  in a cube. Set via names= keyword argument.
;      		- Image blinking also blinks the titles.
;      		- Added SQRT, ASINH stretches. 
;      			See Lupton et al. AJ 118, 1406-1410	for more info on ASINH
;      			stretch
;      			WARNING - ASINH implementation is different from the ATV 2.0b4
;      			version
;   	    - Added code to preserve user's device decomposition setting.
;   	    - Additional / modified keyboard shortcuts
;   	    - keyboard shortcuts work while mouse button down (i.e. during
;   	       vector mode) 
;   	    - Can save image to IDL main-level variable (code due to D. Fanning)
;
;       For the most current version, revision history, instructions,
;       list of known bugs, and further information, go to:
;              http://www.physics.uci.edu/~barth/atv
; 
;-
;----------------------------------------------------------------------
;        atv startup, initialization, and shutdown routines
;----------------------------------------------------------------------

forward_function atv_wcsstring

pro atv_initcommon

; Routine to initialize the atv common blocks.  Use common blocks so
; that other IDL programs can access the atv internal data easily.

common atv_state, state
common atv_color, r_vector, g_vector, b_vector, user_r, user_g, user_b
common atv_pdata, nplot, maxplot, plot_ptr
common atv_images, $
   main_image, $
   main_image_stack, $
   display_image, $
   scaled_image, $
   blink_image1, $
   blink_image2, $
   blink_image3, $
   unblink_image, $  
   pan_image, $
   header_array, $  ; for multiple-file image cube mode
   image_names ; for image cube mode
 

state = {                   $
        version: '2.0b4-MP', $              ; version # of this release
        head_ptr: ptr_new(), $  ; pointer to image header
        astr_ptr: ptr_new(), $  ; pointer to astrometry info structure
        firstimage: 1, $        ; is this the first image?
        block: 0, $             ; are we in blocking mode?
        wcstype: 'none', $      ; coord info type (none/angle/lambda)
        equinox: 'J2000', $     ; equinox of coord system
        display_coord_sys: 'RA--', $ ; coord system displayed
        display_equinox: 'J2000', $ ; equinox of displayed coords
        display_base60: 1B, $   ; Display RA,dec in base 60?
        cunit: '', $            ; wavelength units
        imagename: '', $        ; image file name
        bitdepth: 8, $          ; 8 or 24 bit color mode?
        user_decomposed: 1, $			 ; User's setting of device,/decomposed (outside ATV)
        screen_xsize: 1000, $   ; horizontal size of screen
        screen_ysize: 1000, $   ; vertical size of screen
        base_id: 0L, $          ; id of top-level base
        base_min_size: [512L, 300L], $ ; min size for top-level base
        draw_base_id: 0L, $     ; id of base holding draw window
        draw_window_id: 0L, $   ; window id of draw window
        draw_widget_id: 0L, $   ; widget id of draw widget
        mousemode: "color", $   ; color, blink, zoom, or imexam
        mode_droplist_id: 0L, $ ; id of mode droplist widget
        track_window_id: 0L, $  ; widget id of tracking window
        pan_widget_id: 0L, $    ; widget id of pan window
        pan_window_id: 0L, $    ; window id of pan window
        active_window_id: 0L, $ ; user's active window outside atv
        active_window_pmulti: lonarr(5), $ ; user's active window p.multi
        info_base_id: 0L, $     ; id of base holding info bars
        location_bar_id: 0L, $  ; id of (x,y,value) label
        wcs_bar_id: 0L, $       ; id of WCS label widget
        min_text_id: 0L,  $     ; id of min= widget
        max_text_id: 0L, $      ; id of max= widget
        nonlin_text_id: 0L, $      ; id of text for asinh's nonlinearity
	    nonlin_base_id: 0L, $				; id of base for asinh's nonlinearity
        menu_ids: lonarr(35), $ ; list of top menu items
        colorbar_base_id: 0L, $ ; id of colorbar base widget
        colorbar_widget_id: 0L, $ ; widget id of colorbar draw widget
        colorbar_window_id: 0L, $ ; window id of colorbar
        colorbar_height: 6L, $  ; height of colorbar in pixels
        ncolors: 0B, $          ; image colors (!d.table_size - 9)
        box_color: 2, $         ; color for pan box and zoom x
        brightness: 0.5, $      ; initial brightness setting
        contrast: 0.5, $        ; initial contrast setting
        image_min: 0.0, $       ; min(main_image)
        image_max: 0.0, $       ; max(main_image)
        min_value: 0.0, $       ; min data value mapped to colors
        max_value: 0.0, $       ; max data value mapped to colors
        skymode: 0.0, $                ; sky mode value
        skysig: 0.0, $                 ; sky sigma value
        draw_window_size: [512L, 512L], $ ; size of main draw window
        track_window_size: 121L, $     ; size of tracking window
        pan_window_size: 121L, $; size of pan window
        pan_scale: 0.0, $       ; magnification of pan image
        image_size: [0L,0L,0L], $      ; [0:1] gives size of main_image
                                ; [0:2] gives size of main_image_stack
        cur_image_num: 0L, $    ; gives current image number in
                                ; main_image_stack
        curimnum_base_id: 0L, $ ; id of cur_image_num base widget
        curimnum_text_id: 0L, $ ; id of cur_image_num textbox widget
        curimnum_slidebar_id: 0L, $    ; id of cur_image_num slider widget
        scale_mode_droplist_id: 0L, $  ; id of scale droplist widget
        curimnum_minmaxmode: 'Constant', $  ; mode for determining min/max
                                ; of display when changing curimnum
        invert_colormap: 0L, $  ; 0=normal, 1=inverted
        coord: [0L, 0L],  $     ; cursor position in image coords
        scaling: 3L, $          ; 0=linear, 1=log, 2=histeq, 3=asinh, 4=sqrt
        nonlinearity_value: 5000.0, $     ; "nonlinearity" parameter for asinh stretch
        asinh_beta: 0.1, $      ; asinh nonlinearity parameter
        offset: [0L, 0L], $     ; offset to viewport coords
        base_pad: [0L, 0L], $   ; padding around draw base
        zoom_level: 0L, $       ; integer zoom level, 0=normal
        zoom_factor: 1.0, $     ; magnification factor = 2^zoom_level
        rot_angle: 0.0, $       ; current image rotation angle
        invert_image: 'none', $ ; 'none', 'x', 'y', 'xy' image invert
        centerpix: [0L, 0L], $  ; pixel at center of viewport
        cstretch: 0B, $         ; flag = 1 while stretching colors
        pan_offset: [0L, 0L], $ ; image offset in pan window
        frame: 1L, $            ; put frame around ps output?
        framethick: 6, $        ; thickness of frame
        plot_coord: [0L, 0L], $ ; cursor position for a plot
        vector_coord1: [0L, 0L], $ ; 1st cursor position in vector plot  
        vector_coord2: [0L, 0L], $ ; 2nd cursor position in vector plot
        vector_pixmap_id: 0L, $ ; id for vector pixmap 
        vectorpress: 0L, $      ; are we plotting a vector?
        vectorstart: [0L,0L], $ ; device x,y of vector start
        plot_type:'', $         ; plot type for plot window
        lineplot_widget_id: 0L, $ ; id of lineplot widget
        lineplot_window_id: 0L, $ ; id of lineplot window
        lineplot_base_id: 0L, $ ; id of lineplot top-level base
        lineplot_size: [600L, 500L], $ ; size of lineplot window
        lineplot_min_size: [100L, 0L], $ ; min size of lineplot window
        lineplot_pad: [0L, 0L], $ ; padding around lineplot window
        lineplot_xmin_id: 0L, $ ; id of xmin for lineplot windows
        lineplot_xmax_id: 0L, $ ; id of xmax for lineplot windows
        lineplot_ymin_id: 0L, $ ; id of ymin for lineplot windows
        lineplot_ymax_id: 0L, $ ; id of ymax for lineplot windows
        lineplot_xmin: 0.0, $   ; xmin for lineplot windows
        lineplot_xmax: 0.0, $   ; xmax for lineplot windows
        lineplot_ymin: 0.0, $   ; ymin for lineplot windows
        lineplot_ymax: 0.0, $   ; ymax for lineplot windows
        lineplot_xmin_orig: 0.0, $ ; original xmin saved from histplot
        lineplot_xmax_orig: 0.0, $ ; original xmax saved from histplot
        holdrange_base_id: 0L, $ ; base id for 'Hold Range' button
        holdrange_button_id: 0L, $ ; button id for 'Hold Range' button
        holdrange_value: 0, $   ; 0=HoldRange Off, 1=HoldRange On
        histbutton_base_id: 0L, $ ; id of histogram button base
        histplot_binsize_id: 0L, $ ; id of binsize for histogram plot
        x1_pix_id: 0L, $        ; id of x1 pixel for histogram plot
        x2_pix_id: 0L, $        ; id of x2 pixel for histogram plot
        y1_pix_id: 0L, $        ; id of y1 pixel for histogram plot
        y2_pix_id: 0L, $        ; id of y2 pixel for histogram plot
        binsize: 0.0, $         ; binsize for histogram plots
        regionform_id: 0L, $    ; id of region form widget
        reg_ids_ptr: ptr_new(), $ ; ids for region form widget
;        writeimage_ids_ptr: ptr_new(),$ ; ids for writeimage form widget
;        writeformat: 'PNG', $   ; default format for WriteImage
        cursorpos: lonarr(2), $ ; cursor x,y for photometry & stats
        centerpos: fltarr(2), $ ; centered x,y for photometry
        cursorpos_id: 0L, $     ; id of cursorpos widget
        centerpos_id: 0L, $     ; id of centerpos widget
        centerbox_id: 0L, $     ; id of centeringboxsize widget
        radius_id: 0L, $        ; id of radius widget
        innersky_id: 0L, $      ; id of inner sky widget
        outersky_id: 0L, $      ; id of outer sky widget
        magunits: 0, $          ; 0=counts, 1=magnitudes
        skytype: 0, $           ; 0=idlphot,1=median,2=no sky subtract
        photzpt: 25.0, $        ; magnitude zeropoint
		  photplotmode:0, $					 ; photometry plot type: 0=linear 1=log
        skyresult_id: 0L, $     ; id of sky widget
        photresult_id: 0L, $    ; id of photometry result widget
        photerror_id: 0L, $     ; id of photometry error widget
        fwhm_id: 0L, $          ; id of fwhm widget
        radplot_widget_id: 0L, $ ; id of radial profile widget
        radplot_window_id: 0L, $ ; id of radial profile window
        photzoom_window_id: 0L, $ ; id of photometry zoom window
        photzoom_size: 190L, $  ; size in pixels of photzoom window
        showradplot_id: 0L, $   ; id of button to show/hide radplot
        photwarning_id: 0L, $   ; id of photometry warning widget
        photwarning: ' ', $     ; photometry warning text
        centerboxsize: 9L, $    ; centering box size
        aprad: 20.0, $           ; aperture photometry radius
        innersky: 40L, $        ; inner sky radius
        outersky: 60L, $        ; outer sky radius
        headinfo_base_id: 0L, $ ; headinfo base widget id
        pixtable_base_id: 0L, $ ; pixel table base widget id
        pixtable_tbl_id: 0L, $  ; pixel table widget_table id
        stats_base_id: 0L, $    ; base widget for image stats
        statboxsize: 11L, $     ; box size for computing statistics
        statbox_id: 0L, $       ; widget id for stat box size 
        stat_npix_id: 0L, $     ; widget id for # pixels in stats box
        statxcenter_id: 0L, $   ; widget id for stat box x center
        statycenter_id: 0L, $   ; widget id for stat box y center
        statbox_min_id: 0L, $   ; widget id for stat min box
        statbox_max_id: 0L, $   ; widget id for stat max box
        statbox_mean_id: 0L, $  ; widget id for stat mean box
        statbox_median_id: 0L, $; widget id for stat median box
        statbox_stdev_id: 0L, $ ; widget id for stat stdev box
          statzoom_size: 300, $   ; size of statzoom window
          statzoom_widget_id: 0L, $      ; widget id for stat zoom window
          statzoom_window_id: 0L, $      ; window id for stat zoom window
          showstatzoom_id: 0L, $  ; widget id for show/hide button
          pan_pixmap: 0L, $       ; window id of pan pixmap
          current_dir: '', $      ; current readfits directory
          graphicsdevice: '', $   ; screen device
          ispsformon: 0, $        ; is cmps_form running?
          newrefresh: 0, $        ; refresh since last blink?
          window_title: 'atv:', $ ; string for top level title
          title_blink1: '', $     ; window title for 1st blink image
          title_blink2: '', $     ; window title for 2nd blink image
          title_blink3: '', $     ; window title for 3rd blink image
          title_extras: '', $     ; extras for image title
          blinks: 0B, $           ; remembers which images are blinked
          activator: 0, $         ; is "activator" mode on?
          delimiter: '/', $       ; filesystem level delimiter 
          default_align: 1, $     ; align next image by default?
          default_autoscale: 1, $ ; autoscale images by default?
          default_stretch: 0 ,$    ; use previous minmax for new image?
		  autozoom: 1 ,$				 ; zoom images to fit window on open?
		  polarim_lowthresh: 0.0,$		 ; Polarimetry low threshhold
		  polarim_highthresh: 1.0,$		 ; Polarimetry high threshhold
		  polarim_display: 1, $			 ; overplot polarimetry vectors?
		  polarim_plotindex: 0,$		 ; Which plot structure is the polarimetry?
		  polarim_present: 0, $			 ; Do we have polarimetry data?		
		  polarim_lowth_id: 0L, $		 ; widget ID of polarimetry low thresh
		  polarim_highth_id: 0L, $		 ; widget ID of polarimetry high thresh
		  polarim_mag_id: 0L, $		 ; widget ID of polarimetry high thresh
		  polarim_offset_id: 0L, $		 ; widget ID of polarimetry high thresh
		  polarim_display_id: 0L $		 ; widget ID of polarimetry display flag
        }

nplot = 0
maxplot = 5000
plot_ptr = ptrarr(maxplot+1)  ; The 0th element isn't used.

blink_image1 = 0
blink_image2 = 0
blink_image3 = 0

end

;---------------------------------------------------------------------

pro atv_startup

; This routine initializes the atv internal variables, and creates and
; realizes the window widgets.  It is only called by the atv main
; program level, when there is no previously existing atv window.

common atv_state
common atv_color

; common block atvcompileoptions is only used by my script to compile
; a version of atv for the idl virtual machine.

; save the user color table and pmulti first
tvlct, user_r, user_g, user_b, /get

; Read in a color table to initialize !d.table_size
; As a bare minimum, we need the 8 basic colors used by ATV_ICOLOR(),
; plus 2 more for a color map.

;loadct, 0, /silent
if (!d.table_size LT 12) then begin
    message, 'Too few colors available for color table'
    tvlct, user_r, user_g, user_b
    atv_shutdown
endif

; Initialize the common blocks
atv_initcommon

state.active_window_pmulti = !p.multi
!p.multi = 0

osfamily = strupcase(!version.os_family)
case osfamily of
    'UNIX': state.delimiter = '/'
    'WINDOWS': state.delimiter = '\'
    else:
endcase

state.ncolors = !d.table_size - 9

; If compiling atv to make a sav file for the atv virtual machine,
; always do it for 24-bit color with retain & decomposed set.
; Uncomment this block to compile atv for idl vm.  For some reason,
; idl vm gets !d.n_colors=256 even on a 24-bit display, so we need
; this to work around it to force 24-bit mode.
;device, true_color=24
;device, decomposed=0
;device, retain=2
;state.bitdepth=24

; For normal idl operation, use the following.  Comment this block out
; if compiling atv for idl vm.
if (!d.n_colors LE 256) then begin
    state.bitdepth = 8 
endif else begin
    state.bitdepth = 24
;    device, decomposed=0  ; Now handled in setwindow/resetwindow
endelse


state.graphicsdevice = !d.name

state.screen_xsize = (get_screen_size())[0]
state.screen_ysize = (get_screen_size())[1]

; Get the current window id and color table
atv_getwindow


; Define the widgets.  For the widgets that need to be modified later
; on, save their widget ids in state variables

base = widget_base(title = 'atv', $
                   /column, /base_align_right, $
                   app_mbar = top_menu, $
                   uvalue = 'atv_base', $
                   /tlb_size_events)
state.base_id = base

tmp_struct = {cw_pdmenu_s, flags:0, name:''}

top_menu_desc = [ $
                {cw_pdmenu_s, 1, 'File'}, $ ; file menu
                {cw_pdmenu_s, 0, 'ReadFits'}, $
                {cw_pdmenu_s, 0, 'WriteFits'}, $
                {cw_pdmenu_s, 0, 'WritePS'},  $
                {cw_pdmenu_s, 0, 'WriteImage'}, $
                {cw_pdmenu_s, 1, 'Save To IDL variable...'}, $
                {cw_pdmenu_s, 0, 'Save Image to IDL variable'}, $
 				{cw_pdmenu_s, 0, 'Save Image Cube to IDL variable'}, $
             	{cw_pdmenu_s, 2, 'Save FITS Header to IDL variable'}, $
                {cw_pdmenu_s, 0, '--------------'}, $
                {cw_pdmenu_s, 0, 'GetImage:'}, $
                {cw_pdmenu_s, 0, ' DSS'}, $
                {cw_pdmenu_s, 0, ' FIRST'}, $
                {cw_pdmenu_s, 0, '--------------'}, $
                {cw_pdmenu_s, 2, 'Quit'}, $
                {cw_pdmenu_s, 1, 'ColorMap'}, $ ; color menu
                {cw_pdmenu_s, 0, 'Grayscale'}, $
                {cw_pdmenu_s, 0, 'Blue-White'}, $
                {cw_pdmenu_s, 0, 'Red-Orange'}, $
                {cw_pdmenu_s, 0, 'Green-White'}, $
                {cw_pdmenu_s, 0, 'Red-Purple'}, $
                {cw_pdmenu_s, 0, 'Blue-Red'}, $
                {cw_pdmenu_s, 0, 'Rainbow'}, $
                {cw_pdmenu_s, 0, 'Rainbow18'}, $
                {cw_pdmenu_s, 0, 'BGRY'}, $
                {cw_pdmenu_s, 0, 'GRBW'}, $
                {cw_pdmenu_s, 0, 'Standard Gamma-II'}, $
                {cw_pdmenu_s, 0, 'Prism'}, $
                {cw_pdmenu_s, 0, '16 Level'}, $
                {cw_pdmenu_s, 0, 'Stern Special'}, $
                {cw_pdmenu_s, 0, 'Haze'}, $
                {cw_pdmenu_s, 0, 'Blue-Pastel-Red'}, $
                {cw_pdmenu_s, 0, 'Mac'}, $
                {cw_pdmenu_s, 0, 'Blue-Red 2'}, $
                {cw_pdmenu_s, 2, 'ATV Special'}, $
                {cw_pdmenu_s, 1, 'Scaling'}, $ ; scaling menu
                {cw_pdmenu_s, 0, 'Asinh'}, $
                {cw_pdmenu_s, 0, 'Log'}, $
                {cw_pdmenu_s, 0, 'Linear'}, $
                {cw_pdmenu_s, 0, 'HistEq'}, $
                {cw_pdmenu_s, 0, 'Square Root'}, $
                {cw_pdmenu_s, 6, 'Asinh Settings'}, $
                {cw_pdmenu_s, 1, 'Labels'}, $ ; labels menu
                {cw_pdmenu_s, 0, 'TextLabel'}, $
                {cw_pdmenu_s, 0, 'Arrow'}, $
                {cw_pdmenu_s, 0, 'Contour'}, $
                {cw_pdmenu_s, 0, 'Compass'}, $
                {cw_pdmenu_s, 0, 'ScaleBar'}, $
                {cw_pdmenu_s, 0, 'Polarimetry'}, $
                {cw_pdmenu_s, 0, 'Region'}, $
                {cw_pdmenu_s, 0, 'WCS Grid'}, $
                {cw_pdmenu_s, 0, 'EraseLast'}, $
                {cw_pdmenu_s, 0, 'EraseAll'}, $
                {cw_pdmenu_s, 4, 'LoadRegions'}, $
                {cw_pdmenu_s, 2, 'SaveRegions'}, $
                {cw_pdmenu_s, 1, 'Blink'}, $
                {cw_pdmenu_s, 0, 'SetBlink1'}, $
                {cw_pdmenu_s, 0, 'SetBlink2'}, $
                {cw_pdmenu_s, 2, 'SetBlink3'}, $
                {cw_pdmenu_s, 1, 'Rotate/Zoom'}, $
                {cw_pdmenu_s, 0, 'Rotate'}, $
                {cw_pdmenu_s, 0, '90 deg'}, $
                {cw_pdmenu_s, 0, '180 deg'}, $
                {cw_pdmenu_s, 0, '270 deg'}, $
                {cw_pdmenu_s, 0, '--------------'}, $
                {cw_pdmenu_s, 0, 'Invert X'}, $
                {cw_pdmenu_s, 0, 'Invert Y'}, $
                {cw_pdmenu_s, 0, 'Invert XY'}, $
                {cw_pdmenu_s, 0, '--------------'}, $
                {cw_pdmenu_s, 0, '1/16x'}, $
                {cw_pdmenu_s, 0, '1/8x'}, $
                {cw_pdmenu_s, 0, '1/4x'}, $
                {cw_pdmenu_s, 0, '1/2x'}, $
                {cw_pdmenu_s, 0, '1x'}, $
                {cw_pdmenu_s, 0, '2x'}, $
                {cw_pdmenu_s, 0, '4x'}, $
                {cw_pdmenu_s, 0, '8x'}, $
                {cw_pdmenu_s, 2, '16x'}, $
                {cw_pdmenu_s, 1, 'ImageInfo'}, $ ;info menu
                {cw_pdmenu_s, 0, 'ImageHeader'}, $
                {cw_pdmenu_s, 0, 'Photometry'}, $
                {cw_pdmenu_s, 0, 'Statistics'}, $
                {cw_pdmenu_s, 0, 'ImageHeader'}, $
                {cw_pdmenu_s, 0, 'Print Filename'}, $
                {cw_pdmenu_s, 0, 'Pixel Table'}, $
                {cw_pdmenu_s, 0, 'Load Regions'}, $
                {cw_pdmenu_s, 0, 'Save Regions'}, $
                {cw_pdmenu_s, 0, 'Archive Image'}, $
                {cw_pdmenu_s, 0, '--------------'}, $
                {cw_pdmenu_s, 0, 'RA,dec (J2000)'}, $
                {cw_pdmenu_s, 0, 'RA,dec (B1950)'}, $
                {cw_pdmenu_s, 0, '--------------'}, $
                {cw_pdmenu_s, 0, 'RA,dec (J2000) deg'}, $
                {cw_pdmenu_s, 0, 'Galactic'}, $
                {cw_pdmenu_s, 0, 'Ecliptic (J2000)'}, $
                {cw_pdmenu_s, 2, 'Native'}, $
                {cw_pdmenu_s, 1, 'Help'}, $ ; help menu
                {cw_pdmenu_s, 2, 'ATV Help'} $
                ]

top_menu = cw_pdmenu(top_menu, top_menu_desc, $
                     ids = state.menu_ids, $
                     /mbar, $
                     /help, $
                     /return_name, $
                     uvalue = 'top_menu')

track_base = widget_base(base, /row)

state.info_base_id = widget_base(track_base, /column, /base_align_right,xpad=0,ypad=0)
buttonbar_base = widget_base(base, column=2,xpad=0,ypad=0)

state.curimnum_base_id = widget_base(base, $
                                     /base_align_right, column=3, $
                                     frame = 1, xsize=1, ysize=1, map=0)

state.draw_base_id = widget_base(base, $
                                 /column, /base_align_left, $
                                 uvalue = 'draw_base', $
                                 frame = 2)

state.colorbar_base_id = widget_base(base, $
                                     uvalue = 'cqolorbar_base', $
                                     /column, /base_align_left, $
                                     frame = 2)

state.min_text_id = cw_field(state.info_base_id, $
                             uvalue = 'min_text', $
                             /floating,  $
                             title = 'Min=', $
                             value = state.min_value,  $
                             /return_events, $
                             xsize = 12)

state.max_text_id = cw_field(state.info_base_id, $
                             uvalue = 'max_text', $
                             /floating,  $
                             title = 'Max=', $
                             value = state.max_value, $
                             /return_events, $
                             xsize = 12)

state.nonlin_base_id =    widget_base(state.info_base_id, /column,xpad=0,ypad=0,map=0,/base_align_right,xsize=1,ysize=1)

state.nonlin_text_id = cw_field(state.nonlin_base_id, $
                             uvalue = 'nonlin_text', $
                             /floating,  $
                             title = 'Nonlinearity=', $
                             value = state.nonlinearity_value, $
                             /return_events, $
                             xsize = 12)


tmp_string = string(1000, 1000, 1.0e-10, $
                    format = '("(",i5,",",i5,") ",g12.5)' )

state.location_bar_id = widget_label (state.info_base_id, $
                                      value = tmp_string,  $
                                      uvalue = 'location_bar',  frame = 1)

tmp_string = string(12, 12, 12.001, -60, 60, 60.01, ' J2000', $
        format = '(i2,":",i2,":",f6.3,"  ",i3,":",i2,":",f5.2," ",a6)' )
    
state.wcs_bar_id = widget_label (state.info_base_id, $
                                 value = tmp_string,  $
                                 uvalue = 'wcs_bar',  frame = 1)

state.pan_widget_id = widget_draw(track_base, $
                                  xsize = state.pan_window_size, $
                                  ysize = state.pan_window_size, $
                                  frame = 2, uvalue = 'pan_window', $
                                  /button_events, /motion_events)

track_window = widget_draw(track_base, $
                           xsize=state.track_window_size, $
                           ysize=state.track_window_size, $
                           frame=2, uvalue='track_window')

modebase = widget_base(buttonbar_base, /row, /base_align_center)
modelist = ['Color', 'Zoom', 'Blink', 'ImExam','Measure','Region', 'Vector','Print']
mode_droplist_id = widget_droplist(modebase, $
                                   frame = 1, $
                                   title = 'MouseMode:', $
                                   uvalue = 'mode', $
                                   value = modelist)
state.mode_droplist_id = mode_droplist_id


button_base = widget_base(buttonbar_base, row=2)

invert_button = widget_button(button_base, $
                              value = 'Invert', $
                              uvalue = 'invert')

restretch_button = widget_button(button_base, $
                             value = 'Restretch', $
                             uvalue = 'restretch_button')

autoscale_button = widget_button(button_base, $
                                 uvalue = 'autoscale_button', $
                                 value = 'AutoScale')

fullrange_button = widget_button(button_base, $
                                 uvalue = 'full_range', $
                                 value = 'FullRange')

dummy_spacing_widget = widget_label(button_base,value='')

zoomin_button = widget_button(button_base, $
                              value = 'ZoomIn', $
                              uvalue = 'zoom_in')

zoomout_button = widget_button(button_base, $
                               value = 'ZoomOut', $
                               uvalue = 'zoom_out')

zoomone_button = widget_button(button_base, $
                               value = 'Zoom1', $
                               uvalue = 'zoom_one')

fullview_button = widget_button(button_base, $
                                value = 'FullView', $
                                uvalue = 'fullview')

center_button = widget_button(button_base, $
                              value = 'Center', $
                              uvalue = 'center')

;done_button = widget_button(button_base, $
;                            value = 'Done', $
;                            uvalue = 'done')

state.curimnum_text_id = cw_field(state.curimnum_base_id, $
                             uvalue = 'curimnum_text', $
                             /long,  $
                             /row, $
                             title = 'Image #=', $
                             value = state.cur_image_num, $
                             /return_events, $
                             xsize = 5)

state.curimnum_slidebar_id = widget_slider(state.curimnum_base_id, $
                            /drag, $ 
                            max = 1, $
                            min = 0, $
                            scr_xsize = 290L, $
                            sensitive = 1, $
                            scroll = 1L, $
                            /suppress_value, $
                            uvalue = 'curimnum_slidebar', $
                            value = 0, $
                            vertical = 0)

modelist = ['Constant', 'AutoScale', 'Min/Max', 'Zero/Max']
state.scale_mode_droplist_id = widget_droplist(state.curimnum_base_id, $
                                   uvalue = 'curimnum_minmaxmode', $
                                   value = modelist)

; Set widget y size for small screens
state.draw_window_size[1] = state.draw_window_size[1] < $
  (state.screen_ysize - 300)

state.draw_widget_id = widget_draw(state.draw_base_id, $
                                   uvalue = 'draw_window', $
                                   /motion_events,  /button_events, $
                                   keyboard_events=2, $
                                   scr_xsize = state.draw_window_size[0], $
                                   scr_ysize = state.draw_window_size[1]) 

state.colorbar_widget_id = widget_draw(state.colorbar_base_id, $
                                       uvalue = 'colorbar', $
                                       scr_xsize = state.draw_window_size[0], $
                                       scr_ysize = state.colorbar_height)

; Create the widgets on screen

widget_control, base, /realize
widget_control, state.pan_widget_id, draw_motion_events = 0

; get the window ids for the draw widgets

widget_control, track_window, get_value = tmp_value
state.track_window_id = tmp_value
widget_control, state.draw_widget_id, get_value = tmp_value
state.draw_window_id = tmp_value
widget_control, state.pan_widget_id, get_value = tmp_value
state.pan_window_id = tmp_value
widget_control, state.colorbar_widget_id, get_value = tmp_value
state.colorbar_window_id = tmp_value

; set the event handlers

widget_control, top_menu, event_pro = 'atv_topmenu_event'
widget_control, state.draw_widget_id, event_pro = 'atv_draw_event'
widget_control, state.pan_widget_id, event_pro = 'atv_pan_event'

; Find window padding sizes needed for resizing routines.
; Add extra padding for menu bar, since this isn't included in 
; the geometry returned by widget_info.
; Also add extra padding for margin (frame) in draw base.

basegeom = widget_info(state.base_id, /geometry)
drawbasegeom = widget_info(state.draw_base_id, /geometry)


; Initialize the vectors that hold the current color table.
; See the routine atv_stretchct to see why we do it this way.

r_vector = bytarr(state.ncolors)
g_vector = bytarr(state.ncolors)
b_vector = bytarr(state.ncolors)

atv_getct, 0
state.invert_colormap = 0

; Create a pixmap window to hold the pan image
window, /free, xsize=state.pan_window_size, ysize=state.pan_window_size, $
  /pixmap
state.pan_pixmap = !d.window
atv_resetwindow

atv_colorbar

widget_control, state.base_id, tlb_get_size=tmp_event
state.base_pad = tmp_event - state.draw_window_size


end

;--------------------------------------------------------------------

pro atv_colorbar

; Routine to tv the colorbar at the bottom of the atv window

common atv_state

atv_setwindow, state.colorbar_window_id

xsize = (widget_info(state.colorbar_widget_id, /geometry)).xsize

b = congrid( findgen(state.ncolors), xsize) + 8
c = replicate(1, state.colorbar_height)
a = b # c

tv, a

atv_resetwindow

end

;-------------------------------------------------------------------

pro atvclear

; displays a small blank image, useful for clearing memory if atv is
; displaying a huge image.

atv, fltarr(10,10)

end

;--------------------------------------------------------------------

pro atv_shutdown, windowid

; routine to kill the atv window(s) and clear variables to conserve
; memory when quitting atv.  The windowid parameter is used when
; atv_shutdown is called automatically by the xmanager, if atv is
; killed by the window manager.

common atv_images
common atv_state
common atv_color
common atv_pdata


; reset color table and pmulti to user values
tvlct, user_r, user_g, user_b
!p.multi = state.active_window_pmulti

; Kill top-level base if it still exists
if (xregistered ('atv')) then widget_control, state.base_id, /destroy

; Destroy all pointers to plots and their heap variables: this runs
; ptr_free on any existing plot pointers
if (nplot GT 0) then begin
    atverase, /norefresh
endif

if (size(state, /tname) EQ 'STRUCT') then begin
    if (!d.name EQ state.graphicsdevice) then wdelete, state.pan_pixmap
    if (ptr_valid(state.head_ptr)) then ptr_free, state.head_ptr
    if (ptr_valid(state.astr_ptr)) then ptr_free, state.astr_ptr
endif

; Clean up saved variables in common blocks to conserve memory.
; Previously this was done using delvarx, but since delvarx uses an
; execute function, it's incompatible with IDL virtual machine.  So,
; just set these variables to zero.

plot_ptr=0
maxplot=0
main_image=0
display_image=0
scaled_image=0
blink_image1=0
blink_image2=0
blink_image3=0
unlink_image=0
pan_image=0
r_vector=0
g_vector=0
b_vector=0
user_r=0
user_g=0
user_b=0
state=0


return
end

;--------------------------------------------------------------------
;                  main atv event loops
;--------------------------------------------------------------------

pro atv_topmenu_event, event

; Event handler for top menu

common atv_state
common atv_images

widget_control, event.id, get_uvalue = event_name

if (!d.name NE state.graphicsdevice and event_name NE 'Quit') then return
if (state.bitdepth EQ 24) then true = 1 else true = 0

; Need to get active window here in case mouse goes to menu from top
; of atv window without entering the main base
atv_getwindow


case event_name of
    
; File menu options:
    'ReadFits': begin
        atv_readmultifits, newimage=newimage
        if (newimage EQ 1) then begin
            atv_getstats, align=state.default_align
            atv_settitle
            if (state.default_align EQ 0) then begin
                state.zoom_level =  0
                state.zoom_factor = 1.0
            endif
            if (state.default_stretch EQ 0 AND $
                state.default_autoscale EQ 1) then atv_autoscale
            if (state.firstimage EQ 1) then atv_autoscale
            atv_set_minmax
            atv_displayall
			if state.autozoom then atv_autozoom  ; MDP 2007-10-23. suddenly this was crashing on opening new NIRC2 images until I moved this line down here after the displayall. ???
            state.firstimage = 0
        endif
    end
    'WriteFits': atv_writefits
    'WritePS' : atv_writeps
;    'WriteImage': atv_writeimage
    'PNG': atv_writeimage, 'png'
    'JPEG': atv_writeimage, 'jpg'
    'TIFF': atv_writeimage, 'tiff'
    'Save Image to IDL variable': atv_SaveToVariable, "Image"
    'Save Image Cube to IDL variable': atv_SaveToVariable, "Cube"
    'Save FITS Header to IDL variable': atv_SaveToVariable,  "Header"
    'GetImage:':
    ' DSS': atv_getdss
    ' FIRST': atv_getfirst
    'LoadRegions': atv_loadregion
    'SaveRegions': atv_saveregion
    'Quit':     if (state.activator EQ 0) then atv_shutdown $
      else state.activator = 0
; ColorMap menu options:            
    'Grayscale': atv_getct, 0
    'Blue-White': atv_getct, 1
    'GRBW': atv_getct, 2
    'Red-Orange': atv_getct, 3
    'BGRY': atv_getct, 4
    'Standard Gamma-II': atv_getct, 5
    'Prism': atv_getct, 6
    'Red-Purple': atv_getct, 7
    'Green-White': atv_getct, 8
    'Blue-Red': atv_getct, 11
    '16 Level': atv_getct, 12
    'Rainbow': atv_getct, 13
    'Stern Special': atv_getct, 15
    'Haze' : atv_getct, 16
    'Blue-Pastel-Red': atv_getct, 17
    'Mac': atv_getct, 25
    'Plasma': atv_getct, 32
    'Blue-Red 2': atv_getct, 33
    'Rainbow18': atv_getct, 38
    'ATV Special': atv_makect, event_name
    'Velocity1': atv_makect, event_name
    'Velocity2': atv_makect, event_name
; Scaling options:
    'Linear': atv_setscaling, 0
    'Log': atv_setscaling, 1
    'HistEq': atv_setscaling, 2
    'Asinh': atv_setscaling, 3
    'Square Root': atv_setscaling, 4
    'Asinh Settings': begin
        atv_setasinh
    end

; Label options:
    'TextLabel': atv_textlabel
    'Arrow': atv_setarrow
    'Contour': atv_oplotcontour
    'Compass': atv_setcompass
    'ScaleBar': atv_setscalebar
	'Polarimetry': atv_polarim
    'Region': atv_setregion
    'WCS Grid': atv_wcsgridlabel
    'EraseLast': atverase, 1
    'EraseAll': atverase

; Blink options:
    'SetBlink1': begin   
        atv_setwindow, state.draw_window_id
        blink_image1 = tvrd(true = true) 
        state.title_blink1 = state.window_title
    end
    'SetBlink2': begin   
        atv_setwindow, state.draw_window_id
        blink_image2 = tvrd(true = true)
        state.title_blink2 = state.window_title
    end
    'SetBlink3': begin   
        atv_setwindow, state.draw_window_id
        blink_image3 = tvrd(true = true)
        state.title_blink3 = state.window_title
    end

    'MakeRGB' : atv_makergb

; Zoom/Rotate options
    '1/16x': atv_zoom, 'onesixteenth'
    '1/8x': atv_zoom, 'oneeighth'
    '1/4x': atv_zoom, 'onefourth'
    '1/2x': atv_zoom, 'onehalf'
    '1x': atv_zoom, 'one'
    '2x': atv_zoom, 'two'
    '4x': atv_zoom, 'four'
    '8x': atv_zoom, 'eight'
    '16x': atv_zoom, 'sixteen'
    'Zoom In': atv_zoom, 'in'
    'Zoom Out': atv_zoom, 'out'
    'Center': begin
        state.centerpix = round(state.image_size[0:1] / 2.)
        atv_refresh
    end
    'None': atv_invert, 'none'
    'Invert X': atv_invert, 'x'
    'Invert Y': atv_invert, 'y'
    'Invert XY': atv_invert, 'xy'
    'Rotate': atv_rotate, '0', /get_angle
    '0 deg': atv_rotate, '0'
    '90 deg': atv_rotate, '90'
    '180 deg': atv_rotate, '180'
    '270 deg': atv_rotate, '270'


; Info options:
    'Photometry': atv_apphot
    'ImageHeader': atv_headinfo
    'Print Filename': print, state.imagename
    'Statistics': atv_showstats
    'Pixel Table': atv_pixtable
    'Load Regions': atv_regionfilelabel
    'Save Regions': atv_saveregion
    'Archive Image': atv_getimage

; Coordinate system options:
    '--------------':
    'RA,dec (J2000)': BEGIN 
       state.display_coord_sys = 'RA--'
       state.display_equinox = 'J2000'
       state.display_base60 = 1B
       heuler, *(state.head_ptr), /celestial
       atv_gettrack             ; refresh coordinate window
       atv_displayall           ; replot all overplots, including WCS grid
    END 
    'RA,dec (B1950)': BEGIN 
       state.display_coord_sys = 'RA--'
       state.display_equinox = 'B1950'
       state.display_base60 = 1B
	   	; FIXME will crash if no head_ptr present!!
       heuler, *(state.head_ptr), /celestial
       atv_gettrack             ; refresh coordinate window
       atv_displayall           ; replot all overplots, including WCS grid
    END
    'RA,dec (J2000) deg': BEGIN 
       state.display_coord_sys = 'RA--'
       state.display_equinox = 'J2000'
       state.display_base60 = 0B
       heuler, *(state.head_ptr), /celestial
       atv_gettrack             ; refresh coordinate window
       atv_displayall           ; replot all overplots, including WCS grid
    END 
    'Galactic': BEGIN 
       state.display_coord_sys = 'GLON'
       heuler, *(state.head_ptr), /galactic
       atv_gettrack             ; refresh coordinate window
       atv_displayall           ; replot all overplots, including WCS grid
    END 
    'Ecliptic (J2000)': BEGIN 
       state.display_coord_sys = 'ELON'
       state.display_equinox = 'J2000'
       heuler, *(state.head_ptr), /ecliptic
       atv_gettrack             ; refresh coordinate window
       atv_displayall           ; replot all overplots, including WCS grid
    END 
    'Native': BEGIN 
       IF (state.wcstype EQ 'angle') THEN BEGIN 
          state.display_coord_sys = strmid((*state.astr_ptr).ctype[0], 0, 4)
          state.display_equinox = state.equinox
          atv_gettrack          ; refresh coordinate window
       ENDIF 
    END 

    
; Help options:            
    'ATV Help': atv_help
    
    else: print, 'Unknown event in file menu!'
endcase

; Need to test whether atv is still alive, since the quit option
; might have been selected.
if (xregistered('atv', /noshow)) then atv_resetwindow


end

;--------------------------------------------------------------------

pro atv_draw_event, event

; top-level event handler for draw widget events

common atv_state

if (!d.name NE state.graphicsdevice) then return

if (event.type EQ 0 or event.type EQ 1 or event.type EQ 2) then begin
    case state.mousemode of
        'color':  atv_draw_color_event, event
        'zoom':   atv_draw_zoom_event, event
        'blink':  atv_draw_blink_event, event
        'imexam': atv_draw_phot_event, event
        'measure': atv_draw_vector_event, event, /measure
        'region': atv_draw_region_event, event
        'vector': atv_draw_vector_event, event
        'print': atv_draw_print_event, event
    endcase
endif

if (event.type EQ 5 or event.type EQ 6) then $
  atv_draw_keyboard_event, event

if (xregistered('atv', /noshow)) then $
  widget_control, state.draw_widget_id, /sensitive, /input_focus

end

;--------------------------------------------------------------------

pro atv_draw_color_event, event

; Event handler for color mode

common atv_state
common atv_images

;if (!d.name NE state.graphicsdevice) then return

case event.type of
    0: begin           ; button press
        if (event.press EQ 1) then begin
            state.cstretch = 1
			atv_setwindow,state.draw_window_id,/nostretch
            atv_stretchct, event.x, event.y, /getcursor
            atv_resetwindow
            atv_colorbar
        endif else begin
            atv_zoom, 'none', /recenter
        endelse
    end
    1: begin
        state.cstretch = 0  ; button release
        if (state.bitdepth EQ 24) then atv_refresh
        atv_draw_motion_event, event
    end
    2: begin                ; motion event
        if (state.cstretch EQ 1) then begin
			atv_setwindow,state.draw_window_id,/nostretch
            atv_stretchct, event.x, event.y, /getcursor
			atv_resetwindow
            if (state.bitdepth EQ 24) then atv_refresh, /fast
        endif else begin 
            atv_draw_motion_event, event
        endelse
    end 
endcase

widget_control, state.draw_widget_id, /sensitive, /input_focus

end

;--------------------------------------------------------------------

pro atv_draw_keyboard_event, event


common atv_state
common atv_images
common atv_color

forward_function atv_wcsstring

; Only want to look for key presses, not key releases.
if (event.release EQ 1) then return

if (event.type EQ 5) then begin

    eventchar = string(event.ch)

    if (!d.name NE state.graphicsdevice and eventchar NE 'q') then return
    if (state.bitdepth EQ 24) then true = 1 else true = 0
    
    case eventchar of
        '1': atv_move_cursor, eventchar
        '2': atv_move_cursor, eventchar
        '3': atv_move_cursor, eventchar
        '4': atv_move_cursor, eventchar
        '6': atv_move_cursor, eventchar
        '7': atv_move_cursor, eventchar
        '8': atv_move_cursor, eventchar
        '9': atv_move_cursor, eventchar
        'r': atv_rowplot, /newcoord
        'c': atv_colplot, /newcoord
        's': atv_surfplot, /newcoord
        't': atv_contourplot, /newcoord
        'h': atv_histplot, /newcoord
        'p': atv_apphot
        'i': atv_showstats
        'm': atv_changemode
		'b': atv_changeimage,/previous
		'n': atv_changeimage,/next
		'-': atv_zoom, 'out'
		'+': atv_zoom, 'in'
		'=': atv_zoom, 'in'
		'P': begin
			print,state.coord[0],state.coord[1]
  			if (ptr_valid(state.astr_ptr)) then if (state.wcstype EQ 'angle') then begin
		      xy2ad, state.coord[0], state.coord[1], *(state.astr_ptr), lon, lat
		
		      wcsstring = atv_wcsstring(lon, lat, (*state.astr_ptr).ctype,  $
		                                state.equinox, state.display_coord_sys, $
		                                state.display_equinox, state.display_base60)
			  print, wcsstring
			endif
		
			end
        '!': begin  
            atv_setwindow, state.draw_window_id
            blink_image1 = tvrd(true = true) 
            atv_resetwindow
        end
        '@': begin  
            atv_setwindow, state.draw_window_id
            blink_image2 = tvrd(true = true) 
            atv_resetwindow
        end
        '#': begin  
            atv_setwindow, state.draw_window_id
            blink_image3 = tvrd(true = true) 
            atv_resetwindow
        end
        'q': if (state.activator EQ 0) then atv_shutdown $
        else state.activator = 0
        'Q': if (state.activator EQ 0) then atv_shutdown $
        else state.activator = 0
 
		'l': atv_setscaling, 0
		'g': atv_setscaling, 1
		'h': atv_setscaling, 4
		'a': begin 
				atv_autoscale
	        	atv_displayall
			end
		'z': begin
				state.max_value = state.image_max
				state.min_value = 0.0
				state.scaling=1
	        	atv_set_minmax
				atv_displayall
			end
		'f': begin 
	        state.min_value = state.image_min
	        state.max_value = state.image_max
	        if state.min_value GE state.max_value then begin
	            state.min_value = state.max_value - 1
	            state.max_value = state.max_value + 1
	        endif
	        atv_set_minmax
	        atv_displayall
			end
    'G': atv_setregion
    'H': atv_histplot
    'J': atv_gaussrowplot
    'K': atv_gausscolplot
    'd': if (state.image_size[2] gt 1) then begin 
           atv_slice3dplot 
         endif else begin
           atv_message, 'Image must be 3D for pixel slice', $
             msgtype='error', /window
           return
         endelse
    'Z': atv_pixtable
	'E': atverase ; Erase
	'B': begin ; Blink
		; FiXME this doesn't work??
		   atv_setwindow, state.draw_window_id
       	   blink_image1 = tvrd(true = true) 
           state.title_blink1 = state.window_title
		   ;state.mousemode = 'blink'
		   	; need to update the mouse mode selector bar here

		 end

    	else:  print,eventchar ;any other key press does nothing
    endcase
 
endif

; Starting with IDL 6.0, can generate events on arrow keys:
if (event.type EQ 6) then begin
    case event.key of
        5: atv_move_cursor, '4'
        6: atv_move_cursor, '6'
        7: atv_move_cursor, '8'
        8: atv_move_cursor, '2'
        else:
    endcase
endif 
   
if (xregistered('atv', /noshow)) then $
  widget_control, state.draw_widget_id, /sensitive, /input_focus

end

;-------------------------------------------------------------------

pro atv_activate

; This routine is a workaround to use when you hit an error message or
; a "stop" command in another program while running atv.  If you want
; atv to become active again without typing "retall" and losing your
; current session variables, type "atv_activate" to temporarily
; activate atv again.  This will de-activate the command line but
; allow atv to be used until you hit "q" or click "done" in atv.

; Also, if you need to call atv from a command-line idl program and
; have that program wait until you're done looking at an image in atv
; before moving on to its next step, you can call atv_activate after
; sending your image to atv.  This will make your external program
; stop until you quit out of atv_activate mode.

common atv_state

if (not(xregistered('atv', /noshow))) then begin
    print, 'No ATV window currently exists.'
    return
endif

state.activator = 1
activator = 1

while (activator EQ 1) do begin

    wait, 0.01
    void = widget_event(/nowait)
    
; If atv is killed by the window manager, then by the time we get here
; the state structure has already been destroyed by atv_shutdown.
    if (size(state, /type) NE 8) then begin
        activator = 0
    endif else begin
        activator = state.activator
    endelse
    
endwhile

widget_control, /hourglass

end

;-------------------------------------------------------------------

pro atv_changemode

; Use 'm' keypress to cycle through mouse modes

common atv_state

case state.mousemode of
    'color': begin
        state.mousemode = 'zoom'
        widget_control, state.mode_droplist_id, set_droplist_select=1
    end
    'zoom': begin
        state.mousemode = 'blink'
        widget_control, state.mode_droplist_id, set_droplist_select=2
    end
    'blink': begin
        state.mousemode = 'imexam'
        widget_control, state.mode_droplist_id, set_droplist_select=3
    end
    'imexam': begin
        state.mousemode = 'vector'
        widget_control, state.mode_droplist_id, set_droplist_select=4
    end
    'vector': begin
        state.mousemode = 'color'
        widget_control, state.mode_droplist_id, set_droplist_select=0
    end
endcase


end

;------------------------------------------------------------------

pro atv_draw_zoom_event, event

; Event handler for zoom mode

common atv_state
 
if (!d.name NE state.graphicsdevice) then return

if (event.type EQ 0) then begin 
    case event.press of
        1: atv_zoom, 'in', /recenter
        2: atv_zoom, 'none', /recenter
        4: atv_zoom, 'out', /recenter
    endcase
endif

if (event.type EQ 2) then atv_draw_motion_event, event

if (xregistered('atv', /noshow)) then $
  widget_control, state.draw_widget_id, /sensitive, /input_focus

end

;---------------------------------------------------------------------

pro atv_draw_blink_event, event

; Event handler for blink mode

common atv_state
common atv_images

if (!d.name NE state.graphicsdevice) then return
if (state.bitdepth EQ 24) then true = 1 else true = 0

case event.type of
    0: begin                    ; button press
        atv_setwindow, state.draw_window_id
                                ; define the unblink image if needed
        if ((state.newrefresh EQ 1) AND (state.blinks EQ 0)) then begin
            unblink_image = tvrd(true = true)
            state.newrefresh = 0
        endif
        
        case event.press of
            1: if n_elements(blink_image1) GT 1 then begin
              tv, blink_image1, true = true
              widget_control, state.base_id, $
                tlb_set_title = state.title_blink1
               endif
            2: if n_elements(blink_image2) GT 1 then begin
              tv, blink_image2, true = true
              widget_control, state.base_id, $
                tlb_set_title = state.title_blink2              
               endif
            4: if n_elements(blink_image3) GT 1 then begin
              tv, blink_image3, true = true  
              widget_control, state.base_id, $
                tlb_set_title = state.title_blink3              
               endif
            else: event.press = 0 ; in case of errors
        endcase
        state.blinks = (state.blinks + event.press) < 7
    end
    
    1: begin                    ; button release
        if (n_elements(unblink_image) EQ 0) then return ; just in case
        atv_setwindow, state.draw_window_id
        state.blinks = (state.blinks - event.release) > 0
        case state.blinks of
            0: begin 
               tv, unblink_image, true = true
               widget_control, state.base_id, $
                 tlb_set_title = state.window_title
            end
            1: if n_elements(blink_image1) GT 1 then begin
                 tv, blink_image1, true = true
                 widget_control, state.base_id, $
                   tlb_set_title = state.title_blink1
               endif else begin
                 tv, unblink_image, true = true
               endelse
            2: if n_elements(blink_image2) GT 1 then begin
                 tv, blink_image2, true = true
                 widget_control, state.base_id, $
                   tlb_set_title = state.title_blink2
               endif else begin
                 tv, unblink_image, true = true
               endelse
           3: if n_elements(blink_image1) GT 1 then begin
                tv, blink_image1, true = true
                widget_control, state.base_id, $
                  tlb_set_title = state.title_blink1
              endif else if n_elements(blink_image2) GT 1 then begin
                tv, blink_image2, true = true
              endif else begin
                tv, unblink_image, true = true
            endelse
            4: if n_elements(blink_image3) GT 1 then begin
                 tv, blink_image3, true = true
                 widget_control, state.base_id, $
                   tlb_set_title = state.window_title
               endif else begin
                 tv, unblink_image, true = true
               endelse
            5: if n_elements(blink_image1) GT 1 then begin
                tv, blink_image1, true = true
                widget_control, state.base_id, $
                  tlb_set_title = state.title_blink1
            endif else if n_elements(blink_image3) GT 1 then begin
                tv, blink_image3, true = true
            endif else begin
                tv, unblink_image, true = true
            endelse
            6: if n_elements(blink_image2) GT 1 then begin
                tv, blink_image2, true = true
                widget_control, state.base_id, $
                  tlb_set_title = state.title_blink2
            endif else if n_elements(blink_image4) GT 1 then begin
                tv, blink_image4, true = true
                widget_control, state.base_id, $
                  tlb_set_title = state.title_blink3
            endif else begin
                tv, unblink_image, true = true
            endelse
            else: begin         ; check for errors
                state.blinks = 0
                tv, unblink_image, true = true
            end
        endcase
    end
    2: atv_draw_motion_event, event ; motion event
endcase

widget_control, state.draw_widget_id, /sensitive, /input_focus
atv_resetwindow

end

;-------------------------------------------------------------------

pro atv_draw_phot_event, event

; Event handler for ImExam mode

common atv_state
common atv_images

if (!d.name NE state.graphicsdevice) then return

if (event.type EQ 0) then begin
    case event.press of
        1: atv_apphot
        2: atv_zoom, 'none', /recenter
        4: atv_showstats
        else: 
    endcase
endif

if (event.type EQ 2) then atv_draw_motion_event, event

widget_control, state.draw_widget_id, /sensitive, /input_focus


end

;--------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro atv_draw_region_event, event

; Event handler for Region mode
common atv_state
common atv_images

if (!d.name NE state.graphicsdevice) then return

case (event.type) of
    0: begin
      case event.press of
        1: atv_setregion
        2: 
        4: 
        else: 
      endcase
   end
   2: atv_draw_motion_event, event
   5: atv_draw_keyboard_event, event
else:
endcase

widget_control, state.draw_widget_id, /input_focus

widget_control, state.draw_widget_id, /clear_events
;widget_control, state.keyboard_text_id, /sensitive, /input_focus

end

;--------------------------------------------------------------------

pro atv_draw_print_event, event
;  Print out the pixel coords clicked on.

common atv_state
common atv_images

if (!d.name NE state.graphicsdevice) then return

atv_setwindow, state.draw_window_id

if (event.type EQ 0) then begin
	print,state.coord[0],state.coord[1]
endif

if (event.type EQ 2) then atv_draw_motion_event, event

widget_control, state.draw_widget_id, /input_focus

widget_control, state.draw_widget_id, /clear_events
;widget_control, state.keyboard_text_id, /sensitive, /input_focus

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;--------------------------------------------------------------------

pro atv_draw_motion_event, event

; Event handler for motion events in draw window

common atv_state

if (!d.name NE state.graphicsdevice) then return

tmp_event = [event.x, event.y]
state.coord = $
  round( (0.5 >  ((tmp_event / state.zoom_factor) + state.offset) $
          < (state.image_size[0:1] - 0.5) ) - 0.5)
atv_gettrack

widget_control, state.draw_widget_id, /sensitive, /input_focus

;if atv_pixtable on, then create a 5x5 array of pixel values and the 
;X & Y location strings that are fed to the pixel table

if (xregistered('atv_pixtable', /noshow)) then atv_pixtable_update

end
;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
;--------------------------------------------------------------------

pro atv_draw_vector_event, event, measure=measure

; Check for left button press/depress, then get coords at point 1 and 
; point 2.  Call atv_lineplot.  Calculate vector distance between
; endpoints and plot Vector Distance vs. Pixel Value with atv_vectorplot
; 
; Measure mode added. M. Perrin

common atv_state
common atv_images

if (!d.name NE state.graphicsdevice) then return

;atv_setwindow, state.draw_window_id

case event.type of
    0: begin           ; button press
        if (event.press EQ 1) then begin  ; left button press
          state.vector_coord1[0] = state.coord[0]
          state.vector_coord1[1] = state.coord[1]
          state.vectorstart = [event.x, event.y]
          atv_drawvector, event, measure=measure
          state.vectorpress = 1
        endif
    end
    1: begin           ; button release
        if (event.release EQ 1) then begin  ; left button release
            state.vectorpress = 0
            state.vector_coord2[0] = state.coord[0]
            state.vector_coord2[1] = state.coord[1]
            atv_drawvector, event, measure=measure
            if ~(keyword_set(measure)) then atv_vectorplot, /newcoord
        endif
    end
    2: begin  ; motion event
        atv_draw_motion_event, event 
        if (state.vectorpress EQ 1) then atv_drawvector, event, measure=measure
    end

    5: atv_draw_keyboard_event, event     ; keyboard event
    else:
endcase

widget_control, state.draw_widget_id, /sensitive, /input_focus

end

;----------------------------------------------------------------------

pro atv_drawvector, event, measure=measure

common atv_state

; button press: create initial pixmap and start drawing vector
if (event.type EQ 0) then begin
    window, /free, xsize = state.draw_window_size[0], $
      ysize = state.draw_window_size[1], /pixmap
    state.vector_pixmap_id = !d.window
    device, copy=[0, 0, state.draw_window_size[0], $
                  state.draw_window_size[1], 0, 0, state.draw_window_id]
    atv_resetwindow
endif

; button release: redisplay initial image
if (event.type EQ 1) then begin
    atv_setwindow, state.draw_window_id
    device, copy=[0, 0, state.draw_window_size[0], $
                  state.draw_window_size[1], 0, 0, state.vector_pixmap_id]
	if keyword_set(measure) then atv_labelmeasure,event
    atv_resetwindow
    wdelete, state.vector_pixmap_id
endif

; motion event: redraw with new vector
if (event.type EQ 2) then begin
    atv_setwindow, state.draw_window_id

    device, copy=[0, 0, state.draw_window_size[0], $
                  state.draw_window_size[1], 0, 0, state.vector_pixmap_id]
    xvector = [state.vectorstart[0], event.x]
    yvector = [state.vectorstart[1], event.y]

    plots, xvector, yvector, /device, color = state.box_color
	if keyword_set(measure) then atv_labelmeasure,event

    atv_resetwindow
endif

end

;----------------------------------------------------------------------

pro atv_labelmeasure, event
	common atv_state

	; for mouse button release events, write a permanent version of the
	; measure vector. For other events, we let drawvector handle drawing the
	; temporary version of the vector.
	if (event.type EQ 1) then begin
		atv_resetwindow ; needed before atvplot to not stomp on device.decomposed
		atvplot,[state.vector_coord1[0],state.coord[0]],[state.vector_coord1[1],state.coord[1]],color = state.box_color
	endif

	atv_resetwindow ; needed before atvxyouts to not stomp on device.decomposed
	if state.wcstype eq 'none' then begin
		distance =sqrt( ((state.vector_coord1[0]-state.coord[0]))^2 $
		               +((state.vector_coord1[1]-state.coord[1]))^2 )
		atvxyouts,(state.vector_coord1[0]+state.coord[0])/2+3,(state.vector_coord1[1]+state.coord[1])/2+3,$
			string(distance,format="(g5.4)")+" pixels",charsize=2
		;print,"distance is "+string(distance)+" pixels"

	endif else begin

		
;;		distance =sqrt( ((state.vector_coord1[0]-state.coord[0])* ((*(state.astr_ptr)).cd[0,0]) * (*(state.astr_ptr)).cdelt[0] )^2 $
;;		               +((state.vector_coord1[1]-state.coord[1])* ((*(state.astr_ptr)).cd[1,1]) * (*(state.astr_ptr)).cdelt[1] )^2 )
;;		distance = distance*60*60 ; convert from degrees to arcsec
;;		dx = ((state.vector_coord1[0]-state.coord[0])* ((*(state.astr_ptr)).cd[0,0]) * (*(state.astr_ptr)).cdelt[0] )
;;		dy = ((state.vector_coord1[1]-state.coord[1])* ((*(state.astr_ptr)).cd[1,1]) * (*(state.astr_ptr)).cdelt[1] )
;;		pa = !radeg* atan(state.vector_coord1[0]-state.coord[0],$
;;			-(state.vector_coord1[1]-state.coord[1]))
;;
		; TODO FIXME
		; it would be better to use xy2ad to convert things here.
		; then we could get more correct results for complicated projections.
		xy2ad, state.vector_coord1[0], state.vector_coord1[1],*(state.astr_ptr) ,  startra, startdec
		xy2ad, state.coord[0], state.coord[1], *(state.astr_ptr) ,  stopra, stopdec
		gcirc, 1, startra/15, startdec, stopra/15, stopdec, distance
		posang, 1, startra/15, startdec, stopra/15, stopdec, pa

		;; reference PA to north
		;getrot, *state.head_ptr, northangle
		;;print, "raw pa", pa, "north", northangle
		;pa -= northangle
		
		if pa lt 0 then pa += 360.0
		if distance lt 100.0 then formatstr = "(g7.4)" else formatstr= "(g6.5)"
		atvxyouts,(state.vector_coord1[0]+state.coord[0])/2+3,(state.vector_coord1[1]+state.coord[1])/2+3,$
			string(distance,format=formatstr)+" arcsec,!C"+string(pa,format=formatstr)+" degr",charsize=2
		;print,"distance is "+string(distance)+" arcseconds"
	endelse
	; For Motion events, don't save annotation 
	if (event.type EQ 2) then atverase,1,/norefresh
	

end
;----------------------------------------------------------------------

pro atv_pan_event, event

; event procedure for moving the box around in the pan window

common atv_state

if (!d.name NE state.graphicsdevice) then return

case event.type of
    0: begin                     ; button press
        widget_control, state.pan_widget_id, draw_motion_events = 1
        atv_pantrack, event
    end
    1: begin                     ; button release
        widget_control, state.pan_widget_id, draw_motion_events = 0
        widget_control, state.pan_widget_id, /clear_events
        atv_pantrack, event
        atv_refresh
    end
    2: begin
        atv_pantrack, event     ; motion event
        widget_control, state.pan_widget_id, /clear_events
    end
    else:
endcase

end

;--------------------------------------------------------------------

pro atv_event, event

; Main event loop for atv top-level base, and for all the buttons.

common atv_state
common atv_images
common atv_color
common atv_pdata

widget_control, event.id, get_uvalue = uvalue

if (!d.name NE state.graphicsdevice and uvalue NE 'done') then return

; Get currently active window
atv_getwindow

case uvalue of

    'atv_base': begin     
        c = where(tag_names(event) EQ 'ENTER', count)
        if (count EQ 0) then begin       ; resize event
            atv_resize
            atv_refresh
        endif
    end

    'mode': case event.index of
        0: state.mousemode = 'color'
        1: state.mousemode = 'zoom'
        2: state.mousemode = 'blink'
        3: state.mousemode = 'imexam'
		4: state.mousemode = 'measure'
		5: state.mousemode = 'region'
		6: state.mousemode = 'vector'
		7: state.mousemode = 'print'
        else: print, 'Unknown mouse mode!'
    endcase

    'invert': begin                  ; invert the color table
        state.invert_colormap = abs(state.invert_colormap - 1)

        r_vector = reverse(r_vector)
        g_vector = reverse(g_vector)
        b_vector = reverse(b_vector)

        atv_setwindow, state.draw_window_id
        atv_stretchct
        atv_resetwindow

; For 24-bit color, need to refresh display after stretching color
; map.  Can refresh in /fast mode if there are no overplots
        if (state.bitdepth EQ 24) then begin
            if ptr_valid(plot_ptr[1]) then begin
                atv_refresh
            endif else begin
                atv_refresh, /fast
            endelse
        endif
    end

    'restretch_button': atv_restretch

    'min_text': begin     ; text entry in 'min = ' box
        atv_get_minmax, uvalue, event.value
        atv_displayall
    end

    'max_text': begin     ; text entry in 'max = ' box
        atv_get_minmax, uvalue, event.value
        atv_displayall
    end

    'nonlin_text': begin     ; text entry in 'nonlin = ' box
        atv_get_minmax, uvalue, event.value
        atv_displayall
    end

    'curimnum_text':begin    ; text entry in 'Image #=' box
		atv_changeimage,event.value
     end

    'curimnum_slidebar':begin    ; slidebar controlling cur_image_num
		atv_changeimage,event.value
    end

    'curimnum_minmaxmode': case event.index of
        3: begin
             state.curimnum_minmaxmode = 'Zero/Max'
             state.min_value = 0
             state.max_value = state.image_max
             atv_set_minmax
             atv_displayall
        end
        2: begin
             state.curimnum_minmaxmode = 'Min/Max'
             state.min_value = state.image_min
             state.max_value = state.image_max
             atv_set_minmax
             atv_displayall
        end
        1: begin
             state.curimnum_minmaxmode = 'AutoScale'
             atv_autoscale
             atv_set_minmax
             atv_displayall
        end
        0: state.curimnum_minmaxmode = 'Constant'
        else: print, 'Unknown Min/Max mode for changing cur_image_num!'
    endcase

    'autoscale_button': begin   ; autoscale the image
        atv_autoscale
        atv_displayall
    end

    'full_range': begin    ; display the full intensity range
        state.min_value = state.image_min
        state.max_value = state.image_max
        if state.min_value GE state.max_value then begin
            state.min_value = state.max_value - 1
            state.max_value = state.max_value + 1
        endif
        atv_set_minmax
        atv_displayall
    end

    'zoom_in':  atv_zoom, 'in'         ; zoom buttons
    'zoom_out': atv_zoom, 'out'
    'zoom_one': atv_zoom, 'one'

    'center': begin   ; center image and preserve current zoom level
        state.centerpix = round(state.image_size[0:1] / 2.)
        atv_refresh
    end

    'fullview': atv_fullview

    'done':  if (state.activator EQ 0) then atv_shutdown $
      else state.activator = 0

    else:  print, 'No match for uvalue....'  ; bad news if this happens

endcase
end

;----------------------------------------------------------------------

pro atv_message, msg_txt, msgtype=msgtype, window=window

; Routine to display an error or warning message.  Message can be
; displayed either to the IDL command line or to a popup window,
; depending on whether /window is set.
; msgtype must be 'warning', 'error', or 'information'.

common atv_state

if (n_elements(window) EQ 0) then window = 0

if (window EQ 1) then begin  ; print message to popup window
    case msgtype of
        'warning': t = dialog_message(msg_txt, dialog_parent = state.base_id)
        'error': t = $
          dialog_message(msg_txt,/error,dialog_parent=state.base_id)
        'information': t = $
          dialog_message(msg_txt,/information,dialog_parent=state.base_id)
        else:
    endcase
endif else begin           ;  print message to IDL console
    message = strcompress(strupcase(msgtype) + ': ' + msg_txt)
    print, message
endelse

end

;-----------------------------------------------------------------------
;      main atv routines for scaling, displaying, cursor tracking...
;-----------------------------------------------------------------------

pro atv_displayall

; Call the routines to scale the image, make the pan image, and
; re-display everything.  Use this if the scaling changes (log/
; linear/ histeq), or if min or max are changed, or if a new image is
; passed to atv.  If the display image has just been moved around or
; zoomed without a change in scaling, then just call atv_refresh
; rather than this routine.

atv_scaleimage
atv_makepan
atv_refresh


end

;---------------------------------------------------------------------

pro atv_refresh, fast = fast

; Make the display image from the scaled_image, and redisplay the pan
; image and tracking image.
; The /fast option skips the steps where the display_image is
; recalculated from the main_image.  The /fast option is used in 24
; bit color mode, when the color map has been stretched but everything
; else stays the same.

common atv_state
common atv_images

atv_getwindow
if (not(keyword_set(fast))) then begin
    atv_getoffset
    atv_getdisplay
    atv_displaymain
    atv_plotall
endif else begin
    atv_displaymain
endelse

; redisplay the pan image and plot the boundary box
atv_setwindow, state.pan_pixmap
erase
tv, pan_image, state.pan_offset[0], state.pan_offset[1]
atv_resetwindow

atv_setwindow, state.pan_window_id
if (not(keyword_set(fast))) then erase
tv, pan_image, state.pan_offset[0], state.pan_offset[1]
atv_resetwindow
atv_drawbox, /norefresh

if (state.bitdepth EQ 24) then atv_colorbar

; redisplay the tracking image
if (not(keyword_set(fast))) then atv_gettrack

atv_resetwindow

state.newrefresh = 1


end

;--------------------------------------------------------------------
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro atv_getdisplay

; make the display image from the scaled image by applying the zoom
; factor and matching to the size of the draw window, and display the
; image.

common atv_state
common atv_images

widget_control, /hourglass   

display_image = bytarr(state.draw_window_size[0], state.draw_window_size[1])

view_min = round(state.centerpix - $
                  (0.5 * state.draw_window_size / state.zoom_factor))
view_max = round(view_min + state.draw_window_size / state.zoom_factor)

view_min = (0 > view_min < (state.image_size[0:1] - 1)) 
view_max = (0 > view_max < (state.image_size[0:1] - 1)) 

newsize = round( (view_max - view_min + 1) * state.zoom_factor) > 1
startpos = abs( round(state.offset * state.zoom_factor) < 0)

; Use interp & center keywords to congrid for zoomfactor < 1 :
; improvement contributed by N. Cunningham, added 4/14/06
if (state.zoom_factor LT 1.0) then begin
    tmp_image = congrid(scaled_image[view_min[0]:view_max[0], $
                                   view_min[1]:view_max[1]], $
                      newsize[0], newsize[1], /center, /interp)
endif else begin
    tmp_image = congrid(scaled_image[view_min[0]:view_max[0], $
                                     view_min[1]:view_max[1]], $
                        newsize[0], newsize[1])
endelse


xmax = newsize[0] < (state.draw_window_size[0] - startpos[0])
ymax = newsize[1] < (state.draw_window_size[1] - startpos[1])

display_image[startpos[0], startpos[1]] = tmp_image[0:xmax-1, 0:ymax-1]
tmp_image = 0

end

;-----------------------------------------------------------------------

pro atv_displaymain

; Display the main image and overplots

common atv_state
common atv_images

atv_setwindow, state.draw_window_id
tv, display_image
atv_resetwindow

end

;--------------------------------------------------------------------

pro atv_getoffset
common atv_state

; Routine to calculate the display offset for the current value of
; state.centerpix, which is the central pixel in the display window.

state.offset = $
  round( state.centerpix - $
         (0.5 * state.draw_window_size / state.zoom_factor) )

end

;----------------------------------------------------------------------


pro atv_makepan

; Make the 'pan' image that shows a miniature version of the full image.

common atv_state
common atv_images

sizeratio = state.image_size[1] / state.image_size[0]

if (sizeratio GE 1) then begin
    state.pan_scale = float(state.pan_window_size) / float(state.image_size[1])
endif else begin
    state.pan_scale = float(state.pan_window_size) / float(state.image_size[0])
endelse

tmp_image = $
  scaled_image[0:state.image_size[0]-1, 0:state.image_size[1]-1]

pan_image = $
  congrid(tmp_image, round(state.pan_scale * state.image_size[0])>1, $
          round(state.pan_scale * state.image_size[1])>1, $
         /center, /interp) 

state.pan_offset[0] = round((state.pan_window_size - (size(pan_image))[1]) / 2)
state.pan_offset[1] = round((state.pan_window_size - (size(pan_image))[2]) / 2)

end

;----------------------------------------------------------------------


pro atv_move_cursor, direction

; Use keypad arrow keys to step cursor one pixel at a time.
; Get the new track image, and update the cursor position.

common atv_state

i = 1L

case direction of
    '2': state.coord[1] = max([state.coord[1] - i, 0])
    '4': state.coord[0] = max([state.coord[0] - i, 0])
    '8': state.coord[1] = min([state.coord[1] + i, state.image_size[1] - i])
    '6': state.coord[0] = min([state.coord[0] + i, state.image_size[0] - i])
    '7': begin
        state.coord[1] = min([state.coord[1] + i, state.image_size[1] - i])
        state.coord[0] = max([state.coord[0] - i, 0])
    end
    '9': begin
        state.coord[1] = min([state.coord[1] + i, state.image_size[1] - i])
        state.coord[0] = min([state.coord[0] + i, state.image_size[0] - i])
    end
    '3': begin
        state.coord[1] = max([state.coord[1] - i, 0])
        state.coord[0] = min([state.coord[0] + i, state.image_size[0] - i])
    end
    '1': begin
        state.coord[1] = max([state.coord[1] - i, 0])
        state.coord[0] = max([state.coord[0] - i, 0])
    end

endcase

newpos = (state.coord - state.offset + 0.5) * state.zoom_factor

atv_setwindow,  state.draw_window_id
tvcrs, newpos[0], newpos[1], /device
atv_resetwindow

atv_gettrack

; If pixel table widget is open, update pixel values and cursor position
if (xregistered('atv_pixtable', /noshow)) then atv_pixtable_update


; Prevent the cursor move from causing a mouse event in the draw window
widget_control, state.draw_widget_id, /clear_events

atv_resetwindow

end

;----------------------------------------------------------------------

pro atv_set_minmax

; Updates the min and max text boxes with new values.

common atv_state

widget_control, state.min_text_id, set_value = string(state.min_value)
widget_control, state.max_text_id, set_value = string(state.max_value)

end

;----------------------------------------------------------------------

pro atv_get_minmax, uvalue, newvalue

; Change the min and max state variables when user inputs new numbers
; in the text boxes. 

common atv_state

case uvalue of

    'min_text': begin
        if (newvalue LT state.max_value) then begin
            state.min_value = newvalue
        endif
    end

    'max_text': begin
        if (newvalue GT state.min_value) then begin
            state.max_value = newvalue
        endif
    end
    'nonlin_text': if (newvalue GT 0) then state.nonlinearity_value = newvalue
endcase

atv_set_minmax

end

;--------------------------------------------------------------------

pro atv_zoom, zchange, recenter = recenter
common atv_state

; Routine to do zoom in/out and recentering of image.  The /recenter
; option sets the new display center to the current cursor position.

case zchange of
    'in':    state.zoom_level = (state.zoom_level + 1) < 6
    'out':   begin
        sizeratio = fix(min(state.image_size[0:1]) / 16.) > 1
        minzoom = -1.*fix(alog(sizeratio)/alog(2.0))
        state.zoom_level = (state.zoom_level - 1) > minzoom
    end
    'onesixteenth': state.zoom_level =  -4
    'oneeighth': state.zoom_level =  -3
    'onefourth': state.zoom_level =  -2
    'onehalf': state.zoom_level =  -1
    'two':   state.zoom_level =  1
    'four':  state.zoom_level =  2
    'eight': state.zoom_level =  3
    'sixteen': state.zoom_level = 4
    'one':   state.zoom_level =  0
    'none':  ; no change to zoom level: recenter on current mouse position
    else:  print,  'problem in atv_zoom!'
endcase

state.zoom_factor = (2.0)^(state.zoom_level)

if (n_elements(recenter) GT 0) then begin
    state.centerpix = state.coord
    atv_getoffset
endif

atv_refresh

if (n_elements(recenter) GT 0) then begin
    newpos = (state.coord - state.offset + 0.5) * state.zoom_factor
    atv_setwindow,  state.draw_window_id
    tvcrs, newpos[0], newpos[1], /device 
    atv_resetwindow
    atv_gettrack
endif

atv_resetwindow

end

;-----------------------------------------------------------------------

pro atv_fullview
common atv_state

; set the zoom level so that the full image fits in the display window

sizeratio = float(state.image_size) / float(state.draw_window_size)
maxratio = (max(sizeratio)) 

state.zoom_level = floor((alog(maxratio) / alog(2.0)) * (-1))
state.zoom_factor = (2.0)^(state.zoom_level)

; recenter
state.centerpix = round(state.image_size / 2.)

atv_refresh

atv_resetwindow

end

;----------------------------------------------------------------------

pro atv_invert, ichange
common atv_state
common atv_images

; Routine to do image axis-inversion (X,Y,X&Y)

case ichange of
    'x': begin
         if ptr_valid(state.head_ptr) then begin
           if (state.invert_image eq 'none') then $
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 1, /silent
           if (state.invert_image eq 'x') then return
           if (state.invert_image eq 'y') then begin
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 2, /silent             
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 1, /silent
           endif
           if (state.invert_image eq 'xy') then $
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 2, /silent   

           head = *(state.head_ptr)
           atv_setheader, head

         endif else begin
           if (state.invert_image eq 'none') then $
             main_image = reverse(main_image)
           if (state.invert_image eq 'x') then return
           if (state.invert_image eq 'y') then begin
             main_image = reverse(main_image,2)
             main_image = reverse(main_image)   
           endif
           if (state.invert_image eq 'xy') then $
             main_image = reverse(main_image,2)
         endelse

         state.invert_image = 'x' 
    end
    
    'y': begin
         if ptr_valid(state.head_ptr) then begin
           if (state.invert_image eq 'none') then $
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 2, /silent
           if (state.invert_image eq 'y') then return
           if (state.invert_image eq 'x') then begin
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 1, /silent             
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 2, /silent
           endif
           if (state.invert_image eq 'xy') then $
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 1, /silent   

           head = *(state.head_ptr)
           atv_setheader, head

         endif else begin
           if (state.invert_image eq 'none') then $
             main_image=reverse(main_image,2)
           if (state.invert_image eq 'x') then begin
             main_image = reverse(main_image)
             main_image=reverse(main_image,2)
           endif
           if (state.invert_image eq 'y') then return
           if (state.invert_image eq 'xy') then $
             main_image = reverse(main_image)
         endelse

         state.invert_image = 'y'
    end
    
    
    'xy': begin
         if ptr_valid(state.head_ptr) then begin
           if (state.invert_image eq 'none') then begin
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 1, /silent
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 2, /silent
           endif
           if (state.invert_image eq 'x') then $
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 2, /silent
           if (state.invert_image eq 'y') then $
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 1, /silent             
           if (state.invert_image eq 'xy') then return

           head = *(state.head_ptr)
           atv_setheader, head

         endif else begin
           if (state.invert_image eq 'none') then begin
             main_image = reverse(main_image)
             main_image = reverse(main_image,2)
           endif
           if (state.invert_image eq 'x') then $
             main_image = reverse(main_image,2)
           if (state.invert_image eq 'y') then $
             main_image = reverse(main_image)
           if (state.invert_image eq 'xy') then return
         endelse

         state.invert_image = 'xy'
    end
    'none': begin ; do not invert; revert to normal (X,Y) axes view
         if ptr_valid(state.head_ptr) then begin
           if (state.invert_image eq 'none') then return
           if (state.invert_image eq 'x') then $
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 1, /silent
           if (state.invert_image eq 'y') then $
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 2, /silent             
           if (state.invert_image eq 'xy') then begin
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 1, /silent
             hreverse, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), 2, /silent
           endif

           head = *(state.head_ptr)
           atv_setheader, head

         endif else begin
           if (state.invert_image eq 'none') then return
           if (state.invert_image eq 'x') then $
             main_image = reverse(main_image)
           if (state.invert_image eq 'y') then $
             main_image = reverse(main_image,2)
           if (state.invert_image eq 'xy') then begin
             main_image = reverse(main_image)
             main_image = reverse(main_image,2)
           endif
         endelse

         state.invert_image = 'none'
    end
    
    else:  print,  'problem in atv_invert!'
endcase

atv_getstats, /align, /noerase

;Redisplay inverted image with current zoom, update pan, and refresh image
atv_displayall

;make sure that the image arrays are updated for line/column plots, etc.
atv_resetwindow

end

;------------------------------------------------------------------

pro atv_rotate, rchange, get_angle=get_angle
common atv_state
common atv_images

; Routine to do image rotation

; If /get_angle set, create widget to enter rotation angle

widget_control, /hourglass

if (keyword_set(get_angle)) then begin

  formdesc = ['0, float, 0.0, label_left=Rotation Angle: ', $
              '1, base, , row', $
              '0, button, Cancel, quit', $
              '0, button, Rotate, quit']    

  textform = cw_form(formdesc, /column, title = 'Rotate')

  if (textform.tag2 EQ 1) then return
  if (textform.tag3 EQ 1) then rchange = textform.tag0

endif

case rchange of
    '0': begin ; do not rotate original - back to 0 degrees rotation
         tmp_rot_angle = (state.rot_angle - 0.)

         if ptr_valid(state.head_ptr) then $
           hrot, main_image, *(state.head_ptr), $
             main_image, *(state.head_ptr), tmp_rot_angle, -1, -1, 1, mis=0 $
         else $
           main_image = rot(main_image,tmp_rot_angle,/interp)

         state.rot_angle = 0.
    end
    '90': begin ; rotate original 90 degrees
          tmp_rot_angle = (state.rot_angle - 90.)

          if ptr_valid(state.head_ptr) then $
            hrot, main_image, *(state.head_ptr), $
              main_image, *(state.head_ptr), tmp_rot_angle, -1, -1, 1, mis=0 $
          else $
            main_image = rot(main_image,tmp_rot_angle,/interp)

          state.rot_angle = 90.
    end
    '180': begin ; rotate original 180 degrees
           tmp_rot_angle = (state.rot_angle - 180.)

           if ptr_valid(state.head_ptr) then $
             hrot, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), tmp_rot_angle, -1, -1, 1, mis=0 $
           else $
             main_image = rot(main_image,tmp_rot_angle,/interp)

           state.rot_angle = 180.
    end
    '270': begin ; rotate original image 270 degrees
           tmp_rot_angle = (state.rot_angle - 270.)

           if ptr_valid(state.head_ptr) then $
             hrot, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), tmp_rot_angle, -1, -1, 1, mis=0 $
           else $
             main_image = rot(main_image,tmp_rot_angle,/interp)

           state.rot_angle = 270.
    end
    else:  begin
           tmp_rot_angle = (state.rot_angle - rchange)

           if ptr_valid(state.head_ptr) then $
             hrot, main_image, *(state.head_ptr), $
               main_image, *(state.head_ptr), tmp_rot_angle, -1, -1, 1, mis=0 $
           else $
             main_image = rot(main_image, tmp_rot_angle,/interp)

           state.rot_angle = rchange
    end

endcase

;Update header information after rotation if header is present
if ptr_valid(state.head_ptr) then begin
  head = *(state.head_ptr)
  atv_setheader, head
endif

atv_getstats, /align, /noerase

;Redisplay image with current zoom, update pan, and refresh image
atv_displayall

;make sure that the image arrays are updated for line/column plots, etc.
atv_resetwindow

end
;################################################################################
;------------------------------------------------------------------

pro atv_getimage

common atv_state
common atv_images

; Retrieve DSS, 2MASS, or IRAS image from STSCI/ESO/IRSA archives and 
; load into ATV.

formdesc = ['0, text, , label_left=Object Name: , width=15', $
            '0, label, OR, CENTER', $
            '0, text, , label_left=RA (Deg J2000): , width=15', $
            '0, text, , label_left=DEC (Deg J2000): , width=15', $
            '0, float, 10.0, label_left=Imsize (Arcminutes): ', $
            '0, droplist, DSS-STSCI|DSS-ESO|2MASS-IRSA|IRAS-IRSA, label_left=Archive:, set_value=0 ', $
            '0, droplist, 1st Generation|2nd Generation Blue|2nd Generation Red|2nd Generation Near-IR|J|H|K_s|12um|25um|60um|100um, label_left=Band:, set_value=0 ', $
            '0, button, SIMBAD|NED, set_value=0, exclusive', $
            '1, base, , row', $
            '0, button, Cancel, quit', $
            '0, button, Retrieve, quit']    

archiveform = cw_form(formdesc, /column, title = 'Get Archive Image')

if (archiveform.tag9 eq 1) then return

if (archiveform.tag10 eq 1) then begin

; First do error checking so that archive and band match

  if (strcompress(archiveform.tag0,/remove_all) eq '' AND $
      strcompress(archiveform.tag2,/remove_all) eq '' AND $
      strcompress(archiveform.tag3,/remove_all) eq '') then begin
      atv_message,'Enter Target or Coordinates', msgtype='error', /window
    return
  endif

  if (archiveform.tag5 eq 0 OR $
      archiveform.tag5 eq 1 AND $
      archiveform.tag6 ne 0 AND $
      archiveform.tag6 ne 1 AND $
      archiveform.tag6 ne 2 AND $
      archiveform.tag6 ne 3) then begin
      atv_message,'Select Appropriate Band for DSS', msgtype='error',/window
    return
  endif

  if (archiveform.tag5 eq 2 AND $
      archiveform.tag6 ne 4 AND $
      archiveform.tag6 ne 5 AND $
      archiveform.tag6 ne 6)then begin
      atv_message,'Select Appropriate Band for 2MASS', msgtype='error',/window
    return
  endif

  if (archiveform.tag5 eq 3 AND $
      archiveform.tag6 ne 7 AND $
      archiveform.tag6 ne 8 AND $
      archiveform.tag6 ne 9 AND $
      archiveform.tag6 ne 10) then begin
      atv_message,'Select Appropriate Band for IRAS', msgtype='error',/window
    return
  endif

  if (archiveform.tag4 lt 0.0) then begin
    atv_message, 'Image Size must be > 0', msgtype='error', /window
    return
  endif
 
; Set image size defaults.  For IRAS ISSA images, imsize must be 0.5,
; 1.0, 2.5, 5.0, or 12.5

  if (strcompress(archiveform.tag4, /remove_all) ne '') then $
    imsize = float(strcompress(archiveform.tag4, /remove_all)) $
  else $
    imsize = 10.0

  if (archiveform.tag5 eq 3) then begin
    if (strcompress(archiveform.tag4, /remove_all) ne '') then begin
      imsize = float(strcompress(archiveform.tag4, /remove_all)) 
      imsize = imsize / 60.
      diff_halfdeg = abs(0.5 - imsize)
      diff_deg = abs(1.0 - imsize)
      diff_2halfdeg = abs(2.5 - imsize)
      diff_5deg = abs(5.0 - imsize)
      diff_12halfdeg = abs(12.5 - imsize)

      diff_arr = [diff_halfdeg, diff_deg, diff_2halfdeg, diff_5deg, $
                  diff_12halfdeg]

      imsize_iras = [0.5, 1.0, 2.5, 5.0, 12.5]
      index_min = where(diff_arr eq min(diff_arr))
      imsize = imsize_iras[index_min]
    endif else begin
      imsize = 1.0
    endelse
  endif

  if (archiveform.tag5 eq 0 OR archiveform.tag5 eq 1) then begin
    if (archiveform.tag4 gt 60.0) then begin
      atv_message, 'DSS Image Size must be <= 60.0 Arcminutes', $
        msgtype='error', /window
      return
    endif
  endif

  widget_control, /hourglass
  image_str = ''

  if (strcompress(archiveform.tag0, /remove_all) ne '') then $
    image_str=strcompress(archiveform.tag0, /remove_all)

  if (strcompress(archiveform.tag2, /remove_all) ne '') then $
    ra_tmp=double(strcompress(archiveform.tag2, /remove_all))

  if (strcompress(archiveform.tag3, /remove_all) ne '') then $
    dec_tmp=double(strcompress(archiveform.tag3, /remove_all))

  if (strcompress(archiveform.tag0, /remove_all) ne '') then $
    target=image_str $
  else $
    target=[ra_tmp,dec_tmp]

  case archiveform.tag6 of 

  0: band='1'
  1: band='2b'
  2: band='2r'
  3: band='2i'
  4: band='j'
  5: band='h'
  6: band='k'
  7: band='12'
  8: band='25'
  9: band='60'
 10: band='100'
else:
  endcase

  case archiveform.tag5 of 

  0: begin
    if (archiveform.tag7 eq 0) then $
      querydss, target, tmpim, hdr, imsize=imsize, survey=band, /stsci $
    else $
      querydss, target, tmpim, hdr, imsize=imsize, survey=band, /stsci, /ned
  end

  1: begin
    if (archiveform.tag7 eq 0) then $
      querydss, target, tmpim, hdr, imsize=imsize, survey=band, /eso $
    else $
      querydss, target, tmpim, hdr, imsize=imsize, survey=band, /eso, /ned
  end

  2: begin
    if (archiveform.tag7 eq 0) then $
      query2mass, target, tmpim, hdr, imsize=imsize, band=band $
    else $
      query2mass, target, tmpim, hdr, imsize=imsize, band=band, /ned
  end

  3: begin
    if (archiveform.tag7 eq 0) then $
      queryiras, target, tmpim, hdr, imsize=imsize, band=band $
    else $
      queryiras, target, tmpim, hdr, imsize=imsize, band=band, /ned
  end
     else:
  endcase

  atv,tmpim,head=hdr
endif


;Reset image rotation angle to 0 and inversion to none
state.rot_angle = 0.
state.invert_image = 'none'

;Make pan image and set image to current zoom/stretch levels
atv_makepan
atv_refresh

;make sure that the image arrays are updated for line/column plots, etc.
atv_resetwindow

end

;------------------------------------------------------------------


pro atv_pixtable

; Create a table widget that will show a 5x5 array of pixel values
; around the current cursor position

if (not(xregistered('atv_pixtable', /noshow))) then begin

  common atv_state
  common atv_images

  state.pixtable_base_id = $
    widget_base(/base_align_right, $
                 group_leader = state.base_id, $
                 /column, $
                 title = 'atv pixel table')

  state.pixtable_tbl_id = widget_table(state.pixtable_base_id,   $
                   value=[0,0], xsize=5, ysize=5, row_labels='', $ 
                   column_labels='', alignment=2, /resizeable_columns)

  pixtable_done = widget_button(state.pixtable_base_id, $
                                value = 'Done', $
                                uvalue = 'pixtable_done')

  widget_control, state.pixtable_base_id, /realize
  xmanager, 'atv_pixtable', state.pixtable_base_id, /no_block

endif

end

;---------------------------------------------------------------------

pro atv_pixtable_event, event

common atv_state

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'pixtable_done': widget_control, event.top, /destroy
    else:
endcase

end

;--------------------------------------------------------------------

pro atv_pixtable_update

  common atv_state
  common atv_images

  zcenter = (0 > state.coord < state.image_size[0:1])

;Check and adjust the zcenter if the cursor is near the edges of the image

  if (zcenter[0] le 2) then zcenter[0] = 2
  if (zcenter[0] gt (state.image_size[0]-3)) then $
      zcenter[0] =  state.image_size[0] - 3

  if (zcenter[1] le 2) then zcenter[1] = 2
  if (zcenter[1] gt (state.image_size[1]-3)) then $
      zcenter[1] = state.image_size[1] - 3 

  pix_values = dblarr(5,5)
  row_labels = strarr(5)
  column_labels = strarr(5)
  boxsize=2

  xmin = 0 > (zcenter[0] - boxsize)
  xmax = (zcenter[0] + boxsize) < (state.image_size[0] - 1) 
  ymin = 0 > (zcenter[1] - boxsize) 
  ymax = (zcenter[1] + boxsize) < (state.image_size[1] - 1)

  row_labels = [strcompress(string(ymax),/remove_all),   $
                strcompress(string(ymin+3),/remove_all), $
                strcompress(string(ymin+2),/remove_all), $
                strcompress(string(ymin+1),/remove_all), $
                strcompress(string(ymin),/remove_all)]

  column_labels = [strcompress(string(xmin),/remove_all),   $
                   strcompress(string(xmin+1),/remove_all), $
                   strcompress(string(xmin+2),/remove_all), $
                   strcompress(string(xmin+3),/remove_all), $
                   strcompress(string(xmax),/remove_all)]

  pix_values = main_image[xmin:xmax, ymin:ymax]
  pix_values = reverse(pix_values, 2, /overwrite)

  widget_control, state.pixtable_tbl_id, set_value = pix_values, $
          column_labels=column_labels, row_labels=row_labels

end

;--------------------------------------------------------------------



pro atv_changeimage,imagenum,next=next,previous=previous

	common atv_state
	common atv_images

	; do nothing if we only have one image. 
	if state.image_size[2] lt 2 then return 
		
	; if we've got a 3d image stack this lets us move between
	; images of the stack. Check to make sure we don't exit the boundaries
	
	if keyword_set(next) then imagenum =  (state.cur_image_num +1) < (state.image_size[2]-1)
	if keyword_set(previous) then imagenum =  (state.cur_image_num -1) > 0

	; check imagenum is withing bounds
	if (imagenum lt 0) or (imagenum gt (state.image_size[2]-1)) then begin

        state.cur_image_num = state.cur_image_num
        widget_control, state.curimnum_text_id, $
                                set_value = string(state.cur_image_num)
        text_warn = 'Please enter a value between 0 and ' + $
                     strcompress(string(state.image_size[2] -1))
        atv_message, text_warn, msgtype='error', /window
		return
	endif


    widget_control, state.curimnum_text_id, $
                            set_value = string(imagenum)
    widget_control, state.curimnum_slidebar_id, $
                            set_value = imagenum

	; don't do anything if the image number doesn't change.
	if (imagenum eq state.cur_image_num) then return
    state.cur_image_num = imagenum
	if n_elements(image_names) gt 1 then begin
		state.title_extras = image_names[state.cur_image_num]
	endif else begin
       	state.title_extras = $
       	  strcompress('Plane ' + string(state.cur_image_num))
	endelse
   
	atv_settitle 
	
    main_image = main_image_stack[*, *, state.cur_image_num]
    atv_getstats,/align,/noerase
    case state.curimnum_minmaxmode of
		'Zero/Max': begin
           state.min_value = 0
           state.max_value = state.image_max
		end
         'Min/Max': begin
           state.min_value = state.image_min
           state.max_value = state.image_max
         end
         'AutoScale': atv_autoscale
         'Constant': donothingvariable = 0
         else: print, 'Unknown Min/Max mode for changing cur_image_num!'
    endcase
    atv_set_minmax
    atv_displayall

end 


;-----------------------------------------------------------------------

pro atv_autoscale

; Routine to auto-scale the image.  

common atv_state 
common atv_images

widget_control, /hourglass

if (n_elements(main_image) LT 5.e5) then begin
    med = median(main_image)
    sig = stddev(main_image,/NaN)
endif else begin   ; resample big images before taking median, to save memory
    boxsize = 10
    rx = state.image_size[0] mod boxsize
    ry = state.image_size[1] mod boxsize
    nx = state.image_size[0] - rx
    ny = state.image_size[1] - ry
    tmp_img = rebin(main_image[0: nx-1, 0: ny-1], $
                    nx/boxsize, ny/boxsize, /sample)
    med = median(tmp_img)
    sig = stddev(temporary(tmp_img),/NaN)
endelse

nhigh = 10
nlow = 2

state.max_value = (med + (nhigh * sig)) < state.image_max
state.min_value = (med - (nlow * sig))  > state.image_min

if (finite(state.min_value) EQ 0) then state.min_value = state.image_min
if (finite(state.max_value) EQ 0) then state.max_value = state.image_max

if (state.min_value GE state.max_value) then begin
    state.min_value = state.min_value - 1
    state.max_value = state.max_value + 1
endif

atv_set_minmax

end  

;--------------------------------------------------------------------
pro atv_autozoom
; Routine to auto-zoomthe image to fit in the current draw window.

common atv_state 
common atv_images

    if (state.image_size[0]*2 le state.draw_window_size[0]) AND $
           (state.image_size[1]*2 le state.draw_window_size[1]) then $
           for tempvar = 1, $
           floor(sqrt(min(state.draw_window_size/state.image_size[0:1],/NaN))) do $
           atv_zoom, 'in'
    if (state.image_size[0] ge state.draw_window_size[0]*2) AND $
           (state.image_size[1] ge state.draw_window_size[1]*2) then $
           for tempvar = 1, $
           floor(sqrt(min(state.image_size[0:1]/state.draw_window_size))) do $
           atv_zoom, 'out'
end

;--------------------------------------------------------------------

pro atv_restretch

; Routine to restretch the min and max to preserve the display
; visually but use the full color map linearly.  Written by DF, and
; tweaked and debugged by AJB.  It doesn't always work exactly the way
; you expect (especially in log-scaling mode), but mostly it works fine.

common atv_state

sx = state.brightness
sy = state.contrast

if (state.scaling EQ 2) then return ; do nothing for hist-eq mode

if (state.scaling EQ 0) then begin
    sfac = (state.max_value-state.min_value)
    state.max_value = sfac*(sx+sy)+state.min_value
    state.min_value = sfac*(sx-sy)+state.min_value
endif

if (state.scaling EQ 1) then begin

    offset = state.min_value - $
      (state.max_value - state.min_value) * 0.01

    sfac = alog10((state.max_value - offset) / (state.min_value - offset))
    state.max_value = 10.^(sfac*(sx+sy)+alog10(state.min_value - offset)) $
      + offset
    state.min_value = 10.^(sfac*(sx-sy)+alog10(state.min_value - offset)) $
      + offset
    
endif


; Try different behavior in asinh mode: usually want to keep the min
; value the same and just adjust the max value.  Seems to work ok.
if (state.scaling EQ 3) then begin
    sfac = asinh(state.max_value / state.asinh_beta) - $
      asinh(state.min_value / state.asinh_beta)

    state.max_value = sinh(sfac*(sx+sy) + $
          asinh(state.min_value/state.asinh_beta))*state.asinh_beta 
endif

; do this differently for 8 or 24 bit color, to prevent flashing
atv_setwindow, state.draw_window_id
if (state.bitdepth EQ 8) then begin
    atv_set_minmax
    atv_displayall
    state.brightness = 0.5      ; reset these
    state.contrast = 0.5
    atv_stretchct
endif else begin
    state.brightness = 0.5      ; reset these
    state.contrast = 0.5
    atv_stretchct
    atv_set_minmax
    atv_displayall
endelse
atv_resetwindow

end

;---------------------------------------------------------------------

function atv_wcsstring, lon, lat, ctype, equinox, disp_type, disp_equinox, $
            disp_base60

; Routine to return a string which displays cursor coordinates.
; Allows choice of various coordinate systems.
; Contributed by D. Finkbeiner, April 2000.
; 29 Sep 2000 - added degree (RA,dec) option DPF
; Apr 2007: AJB added additional error checking to prevent crashes

; ctype - coord system in header
; disp_type - type of coords to display

headtype = strmid(ctype[0], 0, 4)

; need numerical equinox values
IF (equinox EQ 'J2000') THEN num_equinox = 2000.0 ELSE $
  IF (equinox EQ 'B1950') THEN num_equinox = 1950.0 ELSE $
  num_equinox = float(equinox)

IF (disp_equinox EQ 'J2000') THEN num_disp_equinox = 2000.0 ELSE $
  IF (disp_equinox EQ 'B1950') THEN num_disp_equinox = 1950.0 ELSE $
  num_disp_equinox = float(equinox)

; first convert lon,lat to RA,dec (J2000)
CASE headtype OF 
    'GLON': euler, lon, lat, ra, dec, 2 ; J2000
    'ELON': BEGIN 
        euler, lon, lat, ra, dec, 4 ; J2000
        IF num_equinox NE 2000.0 THEN precess, ra, dec, num_equinox, 2000.0
    END 
    'RA--': BEGIN    
        ra = lon
        dec = lat
        IF num_equinox NE 2000.0 THEN precess, ra, dec, num_equinox, 2000.0
    END 
    'DEC-': BEGIN       ; for SDSS images
        ra = lon
        dec = lat
        IF num_equinox NE 2000.0 THEN precess, ra, dec, num_equinox, 2000.0
    END
    else: begin
        wcsstring = '---No WCS Info---'
		if keyword_set(state) then begin
			widget_control, state.wcs_bar_id, set_value = wcsstring
			state.wcstype = 'none'
		endif
        return, wcsstring
    end
ENDCASE  

; Now convert RA,dec (J2000) to desired display coordinates:  

IF (disp_type[0] EQ 'RA--' or disp_type[0] EQ 'DEC-') THEN BEGIN 
; generate (RA,dec) string 
   disp_ra  = ra
   disp_dec = dec
   IF num_disp_equinox NE 2000.0 THEN precess, disp_ra, disp_dec, $
     2000.0, num_disp_equinox

   IF disp_base60 THEN BEGIN ; (hh:mm:ss) format
      
      neg_dec  = disp_dec LT 0
      radec, disp_ra, abs(disp_dec), ihr, imin, xsec, ideg, imn, xsc
      wcsstring = string(ihr, imin, xsec, ideg, imn, xsc, disp_equinox, $
         format = '(i2.2,":",i2.2,":",f6.3,"   ",i2.2,":",i2.2,":",f5.2," ",a6)' )
      if (strmid(wcsstring, 6, 1) EQ ' ') then $
        strput, wcsstring, '0', 6
      if (strmid(wcsstring, 21, 1) EQ ' ') then $
        strput, wcsstring, '0', 21
      IF neg_dec THEN strput, wcsstring, '-', 14

   ENDIF ELSE BEGIN ; decimal degree format

      wcsstring = string(disp_ra, disp_dec, disp_equinox, $
                         format='("Deg ",F9.5,",",F9.5,a6)')
   ENDELSE 
ENDIF 
     

IF disp_type[0] EQ 'GLON' THEN BEGIN ; generate (l,b) string
    euler, ra, dec, l, b, 1
    
    wcsstring = string(l, b, format='("Galactic (",F9.5,",",F9.5,")")')
ENDIF 

IF disp_type[0] EQ 'ELON' THEN BEGIN ; generate (l,b) string
    
    disp_ra = ra
    disp_dec = dec
    IF num_disp_equinox NE 2000.0 THEN precess, disp_ra, disp_dec, $
      2000.0, num_disp_equinox
    euler, disp_ra, disp_dec, lam, bet, 3
    
    wcsstring = string(lam, bet, format='("Ecliptic (",F9.5,",",F9.5,")")')
ENDIF 

return, wcsstring
END

;----------------------------------------------------------------------

function atv_wavestring

; function to return string with wavelength info for spectral images

common atv_state

cd = (*state.astr_ptr).cd[0,0]
crpix = (*state.astr_ptr).crpix[0]
crval = (*state.astr_ptr).crval[0]

cunit = sxpar(*state.head_ptr, 'cunit1')
cunit = strcompress(string(cunit), /remove_all)
if (cunit NE '0') then begin
    cunit = strcompress(strupcase(strmid(cunit,0,1)) + strmid(cunit,1), $
                        /remove_all)
endif else begin
    cunit = ''
endelse

shifta = float(sxpar(*state.head_ptr, 'SHIFTA1'))

wavelength = crval + ((state.coord[0] - crpix) * cd) + (shifta * cd)
wstring = string(wavelength, format='(F8.2)')

wavestring = strcompress('Wavelength:  ' + wstring + ' ' + state.cunit)

return, wavestring

end

;--------------------------------------------------------------------


pro atv_gettrack

; Create the image to display in the track window that tracks
; cursor movements.  Also update the coordinate display and the
; (x,y) and pixel value.

common atv_state
common atv_images

; Get x and y for center of track window

zcenter = (0 > state.coord < state.image_size[0:1])

track = bytarr(11,11)
boxsize=5
xmin = 0 > (zcenter[0] - boxsize)
xmax = (zcenter[0] + boxsize) < (state.image_size[0] - 1) 
ymin = 0 > (zcenter[1] - boxsize) 
ymax = (zcenter[1] + boxsize) < (state.image_size[1] - 1)

startx = abs( (zcenter[0] - boxsize) < 0 )
starty = abs( (zcenter[1] - boxsize) < 0 ) 

track[startx,starty] = scaled_image[xmin:xmax,ymin:ymax]
track_image = rebin(track, $
                    state.track_window_size, state.track_window_size, $
                    /sample)

atv_setwindow, state.track_window_id
tv, track_image

; Overplot an X on the central pixel in the track window, to show the
; current mouse position

; Changed central x to be green always
plots, [0.46, 0.54], [0.46, 0.54], /normal, color = state.box_color, psym=0
plots, [0.46, 0.54], [0.54, 0.46], /normal, color = state.box_color, psym=0

; update location bar with x, y, and pixel value

loc_string = $
  string(state.coord[0], $
         state.coord[1], $
         main_image[state.coord[0], $
                    state.coord[1]], $
         format = '("(",i5,",",i5,") ",g12.6)') 

widget_control, state.location_bar_id, set_value = loc_string

; Update coordinate display

if (state.wcstype EQ 'angle') then begin
    xy2ad, state.coord[0], state.coord[1], *(state.astr_ptr), lon, lat

    wcsstring = atv_wcsstring(lon, lat, (*state.astr_ptr).ctype,  $
                              state.equinox, state.display_coord_sys, $
                              state.display_equinox, state.display_base60)

    widget_control, state.wcs_bar_id, set_value = wcsstring

endif    

if (state.wcstype EQ 'lambda') then begin
    wavestring = atv_wavestring()
    widget_control, state.wcs_bar_id, set_value = wavestring
endif

atv_resetwindow

end

;----------------------------------------------------------------------

pro atv_drawbox, norefresh=norefresh

; routine to draw the box on the pan window, given the current center
; of the display image.

common atv_state
common atv_images

atv_setwindow, state.pan_window_id

view_min = round(state.centerpix - $
        (0.5 * state.draw_window_size / state.zoom_factor)) 
view_max = round(view_min + state.draw_window_size / state.zoom_factor) - 1

; Create the vectors which contain the box coordinates

box_x = float((([view_min[0], $
                 view_max[0], $
                 view_max[0], $
                 view_min[0], $
                 view_min[0]]) * state.pan_scale) + state.pan_offset[0]) 

box_y = float((([view_min[1], $
                 view_min[1], $
                 view_max[1], $
                 view_max[1], $
                 view_min[1]]) * state.pan_scale) + state.pan_offset[1]) 

; set limits on box to make sure all sides always appear
w = where(box_x LT 0, count)
if (count GT 0) then box_x[w] = 0
w = where(box_y LT 0, count)
if (count GT 0) then box_y[w] = 0
w = where(box_x GT state.pan_window_size-1, count)
if (count GT 0) then box_x[w] = state.pan_window_size-1
w = where(box_y GT state.pan_window_size-1, count)
if (count GT 0) then box_y[w] = state.pan_window_size-1

; Redraw the pan image and overplot the box
if (not(keyword_set(norefresh))) then $
    device, copy=[0,0,state.pan_window_size, state.pan_window_size, 0, 0, $
                  state.pan_pixmap]

plots, box_x, box_y, /device, color = state.box_color, psym=0

atv_resetwindow

end

;----------------------------------------------------------------------

pro atv_pantrack, event

; routine to track the view box in the pan window during cursor motion

common atv_state

; get the new box coords and draw the new box

tmp_event = [event.x, event.y] 

newpos = state.pan_offset > tmp_event < $
  (state.pan_offset + (state.image_size[0:1] * state.pan_scale))

state.centerpix = round( (newpos - state.pan_offset ) / state.pan_scale)

atv_drawbox
atv_getoffset

end

;----------------------------------------------------------------------

pro atv_resize

; Routine to resize the draw window when a top-level resize event
; occurs.  

common atv_state

widget_control, state.base_id, tlb_get_size=tmp_event

window = (state.base_min_size > tmp_event)

newbase = window - state.base_pad

newxsize = (tmp_event[0] - state.base_pad[0]) > $
  (state.base_min_size[0] - state.base_pad[0]) 
newysize = (tmp_event[1] - state.base_pad[1]) > $
  (state.base_min_size[1] - state.base_pad[1])

widget_control, state.draw_widget_id, $
  scr_xsize = newxsize, scr_ysize = newysize
widget_control, state.colorbar_widget_id, $
  scr_xsize = newxsize, scr_ysize = state.colorbar_height

state.draw_window_size = [newxsize, newysize]

atv_colorbar

widget_control, state.base_id, /clear_events
widget_control, state.draw_base_id, /sensitive, /input_focus

end

;----------------------------------------------------------------------

PRO atv_setscaling,scaling
	common atv_state

	state.scaling = scaling

	; for asinh, disply the nonlinearity parameter
	; otherwise, don't
	if state.scaling eq 3 then begin
    	widget_control, state.nonlin_base_id,map=1,xsize=200,ysize=35
	endif else begin
    	widget_control, state.nonlin_base_id,map=0,xsize=1,ysize=1
	endelse

	atv_displayall

end
;----------------------------------------------------------------------

pro atv_scaleimage

; Create a byte-scaled copy of the image, scaled according to
; the state.scaling parameter.  Add a padding of 5 pixels around the
; image boundary, so that the tracking window can always remain
; centered on an image pixel even if that pixel is at the edge of the
; image.    
;
; We add 8 to the value returned from bytscl to get above the 8 primary
; colors which are used for overplots and annotations. We use a mask to
; only do this for non NAN pixels, so all NANs will always remain black, 
; no matter what color the bottom of the color map is.

common atv_state
common atv_images

; Since this can take some time for a big image, set the cursor 
; to an hourglass until control returns to the event loop.

widget_control, /hourglass

scaled_image=0

nan_mask = main_image eq main_image ; mask out NAN pixels
case state.scaling of
    0: scaled_image = $                 ; linear stretch
      bytscl(main_image, $
             /nan, $
             min=state.min_value, $
             max=state.max_value, $
             top = state.ncolors - 1) + 8*nan_mask
    
    1: begin                            ; log stretch
        offset = state.min_value - $
          (state.max_value - state.min_value) * 0.01

        scaled_image = $        
          bytscl( alog10(main_image - offset), $
                  min=alog10(state.min_value - offset), /nan, $
                  max=alog10(state.max_value - offset),  $
                  top=state.ncolors - 1) + 8*nan_mask
    end
    

    2: scaled_image = $                 ; histogram equalization
      bytscl(hist_equal(main_image, $
                        minv = state.min_value, $    
                        maxv = state.max_value), $
             /nan, top = state.ncolors - 1) + 8*nan_mask
	3: begin				; asinh stretch. requires Dave Fanning's ASINHSCL.PRO
		scaled_image = asinhscl( main_image, $
				min = state.min_value, $
				max = state.max_value,$
				omax = state.ncolors - 1) + 8*nan_mask
			;	nonlinearity = state.nonlinearity_value) + 8*nan_mask

	end
    4: begin                            ; square root stretch
        scaled_image = $        
          bytscl( sqrt(main_image), $
                  min=sqrt(state.min_value), /nan, $
                  max=sqrt(state.max_value),  $
                  top=state.ncolors - 1) + 8*nan_mask
    end
	
    
endcase

end

;----------------------------------------------------------------------

pro atv_setasinh

; get the asinh beta parameter

common atv_state

b = string(state.asinh_beta)

formline = strcompress('0,float,' + b + $
                       ',label_left=Asinh beta parameter: ,width=10')

formdesc = [formline, $
           '0, button, Set beta, quit', $
           '0, button, Cancel, quit']

textform = cw_form(formdesc, ids=ids, /column, $
                 title = 'atv asinh stretch settings')

if (textform.tag2 EQ 1) then return

state.asinh_beta = float(textform.tag0)

atv_displayall

end

;----------------------------------------------------------------------

pro atv_getstats, align=align, noerase=noerase

; Get basic image stats: min and max, and size.
; set align keyword to preserve alignment of previous image

common atv_state
common atv_images

; this routine operates on main_image, which is in the
; atv_images common block

widget_control, /hourglass

oldimagesize = state.image_size

state.image_size[0:1] = (size(main_image))[1:2]

if ((oldimagesize[0] NE state.image_size[0]) OR $
    (oldimagesize[1] NE state.image_size[1])) then align = 0

state.image_min = min(main_image, max=maxx, /nan)
state.image_max = maxx

; Get sky value for autoscaling and asinh stretch.  Eliminate
; zero-valued and NaN pixels from sky calculation, i.e. for HST ACS
; drizzled images, WFPC2 mosaics, or Spitzer images.
w = where(finite(main_image) AND (main_image NE 0.0), goodcount)
if (goodcount GT 25) then begin
    sky, main_image[w], skymode, skysig, /silent
endif else if (goodcount GT 5 AND goodcount LE 25) then begin
    skysig = stddev(main_image[w])
    skymode = median(main_image[w])
endif else if (goodcount LE 5) then begin ; really pathological images
    skysig = 1.
    skymode = 0.
endif

; error checking- in case sky.pro returns a zero or negative sigma
if (skysig LE 0.0) then skysig = stddev(main_image)
if (skysig LE 0.0) then skysig = 1.0

state.skymode = skymode
state.skysig = skysig
state.asinh_beta = state.skysig

if (state.min_value GE state.max_value) then begin
    state.min_value = state.min_value - 1
    state.max_value = state.max_value + 1
endif

; zero the current display position on the center of the image,
; unless user selected /align keyword

state.coord = round(state.image_size[0:1] / 2.)
IF (NOT keyword_set(align) OR (state.firstimage EQ 1)) THEN $
  state.centerpix = round(state.image_size[0:1] / 2.)
atv_getoffset

; Clear all plot annotations
if (not(keyword_set(noerase))) then atverase, /norefresh  

end

;-------------------------------------------------------------------

pro atv_setwindow, windowid,nostretch=nostretch

; replacement for wset.  Reads the current active window first.
; This should be used when the currently active window is an external
; (i.e. non-atv) idl window.  Use atv_setwindow to set the window to
; one of the atv window, then display something to that window, then
; use atv_resetwindow to set the current window back to the currently
; active external window.  Make sure that device is not set to
; postscript, because if it is we can't display anything.
;
; atv_setwindow will now automatically re-stretch the atv color table
; in case the user has changed the color table elsewhere. Set 
; /nostretch to disable this behavior. (This can be useful if you're
; about to call atv_stretchct anyway with different brightness & contrast,
; as it prevents uselessly calling the function twice in a row.)
;
;
; 2004-05-05 This also now stores the user's !P.MULTI setting.
; 2005-12-09 And also the user's device,decomposed setting for 24-bit displays.

common atv_state
common atv_color


state.active_window_pmulti = !p.multi
!p.multi = 0

tvlct, user_r, user_g, user_b, /get

; regenerate atv color table
atv_initcolors
if state.bitdepth gt 8 then begin
	device,get_decomposed=decomp
	state.user_decomposed=decomp
	device,decomposed=0
endif
if ~(keyword_set(nostretch)) then atv_stretchct

if (!d.name NE 'PS') then begin
    state.active_window_id = !d.window
    wset, windowid
endif

; use for debugging
; print, 'atv_setwindow', state.active_window_id


end

;---------------------------------------------------------------------

pro atv_resetwindow

; reset to current active window

common atv_state
common atv_color


; The empty command used below is put there to make sure that all
; graphics to the previous atv window actually get displayed to screen
; before we wset to a different window.  Without it, some line
; graphics would not actually appear on screen.
; Also reset to user's external color map and p.multi.

; use for debugging
; print, 'atv_resetwindow', state.active_window_id

if (!d.name NE 'PS') then begin
    empty
    wset, state.active_window_id
    tvlct, user_r, user_g, user_b
	if state.bitdepth gt 8 then device,decomposed=state.user_decomposed
endif

!p.multi = state.active_window_pmulti


end

;------------------------------------------------------------------

pro atv_getwindow

; get currently active window id

common atv_state

if (!d.name NE 'PS') then begin
    state.active_window_id = !d.window
endif
if state.bitdepth gt 8 then tvlct,/GET,user_r,user_g,user_b
end

;--------------------------------------------------------------------
;    Fits file reading routines
;--------------------------------------------------------------------

; read in multiple fits files. This calls atv readfits some number
; of times as appropriate.
;
; atv_readfits saves things into main_image; 
; atv_readmultifits will concatenate that into main_image_stack as 
;    needed.
pro atv_readmultifits, fitsfilename=fitsfilename, newimage=newimage

common atv_state
common atv_images

newimage = 0
cancelled = 0
; Read in a new image when user goes to the File->ReadFits menu.
; Do a reasonable amount of error-checking first, to prevent unwanted
; crashes. 
if (n_elements(fitsfilename) EQ 0) then window = 1 else window = 0


; Added list of filters for SDSS suffixes, etc.
filterlist = ['*.fits;*.fit;*.fit.gz;*.fits.gz']

; If fitsfilename hasn't been passed to this routine, get filename
; from dialog_pickfile.
if (n_elements(fitsfilename) EQ 0) then begin
    fitsfile = $
      dialog_pickfile(filter = filterlist, $
                      group = state.base_id, $
                      /must_exist, $
					  /multiple_files,$
                      /read, $
                      path = state.current_dir, $
                      get_path = tmp_dir, $
                      title = 'Select Fits Image')        
    if (tmp_dir NE '') then state.current_dir = tmp_dir
	if (n_elements(fitsfile) eq 1) then $
	    if (fitsfile EQ '') then return ; 'cancel' button returns empty string
endif else begin
    fitsfile = fitsfilename
endelse
image_names=fitsfile

if (n_elements(fitsfile) gt 1) then begin
	; read in multiple images into a data cube
	filenames=fitsfile
	message,/info,"Loading multiple fits files - "
	message,/info," they had better all be the same dimensions!"
	; read in first file
	print,filenames[0]
	atv_readfits, fitsfilename=filenames[0], newimage=newimage,head=head
	sz_first = size(main_image)
	if sz_first[0] eq 3 then nz = sz_first[3] else nz=1 ; images PER FILE
	; allocate memory. All images must be the same size.
	main_image_stack = fltarr(sz_first[1],sz_first[2],nz*n_elements(filenames))
	image_names = strarr(nz*n_elements(filenames))
		image_names[0:nz*1-1] = filenames[0]
	main_image_stack[*,*,0] = main_image
	header_array = ptrarr(n_elements(filenames))
	for i=1,n_elements(filenames)-1 do begin

		print,i,"  ",filenames[i]
		atv_readfits, fitsfilename=filenames[i], newimage=newimage

		sz = size(main_image)
	    if sz[0] eq 3 then nzi = sz[3] else nzi=1 ; images PER FILE
		if sz[1] ne sz_first[1] or sz[2] ne sz_first[2] or nzi ne nz then begin
		    atv_message, 'When loading multiple FITS files, they all must have exactly the same shape. File '+filenames[i]+" does not match!", $
      		window = window, msgtype = 'error'
			main_image_stack=fltarr(500,500)
			main_image=fltarr(500,500)
			break
		endif
		main_image_stack[*,*,nz*i]=main_image
		; TODO make this handle multiple planes per image correctly.
		image_names[nz*i:nz*(i+1)-1] = filenames[i]
		
	endfor
	
    main_image = main_image_stack[*, *, 0]

	;if n_elements(image_names) gt 1 then begin
		state.imagename = ""
		state.title_extras = image_names[state.cur_image_num]
	;endif else begin
		;state.imagename = ""
     	;state.title_extras = $
		   	;strcompress('Plane ' + string(state.cur_image_num))
	;endelse

	atv_setheader, head
endif else if (n_elements(fitsfile) eq 1) then begin
	; read in a single image.
	fitsfile=fitsfile[0]
	atv_readfits, fitsfilename=fitsfile, newimage=newimage,head=head
	main_image_stack=main_image
	state.imagename = fitsfile
	atv_setheader, head
	;main_image_stack=main_image
		;stop
endif

state.cur_image_num = 0
if (size(main_image_stack))[0] eq 2 then begin
	; 2D image
	
    widget_control, state.curimnum_base_id,map=0,xsize=1,ysize=1
	state.image_size = [(size(main_image_stack))[1:2],1]
endif else begin 
	; 3D data cube. Either from 1 file, or from multiple files concatenated.
    state.image_size = (size(main_image_stack))[1:3]
    widget_control, state.curimnum_text_id, sensitive = 1, $
             set_value = 0
    widget_control, state.curimnum_slidebar_id, sensitive = 1, $
             set_value = 0, set_slider_min = 0, $
             set_slider_max = state.image_size[2]-1
       widget_control,state.curimnum_base_id,map=1, $
          xsize=state.draw_window_size[0],ysize=45

endelse

widget_control, /hourglass

newimage = 1

;Reset image rotation angle to 0 and inversion to none
state.rot_angle = 0.
state.invert_image = 'none'

;if ( (size(main_image_stack))[0] EQ 3 ) then begin ; case of 3-d stack of images [x,y,n]
;endif else begin
;    widget_control, state.curimnum_text_id, sensitive = 1, $
;             set_value = 0
;    widget_control, state.curimnum_slidebar_id, sensitive = 1, $
;             set_value = 0, set_slider_min = 0, $
;             set_slider_max = state.image_size[2]-1
;
;endelse
;



end


;--------------------------------------------------------------------

pro atv_readfits, fitsfilename=fitsfilename, newimage=newimage,head=head
; this reads in one fits file.
; The file is stored in the variable main_image, regardless of its size.
; It is up to readmultifits to determine what to do with it, in terms of
; concatenating files into a data cube, setting display parameters, etc.

common atv_state
common atv_images


fitsfile=fitsfilename
; Get fits header so we know what kind of image this is.
head = headfits(fitsfile)

; Check validity of fits file header 
if (n_elements(strcompress(head, /remove_all)) LT 2) then begin
    atv_message, 'File does not appear to be a valid FITS image!', $
      window = window, msgtype = 'error'
    return
endif
if (!ERR EQ -1) then begin
    atv_message, $
      'Selected file does not appear to be a valid FITS image!', $
      msgtype = 'error', window = window
    return
endif

; Find out if this is a fits extension file, and how many extensions
; New: use fits_open rather than fits_info
fits_open, fitsfile, fcb, message = message
if (message NE '') then begin
    atv_message, message, msgtype='error', /window
    return
end
numext = fcb.nextend
has_primary_image = (fcb.naxis[0] ge 2)
fits_close, fcb

instrume = strcompress(string(sxpar(head, 'INSTRUME')), /remove_all)
origin = strcompress(sxpar(head, 'ORIGIN'), /remove_all)
naxis = sxpar(head, 'NAXIS')

; check for OSIRIS or NIRC2, which don't use standard keywords
if instrume eq '0' then instrume = strcompress(string(sxpar(head, 'CURRINST')), /remove_all)

; Make sure it's not a 1-d spectrum
if (numext EQ 0 AND naxis LT 2) then begin
    atv_message, 'Selected file is not a 2-d or 3-d FITS image!', $
      window = window, msgtype = 'error'
    return
endif

state.title_extras = ''

; Now call the subroutine that knows how to read in this particular
; data format:
cancelled=0

if ((instrume EQ 'OSIRIS')  AND (naxis EQ 3)) then begin
    atv_osiris_read, fitsfile, head, cancelled
endif else if ((numext GT 0) AND (instrume NE 'WFPC2')) then begin
    atv_fitsext_read, fitsfile, numext, head, cancelled, has_primary_image=has_primary_image
endif else if ((instrume EQ 'WFPC2') AND (naxis EQ 3)) then begin
    atv_wfpc2_read, fitsfile, head, cancelled
endif else if ((naxis EQ 3) AND (origin EQ '2MASS')) then begin
    atv_2mass_read, fitsfile, head, cancelled
endif else begin
    atv_plainfits_read, fitsfile, head, cancelled
endelse

if (cancelled EQ 1) then return

; Make sure it's a 2-d or 3-d image
    if ( (size(main_image))[0] NE 2 AND (size(main_image))[0] NE 3) then begin
        atv_message, 'Selected file is not a 2-D or 3-D fits image!', $
          msgtype = 'error', window = window
          print, 'ERROR: Input data must be a 2-d or 3-d array!'    
		  if (size(main_image))[0] eq 4 then begin
			  print, 'Reforming a 4D array to 3D....'
			  sz = size(main_image)
			  main_image = reform(main_image, sz[1], sz[2], sz[3]*sz[4])
			  stop
		  endif else begin
	        main_image = fltarr(512, 512)
          	newimage = 1
          	return
		  endelse
	endif


end

;----------------------------------------------------------
;  Subroutines for reading specific data formats
;---------------------------------------------------------------

pro atv_fitsext_read, fitsfile, numext, head, cancelled, has_primary_image= has_primary_image

; Fits reader for fits extension files

common atv_state
common atv_images


if numext eq 1 and n_elements(has_primary_image) gt 0 and (has_primary_image eq 0) then begin
	; we are asked to open an image that only contains one image extension and
	; no primary HDU. No need to query, just open that one. 
	extension=1
	
endif else begin
	; it's not clear what extension to open, so ask the user
	; with a dialog box. 



	numlist = ''
	for i = 1, numext do begin
		numlist = strcompress(numlist + string(i) + '|', /remove_all)
	endfor

	numlist = strmid(numlist, 0, strlen(numlist)-1)

	droptext = strcompress('0, droplist, ' + numlist + $
						   ', label_left=Select Extension:, set_value=0')

	formdesc = ['0, button, Read Primary Image, quit', $
				'0, label, OR:', $
				droptext, $
				'0, button, Read Fits Extension, quit', $
				'0, button, Cancel, quit']

	textform = cw_form(formdesc, /column, $
					   title = 'Fits Extension Selector')

	if (textform.tag4 EQ 1) then begin  ; cancelled 
		cancelled = 1
		return                         
	endif

	if (textform.tag3 EQ 1) then begin   ;extension selected
		extension = long(textform.tag2) + 1
	endif else begin
		extension = 0               ; primary image selected
	endelse
endelse

; Make sure it's not a fits table: this would make mrdfits crash
head = headfits(fitsfile, exten=extension)
xten = strcompress(sxpar(head, 'XTENSION'), /remove_all)
if (xten EQ 'BINTABLE') then begin
    atv_message, 'File appears to be a FITS table, not an image.', $
      msgtype='error', /window
    cancelled = 1
    return
endif

if (extension GE 1) then begin
    state.title_extras = strcompress('Extension ' + string(extension))
endif else begin
    state.title_extras = 'Primary Image'
endelse


; Read in the image
main_image=0

; use fits_read so that extension headers will inherit primary header
; keywords.  Needed for HST ACS images.
fits_read, fitsfile, main_image, head, exten_no = extension


end

;----------------------------------------------------------------

pro atv_plainfits_read, fitsfile, head, cancelled

common atv_images

; Fits reader for plain fits files, no extensions.

main_image=0
fits_read, fitsfile, main_image, head

end

;------------------------------------------------------------------

pro atv_wfpc2_read, fitsfile, head, cancelled
    
; Fits reader for 4-panel HST WFPC2 images

common atv_state
common atv_images

droptext = strcompress('0, droplist,PC|WF2|WF3|WF4|Mosaic,' + $
                       'label_left = Select WFPC2 CCD:, set_value=0')

formdesc = [droptext, $
            '0, button, Read WFPC2 Image, quit', $
            '0, button, Cancel, quit']

textform = cw_form(formdesc, /column, title = 'WFPC2 CCD Selector')

if (textform.tag2 EQ 1) then begin ; cancelled
    cancelled = 1
    return                      
endif

ccd = long(textform.tag0) + 1

widget_control, /hourglass
if (ccd LE 4) then begin
    main_image=0
    wfpc2_read, fitsfile, main_image, head, num_chip = ccd
endif

if (ccd EQ 5) then begin
    main_image=0
    wfpc2_read, fitsfile, main_image, head, /batwing
endif
        
case ccd of
    1: state.title_extras = 'PC1'
    2: state.title_extras = 'WF2'
    3: state.title_extras = 'WF3'
    4: state.title_extras = 'WF4'
    5: state.title_extras = 'Mosaic'
    else: state.title_extras = ''
endcase

end

;------------------------------------------------------------------
pro atv_osiris_read, fitsfile, head, cancelled
    
; Fits reader for Keck OSIRIS images
;  Added by Marshall Perrin, 2005-10-13
; 
common atv_state
common atv_images

extension = 0 ; always read the primary HDU

; Make sure it's not a fits table: this would make mrdfits crash
head = headfits(fitsfile, exten=extension)
xten = strcompress(sxpar(head, 'XTENSION'), /remove_all)
if (xten EQ 'BINTABLE') then begin
    atv_message, 'File appears to be a FITS table, not an image.', $
      msgtype='error', /window
    cancelled = 1
    return
endif

if (extension GE 1) then begin
    state.title_extras = strcompress('Extension ' + string(extension))
endif else begin
    state.title_extras = 'Primary Image'
endelse

; Read in the image
delvarx, main_image
main_image = mrdfits(fitsfile, extension, head, /silent, /fscale) 

; OSIRIS writes its FITS files with a very strange axis order.
; at least, the ORP does.
main_image = transpose(main_image,[2,1,0])

; and there are tons of crap pixels
badmask = where(abs(main_image) gt 1e3,badct)
if badct gt 0 then main_image[badmask]=!values.f_nan

end

;----------------------------------------------------------------------

pro atv_2mass_read, fitsfile, head, cancelled
    
; Fits reader for 3-plane 2MASS Extended Source J/H/Ks data cube.
common atv_state
common atv_images

droptext = strcompress('0, droplist,J|H|Ks,' + $
                       'label_left = Select 2MASS Band:, set_value=0')

formdesc = [droptext, $
            '0, button, Read 2MASS Image, quit', $
            '0, button, Cancel, quit']

textform = cw_form(formdesc, /column, title = '2MASS Band Selector')

if (textform.tag2 EQ 1) then begin ; cancelled
    cancelled = 1
    return                     
endif

main_image=0
main_image = mrdfits(fitsfile, 0, head, /silent, /fscale) 

band = long(textform.tag0) 
main_image = main_image[*,*,band]    ; fixed 11/28/2000

case textform.tag0 of
    0: state.title_extras = 'J Band'
    1: state.title_extras = 'H Band'
    2: state.title_extras = 'Ks Band'
    else: state.title_extras = ''
endcase

; fix ctype2 in header to prevent crashes when running xy2ad routine:
if (strcompress(sxpar(head, 'CTYPE2'), /remove_all) EQ 'DEC---SIN') then $
  sxaddpar, head, 'CTYPE2', 'DEC--SIN'

end

;----------------------------------------------------------------------

pro atv_getdss

common atv_state
common atv_images

formdesc = ['0, text, , label_left=Object Name: , width=15, tag=objname', $
            '0, button, NED|SIMBAD, set_value=0, label_left=Object Lookup:, exclusive, tag=lookupsource', $
            '0, label, Or enter J2000 Coordinates:, CENTER', $
            '0, text, , label_left=RA   (hh:mm:ss.ss): , width=15, tag=ra', $
            '0, text, , label_left=Dec (+dd:mm:ss.ss): , width=15, tag=dec', $
            '0, droplist, 1st Generation|2nd Generation Blue|2nd Generation Red|2nd Generation Near-IR, label_left=Band:, set_value=0,tag=band ', $
            '0, float, 10.0, label_left=Image Size (arcmin; max=60): ,tag=imsize', $
            '1, base, , row', $
            '0, button, GetImage, tag=getimage, quit', $
            '0, button, Cancel, tag=cancel, quit']    

archiveform = cw_form(formdesc, /column, title = 'atv: Get DSS Image')

if (archiveform.cancel EQ 1) then return

if (archiveform.imsize LE 0.0 OR archiveform.imsize GT 60.0) then begin
    atv_message, 'Image size must be between 0 and 60 arcmin.', $
      msgtype='error', /window
    return
endif

case archiveform.band of
    0: band = '1'
    1: band = '2b'
    2: band = '2r'
    3: band = '2i'
    else: print, 'error in atv_getdss!'
endcase
    
case archiveform.lookupsource of
    0: ned = 1
    1: ned = 0  ; simbad lookup
endcase

widget_control, /hourglass
if (archiveform.objname NE '') then begin
    ; user entered object name
    querysimbad, archiveform.objname, ra, dec, found=found, ned=ned, $
      errmsg=errmsg
    if (found EQ 0) then begin
        atv_message, errmsg, msgtype='error', /window
        return
    endif
endif else begin
    ;  user entered ra, dec
    rastring = archiveform.ra
    decstring = archiveform.dec
    atv_getradec, rastring, decstring, ra, dec
endelse

; as of nov 2006, stsci server doesn't seem to recognize '2i'
; band in the way it used to.  Use eso server for 2i.
if (band NE '2i') then begin
    querydss, [ra, dec], tmpimg, tmphdr, imsize=archiveform.imsize, $
      survey=band
endif else begin
    querydss, [ra, dec], tmpimg, tmphdr, imsize=archiveform.imsize, $
      survey=band, /eso
endelse

atv, temporary(tmpimg), header=temporary(tmphdr)

end

;-----------------------------------------------------------------

pro atv_getfirst

common atv_state
common atv_images

formdesc = ['0, text, , label_left=Object Name: , width=15, tag=objname', $
            '0, button, NED|SIMBAD, set_value=0, label_left=Object Lookup:, exclusive, tag=lookupsource', $
            '0, label, Or enter J2000 Coordinates:, CENTER', $
            '0, text, , label_left=RA   (hh:mm:ss.ss): , width=15, tag=ra', $
            '0, text, , label_left=Dec (+dd:mm:ss.ss): , width=15, tag=dec', $
            '0, float, 10.0, label_left=Image Size (arcmin; max=30): ,tag=imsize', $
            '1, base, , row', $
            '0, button, GetImage, tag=getimage, quit', $
            '0, button, Cancel, tag=cancel, quit']    

archiveform = cw_form(formdesc, /column, title = 'atv: Get FIRST Image')

if (archiveform.cancel EQ 1) then return

if (archiveform.imsize LE 0.0 OR archiveform.imsize GT 30.0) then begin
    atv_message, 'Image size must be between 0 and 30 arcmin.', $
      msgtype='error', /window
    return
endif

imsize = string(round(archiveform.imsize))

case archiveform.lookupsource of
    0: ned = 1
    1: ned = 0  ; simbad lookup
endcase

widget_control, /hourglass
if (archiveform.objname NE '') then begin
    ; user entered object name
    querysimbad, archiveform.objname, ra, dec, found=found, ned=ned, $
      errmsg=errmsg
    if (found EQ 0) then begin
        atv_message, errmsg, msgtype='error', /window
        return
    endif
    
; convert decimal ra, dec to hms, dms
    sra = sixty(ra/15.0)
    rahour = string(round(sra[0]))
    ramin = string(round(sra[1]))
    rasec = string(sra[2])

    if (dec LT 0) then begin
        decsign = '-'
    endif else begin
        decsign = '+'
    endelse
    sdec = sixty(abs(dec))
 
    decdeg = strcompress(decsign + string(round(sdec[0])), /remove_all)
    decmin = string(round(sdec[1]))
    decsec = string(sdec[2])

endif else begin
    ;  user entered ra, dec
    rastring = archiveform.ra
    decstring = archiveform.dec

    rtmp = rastring
    pos = strpos(rtmp, ':')
    if (pos EQ -1) then pos = strlen(rtmp)
    rahour = strmid(rtmp, 0, pos)
    rtmp = strmid(rtmp, pos+1)
    pos = strpos(rtmp, ':')
    if (pos EQ -1) then pos = strlen(rtmp)
    ramin = strmid(rtmp, 0, pos)
    rtmp = strmid(rtmp, pos+1)
    rasec = rtmp
        
    dtmp = decstring
    pos = strpos(dtmp, ':')
    if (pos EQ -1) then pos = strlen(dtmp)
    decdeg = strmid(dtmp, 0, pos)
    dtmp = strmid(dtmp, pos+1)
    pos = strpos(dtmp, ':')
    if (pos EQ -1) then pos = strlen(dtmp)
    decmin = strmid(dtmp, 0, pos)
    dtmp = strmid(dtmp, pos+1)
    decsec = dtmp

endelse

; build the url to get image

url = 'http://third.ucllnl.org/cgi-bin/firstimage'
url = strcompress(url + '?RA=' + rahour + '%20' + ramin + '%20' + rasec, $
                  /remove_all)
url = strcompress(url + '&Dec=' + decdeg + '%20' + decmin + '%20' + $
                  decsec, /remove_all)
url = strcompress(url + '&Equinox=J2000&ImageSize=' + imsize + $
                  '&MaxInt=10&FITS=1&Download=1', /remove_all)

; now use webget to get the image
result = webget(url)


if (n_elements(result.image) LE 1) then begin
    atv_message, result.text, msgtype='error', /window
    return
endif else begin  ; valid image
    atv, result.image, header=result.imageheader
    result.header = ''
    result.text =  ''
    result.imageheader = ''
    result.image = ''
endelse

end

;-----------------------------------------------------------------

pro atv_getradec, rastring, decstring, ra, dec

; converts ra and dec strings in hh:mm:ss and dd:mm:ss to decimal degrees

rtmp = rastring
pos = strpos(rtmp, ':')
if (pos EQ -1) then pos = strlen(rtmp)
rahour = strmid(rtmp, 0, pos)
rtmp = strmid(rtmp, pos+1)
pos = strpos(rtmp, ':')
if (pos EQ -1) then pos = strlen(rtmp)
ramin = strmid(rtmp, 0, pos)
rtmp = strmid(rtmp, pos+1)
rasec = rtmp


dtmp = decstring
pos = strpos(dtmp, ':')
if (pos EQ -1) then pos = strlen(dtmp)
decdeg = strmid(dtmp, 0, pos)
dtmp = strmid(dtmp, pos+1)
pos = strpos(dtmp, ':')
if (pos EQ -1) then pos = strlen(dtmp)
decmin = strmid(dtmp, 0, pos)
dtmp = strmid(dtmp, pos+1)
decsec = dtmp

ra = 15.0 * ten([rahour, ramin, rasec])
dec = ten([decdeg, decmin, decsec])


end



;-----------------------------------------------------------------------
;     Routines for creating output graphics
;----------------------------------------------------------------------


pro atv_writefits

; Writes image to a FITS file
; If a 3D image, the option to save either the current 2D display or
; the entire cube is possible 

common atv_state
common atv_images

; Get filename to save image

filename = dialog_pickfile(filter = '*.fits', $ 
                           file = 'atv.fits', $
                           dialog_parent =  state.base_id, $
                           path = state.current_dir, $
                           get_path = tmp_dir, $
                           /write)

if (tmp_dir NE '') then state.current_dir = tmp_dir

if (strcompress(filename, /remove_all) EQ '') then return   ; cancel

if (filename EQ state.current_dir) then begin
  atv_message, 'Must indicate filename to save.', msgtype = 'error', /window
  return
endif

IF (state.image_size[2] eq 1) THEN BEGIN

  IF (ptr_valid(state.head_ptr)) THEN $
    writefits, filename, main_image, (*state.head_ptr) $
  ELSE $
    writefits, filename, main_image

ENDIF ELSE BEGIN

  formdesc = ['0, button, Write Current Image Plane, quit', $
              '0, button, Write All Image Planes, quit', $
              '0, button, Cancel, quit']

  textform = cw_form(formdesc, /column, $
                     title = 'Select Image to Write')

  IF (textform.tag0 eq 1) THEN writefits, filename, main_image
  IF (textform.tag1 eq 1) THEN BEGIN
    IF (ptr_valid(state.head_ptr)) THEN $  
      writefits, filename, main_image_stack, (*state.head_ptr) $
    ELSE $
      writefits, filename, main_image_stack
  ENDIF
 
  IF (textform.tag2 eq 1) THEN return

ENDELSE

end


;-----------------------------------------------------------------------


pro atv_writeimage

; Front-end widget to write display image to output

common atv_state
common atv_images

writeimagebase = widget_base(/row)

formdesc = ['0, droplist, PNG|JPG|TIFF,label_left=Format:,set_value=0, TAG=format ', $
            '1, base, , row', $
            '0, button, Filename:, TAG=file ', $
            '2, text, atv.png, width=30', $
            '1, base, , row', $
            '0, button, Cancel, quit, TAG=quit ', $
            '0, button, WriteImage, quit, TAG=write']

writeimageform = cw_form(writeimagebase,formdesc,/column, $
      title='ATV WriteImage', IDS=writeimage_ids_ptr) 

widget_control, writeimagebase, /realize
xmanager, 'atv_writeimage', writeimagebase, /no_block

writeimage_ids_ptr = $
  writeimage_ids_ptr(where(widget_info(writeimage_ids_ptr,/type) eq 3 OR $
  widget_info(writeimage_ids_ptr,/type) eq 8 OR $
  widget_info(writeimage_ids_ptr,/type) eq 1))

if ptr_valid(state.writeimage_ids_ptr) then ptr_free,state.writeimage_ids_ptr
state.writeimage_ids_ptr = ptr_new(writeimage_ids_ptr)


end

;----------------------------------------------------------------------

pro atv_writeimage_event, event

common atv_state


CASE event.tag OF
    
    'FORMAT': BEGIN
        
        widget_control,(*state.writeimage_ids_ptr)[2], get_value=filename
        tagpos = lonarr(5)
        tagpos[0] = strpos(filename, '.jpg', /reverse_search)
        tagpos[1] = strpos(filename, '.jpeg', /reverse_search)
        tagpos[2] = strpos(filename, '.png', /reverse_search)
        tagpos[3] = strpos(filename, '.tiff', /reverse_search)
        tagpos[4] = strpos(filename, '.tif', /reverse_search)
        maxtagpos = max(tagpos)
        if (maxtagpos GT 0) then filename = strmid(filename,0,maxtagpos)
        
        CASE event.value OF
            
            '0': BEGIN
                if (strcompress(filename, /remove_all) EQ '') then $
                  filename = 'atv'
                filename = filename + '.png'
                widget_control,(*state.writeimage_ids_ptr)[2], $
                  set_value=filename
                state.writeformat = 'PNG'
            END
            
            '1': BEGIN
                if (strcompress(filename, /remove_all) EQ '') then $
                  filename = 'atv'
                filename = filename + '.jpg'
                widget_control,(*state.writeimage_ids_ptr)[2], $
                  set_value=filename
                state.writeformat = 'JPG'
            END
            
            '2': BEGIN
                if (strcompress(filename, /remove_all) EQ '') then $
                  filename = 'atv'
                filename = filename + '.tiff'
                widget_control,(*state.writeimage_ids_ptr)[2], $
                  set_value=filename
                state.writeformat = 'TIFF'
            END
            
            ELSE:
        ENDCASE
    END
    
    'FILE': BEGIN
        filename = dialog_pickfile(/write)
        if strcompress(filename, /remove_all) EQ '' then $
          filename = strcompress('atv' + '.' + strlowcase(state.writeformat), $
                                 /remove_all)
        widget_control,(*state.writeimage_ids_ptr)[2], set_value=filename
    END
    
    'WRITE' : BEGIN
        
        atv_setwindow, state.draw_window_id
        
        widget_control,(*state.writeimage_ids_ptr)[0]
        widget_control,(*state.writeimage_ids_ptr)[2], get_value=filename
        filename = filename[0]
        
; set extension if not already set
        tagpos = lonarr(5)
        tagpos[0] = strpos(filename, '.jpg', /reverse_search)
        tagpos[1] = strpos(filename, '.jpeg', /reverse_search)
        tagpos[2] = strpos(filename, '.png', /reverse_search)
        tagpos[3] = strpos(filename, '.tiff', /reverse_search)
        tagpos[4] = strpos(filename, '.tif', /reverse_search)
        maxtagpos = max(tagpos)
        if (maxtagpos EQ -1) then begin
            filename = strcompress(filename + '.' + $
                                   strlowcase(state.writeformat), /remove_all)
            widget_control,(*state.writeimage_ids_ptr)[2], set_value=filename
        endif
        
; check for pre-existing file
        tmp_result = findfile(filename, count = nfiles)
        result = ''
        if (nfiles GT 0) then begin
            mesg = strarr(2)
            mesg[0] = 'Overwrite existing file:'
            mesg[1] = strcompress(filename + '?', /remove_all)
            result =  dialog_message(mesg, $
                                     /default_no, $
                                     dialog_parent = state.base_id, $
                                     /question)                 
        endif
        
        if (strupcase(result) EQ 'NO') then return
        
        CASE state.writeformat OF
            
            'JPG': atv_saveimage, filename, /jpg
            'TIFF': atv_saveimage, filename, /tiff
            'PNG': atv_saveimage, filename, /png
            ELSE:
        ENDCASE
        
        atv_resetwindow
        state.writeformat = 'PNG'  
        
        if ptr_valid(state.writeimage_ids_ptr) then $
          ptr_free, state.writeimage_ids_ptr
        widget_control, event.top, /destroy
    END
    
    'QUIT': BEGIN
        state.writeformat = 'PNG'
        if ptr_valid(state.writeimage_ids_ptr) then $
          ptr_free, state.writeimage_ids_ptr
        widget_control, event.top, /destroy
    END
    
    ELSE:
    
ENDCASE

end

;----------------------------------------------------------------------
pro atv_SaveToVariable, thingname
; Added by Marshall Perrin, based on code from Dave Fanning's 
; XSTRETCH.PRO

	common atv_state
	common atv_images

    case thingname of
        "Image": thing= main_image
        "Cube": thing =  main_image_stack
        "Header": begin
       		if ptr_valid(state.head_ptr) then thing = *(state.head_ptr) $
			else begin
				atv_message,'No FITS header to save! There is no FITS header currently loaded in ATV.', msgtype='error', /window
				return
			endelse
		end
    endcase

            varname = TextBox(Title='Provide Main-Level Variable Name...', Group_Leader=state.base_id, $
               Label=thingname+' Variable Name: ', Cancel=cancelled, XSize=200, Value=thingname)
            ; Dave Fanning says:
            ;
            ;; The ROUTINE_NAMES function is not documented in IDL,
            ;; so it may not always work. This capability has been
            ;; tested in IDL versions 5.3 through 5.6 and found to work.
            ;; People with IDL 6.1 and higher should use SCOPE_VARFETCH to
            ;; set main-level variables. I use the older, undocumented version
            ;; to stay compatible with more users.

            IF NOT cancelled THEN BEGIN
               dummy = Routine_Names(varname, thing, Store=1)
            ENDIF

end




;----------------------------------------------------------------------

pro atv_saveimage, file,  png=png, jpg=jpg, tiff=tiff, $
  quality=quality, dither=dither, cube=cube, quiet=quiet

common atv_state

; This program is a copy of Liam E. Gumley's SAVEIMAGE program,
; modified for use with ATV.

;- Check keyword
output = 'PNG'                  ; default output type
if keyword_set(jpg) then output = 'JPG'
if keyword_set(tiff) then output = 'TIFF'
if (n_elements(quality) eq 0) then quality = 75 ; for jpeg

;- Check for TVRD capable device
if ((!d.flags and 128)) eq 0 then begin
    atv_message, 'Unsupported graphics device- cannot create image.', $
      msgtype='error', /window
    return
endif

depth = state.bitdepth

;- Handle window devices (other than the Z buffer)
if (!d.flags and 256) ne 0 then begin
    
  ;- Copy the contents of the current display to a pixmap
    current_window = state.draw_window_id
    xsize =  state.draw_window_size[0]
    ysize = state.draw_window_size[1]
    window, /free, /pixmap, xsize=xsize, ysize=ysize, retain=2
    device, copy=[0, 0, xsize, ysize, 0, 0, current_window]
    
  ;- Set decomposed color mode for 24-bit displays
    version = float(!version.release)
    if (depth gt 8) then device, get_decomposed=entry_decomposed
    device, decomposed=1
endif

;- Read the pixmap contents into an array
if (depth gt 8) then begin
    image = tvrd(order=0, true=1)
endif else begin
    image = tvrd(order=0)
endelse

;- Handle window devices (other than the Z buffer)
if (!d.flags and 256) ne 0 then begin

  ;- Restore decomposed color mode for 24-bit displays
    if (depth gt 8) then begin
        device, decomposed=entry_decomposed
    endif

  ;- Delete the pixmap
    wdelete, !d.window
    wset, current_window

endif

;- Get the current color table
tvlct, r, g, b, /get

;- If an 8-bit image was read, reduce the number of colors
if (depth le 8) then begin
    reduce_colors, image, index
    r = r[index]
    g = g[index]
    b = b[index]
endif

; write output file

case 1 of
    (output eq 'PNG') : begin
        write_png, file, image, r, g, b
    end

    (output eq 'JPG') or (output eq 'TIFF') : begin

    ;- Convert 8-bit image to 24-bit
        if (depth le 8) then begin
            info = size(image)
            nx = info[1]
            ny = info[2]
            true = bytarr(3, nx, ny)
            true[0, *, *] = r[image]
            true[1, *, *] = g[image]
            true[2, *, *] = b[image]
            image = temporary(true)
        endif

    ;- If TIFF format output, reverse image top to bottom
        if (output eq 'TIFF') then image = reverse(temporary(image), 3)

    ;- Write the image
        case output of
            'JPG' : write_jpeg, file, image, true=1, quality=quality
            'TIFF' : write_tiff, file, image, 1
            else  :
        endcase
        
    end
    else:
endcase


end

;----------------------------------------------------------------------


pro atv_makergb

; Makes an RGB truecolor png image from the 3 blink channels.
; Can be saved using file->writeimage.
; Note- untested for 8-bit displays.  May not work there.

common atv_state
common atv_images

if (n_elements(blink_image1) EQ 1 OR $
    n_elements(blink_image2) EQ 1 OR $
    n_elements(blink_image3) EQ 1) then begin
    
    atv_message, $
      'You need to set the 3 blink channels first to make an RGB image.', $
      msgtype = 'error', /window
    return
endif

atv_getwindow

window, /free, xsize = state.draw_window_size[0], $ 
  ysize = state.draw_window_size[1], /pixmap
tempwindow = !d.window

tv, blink_image1, /true
rimage = tvrd()
tv, blink_image2, /true
gimage = tvrd()
tv, blink_image3, /true
bimage = tvrd()

tcimage = [[[rimage]], [[gimage]], [[bimage]]]

tv, tcimage, true=3

tvlct, rmap, gmap, bmap, /get
image = tvrd(/true)

wdelete, tempwindow

atv_setwindow, state.draw_window_id
tv, image, /true
atv_resetwindow


end

;----------------------------------------------------------------------

pro atv_writeps

; Writes an encapsulated postscript file of the current display.
; Calls cmps_form to get postscript file parameters.

; Note. cmps_form blocks the command line but doesn't block atv
; menus.  If we have one cmps_form active and invoke another one, it
; would crash.  Use state.ispsformon to keep track of whether we have
; one active already or not.

common atv_state
common atv_images
common atv_color

if (state.ispsformon EQ 1) then return

; cmps_form.pro crashes if atv is in blocking mode.
if (state.block EQ 1) then begin
    atv_message, 'PS output is disabled in blocking mode.', $
      msgtype = 'warning', /window
    return
endif

widget_control, /hourglass

view_min = round(state.centerpix - $
                  (0.5 * state.draw_window_size / state.zoom_factor))
; bug fix from N. Cunningham here- modified 4/14/06 to fix centering
; of overplots on the image by subtracting 1 from the max size
view_max = round(view_min + state.draw_window_size $
           / state.zoom_factor - 1)

xsize = (state.draw_window_size[0] / state.zoom_factor) > $
  (view_max[0] - view_min[0] + 1)
ysize = (state.draw_window_size[1] / state.zoom_factor) > $
  (view_max[1] - view_min[1] + 1)

aspect = float(ysize) / float(xsize)
fname = strcompress(state.current_dir + 'atv.ps', /remove_all)

atv_setwindow, state.draw_window_id
tvlct, rr, gg, bb, 8, /get
atv_resetwindow

; make sure that we don't keep the cmps_form window as the active window
external_window_id = !d.window

canceled = 1 ; in case user forcibly closes window with the button in the corner;
			  ; then cancelled doesn't get set by cmps form!
state.ispsformon = 1
psforminfo = cmps_form(cancel = canceled, create = create, $
                     aspect = aspect, parent = state.base_id, $
                     /preserve_aspect, $
                     xsize = 6.0, ysize = 6.0 * aspect, $
                     /color, /encapsulated, $
                     /nocommon, papersize='Letter', $
                     bits_per_pixel=8, $
                     filename = fname, $
                     button_names = ['Create PS File'])
atv_setwindow, external_window_id

state.ispsformon = 0
if (canceled) then return
if (psforminfo.filename EQ '') then return
tvlct, rr, gg, bb, 8


tmp_result = findfile(psforminfo.filename, count = nfiles)

result = ''
if (nfiles GT 0) then begin
    mesg = strarr(2)
    mesg[0] = 'Overwrite existing file:'
    tmp_string = $
      strmid(psforminfo.filename, $
             strpos(psforminfo.filename, state.delimiter, /reverse_search) + 1)
    mesg[1] = strcompress(tmp_string + '?', /remove_all)
    result =  dialog_message(mesg, $
                             /default_no, $
                             dialog_parent = state.base_id, $
                             /question)                 
endif

if (strupcase(result) EQ 'NO') then return
    
widget_control, /hourglass

screen_device = !d.name

; In 8-bit mode, the screen color table will have fewer than 256
; colors.  Stretch out the existing color table to 256 colors for the
; postscript plot.

set_plot, 'ps'

device, _extra = psforminfo

tvlct, rr, gg, bb, 8, /get

rn = congrid(rr, 248)
gn = congrid(gg, 248)
bn = congrid(bb, 248)

tvlct, temporary(rn), temporary(gn), temporary(bn), 8

; Make a full-resolution version of the display image, accounting for
; scalable pixels in the postscript output

newdisplay = bytarr(xsize, ysize)

startpos = abs(round(state.offset) < 0)

view_min = (0 > view_min < (state.image_size[0:1] - 1)) 
view_max = (0 > view_max < (state.image_size[0:1] - 1)) 

dimage = bytscl(scaled_image[view_min[0]:view_max[0], $
                                 view_min[1]:view_max[1]], $
                    top = 247, min=8, max=(!d.table_size-1)) + 8


newdisplay[startpos[0], startpos[1]] = temporary(dimage)

; if there's blank space around the image border, keep it black


tv, newdisplay
atv_plotall


if (state.frame EQ 1) then begin    ; put frame around image
    plot, [0], [0], /nodata, position=[0,0,1,1], $
      xrange=[0,1], yrange=[0,1], xstyle=5, ystyle=5, /noerase
    boxx = [0,0,1,1,0,0]
    boxy = [0,1,1,0,0,1]
    oplot, boxx, boxy, color=0, thick=state.framethick
endif

tvlct, temporary(rr), temporary(gg), temporary(bb), 8


device, /close
set_plot, screen_device


end

;----------------------------------------------------------------------
;       routines for defining the color maps
;----------------------------------------------------------------------

pro atv_stretchct, brightness, contrast,  getcursor = getcursor

; routine to change color stretch for given values of brightness and contrast.
; Complete rewrite 2000-Sep-21 - Doug Finkbeiner
; Updated 12/2006 to allow for brightness,contrast param input
; without changing the state.brightness and state.contrast values.
; Better for surface plots in plot window.

common atv_state
common atv_color

; if GETCURSOR then assume mouse position passed and save as
; state.brightness and state.contrast.  If no params passed, then use
; the current state.brightness and state.contrast.  If b, c passed
; without /getcursor, then make a new color table stretch for that
; brightness and contrast but don't modify the current
; state.brightness and state.contrast

; New in 2.0: scale the contrast by 0.75- gives better contrast by
; default when first starting up, and better in general with asinh
; scaling 

contrastscale=0.75

if (keyword_set(getcursor)) then begin 
    state.brightness = brightness/float(state.draw_window_size[0])
    state.contrast = contrast/float(state.draw_window_size[1])
    x = state.brightness*(state.ncolors-1)
    y = state.contrast*(state.ncolors-1)*contrastscale > 2 
endif else begin
    if (n_elements(brightness) EQ 0 OR n_elements(contrast) EQ 0) then begin
        x = state.brightness*(state.ncolors-1)
        y = state.contrast*(state.ncolors-1)*contrastscale > 2 
    endif else begin
        x = brightness*(state.ncolors-1)
        y = contrast*(state.ncolors-1)*contrastscale > 2
    endelse
endelse

;old version
;x = state.brightness*(state.ncolors-1)
;y = state.contrast*(state.ncolors-1) > 2 ; Minor change by AJB 

high = x+y & low = x-y
diff = (high-low) > 1

slope = float(state.ncolors-1)/diff ;Scale to range of 0 : nc-1
intercept = -slope*low
p = long(findgen(state.ncolors)*slope+intercept) ;subscripts to select
tvlct, r_vector[p], g_vector[p], b_vector[p], 8

end

;------------------------------------------------------------------

pro atv_initcolors

; Load a simple color table with the basic 8 colors in the lowest 
; 8 entries of the color table.  Also set top color to white.

common atv_state

rtiny   = [0, 1, 0, 0, 0, 1, 1, 1]
gtiny = [0, 0, 1, 0, 1, 0, 1, 1]
btiny  = [0, 0, 0, 1, 1, 1, 0, 1]
tvlct, 255*rtiny, 255*gtiny, 255*btiny

tvlct, [255],[255],[255], !d.table_size-1

end

;--------------------------------------------------------------------

pro atv_getct, tablenum

; Read in a pre-defined color table, and invert if necessary.

common atv_color
common atv_state
common atv_images


atv_setwindow, state.draw_window_id
loadct, tablenum, /silent, bottom=8
tvlct, r, g, b, 8, /get

atv_initcolors

r = r[0:state.ncolors-2]
g = g[0:state.ncolors-2]
b = b[0:state.ncolors-2]

if (state.invert_colormap EQ 1) then begin
r = reverse(r)
g = reverse(g)
b = reverse(b)
endif

r_vector = r
g_vector = g
b_vector = b

atv_stretchct

; need this to re-set to external color table
atv_resetwindow

if (state.bitdepth EQ 24 AND (n_elements(pan_image) GT 10) ) then $
  atv_refresh


end

;--------------------------------------------------------------------

pro atv_makect, tablename

; Define new color tables here.  Invert if necessary.

common atv_state
common atv_color

case tablename of
    'ATV Special': begin

       w = findgen(256)
       
       sigma = 60.
       center = 140
       r = 255.* exp(-1.*(w - center)^2 / (2.*sigma^2))
       r[center:255] = 255.
       
       sigma = 60
       center = 255
       g = 255. * exp(-1.*(w - center)^2 / (2.*sigma^2))
       
       sigma = 60
       center = 40
       b = 255. * exp(-1.*(w - center)^2 / (2.*sigma^2))
       
       b[0:center-1] = findgen(center)^0.5 / center^0.5 * 255.
       center = 30
       b[(255-center+1):255] = findgen(center)^2 / center^2 *255.
       
   end

    'Velocity2': begin
       r = fltarr(256)
       r[0:127] = 128. - findgen(128)
       r[128:255] = 255
       
       g = fltarr(256)
       g[0:127] = findgen(128)^1.5
       g[128:255] = reverse(g[0:127])
       g = g / max(g) * 255.
       
       b = 255. - findgen(256)
       b[128:255] = findgen(128)^3 / 128.^2
      
    end

    'Velocity1': begin
       w = findgen(256)

       sigma = 25.
       center = 170
       r = 255.* exp(-1.*(w - center)^2 / (2.*sigma^2))
       r[center:255] = 255.
       sigma = 30.
       center = 0.
       r = r + 100.* exp(-1.*(w - center)^2 / (2.*sigma^2))
       
       sigma = 30.
       center1 = 100.
       g = fltarr(256)
       g[0:center1] = 255. * exp(-1.*(w[0:center1] - center1)^2 / (2.*sigma^2))
       sigma = 60.
       center2 = 140.
       g[center1:center2] = 255.
       g[center2:255] = $
          255. * exp(-1.*(w[center2:255] - center2)^2 / (2.*sigma^2))
       
       sigma = 40.
       center = 70
       b = 255.* exp(-1.*(w - center)^2 / (2.*sigma^2))
       b[0:center] = 255.
           
    end
    
; add more color table definitions here as needed...
    else: return

endcase

r = congrid(r, state.ncolors)
g = congrid(g, state.ncolors)
b = congrid(b, state.ncolors)


if (state.invert_colormap EQ 1) then begin
r = reverse(r)
g = reverse(g)
b = reverse(b)
endif

r_vector = temporary(r)
g_vector = temporary(g)
b_vector = temporary(b)

atv_stretchct

; need this to preserve external color map
atv_resetwindow

if (state.bitdepth EQ 24) then atv_refresh

end

;----------------------------------------------------------------------

function atv_icolor, color

; Routine to reserve the bottom 8 colors of the color table
; for plot overlays and line plots.

if (n_elements(color) EQ 0) then return, 1

ncolor = N_elements(color)

; If COLOR is a string or array of strings, then convert color names
; to integer values
if (size(color,/tname) EQ 'STRING') then begin ; Test if COLOR is a string
    
; Detemine the default color for the current device
    if (!d.name EQ 'X') then defcolor = 7 $ ; white for X-windows
    else defcolor = 0           ; black otherwise
    
    icolor = 0 * (color EQ 'black') $
      + 1 * (color EQ 'red') $
      + 2 * (color EQ 'green') $
      + 3 * (color EQ 'blue') $
      + 4 * (color EQ 'cyan') $
      + 5 * (color EQ 'magenta') $
      + 6 * (color EQ 'yellow') $
      + 7 * (color EQ 'white') $
      + defcolor * (color EQ 'default')
    
endif else begin
    icolor = long(color)
endelse

return, icolor
end 
 
;---------------------------------------------------------------------
;    routines dealing with image header, title,  and related info
;--------------------------------------------------------------------

pro atv_settitle

; Update title bar with the image file name

common atv_state

if (strlen(state.imagename) EQ 0) then begin
    if (strlen(state.title_extras) EQ 0) then begin
        state.window_title = 'atv'
    endif else begin
        state.window_title = strcompress('atv:  ' + state.title_extras)
    endelse
endif else begin

    slash = strpos(state.imagename, state.delimiter, /reverse_search)
    if (slash NE -1) then name = strmid(state.imagename, slash+1) $
      else name = state.imagename
    state.window_title = strcompress('atv:  '+name + '  ' + state.title_extras)
endelse

widget_control, state.base_id, tlb_set_title = state.window_title

end

;----------------------------------------------------------------------

pro atv_setheader, head

; Routine to keep the image header using a pointer to a 
; heap variable.  If there is no header (i.e. if atv has just been
; passed a data array rather than a filename), then make the
; header pointer a null pointer.  Get astrometry info from the 
; header if available.  If there's no astrometry information, set 
; state.astr_ptr to be a null pointer.

common atv_state

; Kill the header info window when a new image is read in

if (xregistered('atv_headinfo')) then begin
    widget_control, state.headinfo_base_id, /destroy
endif

if (xregistered('atv_stats')) then begin
    widget_control, state.stats_base_id, /destroy
endif

state.cunit = ''

if (n_elements(head) LE 1) then begin
; If there's no image header...
    state.wcstype = 'none'
    ptr_free, state.head_ptr
    state.head_ptr = ptr_new()
    ptr_free, state.astr_ptr
    state.astr_ptr = ptr_new()
    widget_control, state.wcs_bar_id, set_value = '---No WCS Info---'
    return
endif

ptr_free, state.head_ptr
state.head_ptr = ptr_new(head)

; Get astrometry information from header, if it exists
ptr_free, state.astr_ptr        ; kill previous astrometry info
state.astr_ptr = ptr_new()
extast, head, astr, noparams

; No valid astrometry in header
if (noparams EQ -1) then begin 
    widget_control, state.wcs_bar_id, set_value = '---No WCS Info---'
    state.wcstype = 'none'
    return
endif

; Here: add escape clauses for any WCS types that cause crashes.  Add
; more as needed 
checkastr = strcompress(string(astr.ctype[0]), /remove_all)
if ( (checkastr EQ 'PIXEL') OR $
     (checkastr EQ '') OR $
     (checkastr EQ 'COLUMN#') ) then begin
    widget_control, state.wcs_bar_id, set_value = '---No WCS Info---'
    state.wcstype = 'none'
    return
endif

if (checkastr EQ 'RA---TNX') then begin
   widget_control, state.wcs_bar_id, set_value = '---No WCS Info---'
   state.wcstype = 'none'
   print
   print, 'WARNING- WCS info is in unsupported TNX format.'
   return
endif

; Image is a 2-d calibrated spectrum:
; (these keywords work for HST STIS 2-d spectral images)
if (astr.ctype[0] EQ 'LAMBDA' OR astr.ctype[0] EQ 'WAVE') then begin
    state.wcstype = 'lambda'
    state.astr_ptr = ptr_new(astr)
    widget_control, state.wcs_bar_id, set_value = '                 '

    state.cunit = sxpar(*state.head_ptr, 'cunit1')
    state.cunit = strcompress(string(state.cunit), /remove_all)
    if (state.cunit NE '0') then begin
        state.cunit = strcompress(strupcase(strmid(state.cunit,0,1)) + $
                                  strmid(state.cunit,1), $
                            /remove_all)
    endif else begin
        state.cunit = ''
    endelse   
    return
endif

; 2-D wavelength calibrated spectrum from iraf gemini reductions:
if (string(sxpar(head, 'WAT1_001')) EQ $
    'wtype=linear label=Wavelength units=angstroms') then begin
    state.wcstype = 'lambda'
    state.astr_ptr = ptr_new(astr)
    widget_control, state.wcs_bar_id, set_value = '                 '
    state.cunit = 'Angstrom'
    return
endif


; final error check on WCS, in case it's in a format that can't be
; understood by the idlastro routines.
catch, error_status

if (error_status NE 0) then begin
   print
   print, 'Warning: WCS information could not be understood.'
   wcsstring = '---No WCS Info---'
   state.wcstype='none'
   return
endif

; see if coordinates can be extracted without an error
xy2ad, 0, 0, astr, lon, lat

catch, /cancel


; Good astrometry info in header:
state.wcstype = 'angle'
widget_control, state.wcs_bar_id, set_value = '                 '

; Check for GSS type header  
if strmid( astr.ctype[0], 5, 3) EQ 'GSS' then begin
    hdr1 = head
    gsss_STDAST, hdr1
    extast, hdr1, astr, noparams
endif

; Create a pointer to the header info
state.astr_ptr = ptr_new(astr)

; Get the equinox of the coordinate system
equ = get_equinox(head, code)

if (code NE -1) then begin
    if (equ EQ 2000.0) then state.equinox = 'J2000'
    if (equ EQ 1950.0) then state.equinox = 'B1950'
    if (equ NE 2000.0 and equ NE 1950.0) then $
      state.equinox = string(equ, format = '(f6.1)')
endif else begin
    IF (strmid(astr.ctype[0], 0, 4) EQ 'GLON') THEN BEGIN 
        state.equinox = 'J2000' ; (just so it is set)
    ENDIF ELSE BEGIN                          
; HST/ACS images don't have an equinox or epoch in the header! And,
; there's no information in the headers for individual extensions that
; uniquely identifies the data as being from ACS.  If we've gotten to
; this point already, then we know that there is WCS info in the
; header, though.  Assume 2000.0.  This is a temporary fix until I
; think of a better way to deal with it.
        state.equinox = 'J2000'
        print, 'Warning: WCS equinox not given in image header.  Assuming J2000.'

; Old version: this turned off WCS when equinox was missing:
;        ptr_free, state.astr_ptr    ; clear pointer
;        state.astr_ptr = ptr_new()
;        state.equinox = 'J2000'
;        state.wcstype = 'none'
;        widget_control, state.wcs_bar_id, set_value = '---No WCS Info---'
    ENDELSE 
endelse

; Set default display to native system in header
state.display_equinox = state.equinox
state.display_coord_sys = strmid(astr.ctype[0], 0, 4)

end

;---------------------------------------------------------------------


pro atv_headinfo

common atv_state

; If there's no header, kill the headinfo window and exit this
; routine.
if (not(ptr_valid(state.head_ptr))) then begin
    if (xregistered('atv_headinfo')) then begin
        widget_control, state.headinfo_base_id, /destroy
    endif

    atv_message, 'No header information available for this image!', $
      msgtype = 'error', /window
    return
endif


; If the user asks to redisplay the header info, and the header window
; is already up, then kill the old one so we can create the new.
; This is sort of unnecessary, but has the effect of raising it to the top.
   if (xregistered('atv_headinfo')) then begin
        widget_control, state.headinfo_base_id, /destroy
    endif

; If there is header information but not headinfo window,
; create the headinfo window.
if (not(xregistered('atv_headinfo', /noshow))) then begin

    headinfo_base = $
      widget_base(/base_align_right, $
                  group_leader = state.base_id, $
                  /column, $
                  title = 'atv image header information', $
                  uvalue = 'headinfo_base')
    state.headinfo_base_id = headinfo_base

    h = *(state.head_ptr)

    headinfo_text = widget_text(headinfo_base, $
                            /scroll, $
                            value = h, $
                            xsize = 85, $
                            ysize = 24)
    
    headinfo_done = widget_button(headinfo_base, $
                              value = 'Done', $
                              uvalue = 'headinfo_done')

    widget_control, headinfo_base, /realize
    xmanager, 'atv_headinfo', headinfo_base, /no_block

endif


end

;---------------------------------------------------------------------

pro atv_headinfo_event, event

common atv_state

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'headinfo_done': widget_control, event.top, /destroy
    else:
endcase

end

;----------------------------------------------------------------------
;             routines to do plot overlays
;----------------------------------------------------------------------

pro atv_plot1plot, iplot
common atv_pdata
common atv_state

; Plot a point or line overplot on the image

atv_setwindow, state.draw_window_id

widget_control, /hourglass

oplot, [(*(plot_ptr[iplot])).x], [(*(plot_ptr[iplot])).y], $
  _extra = (*(plot_ptr[iplot])).options

atv_resetwindow
state.newrefresh=1
end

;----------------------------------------------------------------------

pro atv_plot1text, iplot
common atv_pdata
common atv_state

; Plot a text overlay on the image
atv_setwindow, state.draw_window_id

widget_control, /hourglass

xyouts, (*(plot_ptr[iplot])).x, (*(plot_ptr[iplot])).y, $
  (*(plot_ptr[iplot])).text, _extra = (*(plot_ptr[iplot])).options

atv_resetwindow
state.newrefresh=1
end

;----------------------------------------------------------------------

pro atv_plot1arrow, iplot
common atv_pdata
common atv_state

; Plot a arrow overlay on the image
atv_setwindow, state.draw_window_id

widget_control, /hourglass

arrow, (*(plot_ptr[iplot])).x1, (*(plot_ptr[iplot])).y1, $
  (*(plot_ptr[iplot])).x2, (*(plot_ptr[iplot])).y2, $
  _extra = (*(plot_ptr[iplot])).options, /data

atv_resetwindow
state.newrefresh=1
end


;--------------------------------------------------------------------------------
pro atvpol, q, u, magnification=magnification, noxmcheck=noxmcheck,$
	polmask=polmask,$
	_extra = options
common atv_pdata
common atv_state

; Routine to read in polarizatin plot data and options, store in a heap
; variable structure, and plot the polarization
if not(keyword_set(noxmcheck)) then $
if (not(xregistered('atv', /noshow))) then begin
    print, 'You need to start ATV first!'
    return
endif

if (N_params() LT 1) then begin
   print, 'Too few parameters for ATVPLOT.'
   return
endif

if (n_elements(options) EQ 0) then options = {color: 'red'}

if (nplot LT maxplot) then begin
   nplot = nplot + 1

;  convert color names to index numbers, and set default=red
   c = where(tag_names(options) EQ 'COLOR', count)
   if (count EQ 0) then options = create_struct(options, 'color', 'red')
   options.color = atv_icolor(options.color)

   c = where(tag_names(options) EQ 'MINVAL', count)
   if (count EQ 0) then options = create_struct(options, 'minval', '0.02')

   c = where(tag_names(options) EQ 'THETAOFFSET', count)
   if (count EQ 0) then options = create_struct(options, 'thetaoffset', '0.00')
	;stop
   c = where(tag_names(options) EQ 'POLMASK', count)
; Don't create a polmask if there isn't one supplied. 
   
;   if not(keyword_set(polmask)) then polmask = byte(q*0)
   if (count EQ 0) then if n_elements(polmask) gt 0 then options = create_struct(options, 'polmask', polmask)


   if not(keyword_set(magnification)) then magnification=1e6
	magnification = float(magnification)
   options = create_struct(options, 'magnification', magnification)


   pstruct = {type: 'polarization',   $     ; points
              q: q,             $     ; q stokes
              u: u,             $     ; u stokes
              options: options  $     ; plot keyword options
             }

   plot_ptr[nplot] = ptr_new(pstruct)
	state.polarim_present=1
	state.polarim_plotindex=nplot

   atv_plotwindow
   atv_plot1pol, nplot

endif else begin
   print, 'Too many calls to ATVPLOT.'
endelse

end



;---------------------------------------------------------------------
pro atv_plot1pol, iplot
; a version of polvect.pro modified for atv. 
;
; This overprints polarization vectors onto the image.
common atv_pdata
common atv_state


if state.polarim_display eq 0 then return

atv_setwindow, state.draw_window_id

widget_control, /hourglass

resample =  6 <  8/state.zoom_factor > 1
lengthscale=4 <  4./state.zoom_factor > 1
lengthscale = 1.0*lengthscale*(*(plot_ptr[iplot])).options.magnification
minmag = (*(plot_ptr[iplot])).options.minval
;print,"zoom: ",state.zoom_factor
;print,"resample: ",resample
;print,"minmag: ",minmag
;print,"lengthscale: ",lengthscale

if not(keyword_set(resample)) then resample=1 ; for use later to simplify the code...
color = (*(plot_ptr[iplot])).options.color
thetaoffset = (*(plot_ptr[iplot])).options.thetaoffset
;print,"thetaoffset: ",thetaoffset

qo = (*(plot_ptr[iplot])).q
uo = (*(plot_ptr[iplot])).u
if resample eq 1 then begin
    q = (*(plot_ptr[iplot])).q
    u = (*(plot_ptr[iplot])).u
endif else begin
    sz = size((*(plot_ptr[iplot])).q)
    q = congrid((*(plot_ptr[iplot])).q,sz[1]/resample,sz[2]/resample,/int)
    if arg_present(xo) then x = congrid(xo,sz[1]/resample)
    u = congrid((*(plot_ptr[iplot])).u,sz[1]/resample,sz[2]/resample,/int)
    if arg_present(xo) then y = congrid(yo,sz[2]/resample)
endelse
sz = size(q)
x = (findgen(sz[1]))*resample
y = (findgen(sz[2]))*resample
mag = sqrt(u^2.+q^2.)             ;magnitude.
nbad = 0                        ;# of missing points
        ;if n_elements(missing) gt 0 then begin
                ;good = where(mag lt missing)
                ;if keyword_set(dots) then bad = where(mag ge missing, nbad)
        ;endif else begin
                szmag = size(mag)



c = where(tag_names((*(plot_ptr[iplot])).options) EQ 'POLMASK', polmaskpresent)

	; try rescaling?
;	mag = mag-minmag
;	minmag=0

	if polmaskpresent ne 0 then begin
		maskresize = congrid( (*(plot_ptr[iplot])).options.polmask, sz[1],sz[2])
		good = where (maskresize ne 0)
	endif else good  =where (mag gt state.polarim_lowthresh and mag lt state.polarim_highthresh)

if n_elements(good) eq 1 then return
;;stop
        ;if n_elements(missing) gt 0 then begin
                ;good = where(mag lt missing)
                ;if keyword_set(dots) then bad = where(mag ge missing, nbad)
        ;endif else begin
				

        ugood = u[good]
        qgood = q[good]
        x0 = min(x,max=x1,/NaN)                     ;get scaling
        y0 = min(y,max=y1,/NaN)
    	x_step=(x1-x0)/(sz[1]-1.0)   ; Convert to float. Integer math
    	y_step=(y1-y0)/(sz[2]-1.0)   ; could result in divide by 0
    	theta = 0.5 * atan(u,q)+thetaoffset*!dtor
		;stop
    	maxmag = .50

    ; remember position angle is defined starting at NORTH so
    ; the usual positions of sin and cosine are swapped.
    deltax = -lengthscale * mag * sin(theta)/2
    deltay = lengthscale * mag * cos(theta)/2
    x_b0=x0-x_step
    x_b1=x1+x_step
    y_b0=y0-y_step
    y_b1=y1+y_step
    if n_elements(clip) eq 0 then $
        clip = [!x.crange[0],!y.crange[0],!x.crange[1],!y.crange[1]]
    r = 0                          ;len of arrow head

	ys = y[good /sz[1]]
	xs = x[good mod sz[1]]
	x0 = xs-deltax[good]
	x1 = xs+deltax[good]
	y0 = ys-deltay[good]
	y1 = ys+deltay[good]
    for i=0l,n_elements(good)-1 do begin     ;Each point

			;if mag[good[i]] lt minmag then continue ; skip if it's too small
                ;x0 = x[good[i] mod (sz[1])]        ;get coords of start & end
                ;dx = deltax[good[i]]
                ;y0 = y[good[i] / (sz[1])]
                ;dy = deltay[good[i]]

                ;if keyword_set(badmask) then $
                        ;if badmask[x0,y0] eq 0 then continue ; don't print for b
                ;if x0 eq 0 or y0 eq 0 then continue ; don't print on edges.
                if deltax[good[i]] gt 20 or deltay[good[i]] gt 20 then continue; don't print huge garbage 
                plots,[x0[i],x1[i]], [y0[i],y1[i]],$
                      color=color,noclip=0
    endfor

;%%%%%%%%%%%%%%%%%%%%%%%%%

atv_resetwindow
state.newrefresh=1
end

;----------------------------------------------------------------------

function atv_degperpix, hdr 
             
; This program calculates the pixel scale (deg/pixel) and returns the value

common atv_state

On_error,2                      ;Return to caller

extast, hdr, bastr, noparams    ;extract astrom params in deg.

a = bastr.crval[0]
d = bastr.crval[1]

factor = 60.0                   ;conversion factor from deg to arcmin
d1 = d + (1/factor)             ;compute x,y of crval + 1 arcmin

proj = strmid(bastr.ctype[0],5,3)

case proj of 
    'GSS': gsssadxy, bastr, [a,a], [d,d1], x, y
    else:  ad2xy, [a,a], [d,d1], bastr, x, y 
endcase

dmin = sqrt( (x[1]-x[0])^2 + (y[1]-y[0])^2 ) ;det. size in pixels of 1 arcmin

; Convert to degrees per pixel and return scale
degperpix = 1. / dmin / 60.

return, degperpix
end

;----------------------------------------------------------------------

function atv_wcs2pix, coords, coord_sys=coord_sys, line=line

common atv_state

; check validity of state.astr_ptr and state.head_ptr before
; proceeding to grab wcs information

if ptr_valid(state.astr_ptr) then begin
    ctype = (*state.astr_ptr).ctype
    equinox = state.equinox
    disp_type = state.display_coord_sys
    disp_equinox = state.display_equinox
    disp_base60 = state.display_base60
    bastr = *(state.astr_ptr)
    
; function to convert an ATV region from wcs coordinates to pixel coordinates
    degperpix = atv_degperpix(*(state.head_ptr))
    
; need numerical equinox values
    IF (equinox EQ 'J2000') THEN num_equinox = 2000.0 ELSE $
      IF (equinox EQ 'B1950') THEN num_equinox = 1950.0 ELSE $
      num_equinox = float(equinox)
    
    headtype = strmid(ctype[0], 0, 4)
    n_coords = n_elements(coords)
endif

case coord_sys of
    
    'j2000': begin
        if (strpos(coords[0], ':')) ne -1 then begin
            ra_arr = strsplit(coords[0],':',/extract)
            dec_arr = strsplit(coords[1],':',/extract)
            ra = ten(float(ra_arr[0]), float(ra_arr[1]), $
                     float(ra_arr[2])) * 15.0
            dec = ten(float(dec_arr[0]), float(dec_arr[1]), $
                      float(dec_arr[2]))
            if (keyword_set(line)) then begin
                ra1_arr = strsplit(coords[2],':',/extract)
                dec1_arr = strsplit(coords[3],':',/extract)
                ra1 = ten(float(ra1_arr[0]), float(ra1_arr[1]), $
                          float(ra1_arr[2])) * 15.0
                dec1 = ten(float(dec1_arr[0]), float(dec1_arr[1]), $
                           float(dec1_arr[2]))
            endif
        endif else begin        ; coordinates in degrees
            ra=float(coords[0])
            dec=float(coords[1])
            if (keyword_set(line)) then begin
                ra1=float(coords[2])
                dec1=float(coords[3])  
            endif
        endelse
        
        if (not keyword_set(line)) then begin
            if (n_coords ne 6) then $
              coords[2:n_coords-2] = $
              strcompress(string(float(coords[2:n_coords-2]) / $
                                 (degperpix * 60.)),/remove_all) $
            else $
              coords[2:n_coords-3] = $
              strcompress(string(float(coords[2:n_coords-3]) / $
                                 (degperpix * 60.)),/remove_all)
        endif
        
    end
    
    'b1950': begin
        if (strpos(coords[0], ':')) ne -1 then begin
            ra_arr = strsplit(coords[0],':',/extract)
            dec_arr = strsplit(coords[1],':',/extract)
            ra = ten(float(ra_arr[0]), float(ra_arr[1]), $
                     float(ra_arr[2])) * 15.0
            dec = ten(float(dec_arr[0]), float(dec_arr[1]), float(dec_arr[2]))
            precess, ra, dec, 1950.0, 2000.0
            if (keyword_set(line)) then begin
                ra1_arr = strsplit(coords[2],':',/extract)
                dec1_arr = strsplit(coords[3],':',/extract)
                ra1 = ten(float(ra1_arr[0]), float(ra1_arr[1]), $
                          float(ra1_arr[2])) * 15.0
                dec1 = ten(float(dec1_arr[0]), float(dec1_arr[1]), $
                           float(dec1_arr[2]))
                precess, ra1, dec1, 1950.0,2000.0
            endif
        endif else begin      ; convert B1950 degrees to J2000 degrees
            ra = float(coords[0])
            dec = float(coords[1]) 
            precess, ra, dec, 1950.0, 2000.0
            if (keyword_set(line)) then begin
                ra1=float(coords[2])
                dec1=float(coords[3])
                precess, ra1, dec1, 1950., 2000.0 
            endif
        endelse
        
        if (not keyword_set(line)) then begin
            if (n_coords ne 6) then $
              coords[2:n_coords-2] = $
              strcompress(string(float(coords[2:n_coords-2]) / $
                                 (degperpix * 60.)),/remove_all) $
    else $
              coords[2:n_coords-3] = $
              strcompress(string(float(coords[2:n_coords-3]) / $
                                 (degperpix * 60.)),/remove_all)
        endif
    end
    
    'galactic': begin           ; convert galactic to J2000 degrees
        euler, float(coords[0]), float(coords[1]), ra, dec, 2
        if (not keyword_set(line)) then begin
            if (n_coords ne 6) then $
              coords[2:n_coords-2] = $
              strcompress(string(float(coords[2:n_coords-2]) / $
                                 (degperpix * 60.)),/remove_all) $
            else $
              coords[2:n_coords-3] = $
              strcompress(string(float(coords[2:n_coords-3]) / $
                                 (degperpix * 60.)),/remove_all)
        endif else begin
            euler, float(coords[2]), float(coords[3]), ra1, dec1, 2
        endelse
    end
    
    'ecliptic': begin           ; convert ecliptic to J2000 degrees
  euler, float(coords[0]), float(coords[1]), ra, dec, 4
  if (not keyword_set(line)) then begin
      if (n_coords ne 6) then $ 
        coords[2:n_coords-2] = $
        strcompress(string(float(coords[2:n_coords-2]) / $
                           (degperpix * 60.)),/remove_all) $
      else $
        coords[2:n_coords-3] = $
        strcompress(string(float(coords[2:n_coords-3]) / $
                           (degperpix * 60.)),/remove_all)
  endif else begin
      euler, float(coords[2]), float(coords[3]), ra1, dec1, 4
  endelse
end

'current': begin
    ra_arr = strsplit(coords[0],':',/extract)
    dec_arr = strsplit(coords[1],':',/extract)
    ra = ten(float(ra_arr[0]), float(ra_arr[1]), float(ra_arr[2])) * 15.0
    dec = ten(float(dec_arr[0]), float(dec_arr[1]), float(dec_arr[2]))
    if (not keyword_set(line)) then begin
        coords[2] = strcompress(string(float(coords[2]) / $
                                       (degperpix * 60.)),/remove_all)
        if (n_coords gt 3) then $
          coords[3] = strcompress(string(float(coords[3]) / $
                                         (degperpix * 60.)),/remove_all)
    endif else begin
      ra1_arr = strsplit(coords[2],':',/extract)
      dec1_arr = strsplit(coords[3],':',/extract)
      ra1 = ten(float(ra1_arr[0]), float(ra1_arr[1]), float(ra1_arr[2])) * 15.0
      dec1 = ten(float(dec1_arr[0]), float(dec1_arr[1]), float(dec1_arr[2]))
  endelse
  
  if (num_equinox ne 2000.) then begin
      precess, ra, dec, num_equinox, 2000.
      if (keyword_set(line)) then precess, ra1, dec1, num_equinox, 2000.
  endif
  
end

'pixel': begin
; Do nothing when pixel.  Will pass pixel coords array back.
end

else: 

endcase

if (ptr_valid(state.astr_ptr) AND coord_sys ne 'pixel') then begin
    
    if (num_equinox ne 2000) then begin
        precess, ra, dec, 2000., num_equinox
        if (keyword_set(line)) then precess, ra1, dec1, 2000., num_equinox
    endif
    
    proj = strmid(ctype[0],5,3)
    
    case proj of 
        'GSS': begin
            gsssadxy, bastr, ra, dec, x, y
            if (keyword_set(line)) then gsssadxy, bastr, ra1, dec1, x1, y1
        end
        else: begin
            ad2xy, ra, dec, bastr, x, y 
            if (keyword_set(line)) then ad2xy, ra1, dec1, bastr, x1, y1 
        end
    endcase
    
    coords[0] = strcompress(string(x),/remove_all)
    coords[1] = strcompress(string(y),/remove_all)
    if (keyword_set(line)) then begin
        coords[2] = strcompress(string(x1),/remove_all)
        coords[3] = strcompress(string(y1),/remove_all)
    endif
endif

return, coords
END

;----------------------------------------------------------------------

function atv_set_spaces, minval, spacing, number, extra, spaces
; Routine to set tick lines about image

; Actually set the number of tick marks to 4 more than what should be drawn to
; have marks at either edge

number = number + extra
if (minval eq -90.) then offset = 0 else if (extra lt 0.) then $
  offset = -extra else offset = float(extra/2.0)
  spaces = (indgen(number) - offset) * spacing + $
           round(minval/spacing) * spacing 
  return, number
end

;----------------------------------------------------------------------

pro atv_find_spacing, range, number, coord_type, spacing

; range : coordinate range of image in degrees
; number : number of lines 
; coord_type : 0= h m s
;              1= decimal degrees
;
; spacing : array containing coordinate spacing of lines

; Determine tick spacing in coordinate units
; Need not to round up if a coordinate range is complete, otherwise, there is
; an offset between coordinate lines which should overlap
tx = range / number

; Need to find best spacing to use (in terms of some reasonable multiples
; of the coordinate units , this is assuming that the units are in degrees

if (tx gt 50.0) then tx = round(tx/5.0) * 5.0 else $
if (tx gt 10.0) then tx = round(tx) * 1.0 else $
if (tx gt 1.0) then tx = round(tx/0.5) * 0.5 else $

if (coord_type eq 0) then begin
  if (tx gt 0.5) then tx = round(tx/(10.0/60.0)) * (10.0/60.0) else $
  if (tx gt 10.0/60.0) then tx = round(tx/(2.0/60.0)) * (2.0/60.0) else $
  if (tx gt 5.0/60.0) then tx = round(tx/(1.0/60.0)) * (1.0/60.0) else $
  if (tx gt 1.0/60.0) then tx = round(tx/(0.5/60.0)) * (0.5/60.0) else $
  if (tx gt 0.5/60.0) then tx = round(tx/(10.0/3600.0)) * (10.0/3600.0) else $
  if (tx gt 10.0/3600.0) then tx = round(tx/(1.0/3600.0)) * (1.0/3600.0) $
  else tx = round(tx/(1.0/36000.0)) * (1.0/36000.0) 
endif else begin
  if (tx gt 0.5) then tx = round(tx/0.2) * 0.2 else $
  if (tx gt 0.1) then tx = round(tx/0.1) * 0.1 else $
  if (tx gt 0.05) then tx = round(tx/0.05) * 0.05 else $
  if (tx gt 0.01) then tx = round(tx/0.01) * 0.01 else $
  if (tx gt 0.005) then tx = round(tx/0.005) * 0.005 else $
  if (tx gt 0.001) then tx = round(tx/0.001) * 0.001 else $
  if (tx gt 0.0005) then tx = round(tx/0.0005) * 0.0005 $
  else tx = round(tx/0.0001) * 0.0001
endelse

spacing = tx
return

end

;----------------------------------------------------------------------

pro atv_plot1region, iplot
common atv_pdata
common atv_state

; Plot a region overlay on the image
atv_setwindow, state.draw_window_id

widget_control, /hourglass

reg_array = (*(plot_ptr[iplot])).reg_array
n_reg = n_elements(reg_array)

for i=0, n_reg-1 do begin
    open_parenth_pos = strpos(reg_array[i],'(')
    close_parenth_pos = strpos(reg_array[i],')')   
    reg_type = strcompress(strmid(reg_array[i],0,open_parenth_pos),/remove_all)
    length = close_parenth_pos - open_parenth_pos
    coords_str = strcompress(strmid(reg_array[i], open_parenth_pos+1, $
                                    length-1),/remove_all)
    coords_arr = strsplit(coords_str,',',/extract) 
    n_coords = n_elements(coords_arr)
    color_begin_pos = strpos(strlowcase(reg_array[i]), 'color')
    text_pos = strpos(strlowcase(reg_array[i]), 'text')

  if (color_begin_pos ne -1) then begin
    color_equal_pos = strpos(reg_array[i], '=', color_begin_pos)
   endif

  text_begin_pos = strpos(reg_array[i], '{')

; Text for region
  if (text_begin_pos ne -1) then begin
    text_end_pos = strpos(reg_array[i], '}')
    text_len = (text_end_pos-1) - (text_begin_pos)
    text_str = strmid(reg_array[i], text_begin_pos+1, text_len)
    color_str = ''

; Color & Text for region
    if (color_begin_pos ne -1) then begin
    ; Compare color_begin_pos to text_begin_pos to tell which is first
      
      case (color_begin_pos lt text_begin_pos) of
        0: begin
             ;text before color
           color_str = strcompress(strmid(reg_array[i], color_equal_pos+1, $
                       strlen(reg_array[i])), /remove_all)
        end
        1: begin
             ;color before text
           len_color = (text_pos-1) - color_equal_pos
           color_str = strcompress(strmid(reg_array[i], color_equal_pos+1, $
                                   len_color), /remove_all)
        end
      else:
      endcase
    endif

  endif else begin

; Color but no text for region
    if (color_begin_pos ne -1) then begin
      color_str = strcompress(strmid(reg_array[i], color_equal_pos+1, $
                  strlen(reg_array[i])), /remove_all)

; Neither color nor text for region
    endif else begin
      color_str = ''
    endelse

    text_str = ''

  endelse

  index_j2000 = where(strlowcase(coords_arr) eq 'j2000')
  index_b1950 = where(strlowcase(coords_arr) eq 'b1950')
  index_galactic = where(strlowcase(coords_arr) eq 'galactic')
  index_ecliptic = where(strlowcase(coords_arr) eq 'ecliptic')

  index_coord_system = where(strlowcase(coords_arr) eq 'j2000') AND $
                       where(strlowcase(coords_arr) eq 'b1950') AND $
                       where(strlowcase(coords_arr) eq 'galactic') AND $
                       where(strlowcase(coords_arr) eq 'ecliptic')

  index_coord_system = index_coord_system[0]

if (index_coord_system ne -1) then begin

; Check that a WCS region is not overplotted on image with no WCS
  if (NOT ptr_valid(state.astr_ptr)) then begin
    atv_message, 'WCS Regions cannot be displayed on image without WCS', $
      msgtype='error', /window
    ;Erase pstruct that was formed for this region.
    atverase, 1
    return
  endif

  case strlowcase(coords_arr[index_coord_system]) of
  'j2000': begin
     if (strlowcase(reg_type) ne 'line') then $
       coords_arr = atv_wcs2pix(coords_arr, coord_sys='j2000') $
     else $
       coords_arr = atv_wcs2pix(coords_arr, coord_sys='j2000', /line) 
   end
  'b1950': begin
     if (strlowcase(reg_type) ne 'line') then $
       coords_arr = atv_wcs2pix(coords_arr, coord_sys='b1950') $
     else $
       coords_arr = atv_wcs2pix(coords_arr, coord_sys='b1950', /line)
   end
  'galactic': begin
     if (strlowcase(reg_type) ne 'line') then $
       coords_arr = atv_wcs2pix(coords_arr, coord_sys='galactic') $
     else $
       coords_arr = atv_wcs2pix(coords_arr, coord_sys='galactic', /line)
   end
  'ecliptic': begin
     if (strlowcase(reg_type) ne 'line') then $
       coords_arr = atv_wcs2pix(coords_arr, coord_sys='ecliptic') $
     else $
       coords_arr = atv_wcs2pix(coords_arr, coord_sys='ecliptic', /line)
   end
  else: 
  endcase
endif else begin

  if (strpos(coords_arr[0], ':')) ne -1 then begin

; Check that a WCS region is not overplotted on image with no WCS
    if (NOT ptr_valid(state.astr_ptr)) then begin
      atv_message, 'WCS Regions cannot be displayed on image without WCS', $
        msgtype='error', /window
      return
    endif

    if (strlowcase(reg_type) ne 'line') then $
      coords_arr = atv_wcs2pix(coords_arr,coord_sys='current') $
    else $
      coords_arr = atv_wcs2pix(coords_arr,coord_sys='current', /line)
  endif else begin
    if (strlowcase(reg_type) ne 'line') then $
      coords_arr = atv_wcs2pix(coords_arr,coord_sys='pixel') $
    else $
      coords_arr = atv_wcs2pix(coords_arr,coord_sys='pixel', /line)
  endelse

endelse

  CASE strlowcase(color_str) OF

  'red':     (*(plot_ptr[iplot])).options.color = '1'
  'black':   (*(plot_ptr[iplot])).options.color = '0'
  'green':   (*(plot_ptr[iplot])).options.color = '2'
  'blue':    (*(plot_ptr[iplot])).options.color = '3'
  'cyan':    (*(plot_ptr[iplot])).options.color = '4'
  'magenta': (*(plot_ptr[iplot])).options.color = '5'
  'yellow':  (*(plot_ptr[iplot])).options.color = '6'
  'white':   (*(plot_ptr[iplot])).options.color = '7'
  ELSE:      (*(plot_ptr[iplot])).options.color = '1'

  ENDCASE

  atv_setwindow,state.draw_window_id
  atv_plotwindow  
  
  case strlowcase(reg_type) of

    'circle': begin
        xcenter = (float(coords_arr[0]) - state.offset[0] + 0.5) * $
                   state.zoom_factor
        ycenter = (float(coords_arr[1]) - state.offset[1] + 0.5) * $
                   state.zoom_factor

        radius = float(coords_arr[2]) * state.zoom_factor
        tvcircle, radius, xcenter, ycenter, /device, $
          _extra = (*(plot_ptr[iplot])).options

        if (text_str ne '') then xyouts, xcenter, ycenter, text_str, $
          alignment=0.5, _extra = (*(plot_ptr[iplot])).options, /device
    end
    'box': begin
        angle = 0 ; initialize angle to 0
        if (n_coords ge 4) then begin
          xcenter = (float(coords_arr[0]) - state.offset[0] + 0.5) * $
                     state.zoom_factor
          ycenter = (float(coords_arr[1]) - state.offset[1] + 0.5) * $
                     state.zoom_factor
          xwidth = float(coords_arr[2]) * state.zoom_factor
          ywidth = float(coords_arr[3]) * state.zoom_factor
          if (n_coords ge 5) then angle = float(coords_arr[4])
        endif
        width_arr = [xwidth,ywidth]  
        ; angle = -angle because tvbox rotates clockwise
        tvbox, width_arr, xcenter, ycenter, angle=-angle, $
          _extra = (*(plot_ptr[iplot])).options

        if (text_str ne '') then xyouts, xcenter, ycenter, text_str, $
          alignment=0.5, _extra = (*(plot_ptr[iplot])).options, /device
    end
    'ellipse': begin
        angle = 0 ; initialize angle to 0
        if (n_coords ge 4) then begin
          xcenter = (float(coords_arr[0]) - state.offset[0] + 0.5) * $
                     state.zoom_factor
          ycenter = (float(coords_arr[1]) - state.offset[1] + 0.5) * $
                     state.zoom_factor
          xradius = float(coords_arr[2]) * state.zoom_factor
          yradius = float(coords_arr[3]) * state.zoom_factor
          if (n_coords ge 5) then angle = float(coords_arr[4])
        endif

       ; Correct angle for default orientation used by tvellipse
        angle=angle+180.

        if (xcenter ge 0.0 and ycenter ge 0.0) then $
          tvellipse, xradius, yradius, xcenter, ycenter, angle, $
            _extra = (*(plot_ptr[iplot])).options

        if (text_str ne '') then xyouts, xcenter, ycenter, text_str, $
          alignment=0.5, _extra = (*(plot_ptr[iplot])).options, /device
    end
    'polygon': begin
       n_vert = n_elements(coords_arr) / 2
       xpoints = fltarr(n_vert)
       ypoints = fltarr(n_vert)
       for vert_i = 0, n_vert - 1 do begin
         xpoints[vert_i] = coords_arr[vert_i*2]
         ypoints[vert_i] = coords_arr[vert_i*2+1]
       endfor

       if (xpoints[0] ne xpoints[n_vert-1] OR $
           ypoints[0] ne ypoints[n_vert-1]) then begin
         xpoints1 = fltarr(n_vert+1)
         ypoints1 = fltarr(n_vert+1)
         xpoints1[0:n_vert-1] = xpoints
         ypoints1[0:n_vert-1] = ypoints
         xpoints1[n_vert] = xpoints[0]
         ypoints1[n_vert] = ypoints[0]
         xpoints = xpoints1
         ypoints = ypoints1
       endif

       xcenter = total(xpoints) / n_elements(xpoints)
       ycenter = total(ypoints) / n_elements(ypoints)

       plots, xpoints, ypoints,  $
           _extra = (*(plot_ptr[iplot])).options         

       if (text_str ne '') then xyouts, xcenter, ycenter, text_str, $
         alignment=0.5, _extra = (*(plot_ptr[iplot])).options, /device
    end
    'line': begin
        x1 = (float(coords_arr[0]) - state.offset[0] + 0.5) * $
                    state.zoom_factor
        y1 = (float(coords_arr[1]) - state.offset[1] + 0.5) * $
                    state.zoom_factor
        x2 = (float(coords_arr[2]) - state.offset[0] + 0.5) * $
                    state.zoom_factor
        y2 = (float(coords_arr[3]) - state.offset[1] + 0.5) * $
                    state.zoom_factor

        xpoints = [x1,x2]
        ypoints = [y1,y2]
        xcenter = total(xpoints) / n_elements(xpoints)
        ycenter = total(ypoints) / n_elements(ypoints)

        plots, xpoints, ypoints, /device, $
          _extra = (*(plot_ptr[iplot])).options

        if (text_str ne '') then xyouts, xcenter, ycenter, text_str, $
          alignment=0.5, _extra = (*(plot_ptr[iplot])).options, /device
    end
    else: begin

    end

    endcase

endfor

atv_resetwindow
state.newrefresh=1
end

;----------------------------------------------------------------------

pro atv_plot1wcsgrid, iplot
common atv_pdata
common atv_state
common atv_images

wcslabelcolor = (*(plot_ptr[iplot])).options.wcslabelcolor
gridcolor = (*(plot_ptr[iplot])).options.gridcolor
charsize = (*(plot_ptr[iplot])).options.charsize
charthick = (*(plot_ptr[iplot])).options.charthick

if (NOT ptr_valid(state.astr_ptr)) then begin
  atverase, 1
  atv_message, 'Cannot Display WCS Grid On Image Without WCS Coordinates.', $
      msgtype = 'error', /window  
  return
endif

headtype = strmid((*state.astr_ptr).ctype[0], 0, 4)

if (state.wcstype EQ 'angle') then begin

; Create local header variable to use for WCS grid.
  hdr = *(state.head_ptr)
 
; need numerical equinox values
  IF (state.equinox EQ 'J2000') THEN num_equinox = 2000.0 ELSE $
    IF (state.equinox EQ 'B1950') THEN num_equinox = 1950.0 ELSE $
    num_equinox = float(state.equinox)

  IF (state.display_equinox EQ 'J2000') THEN num_disp_equinox = 2000.0 ELSE $
    IF (state.display_equinox EQ 'B1950') THEN num_disp_equinox = 1950.0 ELSE $
    num_disp_equinox = float(state.equinox)

; Add EQUINOX to hdr if it does not exist in order to precess to
; display equinox
  year = GET_EQUINOX(hdr, code)    ;YEAR of hdr equinox
  IF code EQ -1 THEN $     
    sxaddpar, hdr, 'EQUINOX', num_equinox

  IF (num_equinox ne 2000.0) THEN $ 
    hprecess, hdr, 2000.0

; Now convert the hdr variable to the display coordinate system
  CASE state.display_coord_sys OF 
    'RA--': heuler, hdr, /celestial
    'GLON': heuler, hdr, /galactic
    'ELON': heuler, hdr, /ecliptic
    ELSE:
  ENDCASE  

; Now precess header to display equinox
  IF (num_equinox ne num_disp_equinox) THEN $
    hprecess, hdr, num_disp_equinox

; Extract an astrometry structure from hdr variable
  extast, hdr, astr

; Now operate on hdr variable to find grid coordinates, spacing, etc.

  x=findgen(n_elements(main_image[*,0]))
  nx = n_elements(main_image[*,0])
  y=findgen(n_elements(main_image[0,*]))
  ny = n_elements(main_image[0,*])
  
  x_bottom = x
  y_bottom = 0.
  x_top = x
  y_top = ny - 1
  x_left = 0.
  y_left = y
  x_right = nx - 1
  y_right = y

  xy2ad, x_bottom, y_bottom, astr, lon_bottom, lat_bottom
  xy2ad, x_top, y_top, astr, lon_top, lat_top
  xy2ad, x_left, y_left, astr, lon_left, lat_left
  xy2ad, x_right, y_right, astr, lon_right, lat_right

; Now create min/max lon/lat arrays
  lon_min = min([lon_bottom,lon_top,lon_left,lon_right])
  lon_max = max([lon_bottom,lon_top,lon_left,lon_right])
  lat_min = min([lat_bottom,lat_top,lat_left,lat_right])
  lat_max = max([lat_bottom,lat_top,lat_left,lat_right])

; Search for the poles for currently displayed coordinate system.
; Get positions of North and South Poles and check if in the image.  
  ad2xy, 0., 90., astr, x_npole, y_npole 
  ad2xy, 0., -90., astr, x_spole, y_spole

  north_diff = abs(90. - lat_max)
  south_diff = abs(-90. - lat_max)

  if (x_npole gt 0. and x_npole lt x_right and $
      y_npole gt 0. and y_npole lt y_top and $
      north_diff lt south_diff) then lat_max = 90.

  if (x_spole gt 0. and x_spole lt x_right and $
      y_spole gt 0. and y_spole lt y_top and $
      north_diff gt south_diff) then lat_min = -90.

; Adjust deltalon, lon_min, lon_max when Meridian in image.

  IF (round(lon_min) eq 0 AND round(lon_max) eq 360 AND $
      lat_min ne -90. AND lat_max ne 90.) THEN BEGIN

    ind_bottom = where(lon_bottom gt 180. AND lon_bottom lt 360.,count_bottom)
    ind_top = where(lon_top gt 180. AND lon_top lt 360.,count_top)
    ind_left = where(lon_left gt 180. AND lon_left lt 360.,count_left)
    ind_right = where(lon_right gt 180. AND lon_right lt 360.,count_right)
    if count_bottom ne 0 then $
      lon_bottom[ind_bottom] = lon_bottom[ind_bottom] - 360.
    if count_top ne 0 then $
      lon_top[ind_top] = lon_top[ind_top] - 360.
    if count_left ne 0 then $
      lon_left[ind_left] = lon_left[ind_left] - 360.
    if count_right ne 0 then $
      lon_right[ind_right] = lon_right[ind_right] - 360.
    lon_min = min([lon_bottom,lon_top,lon_left,lon_right])
    lon_max = max([lon_bottom,lon_top,lon_left,lon_right])

  ENDIF

  deltalon = lon_max - lon_min
  deltalat = lat_max - lat_min

  CASE (state.display_base60) OF
    0: BEGIN
       atv_find_spacing, deltalat, 5, 1, lat_spacing
       lat_tics = atv_set_spaces(lat_min, lat_spacing, 13, 0, lat_spaces)
       atv_find_spacing, deltalon, 5, 5, lon_spacing
       lon_tics = atv_set_spaces(lon_min, lon_spacing, 13, 0, lon_spaces)
    END
    1: BEGIN
       atv_find_spacing, deltalat, 5, 0, lat_spacing
       lat_tics = atv_set_spaces(lat_min, lat_spacing, 13, 0, lat_spaces)
       atv_find_spacing, deltalon, 5, 1, lon_spacing
       lon_tics = atv_set_spaces(lon_min, lon_spacing, 13, 0, lon_spaces)
    END
   ELSE:
  ENDCASE

; Make adjustments when Pole is in image
  if (lat_min eq -90. or lat_max eq 90.) then begin
    lon_spaces=[0.,30.,60.,90.,120.,150.,180.,210.,240.,270.,300.,330.,360.]

    tmp_index = where(lat_spaces gt 90., tmpcnt)
    if tmpcnt ne 0 then $
      lat_spaces[tmp_index] = ((90 - (lat_spaces[tmp_index] - 90.)))
  endif

  v0 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat0 = replicate(lat_spaces[0],n_elements(main_image[*,1]))

  v1 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat1 = replicate(lat_spaces[1],n_elements(main_image[*,1]))

  v2 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat2 = replicate(lat_spaces[2],n_elements(main_image[*,1]))

  v3 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat3 = replicate(lat_spaces[3],n_elements(main_image[*,1]))

  v4 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat4 = replicate(lat_spaces[4],n_elements(main_image[*,1]))

  v5 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat5 = replicate(lat_spaces[5],n_elements(main_image[*,1]))

  v6 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat6 = replicate(lat_spaces[6],n_elements(main_image[*,1]))

  v7 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat7 = replicate(lat_spaces[7],n_elements(main_image[*,1]))

  v8 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat8 = replicate(lat_spaces[8],n_elements(main_image[*,1]))

  v9 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat9 = replicate(lat_spaces[9],n_elements(main_image[*,1]))

  v10 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat10 = replicate(lat_spaces[10],n_elements(main_image[*,1]))

  v11 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat11 = replicate(lat_spaces[11],n_elements(main_image[*,1]))

  v12 = findgen(n_elements(main_image[*,1])) * deltalon / $
       (n_elements(main_image[*,1])-1) + lon_min
  vlat12 = replicate(lat_spaces[12],n_elements(main_image[*,1]))



  vv0 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon0 = replicate(lon_spaces[0],n_elements(main_image[1,*]))

  vv1 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon1 = replicate(lon_spaces[1],n_elements(main_image[1,*]))
 
  vv2 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon2 = replicate(lon_spaces[2],n_elements(main_image[1,*]))

  vv3 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon3 = replicate(lon_spaces[3],n_elements(main_image[1,*]))

  vv4 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon4 = replicate(lon_spaces[4],n_elements(main_image[1,*]))

  vv5 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon5 = replicate(lon_spaces[5],n_elements(main_image[1,*]))

  vv6 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon6 = replicate(lon_spaces[6],n_elements(main_image[1,*]))

  vv7 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon7 = replicate(lon_spaces[7],n_elements(main_image[1,*]))

  vv8 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon8 = replicate(lon_spaces[8],n_elements(main_image[1,*]))

  vv9 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon9 = replicate(lon_spaces[9],n_elements(main_image[1,*]))

  vv10 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon10 = replicate(lon_spaces[10],n_elements(main_image[1,*]))

  vv11 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon11 = replicate(lon_spaces[11],n_elements(main_image[1,*]))

  vv12 = findgen(n_elements(main_image[1,*])) * deltalat / $
       (n_elements(main_image[1,*])-1) + lat_min
  vlon12 = replicate(lon_spaces[12],n_elements(main_image[1,*]))


; When Meridian in image, negative lon_spaces exist--add 360 
    tmp = where(lon_spaces lt 0., tmp_cnt)
    if tmp_cnt ne 0 then lon_spaces[tmp] = lon_spaces[tmp] + 360.

  ad2xy, v0, vlat0, astr, xlat0, ylat0 
  ad2xy, v1, vlat1, astr, xlat1, ylat1 
  ad2xy, v2, vlat2, astr, xlat2, ylat2 
  ad2xy, v3, vlat3, astr, xlat3, ylat3 
  ad2xy, v4, vlat4, astr, xlat4, ylat4 
  ad2xy, v5, vlat5, astr, xlat5, ylat5 
  ad2xy, v6, vlat6, astr, xlat6, ylat6 
  ad2xy, v7, vlat7, astr, xlat7, ylat7 
  ad2xy, v8, vlat8, astr, xlat8, ylat8 
  ad2xy, v9, vlat9, astr, xlat9, ylat9 
  ad2xy, v10, vlat10, astr, xlat10, ylat10 
  ad2xy, v11, vlat11, astr, xlat11, ylat11
  ad2xy, v12, vlat12, astr, xlat12, ylat12

  ad2xy, vlon0, vv0, astr,  xlon0, ylon0 
  ad2xy, vlon1, vv1, astr,  xlon1, ylon1 
  ad2xy, vlon2, vv2, astr,  xlon2, ylon2 
  ad2xy, vlon3, vv3, astr,  xlon3, ylon3 
  ad2xy, vlon4, vv4, astr,  xlon4, ylon4 
  ad2xy, vlon5, vv5, astr,  xlon5, ylon5 
  ad2xy, vlon6, vv6, astr,  xlon6, ylon6 
  ad2xy, vlon7, vv7, astr,  xlon7, ylon7 
  ad2xy, vlon8, vv8, astr,  xlon8, ylon8 
  ad2xy, vlon9, vv9, astr,  xlon9, ylon9 
  ad2xy, vlon10, vv10, astr,  xlon10, ylon10 
  ad2xy, vlon11, vv11, astr,  xlon11, ylon11
  ad2xy, vlon12, vv12, astr,  xlon12, ylon12

  ad2xy, lon_spaces[2], lat_spaces[0], astr, x_latline0, y_latline0
  ad2xy, lon_spaces[2], lat_spaces[1], astr, x_latline1, y_latline1
  ad2xy, lon_spaces[2], lat_spaces[2], astr, x_latline2, y_latline2
  ad2xy, lon_spaces[2], lat_spaces[3], astr, x_latline3, y_latline3
  ad2xy, lon_spaces[2], lat_spaces[4], astr, x_latline4, y_latline4
  ad2xy, lon_spaces[2], lat_spaces[5], astr, x_latline5, y_latline5
  ad2xy, lon_spaces[2], lat_spaces[6], astr, x_latline6, y_latline6
  ad2xy, lon_spaces[2], lat_spaces[7], astr, x_latline7, y_latline7
  ad2xy, lon_spaces[2], lat_spaces[8], astr, x_latline8, y_latline8
  ad2xy, lon_spaces[2], lat_spaces[9], astr, x_latline9, y_latline9
  ad2xy, lon_spaces[2], lat_spaces[10], astr, x_latline10, y_latline10
  ad2xy, lon_spaces[2], lat_spaces[11], astr, x_latline11, y_latline11
  ad2xy, lon_spaces[2], lat_spaces[12], astr, x_latline12, y_latline12

  ad2xy, lon_spaces[0], lat_spaces[2], astr, x_lonline0, y_lonline0
  ad2xy, lon_spaces[1], lat_spaces[2], astr, x_lonline1, y_lonline1
  ad2xy, lon_spaces[2], lat_spaces[2], astr, x_lonline2, y_lonline2
  ad2xy, lon_spaces[3], lat_spaces[2], astr, x_lonline3, y_lonline3
  ad2xy, lon_spaces[4], lat_spaces[2], astr, x_lonline4, y_lonline4
  ad2xy, lon_spaces[5], lat_spaces[2], astr, x_lonline5, y_lonline5
  ad2xy, lon_spaces[6], lat_spaces[2], astr, x_lonline6, y_lonline6
  ad2xy, lon_spaces[7], lat_spaces[2], astr, x_lonline7, y_lonline7
  ad2xy, lon_spaces[8], lat_spaces[2], astr, x_lonline8, y_lonline8
  ad2xy, lon_spaces[9], lat_spaces[2], astr, x_lonline9, y_lonline9
  ad2xy, lon_spaces[10], lat_spaces[2], astr, x_lonline10, y_lonline10
  ad2xy, lon_spaces[11], lat_spaces[2], astr, x_lonline11, y_lonline11
  ad2xy, lon_spaces[12], lat_spaces[2], astr, x_lonline12, y_lonline12

; Determine orientation for labels

  xlat2_diff = abs(xlat2 - x_latline2[0])
  ylat2_diff = abs(ylat2 - y_latline2[0])
  xlat2_diff_index = where(xlat2_diff eq min(xlat2_diff))
  ylat2_diff_index = where(ylat2_diff eq min(ylat2_diff))
  x1 = xlat2[xlat2_diff_index] 
  y1 = ylat2[xlat2_diff_index]
  x2 = xlat2[xlat2_diff_index+2]
  y2 = ylat2[xlat2_diff_index+2]
  deltax = x1 - x2
  deltay = y1 - y2
  dy_dx = deltay / deltax
  latlabel_orientation = (180./!dpi * atan(dy_dx))

  xlon2_diff = abs(xlon2 - x_lonline2[0])
  ylon2_diff = abs(ylon2 - y_lonline2[0])
  xlon2_diff_index = where(xlon2_diff eq min(xlon2_diff))
  ylon2_diff_index = where(ylon2_diff eq min(ylon2_diff))
  x1 = xlon2[xlon2_diff_index] 
  y1 = ylon2[xlon2_diff_index]
  x2 = xlon2[xlon2_diff_index+2]
  y2 = ylon2[xlon2_diff_index+2]
  deltax = x1 - x2
  deltay = y1 - y2
  dy_dx = deltay / deltax
  lonlabel_orientation = (180./!dpi * atan(dy_dx))

; Check for label orientations of -90. where divide by 0 occurs
if (latlabel_orientation[0] eq 0.0) then lonlabel_orientation = -90.
if (lonlabel_orientation[0] eq 0.0) then latlabel_orientation = -90.

  index_vlat0 = where(xlat0 ge 0.0 AND xlat0 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat0 ge 0.0 AND ylat0 le $
                      (n_elements(main_image[0,*])-1),count_vlat0)

  index_vlat1 = where(xlat1 ge 0.0 AND xlat1 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat1 ge 0.0 AND ylat1 le $
                      (n_elements(main_image[0,*])-1), count_vlat1)

  index_vlat2 = where(xlat2 ge 0.0 AND xlat2 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat2 ge 0.0 AND ylat2 le $
                      (n_elements(main_image[0,*])-1), count_vlat2)

  index_vlat3 = where(xlat3 ge 0.0 AND xlat3 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat3 ge 0.0 AND ylat3 le $
                      (n_elements(main_image[0,*])-1), count_vlat3)

  index_vlat4 = where(xlat4 ge 0.0 AND xlat4 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat4 ge 0.0 AND ylat4 le $
                      (n_elements(main_image[0,*])-1), count_vlat4)

  index_vlat5 = where(xlat5 ge 0.0 AND xlat5 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat5 ge 0.0 AND ylat5 le $
                      (n_elements(main_image[0,*])-1), count_vlat5)

  index_vlat6 = where(xlat6 ge 0.0 AND xlat6 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat6 ge 0.0 AND ylat6 le $
                      (n_elements(main_image[0,*])-1), count_vlat6)

  index_vlat7 = where(xlat7 ge 0.0 AND xlat7 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat7 ge 0.0 AND ylat7 le $
                      (n_elements(main_image[0,*])-1), count_vlat7)

  index_vlat8 = where(xlat8 ge 0.0 AND xlat8 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat8 ge 0.0 AND ylat8 le $
                      (n_elements(main_image[0,*])-1), count_vlat8)

  index_vlat9 = where(xlat9 ge 0.0 AND xlat9 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat9 ge 0.0 AND ylat9 le $
                      (n_elements(main_image[0,*])-1), count_vlat9)

  index_vlat10 = where(xlat10 ge 0.0 AND xlat10 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat10 ge 0.0 AND ylat10 le $
                      (n_elements(main_image[0,*])-1), count_vlat10)

  index_vlat11 = where(xlat11 ge 0.0 AND xlat11 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat11 ge 0.0 AND ylat11 le $
                      (n_elements(main_image[0,*])-1), count_vlat11)

  index_vlat12 = where(xlat12 ge 0.0 AND xlat12 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylat12 ge 0.0 AND ylat12 le $
                      (n_elements(main_image[0,*])-1), count_vlat12)


  count_vlat = [count_vlat0, count_vlat1, count_vlat2, count_vlat3, $
                count_vlat4, count_vlat5, count_vlat6, count_vlat7, $
                count_vlat8, count_vlat9, count_vlat10, count_vlat11, $
                count_vlat12]


  index_vlon0 = where(xlon0 ge 0.0 AND xlon0 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon0 ge 0.0 AND ylon0 le $
                      (n_elements(main_image[0,*])-1),count_vlon0)

  index_vlon1 = where(xlon1 ge 0.0 AND xlon1 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon1 ge 0.0 AND ylon1 le $
                      (n_elements(main_image[0,*])-1), count_vlon1)

  index_vlon2 = where(xlon2 ge 0.0 AND xlon2 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon2 ge 0.0 AND ylon2 le $
                      (n_elements(main_image[0,*])-1), count_vlon2)

  index_vlon3 = where(xlon3 ge 0.0 AND xlon3 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon3 ge 0.0 AND ylon3 le $
                      (n_elements(main_image[0,*])-1), count_vlon3)

  index_vlon4 = where(xlon4 ge 0.0 AND xlon4 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon4 ge 0.0 AND ylon4 le $
                      (n_elements(main_image[0,*])-1), count_vlon4)

  index_vlon5 = where(xlon5 ge 0.0 AND xlon5 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon5 ge 0.0 AND ylon5 le $
                      (n_elements(main_image[0,*])-1), count_vlon5)

  index_vlon6 = where(xlon6 ge 0.0 AND xlon6 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon6 ge 0.0 AND ylon6 le $
                      (n_elements(main_image[0,*])-1), count_vlon6)

  index_vlon7 = where(xlon7 ge 0.0 AND xlon7 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon7 ge 0.0 AND ylon7 le $
                      (n_elements(main_image[0,*])-1), count_vlon7)

  index_vlon8 = where(xlon8 ge 0.0 AND xlon8 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon8 ge 0.0 AND ylon8 le $
                      (n_elements(main_image[0,*])-1), count_vlon8)

  index_vlon9 = where(xlon9 ge 0.0 AND xlon9 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon9 ge 0.0 AND ylon9 le $
                      (n_elements(main_image[0,*])-1), count_vlon9)

  index_vlon10 = where(xlon10 ge 0.0 AND xlon10 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon10 ge 0.0 AND ylon10 le $
                      (n_elements(main_image[0,*])-1), count_vlon10)

  index_vlon11 = where(xlon11 ge 0.0 AND xlon11 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon11 ge 0.0 AND ylon11 le $
                      (n_elements(main_image[0,*])-1), count_vlon11)

  index_vlon12 = where(xlon12 ge 0.0 AND xlon12 le $
                      (n_elements(main_image[*,0])-1) AND $
                      ylon12 ge 0.0 AND ylon12 le $
                      (n_elements(main_image[0,*])-1), count_vlon12)

  count_vlon = [count_vlon0, count_vlon1, count_vlon2, count_vlon3, $
                count_vlon4, count_vlon5, count_vlon6, count_vlon7, $
                count_vlon8, count_vlon9, count_vlon10, count_vlon11, $
                count_vlon12]


  atv_setwindow, state.draw_window_id


  if count_vlat0 ne 0. then plots,xlat0[index_vlat0],ylat0[index_vlat0], $
     color=gridcolor
  if count_vlat1 ne 0. then plots,xlat1[index_vlat1],ylat1[index_vlat1], $
     color=gridcolor
  if count_vlat2 ne 0. then plots,xlat2[index_vlat2],ylat2[index_vlat2], $
     color=gridcolor
  if count_vlat3 ne 0. then plots,xlat3[index_vlat3],ylat3[index_vlat3], $
     color=gridcolor
  if count_vlat4 ne 0. then plots,xlat4[index_vlat4],ylat4[index_vlat4], $
     color=gridcolor
  if count_vlat5 ne 0. then plots,xlat5[index_vlat5],ylat5[index_vlat5], $
     color=gridcolor
  if count_vlat6 ne 0. then plots,xlat6[index_vlat6],ylat6[index_vlat6], $
     color=gridcolor
  if count_vlat7 ne 0. then plots,xlat7[index_vlat7],ylat7[index_vlat7], $
     color=gridcolor
  if count_vlat8 ne 0. then plots,xlat8[index_vlat8],ylat8[index_vlat8], $
     color=gridcolor
  if count_vlat9 ne 0. then plots,xlat9[index_vlat9],ylat9[index_vlat9], $
     color=gridcolor
  if count_vlat10 ne 0. then plots,xlat10[index_vlat10],ylat10[index_vlat10], $
     color=gridcolor
  if count_vlat11 ne 0. then plots,xlat11[index_vlat11],ylat11[index_vlat11], $
     color=gridcolor
  if count_vlat12 ne 0. then plots,xlat12[index_vlat12],ylat12[index_vlat12], $
     color=gridcolor

  if count_vlon0 ne 0. then plots,xlon0[index_vlon0],ylon0[index_vlon0], $
     color=gridcolor
  if count_vlon1 ne 0. then plots,xlon1[index_vlon1],ylon1[index_vlon1], $
     color=gridcolor
  if count_vlon2 ne 0. then plots,xlon2[index_vlon2],ylon2[index_vlon2], $
     color=gridcolor
  if count_vlon3 ne 0. then plots,xlon3[index_vlon3],ylon3[index_vlon3], $
     color=gridcolor
  if count_vlon4 ne 0. then plots,xlon4[index_vlon4],ylon4[index_vlon4], $
     color=gridcolor
  if count_vlon5 ne 0. then plots,xlon5[index_vlon5],ylon5[index_vlon5], $
     color=gridcolor
  if count_vlon6 ne 0. then plots,xlon6[index_vlon6],ylon6[index_vlon6], $
     color=gridcolor
  if count_vlon7 ne 0. then plots,xlon7[index_vlon7],ylon7[index_vlon7], $
     color=gridcolor
  if count_vlon8 ne 0. then plots,xlon8[index_vlon8],ylon8[index_vlon8], $
     color=gridcolor
  if count_vlon9 ne 0. then plots,xlon9[index_vlon9],ylon9[index_vlon9], $
     color=gridcolor
  if count_vlon10 ne 0. then plots,xlon10[index_vlon10],ylon10[index_vlon10], $
     color=gridcolor
  if count_vlon11 ne 0. then plots,xlon11[index_vlon11],ylon11[index_vlon11], $
     color=gridcolor
  if count_vlon12 ne 0. then plots,xlon12[index_vlon12],ylon12[index_vlon12], $
     color=gridcolor

; Create label strings for different coordinate systems

  CASE (state.display_coord_sys) OF

    'RA--': BEGIN

      IF (state.display_base60 eq 1) THEN BEGIN

        lon0_arr = sixty(lon_spaces[0]/15.)
        lon1_arr = sixty(lon_spaces[1]/15.)
        lon2_arr = sixty(lon_spaces[2]/15.)
        lon3_arr = sixty(lon_spaces[3]/15.)
        lon4_arr = sixty(lon_spaces[4]/15.)
        lon5_arr = sixty(lon_spaces[5]/15.)
        lon6_arr = sixty(lon_spaces[6]/15.)
        lon7_arr = sixty(lon_spaces[7]/15.)
        lon8_arr = sixty(lon_spaces[8]/15.)
        lon9_arr = sixty(lon_spaces[9]/15.)
        lon10_arr = sixty(lon_spaces[10]/15.)
        lon11_arr = sixty(lon_spaces[11]/15.)
        lon12_arr = sixty(lon_spaces[12]/15.)
        lat0_arr = sixty(lat_spaces[0])
        lat1_arr = sixty(lat_spaces[1])
        lat2_arr = sixty(lat_spaces[2])
        lat3_arr = sixty(lat_spaces[3])
        lat4_arr = sixty(lat_spaces[4])
        lat5_arr = sixty(lat_spaces[5])
        lat6_arr = sixty(lat_spaces[6])
        lat7_arr = sixty(lat_spaces[7])
        lat8_arr = sixty(lat_spaces[8])
        lat9_arr = sixty(lat_spaces[9])
        lat10_arr = sixty(lat_spaces[10])
        lat11_arr = sixty(lat_spaces[11])
        lat12_arr = sixty(lat_spaces[12])

        lon0_hh=strcompress(string(fix(lon0_arr[0])),/remove_all)
        lon0_mm=strcompress(string(fix(lon0_arr[1])),/remove_all)
        lon0_ss=strcompress(string(fix(lon0_arr[2])),/remove_all)
        if (strlen(lon0_hh) lt 2) then lon0_hh = '0' + lon0_hh
        if (strlen(lon0_mm) lt 2) then lon0_mm = '0' + lon0_mm
        if (strlen(lon0_ss) lt 2) then lon0_ss = '0' + lon0_ss
        lon0_str = lon0_hh + ':' + lon0_mm + ':' + lon0_ss

        lon1_hh=strcompress(string(fix(lon1_arr[0])),/remove_all)
        lon1_mm=strcompress(string(fix(lon1_arr[1])),/remove_all)
        lon1_ss=strcompress(string(fix(lon1_arr[2])),/remove_all)
        if (strlen(lon1_hh) lt 2) then lon1_hh = '0' + lon1_hh
        if (strlen(lon1_mm) lt 2) then lon1_mm = '0' + lon1_mm
        if (strlen(lon1_ss) lt 2) then lon1_ss = '0' + lon1_ss
        lon1_str = lon1_hh + ':' + lon1_mm + ':' + lon1_ss

        lon2_hh=strcompress(string(fix(lon2_arr[0])),/remove_all)
        lon2_mm=strcompress(string(fix(lon2_arr[1])),/remove_all)
        lon2_ss=strcompress(string(fix(lon2_arr[2])),/remove_all)
        if (strlen(lon2_hh) lt 2) then lon2_hh = '0' + lon2_hh
        if (strlen(lon2_mm) lt 2) then lon2_mm = '0' + lon2_mm
        if (strlen(lon2_ss) lt 2) then lon2_ss = '0' + lon2_ss
        lon2_str = lon2_hh + ':' + lon2_mm + ':' + lon2_ss

        lon3_hh=strcompress(string(fix(lon3_arr[0])),/remove_all)
        lon3_mm=strcompress(string(fix(lon3_arr[1])),/remove_all)
        lon3_ss=strcompress(string(fix(lon3_arr[2])),/remove_all)
        if (strlen(lon3_hh) lt 2) then lon3_hh = '0' + lon3_hh
        if (strlen(lon3_mm) lt 2) then lon3_mm = '0' + lon3_mm
        if (strlen(lon3_ss) lt 2) then lon3_ss = '0' + lon3_ss
        lon3_str = lon3_hh + ':' + lon3_mm + ':' + lon3_ss

        lon4_hh=strcompress(string(fix(lon4_arr[0])),/remove_all)
        lon4_mm=strcompress(string(fix(lon4_arr[1])),/remove_all)
        lon4_ss=strcompress(string(fix(lon4_arr[2])),/remove_all)
        if (strlen(lon4_hh) lt 2) then lon4_hh = '0' + lon4_hh
        if (strlen(lon4_mm) lt 2) then lon4_mm = '0' + lon4_mm
        if (strlen(lon4_ss) lt 2) then lon4_ss = '0' + lon4_ss
        lon4_str = lon4_hh + ':' + lon4_mm + ':' + lon4_ss

        lon5_hh=strcompress(string(fix(lon5_arr[0])),/remove_all)
        lon5_mm=strcompress(string(fix(lon5_arr[1])),/remove_all)
        lon5_ss=strcompress(string(fix(lon5_arr[2])),/remove_all)
        if (strlen(lon5_hh) lt 2) then lon5_hh = '0' + lon5_hh
        if (strlen(lon5_mm) lt 2) then lon5_mm = '0' + lon5_mm
        if (strlen(lon5_ss) lt 2) then lon5_ss = '0' + lon5_ss
        lon5_str = lon5_hh + ':' + lon5_mm + ':' + lon5_ss

        lon6_hh=strcompress(string(fix(lon6_arr[0])),/remove_all)
        lon6_mm=strcompress(string(fix(lon6_arr[1])),/remove_all)
        lon6_ss=strcompress(string(fix(lon6_arr[2])),/remove_all)
        if (strlen(lon6_hh) lt 2) then lon6_hh = '0' + lon6_hh
        if (strlen(lon6_mm) lt 2) then lon6_mm = '0' + lon6_mm
        if (strlen(lon6_ss) lt 2) then lon6_ss = '0' + lon6_ss
        lon6_str = lon6_hh + ':' + lon6_mm + ':' + lon6_ss

        lon7_hh=strcompress(string(fix(lon7_arr[0])),/remove_all)
        lon7_mm=strcompress(string(fix(lon7_arr[1])),/remove_all)
        lon7_ss=strcompress(string(fix(lon7_arr[2])),/remove_all)
        if (strlen(lon7_hh) lt 2) then lon7_hh = '0' + lon7_hh
        if (strlen(lon7_mm) lt 2) then lon7_mm = '0' + lon7_mm
        if (strlen(lon7_ss) lt 2) then lon7_ss = '0' + lon7_ss
        lon7_str = lon7_hh + ':' + lon7_mm + ':' + lon7_ss

        lon8_hh=strcompress(string(fix(lon8_arr[0])),/remove_all)
        lon8_mm=strcompress(string(fix(lon8_arr[1])),/remove_all)
        lon8_ss=strcompress(string(fix(lon8_arr[2])),/remove_all)
        if (strlen(lon8_hh) lt 2) then lon8_hh = '0' + lon8_hh
        if (strlen(lon8_mm) lt 2) then lon8_mm = '0' + lon8_mm
        if (strlen(lon8_ss) lt 2) then lon8_ss = '0' + lon8_ss
        lon8_str = lon8_hh + ':' + lon8_mm + ':' + lon8_ss

        lon9_hh=strcompress(string(fix(lon9_arr[0])),/remove_all)
        lon9_mm=strcompress(string(fix(lon9_arr[1])),/remove_all)
        lon9_ss=strcompress(string(fix(lon9_arr[2])),/remove_all)
        if (strlen(lon9_hh) lt 2) then lon9_hh = '0' + lon9_hh
        if (strlen(lon9_mm) lt 2) then lon9_mm = '0' + lon9_mm
        if (strlen(lon9_ss) lt 2) then lon9_ss = '0' + lon9_ss
        lon9_str = lon9_hh + ':' + lon9_mm + ':' + lon9_ss

        lon10_hh=strcompress(string(fix(lon10_arr[0])),/remove_all)
        lon10_mm=strcompress(string(fix(lon10_arr[1])),/remove_all)
        lon10_ss=strcompress(string(fix(lon10_arr[2])),/remove_all)
        if (strlen(lon10_hh) lt 2) then lon10_hh = '0' + lon10_hh
        if (strlen(lon10_mm) lt 2) then lon10_mm = '0' + lon10_mm
        if (strlen(lon10_ss) lt 2) then lon10_ss = '0' + lon10_ss
        lon10_str = lon10_hh + ':' + lon10_mm + ':' + lon10_ss

        lon11_hh=strcompress(string(fix(lon11_arr[0])),/remove_all)
        lon11_mm=strcompress(string(fix(lon11_arr[1])),/remove_all)
        lon11_ss=strcompress(string(fix(lon11_arr[2])),/remove_all)
        if (strlen(lon11_hh) lt 2) then lon11_hh = '0' + lon11_hh
        if (strlen(lon11_mm) lt 2) then lon11_mm = '0' + lon11_mm
        if (strlen(lon11_ss) lt 2) then lon11_ss = '0' + lon11_ss
        lon11_str = lon11_hh + ':' + lon11_mm + ':' + lon11_ss

        lon12_hh=strcompress(string(fix(lon12_arr[0])),/remove_all)
        lon12_mm=strcompress(string(fix(lon12_arr[1])),/remove_all)
        lon12_ss=strcompress(string(fix(lon12_arr[2])),/remove_all)
        if (strlen(lon12_hh) lt 2) then lon12_hh = '0' + lon12_hh
        if (strlen(lon12_mm) lt 2) then lon12_mm = '0' + lon12_mm
        if (strlen(lon12_ss) lt 2) then lon12_ss = '0' + lon12_ss
        lon12_str = lon12_hh + ':' + lon12_mm + ':' + lon12_ss


        lat0_dd=strcompress(string(fix(lat0_arr[0])),/remove_all)
        lat0_mm=strcompress(string(fix(lat0_arr[1])),/remove_all)
        lat0_ss=strcompress(string(fix(lat0_arr[2])),/remove_all)
        if (strlen(lat0_dd) lt 2) then lat0_dd = '0' + lat0_dd
        if (strmid(lat0_dd,0,1) eq '-' AND strlen(lat0_dd) lt 3) then begin
          strput, lat0_dd, '0', 0
          lat0_dd = '-' + lat0_dd
        endif
        if (strlen(lat0_mm) lt 2) then lat0_mm = '0' + lat0_mm
        if (strlen(lat0_ss) lt 2) then lat0_ss = '0' + lat0_ss
        lat0_str = lat0_dd + ':' + lat0_mm + ':' + lat0_ss

        lat1_dd=strcompress(string(fix(lat1_arr[0])),/remove_all)
        lat1_mm=strcompress(string(fix(lat1_arr[1])),/remove_all)
        lat1_ss=strcompress(string(fix(lat1_arr[2])),/remove_all)
        if (strlen(lat1_dd) lt 2) then lat1_dd = '0' + lat1_dd
        if (strmid(lat1_dd,0,1) eq '-' AND strlen(lat1_dd) lt 3) then begin
          strput, lat1_dd, '0', 0
          lat1_dd = '-' + lat1_dd
        endif
        if (strlen(lat1_mm) lt 2) then lat1_mm = '0' + lat1_mm
        if (strlen(lat1_ss) lt 2) then lat1_ss = '0' + lat1_ss
        lat1_str = lat1_dd + ':' + lat1_mm + ':' + lat1_ss

        lat2_dd=strcompress(string(fix(lat2_arr[0])),/remove_all)
        lat2_mm=strcompress(string(fix(lat2_arr[1])),/remove_all)
        lat2_ss=strcompress(string(fix(lat2_arr[2])),/remove_all)
        if (strlen(lat2_dd) lt 2) then lat2_dd = '0' + lat2_dd
        if (strmid(lat2_dd,0,1) eq '-' AND strlen(lat2_dd) lt 3) then begin
          strput, lat2_dd, '0', 0
          lat2_dd = '-' + lat2_dd
        endif
        if (strlen(lat2_mm) lt 2) then lat2_mm = '0' + lat2_mm
        if (strlen(lat2_ss) lt 2) then lat2_ss = '0' + lat2_ss
        lat2_str = lat2_dd + ':' + lat2_mm + ':' + lat2_ss

        lat3_dd=strcompress(string(fix(lat3_arr[0])),/remove_all)
        lat3_mm=strcompress(string(fix(lat3_arr[1])),/remove_all)
        lat3_ss=strcompress(string(fix(lat3_arr[2])),/remove_all)
        if (strlen(lat3_dd) lt 2) then lat3_dd = '0' + lat3_dd
        if (strmid(lat3_dd,0,1) eq '-' AND strlen(lat3_dd) lt 3) then begin
          strput, lat3_dd, '0', 0
          lat3_dd = '-' + lat3_dd
        endif
        if (strlen(lat3_mm) lt 2) then lat3_mm = '0' + lat3_mm
        if (strlen(lat3_ss) lt 2) then lat3_ss = '0' + lat3_ss
        lat3_str = lat3_dd + ':' + lat3_mm + ':' + lat3_ss

        lat4_dd=strcompress(string(fix(lat4_arr[0])),/remove_all)
        lat4_mm=strcompress(string(fix(lat4_arr[1])),/remove_all)
        lat4_ss=strcompress(string(fix(lat4_arr[2])),/remove_all)
        if (strlen(lat4_dd) lt 2) then lat4_dd = '0' + lat4_dd
        if (strmid(lat4_dd,0,1) eq '-' AND strlen(lat4_dd) lt 3) then begin
          strput, lat4_dd, '0', 0
          lat4_dd = '-' + lat4_dd
        endif
        if (strlen(lat4_mm) lt 2) then lat4_mm = '0' + lat4_mm
        if (strlen(lat4_ss) lt 2) then lat4_ss = '0' + lat4_ss
        lat4_str = lat4_dd + ':' + lat4_mm + ':' + lat4_ss

        lat5_dd=strcompress(string(fix(lat5_arr[0])),/remove_all)
        lat5_mm=strcompress(string(fix(lat5_arr[1])),/remove_all)
        lat5_ss=strcompress(string(fix(lat5_arr[2])),/remove_all)
        if (strlen(lat5_dd) lt 2) then lat5_dd = '0' + lat5_dd
        if (strmid(lat5_dd,0,1) eq '-' AND strlen(lat5_dd) lt 3) then begin
          strput, lat5_dd, '0', 0
          lat5_dd = '-' + lat5_dd
        endif
        if (strlen(lat5_mm) lt 2) then lat5_mm = '0' + lat5_mm
        if (strlen(lat5_ss) lt 2) then lat5_ss = '0' + lat5_ss
        lat5_str = lat5_dd + ':' + lat5_mm + ':' + lat5_ss

        lat6_dd=strcompress(string(fix(lat6_arr[0])),/remove_all)
        lat6_mm=strcompress(string(fix(lat6_arr[1])),/remove_all)
        lat6_ss=strcompress(string(fix(lat6_arr[2])),/remove_all)
        if (strlen(lat6_dd) lt 2) then lat6_dd = '0' + lat6_dd
        if (strmid(lat6_dd,0,1) eq '-' AND strlen(lat6_dd) lt 3) then begin
          strput, lat6_dd, '0', 0
          lat6_dd = '-' + lat6_dd
        endif
        if (strlen(lat6_mm) lt 2) then lat6_mm = '0' + lat6_mm
        if (strlen(lat6_ss) lt 2) then lat6_ss = '0' + lat6_ss
        lat6_str = lat6_dd + ':' + lat6_mm + ':' + lat6_ss

        lat7_dd=strcompress(string(fix(lat7_arr[0])),/remove_all)
        lat7_mm=strcompress(string(fix(lat7_arr[1])),/remove_all)
        lat7_ss=strcompress(string(fix(lat7_arr[2])),/remove_all)
        if (strlen(lat7_dd) lt 2) then lat7_dd = '0' + lat7_dd
        if (strmid(lat7_dd,0,1) eq '-' AND strlen(lat7_dd) lt 3) then begin
          strput, lat7_dd, '0', 0
          lat7_dd = '-' + lat7_dd
        endif
        if (strlen(lat7_mm) lt 2) then lat7_mm = '0' + lat7_mm
        if (strlen(lat7_ss) lt 2) then lat7_ss = '0' + lat7_ss
        lat7_str = lat7_dd + ':' + lat7_mm + ':' + lat7_ss

        lat8_dd=strcompress(string(fix(lat8_arr[0])),/remove_all)
        lat8_mm=strcompress(string(fix(lat8_arr[1])),/remove_all)
        lat8_ss=strcompress(string(fix(lat8_arr[2])),/remove_all)
        if (strlen(lat8_dd) lt 2) then lat8_dd = '0' + lat8_dd
        if (strmid(lat8_dd,0,1) eq '-' AND strlen(lat8_dd) lt 3) then begin
          strput, lat8_dd, '0', 0
          lat8_dd = '-' + lat8_dd
        endif
        if (strlen(lat8_mm) lt 2) then lat8_mm = '0' + lat8_mm
        if (strlen(lat8_ss) lt 2) then lat8_ss = '0' + lat8_ss
        lat8_str = lat8_dd + ':' + lat8_mm + ':' + lat8_ss

        lat9_dd=strcompress(string(fix(lat9_arr[0])),/remove_all)
        lat9_mm=strcompress(string(fix(lat9_arr[1])),/remove_all)
        lat9_ss=strcompress(string(fix(lat9_arr[2])),/remove_all)
        if (strlen(lat9_dd) lt 2) then lat9_dd = '0' + lat0_dd
        if (strmid(lat9_dd,0,1) eq '-' AND strlen(lat9_dd) lt 3) then begin
          strput, lat9_dd, '0', 0
          lat9_dd = '-' + lat9_dd
        endif
        if (strlen(lat9_mm) lt 2) then lat9_mm = '0' + lat9_mm
        if (strlen(lat9_ss) lt 2) then lat9_ss = '0' + lat9_ss
        lat9_str = lat9_dd + ':' + lat9_mm + ':' + lat9_ss

        lat10_dd=strcompress(string(fix(lat10_arr[0])),/remove_all)
        lat10_mm=strcompress(string(fix(lat10_arr[1])),/remove_all)
        lat10_ss=strcompress(string(fix(lat10_arr[2])),/remove_all)
        if (strlen(lat10_dd) lt 2) then lat10_dd = '0' + lat10_dd
        if (strmid(lat10_dd,0,1) eq '-' AND strlen(lat10_dd) lt 3) then begin
          strput, lat10_dd, '0', 0
          lat10_dd = '-' + lat10_dd
        endif
        if (strlen(lat10_mm) lt 2) then lat10_mm = '0' + lat10_mm
        if (strlen(lat10_ss) lt 2) then lat10_ss = '0' + lat10_ss
        lat10_str = lat10_dd + ':' + lat10_mm + ':' + lat10_ss

        lat11_dd=strcompress(string(fix(lat11_arr[0])),/remove_all)
        lat11_mm=strcompress(string(fix(lat11_arr[1])),/remove_all)
        lat11_ss=strcompress(string(fix(lat11_arr[2])),/remove_all)
        if (strlen(lat11_dd) lt 2) then lat11_dd = '0' + lat11_dd
        if (strmid(lat11_dd,0,1) eq '-' AND strlen(lat11_dd) lt 3) then begin
          strput, lat11_dd, '0', 0
          lat11_dd = '-' + lat11_dd
        endif
        if (strlen(lat11_mm) lt 2) then lat11_mm = '0' + lat11_mm
        if (strlen(lat11_ss) lt 2) then lat11_ss = '0' + lat11_ss
        lat11_str = lat11_dd + ':' + lat11_mm + ':' + lat11_ss

        lat12_dd=strcompress(string(fix(lat12_arr[0])),/remove_all)
        lat12_mm=strcompress(string(fix(lat12_arr[1])),/remove_all)
        lat12_ss=strcompress(string(fix(lat12_arr[2])),/remove_all)
        if (strlen(lat12_dd) lt 2) then lat12_dd = '0' + lat12_dd
        if (strmid(lat12_dd,0,1) eq '-' AND strlen(lat12_dd) lt 3) then begin
          strput, lat12_dd, '0', 0
          lat12_dd = '-' + lat12_dd
        endif
        if (strlen(lat12_mm) lt 2) then lat12_mm = '0' + lat12_mm
        if (strlen(lat12_ss) lt 2) then lat12_ss = '0' + lat12_ss
        lat12_str = lat12_dd + ':' + lat12_mm + ':' + lat12_ss

      ENDIF ELSE BEGIN

        lon0_str = strcompress(string(lon_spaces[0]),/remove_all)
        lon1_str = strcompress(string(lon_spaces[1]),/remove_all)
        lon2_str = strcompress(string(lon_spaces[2]),/remove_all)
        lon3_str = strcompress(string(lon_spaces[3]),/remove_all)
        lon4_str = strcompress(string(lon_spaces[4]),/remove_all)
        lon5_str = strcompress(string(lon_spaces[5]),/remove_all)
        lon6_str = strcompress(string(lon_spaces[6]),/remove_all)
        lon7_str = strcompress(string(lon_spaces[7]),/remove_all)
        lon8_str = strcompress(string(lon_spaces[8]),/remove_all)
        lon9_str = strcompress(string(lon_spaces[9]),/remove_all)
        lon10_str = strcompress(string(lon_spaces[10]),/remove_all)
        lon11_str = strcompress(string(lon_spaces[11]),/remove_all)
        lon12_str = strcompress(string(lon_spaces[12]),/remove_all)
        lat0_str = strcompress(string(lat_spaces[0]),/remove_all)
        lat1_str = strcompress(string(lat_spaces[1]),/remove_all)
        lat2_str = strcompress(string(lat_spaces[2]),/remove_all)
        lat3_str = strcompress(string(lat_spaces[3]),/remove_all)
        lat4_str = strcompress(string(lat_spaces[4]),/remove_all)
        lat5_str = strcompress(string(lat_spaces[5]),/remove_all)
        lat6_str = strcompress(string(lat_spaces[6]),/remove_all)
        lat7_str = strcompress(string(lat_spaces[7]),/remove_all)
        lat8_str = strcompress(string(lat_spaces[8]),/remove_all)
        lat9_str = strcompress(string(lat_spaces[9]),/remove_all)
        lat10_str = strcompress(string(lat_spaces[10]),/remove_all)
        lat11_str = strcompress(string(lat_spaces[11]),/remove_all)
        lat12_str = strcompress(string(lat_spaces[12]),/remove_all)

        pos_lon0 = strpos(lon0_str, '.')
        lon0_str = strmid(lon0_str, 0, pos_lon0+4)
        pos_lon1 = strpos(lon1_str, '.')
        lon1_str = strmid(lon1_str, 0, pos_lon1+4)
        pos_lon2 = strpos(lon2_str, '.')
        lon2_str = strmid(lon2_str, 0, pos_lon2+4)
        pos_lon3 = strpos(lon3_str, '.')
        lon3_str = strmid(lon3_str, 0, pos_lon3+4)
        pos_lon4 = strpos(lon4_str, '.')
        lon4_str = strmid(lon4_str, 0, pos_lon4+4)
        pos_lon5 = strpos(lon5_str, '.')
        lon5_str = strmid(lon5_str, 0, pos_lon5+4)
        pos_lon6 = strpos(lon6_str, '.')
        lon6_str = strmid(lon6_str, 0, pos_lon6+4)
        pos_lon7 = strpos(lon7_str, '.')
        lon7_str = strmid(lon7_str, 0, pos_lon7+4)
        pos_lon8 = strpos(lon8_str, '.')
        lon8_str = strmid(lon8_str, 0, pos_lon8+4)
        pos_lon9 = strpos(lon9_str, '.')
        lon9_str = strmid(lon9_str, 0, pos_lon9+4)
        pos_lon10 = strpos(lon10_str, '.')
        lon10_str = strmid(lon10_str, 0, pos_lon10+4)
        pos_lon11 = strpos(lon11_str, '.')
        lon11_str = strmid(lon11_str, 0, pos_lon11+4)
        pos_lon12 = strpos(lon12_str, '.')
        lon12_str = strmid(lon12_str, 0, pos_lon12+4)

        pos_lat0 = strpos(lat0_str, '.')
        lat0_str = strmid(lat0_str, 0, pos_lat0+4)
        pos_lat1 = strpos(lat1_str, '.')
        lat1_str = strmid(lat1_str, 0, pos_lat1+4)
        pos_lat2 = strpos(lat2_str, '.')
        lat2_str = strmid(lat2_str, 0, pos_lat2+4)
        pos_lat3 = strpos(lat3_str, '.')
        lat3_str = strmid(lat3_str, 0, pos_lat3+4)
        pos_lat4 = strpos(lat4_str, '.')
        lat4_str = strmid(lat4_str, 0, pos_lat4+4)
        pos_lat5 = strpos(lat5_str, '.')
        lat5_str = strmid(lat5_str, 0, pos_lat5+4)
        pos_lat6 = strpos(lat6_str, '.')
        lat6_str = strmid(lat6_str, 0, pos_lat6+4)
        pos_lat7 = strpos(lat7_str, '.')
        lat7_str = strmid(lat7_str, 0, pos_lat7+4)
        pos_lat8 = strpos(lat8_str, '.')
        lat8_str = strmid(lat8_str, 0, pos_lat8+4)
        pos_lat9 = strpos(lat9_str, '.')
        lat9_str = strmid(lat9_str, 0, pos_lat9+4)
        pos_lat10 = strpos(lat10_str, '.')
        lat10_str = strmid(lat10_str, 0, pos_lat10+4)
        pos_lat11 = strpos(lat11_str, '.')
        lat11_str = strmid(lat11_str, 0, pos_lat11+4)
        pos_lat12 = strpos(lat12_str, '.')
        lat12_str = strmid(lat12_str, 0, pos_lat12+4)

      ENDELSE

  END

      ELSE: BEGIN

        lon0_str = strcompress(string(lon_spaces[0]),/remove_all)
        lon1_str = strcompress(string(lon_spaces[1]),/remove_all)
        lon2_str = strcompress(string(lon_spaces[2]),/remove_all)
        lon3_str = strcompress(string(lon_spaces[3]),/remove_all)
        lon4_str = strcompress(string(lon_spaces[4]),/remove_all)
        lon5_str = strcompress(string(lon_spaces[5]),/remove_all)
        lon6_str = strcompress(string(lon_spaces[6]),/remove_all)
        lon7_str = strcompress(string(lon_spaces[7]),/remove_all)
        lon8_str = strcompress(string(lon_spaces[8]),/remove_all)
        lon9_str = strcompress(string(lon_spaces[9]),/remove_all)
        lon10_str = strcompress(string(lon_spaces[10]),/remove_all)
        lon11_str = strcompress(string(lon_spaces[11]),/remove_all)
        lon12_str = strcompress(string(lon_spaces[12]),/remove_all)
        lat0_str = strcompress(string(lat_spaces[0]),/remove_all)
        lat1_str = strcompress(string(lat_spaces[1]),/remove_all)
        lat2_str = strcompress(string(lat_spaces[2]),/remove_all)
        lat3_str = strcompress(string(lat_spaces[3]),/remove_all)
        lat4_str = strcompress(string(lat_spaces[4]),/remove_all)
        lat5_str = strcompress(string(lat_spaces[5]),/remove_all)
        lat6_str = strcompress(string(lat_spaces[6]),/remove_all)
        lat7_str = strcompress(string(lat_spaces[7]),/remove_all)
        lat8_str = strcompress(string(lat_spaces[8]),/remove_all)
        lat9_str = strcompress(string(lat_spaces[9]),/remove_all)
        lat10_str = strcompress(string(lat_spaces[10]),/remove_all)
        lat11_str = strcompress(string(lat_spaces[11]),/remove_all)
        lat12_str = strcompress(string(lat_spaces[12]),/remove_all)

        pos_lon0 = strpos(lon0_str, '.')
        lon0_str = strmid(lon0_str, 0, pos_lon0+4)
        pos_lon1 = strpos(lon1_str, '.')
        lon1_str = strmid(lon1_str, 0, pos_lon1+4)
        pos_lon2 = strpos(lon2_str, '.')
        lon2_str = strmid(lon2_str, 0, pos_lon2+4)
        pos_lon3 = strpos(lon3_str, '.')
        lon3_str = strmid(lon3_str, 0, pos_lon3+4)
        pos_lon4 = strpos(lon4_str, '.')
        lon4_str = strmid(lon4_str, 0, pos_lon4+4)
        pos_lon5 = strpos(lon5_str, '.')
        lon5_str = strmid(lon5_str, 0, pos_lon5+4)
        pos_lon6 = strpos(lon6_str, '.')
        lon6_str = strmid(lon6_str, 0, pos_lon6+4)
        pos_lon7 = strpos(lon7_str, '.')
        lon7_str = strmid(lon7_str, 0, pos_lon7+4)
        pos_lon8 = strpos(lon8_str, '.')
        lon8_str = strmid(lon8_str, 0, pos_lon8+4)
        pos_lon9 = strpos(lon9_str, '.')
        lon9_str = strmid(lon9_str, 0, pos_lon9+4)
        pos_lon10 = strpos(lon10_str, '.')
        lon10_str = strmid(lon10_str, 0, pos_lon10+4)
        pos_lon11 = strpos(lon11_str, '.')
        lon11_str = strmid(lon11_str, 0, pos_lon11+4)
        pos_lon12 = strpos(lon12_str, '.')
        lon12_str = strmid(lon12_str, 0, pos_lon12+4)

        pos_lat0 = strpos(lat0_str, '.')
        lat0_str = strmid(lat0_str, 0, pos_lat0+4)
        pos_lat1 = strpos(lat1_str, '.')
        lat1_str = strmid(lat1_str, 0, pos_lat1+4)
        pos_lat2 = strpos(lat2_str, '.')
        lat2_str = strmid(lat2_str, 0, pos_lat2+4)
        pos_lat3 = strpos(lat3_str, '.')
        lat3_str = strmid(lat3_str, 0, pos_lat3+4)
        pos_lat4 = strpos(lat4_str, '.')
        lat4_str = strmid(lat4_str, 0, pos_lat4+4)
        pos_lat5 = strpos(lat5_str, '.')
        lat5_str = strmid(lat5_str, 0, pos_lat5+4)
        pos_lat6 = strpos(lat6_str, '.')
        lat6_str = strmid(lat6_str, 0, pos_lat6+4)
        pos_lat7 = strpos(lat7_str, '.')
        lat7_str = strmid(lat7_str, 0, pos_lat7+4)
        pos_lat8 = strpos(lat8_str, '.')
        lat8_str = strmid(lat8_str, 0, pos_lat8+4)
        pos_lat9 = strpos(lat9_str, '.')
        lat9_str = strmid(lat9_str, 0, pos_lat9+4)
        pos_lat10 = strpos(lat10_str, '.')
        lat10_str = strmid(lat10_str, 0, pos_lat10+4)
        pos_lat11 = strpos(lat11_str, '.')
        lat11_str = strmid(lat11_str, 0, pos_lat11+4)
        pos_lat12 = strpos(lat12_str, '.')
        lat12_str = strmid(lat12_str, 0, pos_lat12+4)
    END

  ENDCASE


  x_lonline = [x_lonline0, x_lonline1, x_lonline2, x_lonline3, x_lonline4, $
               x_lonline5, x_lonline6, x_lonline7, x_lonline8, x_lonline9, $
               x_lonline10, x_lonline11, x_lonline12]
  y_lonline = [y_lonline0, y_lonline1, y_lonline2, y_lonline3, y_lonline4, $
               y_lonline5, y_lonline6, y_lonline7, y_lonline8, y_lonline9, $
               y_lonline10, y_lonline11, y_lonline12]
  lon_str = [lon0_str,lon1_str,lon2_str,lon3_str,lon4_str,lon5_str,lon6_str,$
             lon7_str,lon8_str,lon9_str,lon10_str,lon11_str,lon12_str]

  lon_str[2] = ''

; If pole in image, make last longitude string (360.000 or 24:00:00) blank
  IF (lat_min eq -90. OR lat_max eq 90.) THEN lon_str[12] = ''

  lonindex = where(count_vlon ne 0)

  xyouts, x_lonline[lonindex], y_lonline[lonindex], lon_str[lonindex], $
    charsize=charsize, orientation=lonlabel_orientation, alignment=0.5, $
    color=wcslabelcolor, charthick=charthick

  x_latline = [x_latline0, x_latline1, x_latline2, x_latline3, x_latline4, $
               x_latline5, x_latline6, x_latline7, x_latline8, x_latline9, $
               x_latline10, x_latline11, x_latline12]
  y_latline = [y_latline0, y_latline1, y_latline2, y_latline3, y_latline4, $
               y_latline5, y_latline6, y_latline7, y_latline8, y_latline9, $
               y_latline10, y_latline11, y_latline12]
  lat_str = [lat0_str,lat1_str,lat2_str,lat3_str,lat4_str,lat5_str,lat6_str,$
             lat7_str,lat8_str,lat9_str,lat10_str,lat11_str,lat12_str]

  latindex = where(count_vlat ne 0)

  xyouts, x_latline[latindex], y_latline[latindex], lat_str[latindex], $
    charsize=charsize, orientation=latlabel_orientation, alignment=0.5, $
    color=wcslabelcolor, charthick=charthick

endif

atv_resetwindow
state.newrefresh=1
end



pro atv_plot1contour, iplot
common atv_pdata
common atv_state

; Overplot contours on the image

atv_setwindow, state.draw_window_id
widget_control, /hourglass

xrange = !x.crange
yrange = !y.crange

; The following allows for 2 conditions, depending upon whether X and Y
; are set

dims = size( (*(plot_ptr[iplot])).z,/dim )

if (size( (*(plot_ptr[iplot])).x,/N_elements ) EQ dims[0] $
    AND size( (*(plot_ptr[iplot])).y,/N_elements) EQ dims[1] ) then begin
    
    contour, (*(plot_ptr[iplot])).z, (*(plot_ptr[iplot])).x, $
      (*(plot_ptr[iplot])).y, $
      position=[0,0,1,1], xrange=xrange, yrange=yrange, $
      xstyle=5, ystyle=5, /noerase, $
      _extra = (*(plot_ptr[iplot])).options
    
endif else begin
    
    contour, (*(plot_ptr[iplot])).z, $
      position=[0,0,1,1], xrange=xrange, yrange=yrange, $
      xstyle=5, ystyle=5, /noerase, $
      _extra = (*(plot_ptr[iplot])).options
          
endelse

atv_resetwindow
state.newrefresh=1
end

;---------------------------------------------------------------------

pro atv_plot1compass, iplot

; Uses idlastro routine arrows to plot compass arrows.

common atv_pdata
common atv_state

atv_setwindow, state.draw_window_id

widget_control, /hourglass

arrows, *(state.head_ptr), $
  (*(plot_ptr[iplot])).x, $
  (*(plot_ptr[iplot])).y, $
  thick = (*(plot_ptr[iplot])).thick, $
  charsize = (*(plot_ptr[iplot])).charsize, $
  arrowlen = (*(plot_ptr[iplot])).arrowlen, $
  color = (*(plot_ptr[iplot])).color, $
  notvertex = (*(plot_ptr[iplot])).notvertex, $
  /data

atv_resetwindow
state.newrefresh=1
end

;---------------------------------------------------------------------

pro atv_plot1scalebar, iplot

; uses modified version of idlastro routine arcbar to plot a scalebar

common atv_pdata
common atv_state

atv_setwindow, state.draw_window_id
widget_control, /hourglass

; routine arcbar doesn't recognize color=0, because it uses 
; keyword_set to check the color.  So we need to set !p.color = 0
; to get black if the user wants color=0

!p.color = 0

atv_arcbar, *(state.head_ptr), $
  (*(plot_ptr[iplot])).arclen, $
  position = (*(plot_ptr[iplot])).position, $
  thick = (*(plot_ptr[iplot])).thick, $
  size = (*(plot_ptr[iplot])).size, $
  color = (*(plot_ptr[iplot])).color, $
  seconds = (*(plot_ptr[iplot])).seconds, $
  /data

atv_resetwindow
state.newrefresh=1
end

;----------------------------------------------------------------------

pro atv_arcbar, hdr, arclen, LABEL = label, SIZE = size, THICK = thick, $
                DATA =data, COLOR = color, POSITION = position, $
                NORMAL = normal, SECONDS=SECONDS

common atv_state

; This is a copy of the IDL Astronomy User's Library routine 'arcbar',
; abbreviated for atv and modified to work with zoomed images.  For
; the revision history of the original arcbar routine, look at
; arcbar.pro in the pro/astro subdirectory of the IDL Astronomy User's
; Library.

; Modifications for atv:
; Modified to work with zoomed ATV images, AJB Jan. 2000 
; Moved text label upwards a bit for better results, AJB Jan. 2000

On_error,2                      ;Return to caller
 
extast, hdr, bastr, noparams    ;extract astrom params in deg.
 
if N_params() LT 2 then arclen = 1 ;default size = 1 arcmin

if not keyword_set( SIZE ) then size = 1.0
if not keyword_set( THICK ) then thick = !P.THICK
if not keyword_set( COLOR ) then color = !P.COLOR

a = bastr.crval[0]
d = bastr.crval[1]
if keyword_set(seconds) then factor = 3600.0d else factor = 60.0
d1 = d + (1/factor)             ;compute x,y of crval + 1 arcmin

proj = strmid(bastr.ctype[0],5,3)

case proj of 
    'GSS': gsssadxy, bastr, [a,a], [d,d1], x, y
    else:  ad2xy, [a,a], [d,d1], bastr, x, y 
endcase

dmin = sqrt( (x[1]-x[0])^2 + (y[1]-y[0])^2 ) ;det. size in pixels of 1 arcmin

if (!D.FLAGS AND 1) EQ 1 then begin ;Device have scalable pixels?
    if !X.s[1] NE 0 then begin
        dmin = convert_coord( dmin, 0, /DATA, /TO_DEVICE) - $ 
          convert_coord(    0, 0, /DATA, /TO_DEVICE) ;Fixed Apr 97
        dmin = dmin[0]
    endif else dmin = dmin/sxpar(hdr, 'NAXIS1' ) ;Fixed Oct. 96
endif else  dmin = dmin * state.zoom_factor    ; added by AJB Jan. '00

dmini2 = round(dmin * arclen)

if keyword_set(NORMAL) then begin
    posn = convert_coord(position,/NORMAL, /TO_DEVICE) 
    xi = posn[0] & yi = posn[1]
endif else if keyword_set(DATA) then begin
    posn = convert_coord(position,/DATA, /TO_DEVICE) 
    xi = posn[0] & yi = posn[1]
endif else begin
    xi = position[0]   & yi = position[1]
endelse         


xf = xi + dmini2
dmini3 = dmini2/10       ;Height of vertical end bars = total length/10.

plots,[xi,xf],[yi,yi], COLOR=color, /DEV, THICK=thick
plots,[xf,xf],[ yi+dmini3, yi-dmini3 ], COLOR=color, /DEV, THICK=thick
plots,[xi,xi],[ yi+dmini3, yi-dmini3 ], COLOR=color, /DEV, THICK=thick

if not keyword_set(Seconds) then begin
    if (!D.NAME EQ 'PS') and (!P.FONT EQ 0) then $ ;Postscript Font?
      arcsym='!9'+string(162B)+'!X' else arcsym = "'" 
endif else begin
    if (!D.NAME EQ 'PS') and (!P.FONT EQ 0) then $ ;Postscript Font?
      arcsym = '!9'+string(178B)+'!X' else arcsym = "''" 
endelse
if not keyword_set( LABEL) then begin
    if (arclen LT 1) then arcstr = string(arclen,format='(f4.2)') $
    else arcstr = string(arclen)
    label = strtrim(arcstr,2) + arcsym 
endif

; AJB modified this to move the numerical label upward a bit: 5/8/2000
xyouts,(xi+xf)/2, (yi+(dmini2/10)), label, SIZE = size, COLOR=color,$
  /DEV, alignment=.5, CHARTHICK=thick

return
end

;----------------------------------------------------------------------

pro atv_plotwindow
common atv_state

atv_setwindow, state.draw_window_id

; Set plot window

; improved version by N. Cunningham- different scaling for postscript
; vs non-postscript output  -- added 4/14/06
if !d.name eq 'PS' then begin
   xrange=[state.offset[0], $
           state.offset[0] + state.draw_window_size[0] $
           / state.zoom_factor] - 0.5
   yrange=[state.offset[1], $
           state.offset[1] + state.draw_window_size[1] $
           / state.zoom_factor] - 0.5
endif else begin
   xrange=[state.offset[0] + 0.5 / state.zoom_factor, $
           state.offset[0] + (state.draw_window_size[0] + 0.5) $
           / state.zoom_factor] - 0.5
   yrange=[state.offset[1] + 0.5 / state.zoom_factor, $
           state.offset[1] + (state.draw_window_size[1] + 0.5) $
           / state.zoom_factor] - 0.5
endelse

plot, [0], [0], /nodata, position=[0,0,1,1], $
 xrange=xrange, yrange=yrange, xstyle=5, ystyle=5, /noerase

atv_resetwindow
end

;----------------------------------------------------------------------

pro atv_plot1ellipse, rmax, rmin, xc, yc, pos_ang, _extra = _extra

; This is a modified version of Wayne Landsman's tvellipse, changed so
; that it won't ask for interactive input under any circumstances.

if N_params() LT 5 then pos_ang = 0. ;Default position angle

npoints = 500                   ;Number of points to connect
phi = 2*!pi*(findgen(npoints)/(npoints-1)) ;Divide circle into Npoints
ang = pos_ang/!RADEG            ;Position angle in radians
cosang = cos(ang)
sinang = sin(ang)

x =  rmax*cos(phi)              ;Parameterized equation of ellipse
y =  rmin*sin(phi)

xprime = xc + x*cosang - y*sinang ;Rotate to desired position angle
yprime = yc + x*sinang + y*cosang

plots, round(xprime), round(yprime), color=color, /device,  $
  _STRICT_Extra = _extra

end

;---------------------------------------------------------------------

pro atv_plotall
common atv_state
common atv_pdata

; Routine to overplot all line, text, and contour plots

if (nplot EQ 0) then return

atv_plotwindow

for iplot = 1, nplot do begin
    case (*(plot_ptr[iplot])).type of
        'points'  : atv_plot1plot, iplot
        'text'    : atv_plot1text, iplot
        'arrow'   : atv_plot1arrow, iplot
        'contour' : atv_plot1contour, iplot
        'compass' : atv_plot1compass, iplot
        'scalebar': atv_plot1scalebar, iplot
        'region'  : atv_plot1region, iplot
        'wcsgrid' : atv_plot1wcsgrid, iplot
        'polarization' : atv_plot1pol, iplot
        else      : print, 'Problem in atv_plotall!'   
    endcase
endfor

end

;----------------------------------------------------------------------

pro atvplot, x, y, _extra = options
common atv_pdata
common atv_state

; Routine to read in line plot data and options, store in a heap
; variable structure, and plot the line plot

if (not(xregistered('atv', /noshow))) then begin
    print, 'You need to start ATV first!'
    return
endif

if (N_params() LT 1) then begin
   print, 'Too few parameters for ATVPLOT.'
   return
endif

if (n_elements(options) EQ 0) then options = {color: 'red'}

if (nplot LT maxplot) then begin
   nplot = nplot + 1

;  convert color names to index numbers, and set default=red
   c = where(tag_names(options) EQ 'COLOR', count)
   if (count EQ 0) then options = create_struct(options, 'color', 'red')
   options.color = atv_icolor(options.color)

   pstruct = {type: 'points',   $     ; points
              x: x,             $     ; x coordinate
              y: y,             $     ; y coordinate
              options: options  $     ; plot keyword options
             }

   plot_ptr[nplot] = ptr_new(pstruct)

   atv_plotwindow
   atv_plot1plot, nplot

endif else begin
   print, 'Too many calls to ATVPLOT.'
endelse

end

;----------------------------------------------------------------------

pro atvxyouts, x, y, text, _extra = options
common atv_pdata
common atv_state

; Routine to read in text overplot string and options, store in a heap
; variable structure, and overplot the text

if (not(xregistered('atv', /noshow))) then begin
    print, 'You need to start ATV first!'
    return
endif

if (N_params() LT 3) then begin
   print, 'Too few parameters for ATVXYOUTS'
   return
endif

if (n_elements(options) EQ 0) then options = {color: 'red'}

if (nplot LT maxplot) then begin
   nplot = nplot + 1

;  convert color names to index numbers, and set default=red
   c = where(tag_names(options) EQ 'COLOR', count)
   if (count EQ 0) then options = create_struct(options, 'color', 'red')
   options.color = atv_icolor(options.color)

;  set default font to 1
   c = where(tag_names(options) EQ 'FONT', count)
   if (count EQ 0) then options = create_struct(options, 'font', 1)

   pstruct = {type: 'text',   $       ; type of plot 
              x: x,             $     ; x coordinate
              y: y,             $     ; y coordinate
              text: text,       $     ; text to plot
              options: options  $     ; plot keyword options
             }

   plot_ptr[nplot] = ptr_new(pstruct)

   atv_plotwindow
   atv_plot1text, nplot

endif else begin
   print, 'Too many calls to ATVPLOT.'
endelse

end

;----------------------------------------------------------------------

pro atvarrow, x1, y1, x2, y2, _extra = options
common atv_pdata
common atv_state

; Routine to read in arrow overplot options, store in a heap
; variable structure, and overplot the arrow

if (not(xregistered('atv', /noshow))) then begin
    print, 'You need to start ATV first!'
    return
endif

if (N_params() LT 4) then begin
   print, 'Too few parameters for ATVARROW'
   return
endif

if (n_elements(options) EQ 0) then options = {color: 'red'}

if (nplot LT maxplot) then begin
   nplot = nplot + 1

;  convert color names to index numbers, and set default=red
   c = where(tag_names(options) EQ 'COLOR', count)
   if (count EQ 0) then options = create_struct(options, 'color', 'red')
   options.color = atv_icolor(options.color)

   pstruct = {type: 'arrow',   $       ; type of plot 
              x1: x1,             $     ; x1 coordinate
              y1: y1,             $     ; y1 coordinate
              x2: x2,             $     ; x2 coordinate
              y2: y2,             $     ; y2 coordinate     
              options: options  $     ; plot keyword options
             }

   plot_ptr[nplot] = ptr_new(pstruct)

   atv_plotwindow
   atv_plot1arrow, nplot

endif else begin
   print, 'Too many calls to ATVPLOT.'
endelse

end

;----------------------------------------------------------------------

pro atvregionfile, region_file
common atv_state
common atv_pdata

; Routine to read in region filename, store in a heap variable
; structure, and overplot the regions

if (not(xregistered('atv', /noshow))) then begin
    print, 'You need to start ATV first!'
    return
endif

if (nplot LT maxplot) then begin
   nplot = nplot + 1

options = {color: 'green'}
options.color = atv_icolor(options.color)

readfmt, region_file, 'a200', reg_array, /silent

pstruct = {type:'region', $            ; type of plot
           reg_array: reg_array, $     ; array of regions to plot
           options: options $          ; plot keyword options
          }

plot_ptr[nplot] = ptr_new(pstruct)

atv_plotwindow
atv_plot1region, nplot

endif else begin
  print, 'Too many calls to ATVPLOT.'
endelse

end

;----------------------------------------------------------------------

pro atv_setregion_event, event

; Event handler for atv_setregion.  Region plot structure created from
; information in form widget.  Plotting routine atv_plot1region is
; then called.

common atv_state
common atv_pdata

  CASE event.tag OF
  
  'REG_OPT' : BEGIN
       CASE event.value OF
         '0' : BEGIN
             widget_control,(*state.reg_ids_ptr)[3],Sensitive=1 
             widget_control,(*state.reg_ids_ptr)[4],Sensitive=1
             widget_control,(*state.reg_ids_ptr)[5],Sensitive=1
             widget_control,(*state.reg_ids_ptr)[6],Sensitive=1
             widget_control,(*state.reg_ids_ptr)[7],Sensitive=0         
             widget_control,(*state.reg_ids_ptr)[8],Sensitive=0
             widget_control,(*state.reg_ids_ptr)[9],Sensitive=0         
             widget_control,(*state.reg_ids_ptr)[10],Sensitive=0  
             widget_control,(*state.reg_ids_ptr)[11],Sensitive=0
         END
         '1' : BEGIN
             widget_control,(*state.reg_ids_ptr)[3],Sensitive=1 
             widget_control,(*state.reg_ids_ptr)[4],Sensitive=1
             widget_control,(*state.reg_ids_ptr)[5],Sensitive=1
             widget_control,(*state.reg_ids_ptr)[6],Sensitive=1
             widget_control,(*state.reg_ids_ptr)[7],Sensitive=0         
             widget_control,(*state.reg_ids_ptr)[8],Sensitive=0
             widget_control,(*state.reg_ids_ptr)[9],Sensitive=0         
             widget_control,(*state.reg_ids_ptr)[10],Sensitive=0
             widget_control,(*state.reg_ids_ptr)[11],Sensitive=1         
         END
         '2' : BEGIN
             widget_control,(*state.reg_ids_ptr)[3],Sensitive=1 
             widget_control,(*state.reg_ids_ptr)[4],Sensitive=1
             widget_control,(*state.reg_ids_ptr)[5],Sensitive=1
             widget_control,(*state.reg_ids_ptr)[6],Sensitive=1
             widget_control,(*state.reg_ids_ptr)[7],Sensitive=0
             widget_control,(*state.reg_ids_ptr)[8],Sensitive=0
             widget_control,(*state.reg_ids_ptr)[9],Sensitive=0         
             widget_control,(*state.reg_ids_ptr)[10],Sensitive=0
             widget_control,(*state.reg_ids_ptr)[11],Sensitive=1
         END
         '3' : BEGIN
             widget_control,(*state.reg_ids_ptr)[3],Sensitive=0 
             widget_control,(*state.reg_ids_ptr)[4],Sensitive=0
             widget_control,(*state.reg_ids_ptr)[5],Sensitive=0
             widget_control,(*state.reg_ids_ptr)[6],Sensitive=0
             widget_control,(*state.reg_ids_ptr)[7],Sensitive=1
             widget_control,(*state.reg_ids_ptr)[8],Sensitive=1
             widget_control,(*state.reg_ids_ptr)[9],Sensitive=1         
             widget_control,(*state.reg_ids_ptr)[10],Sensitive=1  
             widget_control,(*state.reg_ids_ptr)[11],Sensitive=0
         END
           ELSE:
     ENDCASE

 END

   'QUIT': BEGIN
         if ptr_valid(state.reg_ids_ptr) then ptr_free, state.reg_ids_ptr
         widget_control, event.top, /destroy
     END

   'DRAW': BEGIN
       IF (nplot LT maxplot) then begin

         nplot = nplot + 1

         reg_type = ['circle','box','ellipse','line']
         reg_color = ['red','black','green','blue','cyan','magenta', $
                      'yellow','white']
         coords_type = ['Pixel', 'J2000','B1950', $
                        'Galactic','Ecliptic', 'Native']
         reg_index = widget_info((*state.reg_ids_ptr)[0], /droplist_select)
         color_index = widget_info((*state.reg_ids_ptr)[1], /droplist_select)
         coords_index = widget_info((*state.reg_ids_ptr)[2], /droplist_select) 
         widget_control,(*state.reg_ids_ptr)[3],get_value=xcenter 
         widget_control,(*state.reg_ids_ptr)[4],get_value=ycenter           
         widget_control,(*state.reg_ids_ptr)[5],get_value=xwidth
         widget_control,(*state.reg_ids_ptr)[6],get_value=ywidth
         widget_control,(*state.reg_ids_ptr)[7],get_value=x1            
         widget_control,(*state.reg_ids_ptr)[8],get_value=y1
         widget_control,(*state.reg_ids_ptr)[9],get_value=x2       
         widget_control,(*state.reg_ids_ptr)[10],get_value=y2
         widget_control,(*state.reg_ids_ptr)[11],get_value=angle
         widget_control,(*state.reg_ids_ptr)[12],get_value=thick
         widget_control,(*state.reg_ids_ptr)[13],get_value=text_str
         text_str = strcompress(text_str[0],/remove_all)
  
         CASE reg_type[reg_index] OF 

         'circle': BEGIN
           region_str = reg_type[reg_index] + '(' + xcenter + ', ' + $
             ycenter + ', ' + xwidth  
           if (coords_index ne 0 and coords_index ne 5) then $
             region_str = region_str + ', ' + coords_type[coords_index]
           region_str = region_str + ') # color=' + reg_color[color_index]
         END

         'box': BEGIN
           region_str = reg_type[reg_index] + '(' + xcenter + ', ' + $
             ycenter + ', ' + xwidth + ', ' + ywidth + ', ' + angle 
           if (coords_index ne 0 and coords_index ne 5) then $
             region_str = region_str + ', ' + coords_type[coords_index]
           region_str = region_str + ') # color=' + reg_color[color_index]
         END

         'ellipse': BEGIN
           region_str = reg_type[reg_index] + '(' + xcenter + ', ' + $
             ycenter + ', ' + xwidth + ', ' + ywidth + ', ' + angle
           if (coords_index ne 0 and coords_index ne 5) then $
             region_str = region_str + ', ' + coords_type[coords_index]
           region_str = region_str + ') # color=' + reg_color[color_index]
         END

         'line': BEGIN
           region_str = reg_type[reg_index] + '(' + x1 + ', ' + y1 + ', ' + $
             x2 + ', ' + y2
           if (coords_index ne 0 and coords_index ne 5) then $
             region_str = region_str + ', ' + coords_type[coords_index]
           region_str = region_str + ') # color=' + reg_color[color_index]
         END

         ELSE: 
         ENDCASE

         if (text_str ne '') then region_str = region_str + $
            'text={' + text_str + '}'

         options = {color: reg_color[color_index], $
                    thick:thick}
         options.color = atv_icolor(options.color)

         pstruct = {type:'region', $          ;type of plot
                    reg_array:[region_str], $ ;region array to plot
                    options: options $
                    }

         plot_ptr[nplot] = ptr_new(pstruct)

         atv_plotwindow
         atv_plot1region, nplot

       ENDIF ELSE BEGIN
         print, 'Too many calls to ATVPLOT.'
       ENDELSE

;       if ptr_valid(state.reg_ids_ptr) then ptr_free, state.reg_ids_ptr
;       widget_control, event.top, /destroy

     END

  ELSE:
  ENDCASE

end

;----------------------------------------------------------------------

pro atv_wcsgrid, _extra = options

common atv_state
common atv_pdata

; Routine to read in wcs overplot options, store in a heap variable
; structure, and overplot the grid

if (not(xregistered('atv', /noshow))) then begin
    print, 'You need to start ATV first!'
    return
endif

if (nplot LT maxplot) then begin
   nplot = nplot + 1

; set default font to 1
  c = where(tag_names(options) EQ 'FONT', count)
  if (count EQ 0) then options = create_struct(options, 'font', 1)

pstruct = {type:'wcsgrid', $            ; type of plot
           options: options $           ; plot keyword options
          }

plot_ptr[nplot] = ptr_new(pstruct)

atv_plotwindow
atv_plot1wcsgrid, nplot

endif else begin
  print, 'Too many calls to ATVPLOT.'
endelse

end

;----------------------------------------------------------------------

pro atvcontour, z, x, y, _extra = options
common atv_pdata
common atv_state

; Routine to read in contour plot data and options, store in a heap
; variable structure, and overplot the contours.  Data to be contoured
; need not be the same dataset displayed in the atv window, but it
; should have the same x and y dimensions in order to align the
; overplot correctly.

if (not(xregistered('atv', /noshow))) then begin
    print, 'You need to start ATV first!'
    return
endif

if (N_params() LT 1) then begin
   print, 'Too few parameters for ATVCONTOUR.'
   return
endif

if (n_params() EQ 1 OR n_params() EQ 2) then begin
    x = 0
    y = 0
endif

if (n_elements(options) EQ 0) then options = {c_color: 'red'}

if (nplot LT maxplot) then begin
   nplot = nplot + 1

;  convert color names to index numbers, and set default=red
   c = where(tag_names(options) EQ 'C_COLOR', count)
   if (count EQ 0) then options = create_struct(options, 'c_color', 'red')
   options.c_color = atv_icolor(options.c_color)

   pstruct = {type: 'contour',  $     ; type of plot
              z: z,             $     ; z values
              x: x,             $     ; x coordinate
              y: y,             $     ; y coordinate
              options: options  $     ; plot keyword options
             }

   plot_ptr[nplot] = ptr_new(pstruct)

   atv_plotwindow
   atv_plot1contour, nplot

endif else begin
   print, 'Too many calls to ATVCONTOUR.'
endelse

end

;----------------------------------------------------------------------

pro atverase, nerase, norefresh = norefresh
common atv_pdata

; Routine to erase line plots from ATVPLOT, text from ATVXYOUTS,
; arrows from ATVARROW and contours from ATVCONTOUR.

if (n_params() LT 1) then begin
    nerase = nplot
endif else begin
    if (nerase GT nplot) then nerase = nplot
endelse

for iplot = nplot - nerase + 1, nplot do begin
    ptr_free, plot_ptr[iplot]
    plot_ptr[iplot] = ptr_new()
endfor

nplot = nplot - nerase

if (NOT keyword_set(norefresh)) then atv_refresh

end

;----------------------------------------------------------------------

pro atv_textlabel

; widget front end for atvxyouts

formdesc = ['0, text, , label_left=Text: , width=15', $
            '0, integer, 0, label_left=x: ', $
            '0, integer, 0, label_left=y: ', $
            '0, droplist, red|black|green|blue|cyan|magenta|yellow|white,label_left=Color:, set_value=0 ', $
            '0, float, 2.0, label_left=Charsize: ', $
            '0, integer, 1, label_left=Charthick: ', $
            '0, integer, 0, label_left=Orientation: ', $
            '1, base, , row', $
            '0, button, Cancel, quit', $
            '0, button, DrawText, quit']
            
textform = cw_form(formdesc, /column, $
                   title = 'atv text label')

if (textform.tag9 EQ 1) then begin
; switch red and black indices
    case textform.tag3 of
        0: labelcolor = 1
        1: labelcolor = 0
        else: labelcolor = textform.tag3
    endcase

    atvxyouts, textform.tag1, textform.tag2, textform.tag0, $
      color = labelcolor, charsize = textform.tag4, $
      charthick = textform.tag5, orientation = textform.tag6
endif

end

;---------------------------------------------------------------------

pro atv_setarrow

; widget front end for atvarrow

formdesc = ['0, integer, , label_left=Tail x: ', $
            '0, integer, , label_left=Tail y: ', $
            '0, integer, , label_left=Head x: ', $
            '0, integer, , label_left=Head y: ', $
            '0, droplist, red|black|green|blue|cyan|magenta|yellow|white,label_left=Color:, set_value=0 ', $
            '0, float, 1.0, label_left=LineThickness: ', $
            '0, float, 1.0, label_left=HeadThickness: ', $
            '1, base, , row', $
            '0, button, Cancel, quit', $
            '0, button, DrawArrow, quit']
            
textform = cw_form(formdesc, /column, $
                   title = 'atv arrow')

if (textform.tag9 EQ 1) then begin
; switch red and black indices
    case textform.tag4 of
        0: labelcolor = 1
        1: labelcolor = 0
        else: labelcolor = textform.tag4
    endcase

    atvarrow, textform.tag0, textform.tag1, $
      textform.tag2, textform.tag3, $
      color = labelcolor, thick = textform.tag5, $
      hthick = textform.tag6

endif

end

;---------------------------------------------------------------------

pro atv_regionfilelabel

; Routine to load region files into ATV

common atv_state
common atv_images

region_file = dialog_pickfile(/read, filter='*.reg')

;set up an array of strings

if (region_file ne '') then atvregionfile, region_file $
else return

end

;---------------------------------------------------------------------

pro atv_setregion

; Widget front-end for plotting individual regions on image

if (not(xregistered('atv_setregion', /noshow))) then begin
  common atv_state
  common atv_images

  regionbase = widget_base(/row)

  formdesc = ['0, droplist, circle|box|ellipse|line,label_left=Region:, set_value=0, TAG=reg_opt ', $
              '0, droplist, red|black|green|blue|cyan|magenta|yellow|white,label_left=Color:, set_value=0, TAG=color_opt ', $
              '0, droplist, Pixel|RA Dec (J2000)|RA Dec (B1950)|Galactic|Ecliptic|Native,label_left=Coords:, set_value=0, TAG=coord_opt ', $
              '0, text, 0, label_left=xcenter: , width=15', $
              '0, text, 0, label_left=ycenter: , width=15', $
              '0, text, 0, label_left=xwidth (Pix/ArcMin): , width=15', $
              '0, text, 0, label_left=ywidth (Pix/ArcMin): , width=15', $
              '0, text, 0, label_left=x1: , width=15', $
              '0, text, 0, label_left=y1: , width=15', $
              '0, text, 0, label_left=x2: , width=15', $
              '0, text, 0, label_left=y2: , width=15', $
              '0, text, 0.0, label_left=Angle: ', $
              '0, integer, 1, label_left=Thick: ', $
              '0, text,  , label_left=Text: ', $
              '1, base, , row', $
              '0, button, Done, quit, TAG=quit ', $
              '0, button, DrawRegion, quit, TAG=draw']

  regionform = cw_form(regionbase,formdesc, /column,title = 'atv region',$
                 IDS=reg_ids_ptr)

  widget_control, regionbase, /REALIZE

  reg_ids_ptr = reg_ids_ptr(where(widget_info(reg_ids_ptr,/type) eq 3 OR $
     widget_info(reg_ids_ptr,/type) eq 8))

  if ptr_valid(state.reg_ids_ptr) then ptr_free,state.reg_ids_ptr

  state.reg_ids_ptr = ptr_new(reg_ids_ptr)

  widget_control,(*state.reg_ids_ptr)[6],sensitive=0
  widget_control,(*state.reg_ids_ptr)[7],sensitive=0
  widget_control,(*state.reg_ids_ptr)[8],sensitive=0
  widget_control,(*state.reg_ids_ptr)[9],sensitive=0
  widget_control,(*state.reg_ids_ptr)[10],sensitive=0
  widget_control,(*state.reg_ids_ptr)[11],sensitive=0

; Check for WCS.  If WCS exists, then convert to display coordinates.
  if (ptr_valid(state.astr_ptr)) then begin
    ; Convert to display coordinates and change droplist selection.

    if (state.wcstype EQ 'angle') then begin
      xy2ad, state.coord[0], state.coord[1], *(state.astr_ptr), lon, lat

      wcsstring = atv_wcsstring(lon, lat, (*state.astr_ptr).ctype,  $
                                state.equinox, state.display_coord_sys, $
                                state.display_equinox, state.display_base60)

      if (strpos(wcsstring, 'J2000') ne -1) then coord_select = 1
      if (strpos(wcsstring, 'B1950') ne -1) then coord_select = 2
      if (strpos(wcsstring, 'Galactic') ne -1) then coord_select = 3
      if (strpos(wcsstring, 'Ecliptic') ne -1) then coord_select = 4

      if (strpos(wcsstring, 'J2000') eq -1 AND $
          strpos(wcsstring, 'B1950') eq -1 AND $
          strpos(wcsstring, 'Galactic') eq -1 AND $
          strpos(wcsstring, 'Ecliptic') eq -1) then coord_select = 5

      wcsstring = repstr(wcsstring,'J2000','')
      wcsstring = repstr(wcsstring,'B1950','')
      wcsstring = repstr(wcsstring,'Deg','')
      wcsstring = repstr(wcsstring,'Galactic','')
      wcsstring = repstr(wcsstring,'Ecliptic','')
      wcsstring = repstr(wcsstring,'(','')
      wcsstring = repstr(wcsstring,')','')

      xcent = strcompress(gettok(wcsstring,','), /remove_all)
      ycent = strcompress(wcsstring, /remove_all)

      widget_control,(*state.reg_ids_ptr)[3], Set_Value = xcent
      widget_control,(*state.reg_ids_ptr)[4], Set_Value = ycent
      widget_control,(*state.reg_ids_ptr)[7], Set_Value = xcent
      widget_control,(*state.reg_ids_ptr)[8], Set_Value = ycent
      widget_control,(*state.reg_ids_ptr)[2], set_droplist_select=coord_select
    endif    

  endif else begin
    widget_control,(*state.reg_ids_ptr)[3], Set_Value = $
      strcompress(string(state.coord[0]), /remove_all)
    widget_control,(*state.reg_ids_ptr)[4], Set_Value = $
      strcompress(string(state.coord[1]), /remove_all)
    widget_control,(*state.reg_ids_ptr)[7], Set_Value = $
      strcompress(string(state.coord[0]), /remove_all)
    widget_control,(*state.reg_ids_ptr)[8], Set_Value = $
      strcompress(string(state.coord[1]), /remove_all)
  endelse

  xmanager, 'atv_setregion', regionbase

endif else begin

  if (ptr_valid(state.astr_ptr)) then begin
    ; Convert to display coordinates and change droplist selection.

    if (state.wcstype EQ 'angle') then begin
      xy2ad, state.coord[0], state.coord[1], *(state.astr_ptr), lon, lat

      wcsstring = atv_wcsstring(lon, lat, (*state.astr_ptr).ctype,  $
                                state.equinox, state.display_coord_sys, $
                                state.display_equinox, state.display_base60)

      if (strpos(wcsstring, 'J2000') ne -1) then coord_select = 1
      if (strpos(wcsstring, 'B1950') ne -1) then coord_select = 2
      if (strpos(wcsstring, 'Galactic') ne -1) then coord_select = 3
      if (strpos(wcsstring, 'Ecliptic') ne -1) then coord_select = 4

      if (strpos(wcsstring, 'J2000') eq -1 AND $
          strpos(wcsstring, 'B1950') eq -1 AND $
          strpos(wcsstring, 'Galactic') eq -1 AND $
          strpos(wcsstring, 'Ecliptic') eq -1) then coord_select = 5

      wcsstring = repstr(wcsstring,'J2000','')
      wcsstring = repstr(wcsstring,'B1950','')
      wcsstring = repstr(wcsstring,'Deg','')
      wcsstring = repstr(wcsstring,'Galactic','')
      wcsstring = repstr(wcsstring,'Ecliptic','')
      wcsstring = repstr(wcsstring,'(','')
      wcsstring = repstr(wcsstring,')','')

      xcent = strcompress(gettok(wcsstring,','), /remove_all)
      ycent = strcompress(wcsstring, /remove_all)

      widget_control,(*state.reg_ids_ptr)[3], Set_Value = xcent
      widget_control,(*state.reg_ids_ptr)[4], Set_Value = ycent
      widget_control,(*state.reg_ids_ptr)[7], Set_Value = xcent
      widget_control,(*state.reg_ids_ptr)[8], Set_Value = ycent
      widget_control,(*state.reg_ids_ptr)[2], set_droplist_select=coord_select
    endif  

  endif else begin
    widget_control,(*state.reg_ids_ptr)[3], Set_Value = $
      strcompress(string(state.coord[0]), /remove_all)
    widget_control,(*state.reg_ids_ptr)[4], Set_Value = $
      strcompress(string(state.coord[1]), /remove_all)
    widget_control,(*state.reg_ids_ptr)[7], Set_Value = $
      strcompress(string(state.coord[0]), /remove_all)
    widget_control,(*state.reg_ids_ptr)[8], Set_Value = $
      strcompress(string(state.coord[1]), /remove_all)
  endelse

endelse

end

;---------------------------------------------------------------------

pro atv_wcsgridlabel

; Front-end widget for WCS labels

formdesc = ['0, droplist, red|black|green|blue|cyan|magenta|yellow|white,label_left=Grid Color:, set_value=7 ', $
            '0, droplist, red|black|green|blue|cyan|magenta|yellow|white,label_left=Label Color:, set_value=2 ', $ 
            '0, float, 1.0, label_left=Charsize: ', $
            '0, integer, 1, label_left=Charthick: ', $
            '1, base, , row', $
            '0, button, Cancel, quit', $
            '0, button, DrawGrid, quit']

gridform=cw_form(formdesc, /column, title = 'ATV WCS Grid')

gridcolor = gridform.tag0
wcslabelcolor = gridform.tag1

if (gridform.tag6 eq 1) then begin
; switch red and black indices
  case gridform.tag0 of 
    0: gridcolor = 1
    1: gridcolor = 0
    else: gridcolor = gridform.tag0
  endcase

  case gridform.tag1 of
    0: wcslabelcolor = 1
    1: wcslabelcolor = 0
    else: wcslabelcolor = gridform.tag1
  endcase

atv_wcsgrid, gridcolor=gridcolor, wcslabelcolor=wcslabelcolor, $
  charsize=gridform.tag2, charthick=gridform.tag3

endif

end

;---------------------------------------------------------------------

pro atv_oplotcontour

; widget front end for atvcontour

common atv_state
common atv_images

minvalstring = strcompress('0, float, ' + string(state.min_value) + $
                           ', label_left=MinValue: , width=15 ')
maxvalstring = strcompress('0, float, ' + string(state.max_value) + $
                           ', label_left=MaxValue: , width=15')

formdesc = ['0, droplist, red|black|green|blue|cyan|magenta|yellow|white,label_left=Color:, set_value=0 ', $
;            '0, float, 1.0, label_left=Charsize: ', $
;            '0, integer, 1, label_left=Charthick: ', $
            '0, droplist, solid|dotted|dashed|dashdot|dashdotdotdot|longdash, label_left=Linestyle: , set_value=0', $
            '0, integer, 1, label_left=LineThickness: ', $
            minvalstring, $
            maxvalstring, $
            '0, integer, 6, label_left=NLevels: ', $
            '1, base, , row,', $
            '0, button, Cancel, quit', $
            '0, button, DrawContour, quit']
            
cform = cw_form(formdesc, /column, $
                   title = 'atv text label')


if (cform.tag8 EQ 1) then begin
; switch red and black indices
    case cform.tag0 of
        0: labelcolor = 1
        1: labelcolor = 0
        else: labelcolor = cform.tag0
    endcase

    atvcontour, main_image, c_color = labelcolor, $
;      c_charsize = cform.tag1, c_charthick = cform.tag2, $
      c_linestyle = cform.tag1, $
      c_thick = cform.tag2, $
      min_value = cform.tag3, max_value = cform.tag4, $, 
      nlevels = cform.tag5
endif

end

;---------------------------------------------------------------------

pro atv_setcompass

; Routine to prompt user for compass parameters

common atv_state
common atv_images
common atv_pdata

if (nplot GE maxplot) then begin
    atv_message, 'Total allowed number of overplots exceeded.', $
      msgtype = 'error', /window
    return
endif
    

if (state.wcstype NE 'angle') then begin 
    atv_message, 'Cannot get coordinate info for this image!', $
      msgtype = 'error', /window
    return
endif

view_min = round(state.centerpix - $
        (0.5 * state.draw_window_size / state.zoom_factor)) 
view_max = round(view_min + state.draw_window_size / state.zoom_factor) - 1

xpos = string(round(view_min[0] + 0.15 * (view_max[0] - view_min[0])))
ypos = string(round(view_min[1] + 0.15 * (view_max[1] - view_min[1])))

xposstring = strcompress('0,integer,'+xpos+',label_left=XCenter: ')
yposstring = strcompress('0,integer,'+ypos+',label_left=YCenter: ')

formdesc = [ $
             xposstring, $
             yposstring, $
             '0, droplist, Vertex of Compass|Center of Compass, label_left = Coordinates Specify:, set_value=0', $
             '0, droplist, red|black|green|blue|cyan|magenta|yellow|white,label_left=Color:, set_value=0 ', $
             '0, integer, 1, label_left=LineThickness: ', $
             '0, float, 1, label_left=Charsize: ', $
             '0, float, 3.5, label_left=ArrowLength: ', $
             '1, base, , row,', $
             '0, button, Cancel, quit', $
             '0, button, DrawCompass, quit']
            
cform = cw_form(formdesc, /column, $
                   title = 'atv compass properties')

if (cform.tag8 EQ 1) then return

cform.tag0 = 0 > cform.tag0 < (state.image_size[0] - 1)
cform.tag1 = 0 > cform.tag1 < (state.image_size[1] - 1)

; switch red and black indices
case cform.tag3 of
    0: labelcolor = 1
    1: labelcolor = 0
    else: labelcolor = cform.tag3
endcase

pstruct = {type: 'compass',  $  ; type of plot
           x: cform.tag0,         $ 
           y: cform.tag1,         $
           notvertex: cform.tag2, $
           color: labelcolor, $
           thick: cform.tag4, $
           charsize: cform.tag5, $
           arrowlen: cform.tag6 $
          }

nplot = nplot + 1
plot_ptr[nplot] = ptr_new(pstruct)

atv_plotwindow
atv_plot1compass, nplot

end

;---------------------------------------------------------------------

pro atv_setscalebar

; Routine to prompt user for scalebar parameters

common atv_state
common atv_images
common atv_pdata

if (nplot GE maxplot) then begin
    atv_message, 'Total allowed number of overplots exceeded.', $
      msgtype = 'error', /window
    return
endif
    

if (state.wcstype NE 'angle') then begin 
    atv_message, 'Cannot get coordinate info for this image!', $
      msgtype = 'error', /window
    return
endif

view_min = round(state.centerpix - $
        (0.5 * state.draw_window_size / state.zoom_factor)) 
view_max = round(view_min + state.draw_window_size / state.zoom_factor) - 1

xpos = string(round(view_min[0] + 0.75 * (view_max[0] - view_min[0])))
ypos = string(round(view_min[1] + 0.15 * (view_max[1] - view_min[1])))

xposstring = strcompress('0,integer,'+xpos+',label_left=X (left end of bar): ')
yposstring = strcompress('0,integer,'+ypos+',label_left=Y (center of bar): ')

formdesc = [ $
             xposstring, $
             yposstring, $
             '0, float, 5.0, label_left=BarLength: ', $
             '0, droplist, arcsec|arcmin, label_left=Units:,set_value=0', $
             '0, droplist, red|black|green|blue|cyan|magenta|yellow|white,label_left=Color:, set_value=0 ', $
             '0, integer, 1, label_left=LineThickness: ', $
             '0, float, 1, label_left=Charsize: ', $
             '1, base, , row,', $
             '0, button, Cancel, quit', $
             '0, button, DrawScalebar, quit']
            
cform = cw_form(formdesc, /column, $
                   title = 'atv scalebar properties')

if (cform.tag8 EQ 1) then return

; switch red and black indices
case cform.tag4 of
    0: labelcolor = 1
    1: labelcolor = 0
    else: labelcolor = cform.tag4
endcase


cform.tag0 = 0 > cform.tag0 < (state.image_size[0] - 1)
cform.tag1 = 0 > cform.tag1 < (state.image_size[1] - 1)
cform.tag3 = abs(cform.tag3 - 1)  ; set default to be arcseconds

arclen = cform.tag2
if (float(round(arclen)) EQ arclen) then arclen = round(arclen)

pstruct = {type: 'scalebar',  $  ; type of plot
           arclen: arclen, $
           seconds: cform.tag3, $
           position: [cform.tag0,cform.tag1], $ 
           color: labelcolor, $
           thick: cform.tag5, $
           size: cform.tag6 $
          }

nplot = nplot + 1
plot_ptr[nplot] = ptr_new(pstruct)

atv_plotwindow
atv_plot1scalebar, nplot

end

;---------------------------------------------------------------------

pro atv_saveregion

; Save currently displayed regions to a file

common atv_state
common atv_pdata

reg_savefile = dialog_pickfile(file='atv.reg', filter='*.reg', /write) 

if (reg_savefile ne '') then begin 
  openw, lun, reg_savefile, /get_lun

  for iplot = 1, nplot do begin
    if ((*(plot_ptr[iplot])).type eq 'region') then begin
      n_regions = n_elements((*(plot_ptr[iplot])).reg_array)
      for n = 0, n_regions - 1 do begin
        printf, lun, strcompress((*(plot_ptr[iplot])).reg_array[n],/remove_all)
      endfor
    endif
  endfor

  close, lun
  free_lun, lun
endif else begin
  return
endelse

end

;---------------------------------------------------------------------
;          routines for drawing in the lineplot window
;---------------------------------------------------------------------

pro atv_lineplot_init

; This routine creates the window for line plots

common atv_state

state.lineplot_base_id = $
  widget_base(group_leader = state.base_id, $
              /row, $
              /base_align_right, $
              title = 'atv plot', $
              /tlb_size_events, $
              uvalue = 'lineplot_base')

state.lineplot_widget_id = $
  widget_draw(state.lineplot_base_id, $
              frame = 0, $
              scr_xsize = state.lineplot_size[0], $
              scr_ysize = state.lineplot_size[1], $
              uvalue = 'lineplot_window')

lbutton_base = $
  widget_base(state.lineplot_base_id, $
              /base_align_bottom, $
              /column, frame=2)

state.histbutton_base_id = $
  widget_base(lbutton_base, $
              /base_align_bottom, $
              /column, map=1)

state.x1_pix_id = $
    cw_field(state.histbutton_base_id, $
             /return_events, $
             /floating, $
             title = 'X1:', $
             uvalue = 'lineplot_newrange', $
             xsize = 12)

state.x2_pix_id = $
    cw_field(state.histbutton_base_id, $
             /return_events, $
             /floating, $
             title = 'X2:', $
             uvalue = 'lineplot_newrange', $
             xsize = 12)

state.y1_pix_id = $
    cw_field(state.histbutton_base_id, $
             /return_events, $
             /floating, $
             title = 'Y1:', $
             uvalue = 'lineplot_newrange', $
             xsize = 12)

state.y2_pix_id = $
    cw_field(state.histbutton_base_id, $
             /return_events, $
             /floating, $
             title = 'Y2:', $
             uvalue = 'lineplot_newrange', $
             xsize = 12)

state.histplot_binsize_id = $
    cw_field(state.histbutton_base_id, $
             /return_events, $
             /floating, $
             title = 'Bin:', $
             uvalue = 'lineplot_newrange', $
             xsize = 12)

state.lineplot_xmin_id = $
  cw_field(lbutton_base, $
           /return_events, $
           /floating, $
           title = 'XMin:', $
           uvalue = 'lineplot_newrange', $
           xsize = 12)

state.lineplot_xmax_id = $
  cw_field(lbutton_base, $
           /return_events, $
           /floating, $
           title = 'XMax:', $
           uvalue = 'lineplot_newrange', $
           xsize = 12)

state.lineplot_ymin_id = $
  cw_field(lbutton_base, $
           /return_events, $
           /floating, $
           title = 'YMin:', $
           uvalue = 'lineplot_newrange', $
           xsize = 12)

state.lineplot_ymax_id = $
  cw_field(lbutton_base, $
           /return_events, $
           /floating, $
           title = 'YMax:', $
           uvalue = 'lineplot_newrange', $
           xsize = 12)


state.holdrange_base_id = $
  widget_base(lbutton_base, $
              row = 1, $
              /nonexclusive, frame=1)

state.holdrange_button_id = $
  widget_button(state.holdrange_base_id, $
                value = 'Hold Ranges', $
                uvalue = 'lineplot_holdrange')

lineplot_fullrange = $
  widget_button(lbutton_base, $
                value = 'FullRange', $
                uvalue = 'lineplot_fullrange')

lineplot_ps = $
  widget_button(lbutton_base, $
                value = 'Create PS', $
                uvalue = 'lineplot_ps')

lineplot_done = $
  widget_button(lbutton_base, $
                value = 'Done', $
                uvalue = 'lineplot_done')

widget_control, state.lineplot_base_id, /realize
widget_control, state.holdrange_button_id, set_button=state.holdrange_value

widget_control, state.lineplot_widget_id, get_value = tmp_value
state.lineplot_window_id = tmp_value

lbuttgeom = widget_info(lbutton_base, /geometry)
state.lineplot_min_size[1] = lbuttgeom.ysize

basegeom = widget_info(state.lineplot_base_id, /geometry)
drawgeom = widget_info(state.lineplot_widget_id, /geometry)

state.lineplot_pad[0] = basegeom.xsize - drawgeom.xsize
state.lineplot_pad[1] = basegeom.ysize - drawgeom.ysize
    
xmanager, 'atv_lineplot', state.lineplot_base_id, /no_block

atv_resetwindow
end

;--------------------------------------------------------------------

pro atv_rowplot, ps=ps, fullrange=fullrange, newcoord=newcoord

; draws a new row plot in the plot window or to postscript output

common atv_state
common atv_images

if (keyword_set(ps)) then begin
    thick = 3
    color = 0
endif else begin
    thick = 1
    color = 7
endelse

if (keyword_set(newcoord)) then state.plot_coord = state.coord

if (not (keyword_set(ps))) then begin
    newplot = 0
    if (not (xregistered('atv_lineplot', /noshow))) then begin
        atv_lineplot_init
        newplot = 1
    endif 

    widget_control, state.histbutton_base_id, map=0
    widget_control, state.holdrange_button_id, sensitive=1

    widget_control, state.lineplot_xmin_id, get_value=xmin
    widget_control, state.lineplot_xmax_id, get_value=xmax
    widget_control, state.lineplot_ymin_id, get_value=ymin
    widget_control, state.lineplot_ymax_id, get_value=ymax

    if (newplot EQ 1 OR state.plot_type NE 'rowplot' OR $
        keyword_set(fullrange) OR $
        (state.holdrange_value EQ 0 AND keyword_set(newcoord))) then begin
        xmin = 0.0
        xmax = state.image_size[0]
        ymin = min(main_image[*,state.plot_coord[1]])
        ymax = max(main_image[*,state.plot_coord[1]]) 
    endif
   
    widget_control, state.lineplot_xmin_id, set_value=xmin
    widget_control, state.lineplot_xmax_id, set_value=xmax
    widget_control, state.lineplot_ymin_id, set_value=ymin
    widget_control, state.lineplot_ymax_id, set_value=ymax

    state.lineplot_xmin = xmin
    state.lineplot_xmax = xmax
    state.lineplot_ymin = ymin
    state.lineplot_ymax = ymax

    state.plot_type = 'rowplot'
    atv_setwindow, state.lineplot_window_id
    erase
    
endif

plot, main_image[*, state.plot_coord[1]], $
  xst = 3, yst = 3, psym = 10, $
  title = strcompress('Plot of row ' + $
                      string(state.plot_coord[1])), $

  xtitle = 'Column', $
  ytitle = 'Pixel Value', $
  color = color, xmargin=[15,3], $
  xran = [state.lineplot_xmin, state.lineplot_xmax], $
  yran = [state.lineplot_ymin, state.lineplot_ymax], $
  thick = thick, xthick = thick, ythick = thick, charthick = thick


if (not (keyword_set(ps))) then begin 
  widget_control, state.lineplot_base_id, /clear_events
  atv_resetwindow
endif

end


;--------------------------------------------------------------------

pro atv_colplot, ps=ps, fullrange=fullrange, newcoord=newcoord

common atv_state
common atv_images

if (keyword_set(ps)) then begin
    thick = 3
    color = 0
endif else begin
    thick = 1
    color = 7
endelse

if (keyword_set(newcoord)) then state.plot_coord = state.coord

if (not (keyword_set(ps))) then begin
    newplot = 0
    if (not (xregistered('atv_lineplot', /noshow))) then begin
        atv_lineplot_init
        newplot = 1
    endif 

    widget_control, state.histbutton_base_id, map=0
    widget_control, state.holdrange_button_id, sensitive=1

    widget_control, state.lineplot_xmin_id, get_value=xmin
    widget_control, state.lineplot_xmax_id, get_value=xmax
    widget_control, state.lineplot_ymin_id, get_value=ymin
    widget_control, state.lineplot_ymax_id, get_value=ymax

    if (newplot EQ 1 OR state.plot_type NE 'colplot' OR $
        keyword_set(fullrange) OR $
       (state.holdrange_value EQ 0 AND keyword_set(newcoord))) then begin
        xmin = 0.0
        xmax = state.image_size[1]
        ymin = min(main_image[state.plot_coord[0],*])
        ymax = max(main_image[state.plot_coord[0],*]) 
    endif
    
    widget_control, state.lineplot_xmin_id, set_value=xmin
    widget_control, state.lineplot_xmax_id, set_value=xmax
    widget_control, state.lineplot_ymin_id, set_value=ymin
    widget_control, state.lineplot_ymax_id, set_value=ymax

    state.lineplot_xmin = xmin
    state.lineplot_xmax = xmax
    state.lineplot_ymin = ymin
    state.lineplot_ymax = ymax

    state.plot_type = 'colplot'
    atv_setwindow, state.lineplot_window_id
    erase
    
endif


plot, main_image[state.plot_coord[0], *], $
  xst = 3, yst = 3, psym = 10, $
  title = strcompress('Plot of column ' + $
                      string(state.plot_coord[0])), $
  xtitle = 'Row', $
  ytitle = 'Pixel Value', $
  color = color, xmargin=[15,3], $
  xran = [state.lineplot_xmin, state.lineplot_xmax], $
  yran = [state.lineplot_ymin, state.lineplot_ymax], $
  thick = thick, xthick = thick, ythick = thick, charthick = thick


if (not (keyword_set(ps))) then begin 
  widget_control, state.lineplot_base_id, /clear_events
  atv_resetwindow
endif

end


;----------------------------------------------------------------------

pro atv_gaussrowplot, ps=ps, update=update

common atv_state
common atv_images

if (not (keyword_set(ps))) then begin

; Only initialize plot window and plot ranges to the min/max ranges
; when gaussrowplot window is not already present or plot window is present
; but last plot was not a gaussrowplot.  Otherwise, use the values
; currently in the min/max boxes

  if (not (xregistered('atv_lineplot', /noshow))) then begin
    atv_lineplot_init
  endif

  widget_control, state.histbutton_base_id, map=0
  widget_control, state.holdrange_butt_id, sensitive=0

  state.plot_type = 'gaussrowplot'
  atv_setwindow, state.lineplot_window_id
  erase

; must store the coordinates in state structure if you want to make a
; PS plot because state.coord array will change if you move cursor
; before pressing 'Create PS' button

  if (not (keyword_set(update))) then state.plot_coord = state.coord

  x2=long((state.plot_coord[0]+10.) < (state.image_size[0]-1.))
  x1=long((state.plot_coord[0]-10.) > 0.)
  y2=long((state.plot_coord[1]+2.) < (state.image_size[1]-1))
  y1=long((state.plot_coord[1]-2.) > 0.)
  x=fltarr(x2-x1+1)
  y=fltarr(x2-x1+1)

  n_x = x2-x1+1
  n_y = y2-y1+1

  for i=0, n_x - 1 do begin
    x[i]=x1+i
    y[i]=total(main_image[x[i],y1:y2])/(n_y)
  endfor

  x_interp=interpol(x,1000)
  y_interp=interpol(y,1000)
  yfit=gaussfit(x_interp,y_interp,a,nterms=4)
  peak = a[0]
  center = a[1]
  fwhm = a[2] * 2.354
  bkg = min(yfit)

  if (not (keyword_set(update))) then begin

    widget_control,state.lineplot_xmin_id, $
      set_value=x[0]

    state.lineplot_xmin = x[0]

    widget_control,state.lineplot_xmax_id, $
      set_value=x[n_x-1]

    state.lineplot_xmax = x[n_x-1]
  
    widget_control,state.lineplot_ymin_id, $
      set_value=min(y)

    state.lineplot_ymin = min(y)

    widget_control,state.lineplot_ymax_id, $
      set_value=(max(y) > max(yfit))

    state.lineplot_ymax = max(y) > max(yfit)

  endif

  title_str = 'Rows ' + $
              strcompress(string(y1),/remove_all) + $
              '-' + strcompress(string(y2),/remove_all) + $
              '   Center=' + strcompress(string(center,format='(f10.2)'),/remove_all) + $
              '   Peak=' + strcompress(string(peak,format='(f10.2)'),/remove_all) + $
              '   FWHM=' + strcompress(string(fwhm,format='(f10.2)'),/remove_all) + $
              '   Bkg=' + strcompress(string(bkg,format='(f10.2)'),/remove_all)

  plot,x,y,psym=1,/ynozero, title = title_str, xtitle='Column (pixels)', $
    ytitle='Pixel Value', $
    color = 7, xst = 3, yst = 3, xmargin=[15,3], $
    xran = [state.lineplot_xmin, state.lineplot_xmax], $
    yran = [state.lineplot_ymin, state.lineplot_ymax]

  oplot, x_interp, yfit

endif else begin

  x2=long((state.plot_coord[0]+10.) < (state.image_size[0]-1.))
  x1=long((state.plot_coord[0]-10.) > 0.)
  y2=long((state.plot_coord[1]+2.) < (state.image_size[1]-1))
  y1=long((state.plot_coord[1]-2.) > 0.)
  x=fltarr(x2-x1+1)
  y=fltarr(x2-x1+1)

  n_x = x2-x1+1
  n_y = y2-y1+1

  for i=0, n_x - 1 do begin
    x[i]=x1+i
    y[i]=total(main_image[x[i],y1:y2])/(n_y)
  endfor

  x_interp=interpol(x,1000)
  y_interp=interpol(y,1000)
  yfit=gaussfit(x_interp,y_interp,a,nterms=4)
  peak = a[0]
  center = a[1]
  fwhm = a[2] * 2.354
  bkg = min(yfit) 

  title_str = 'Rows ' + $
              strcompress(string(y1),/remove_all) + $
              '-' + strcompress(string(y2),/remove_all) + $
              ' Ctr=' + strcompress(string(center,format='(f10.2)'),/remove_all) + $
              ' Pk=' + strcompress(string(peak,format='(f10.2)'),/remove_all) + $
              ' FWHM=' + strcompress(string(fwhm,format='(f10.2)'),/remove_all) + $
              ' Bkg=' + strcompress(string(bkg,format='(f10.2)'),/remove_all)

  plot,x,y,psym=1,/ynozero, title = title_str, xtitle='Column (pixels)', $
    ytitle='Pixel Value', $
    color = 0, xst = 3, yst = 3, xmargin=[15,3], $
    xran = [state.lineplot_xmin, state.lineplot_xmax], $
    yran = [state.lineplot_ymin, state.lineplot_ymax]

  oplot, x_interp, yfit

endelse

if (not (keyword_set(ps))) then begin 
  widget_control, state.lineplot_base_id, /clear_events
  atv_resetwindow
endif

end

;--------------------------------------------------------------------

pro atv_gausscolplot, ps=ps, update=update

common atv_state
common atv_images

if (not (keyword_set(ps))) then begin

; Only initialize plot window and plot ranges to the min/max ranges
; when gausscolplot window is not already present or plot window is present
; but last plot was not a gausscolplot.  Otherwise, use the values
; currently in the min/max boxes

  if (not (xregistered('atv_lineplot', /noshow))) then begin
    atv_lineplot_init
  endif

  widget_control, state.histbutton_base_id, map=0
  widget_control, state.holdrange_butt_id, sensitive=0

  state.plot_type = 'gausscolplot'
  atv_setwindow, state.lineplot_window_id
  erase

; must store the coordinates in state structure if you want to make a
; PS plot because state.coord array will change if you move cursor
; before pressing 'Create PS' button

  if (not (keyword_set(update))) then state.plot_coord = state.coord

  x2=long((state.plot_coord[1]+10.) < (state.image_size[1]-1.))
  x1=long((state.plot_coord[1]-10.) > 0.)
  y2=long((state.plot_coord[0]+2.) < (state.image_size[0]-1))
  y1=long((state.plot_coord[0]-2.) > 0.)
  x=fltarr(x2-x1+1)
  y=fltarr(x2-x1+1)

  n_x = x2-x1+1
  n_y = y2-y1+1

  for i=0, n_x - 1 do begin
    x[i]=x1+i
    y[i]=total(main_image[y1:y2,x[i]])/(n_y)
  endfor

  x_interp=interpol(x,1000)
  y_interp=interpol(y,1000)
  yfit=gaussfit(x_interp,y_interp,a,nterms=4)
  peak = a[0]
  center = a[1]
  fwhm = a[2] * 2.354
  bkg = min(yfit) 

  if (not (keyword_set(update))) then begin

    widget_control,state.lineplot_xmin_id, $
      set_value=x[0]

    state.lineplot_xmin = x[0]

    widget_control,state.lineplot_xmax_id, $
      set_value=x[n_x-1]

    state.lineplot_xmax = x[n_x-1]

    widget_control,state.lineplot_ymin_id, $
      set_value=min(y)

    state.lineplot_ymin = min(y)

    widget_control,state.lineplot_ymax_id, $
      set_value=(max(y) > max(yfit))

    state.lineplot_ymax = max(y) > max(yfit)

  endif

  title_str = 'Columns ' + $
              strcompress(string(y1),/remove_all) + $
              '-' + strcompress(string(y2),/remove_all) + $
              '   Center=' + strcompress(string(center,format='(f10.2)'),/remove_all) + $
              '   Peak=' + strcompress(string(peak,format='(f10.2)'),/remove_all) + $
              '   FWHM=' + strcompress(string(fwhm,format='(f10.2)'),/remove_all) + $
              '   Bkg=' + strcompress(string(bkg,format='(f10.2)'),/remove_all)

  plot,x,y,psym=1,/ynozero, title = title_str, xtitle='Row (pixels)', $
    ytitle='Pixel Value', $
    color = 7, xst = 3, yst = 3, xmargin=[15,3], $
    xran = [state.lineplot_xmin, state.lineplot_xmax], $
    yran = [state.lineplot_ymin, state.lineplot_ymax]

  oplot, x_interp, yfit

endif else begin

  x2=long((state.plot_coord[1]+10.) < (state.image_size[1]-1.))
  x1=long((state.plot_coord[1]-10.) > 0.)
  y2=long((state.plot_coord[0]+2.) < (state.image_size[0]-1))
  y1=long((state.plot_coord[0]-2.) > 0.)
  x=fltarr(x2-x1+1)
  y=fltarr(x2-x1+1)

  n_x = x2-x1+1
  n_y = y2-y1+1

  for i=0, n_x - 1 do begin
    x[i]=x1+i
    y[i]=total(main_image[y1:y2,x[i]])/(n_y)
  endfor

  x_interp=interpol(x,1000)
  y_interp=interpol(y,1000)
  yfit=gaussfit(x_interp,y_interp,a,nterms=4)
  peak = a[0]
  center = a[1]
  fwhm = a[2] * 2.354
  bkg = min(yfit) 

  title_str = 'Cols ' + $
              strcompress(string(y1),/remove_all) + $
              '-' + strcompress(string(y2),/remove_all) + $
              ' Ctr=' + strcompress(string(center,format='(f10.2)'),/remove_all) + $
              ' Pk=' + strcompress(string(peak,format='(f10.2)'),/remove_all) + $
              ' FWHM=' + strcompress(string(fwhm,format='(f10.2)'),/remove_all) + $
              ' Bkg=' + strcompress(string(bkg,format='(f10.2)'),/remove_all)

  plot,x,y,psym=1,/ynozero, title = title_str, xtitle='Row (pixels)', $
    ytitle='Pixel Value', $
    color = 0, xst = 3, yst = 3, xmargin=[15,3], $
    xran = [state.lineplot_xmin, state.lineplot_xmax], $
    yran = [state.lineplot_ymin, state.lineplot_ymax]

  oplot, x_interp, yfit

endelse

if (not (keyword_set(ps))) then begin 
  widget_control, state.lineplot_base_id, /clear_events
  atv_resetwindow
endif

end


;----------------------------------------------------------------------

pro atv_vectorplot, ps=ps, fullrange=fullrange, newcoord=newcoord

common atv_state
common atv_images

if (keyword_set(ps)) then begin
    thick = 3
    color = 0
endif else begin
    thick = 1
    color = 7
endelse


if (state.vector_coord1[0] eq state.vector_coord2[0]) then begin
  atv_colplot
  return
endif 
if (state.vector_coord1[1] eq state.vector_coord2[1]) then begin
  atv_rowplot
  return
endif

d = sqrt((state.vector_coord1[0]-state.vector_coord2[0])^2 + $
         (state.vector_coord1[1]-state.vector_coord2[1])^2)

v_d = fix(d + 1)
dx = (state.vector_coord2[0]-state.vector_coord1[0]) / float(v_d - 1)
dy = (state.vector_coord2[1]-state.vector_coord1[1]) / float(v_d - 1)

x = fltarr(v_d)
y = fltarr(v_d)
vectdist = indgen(v_d)
pixval = fltarr(v_d)

x[0] = state.vector_coord1[0]
y[0] = state.vector_coord1[1]
for i = 1, n_elements(x) - 1 do begin
    x[i] = state.vector_coord1[0] + dx * i
    y[i] = state.vector_coord1[1] + dy * i
endfor



for j = 0, n_elements(x) - 1 do begin
    col = x[j]
    row = y[j]
    floor_col = floor(col)
    ceil_col = ceil(col)
    floor_row = floor(row)
    ceil_row = ceil(row)
    
    pixval[j] = (total([main_image[floor_col,floor_row], $
                        main_image[floor_col,ceil_row], $
                        main_image[ceil_col,floor_row], $
                        main_image[ceil_col,ceil_row]])) / 4.
    
endfor

if (not (keyword_set(ps))) then begin

    newplot = 0
    if (not (xregistered('atv_lineplot', /noshow))) then begin
        atv_lineplot_init
        newplot = 1
    endif
    
    widget_control, state.histbutton_base_id, map=0
    widget_control, state.holdrange_button_id, sensitive=1

    widget_control, state.lineplot_xmin_id, get_value=xmin
    widget_control, state.lineplot_xmax_id, get_value=xmax
    widget_control, state.lineplot_ymin_id, get_value=ymin
    widget_control, state.lineplot_ymax_id, get_value=ymax

    if (newplot EQ 1 OR state.plot_type NE 'vectorplot' OR $
        keyword_set(fullrange) OR $
       (state.holdrange_value EQ 0 AND keyword_set(newcoord))) then begin
        xmin = 0.0
        xmax = max(vectdist)
        ymin = min(pixval)
        ymax = max(pixval) 
        
    endif 

    widget_control, state.lineplot_xmin_id, set_value=xmin
    widget_control, state.lineplot_xmax_id, set_value=xmax
    widget_control, state.lineplot_ymin_id, set_value=ymin
    widget_control, state.lineplot_ymax_id, set_value=ymax

    state.lineplot_xmin = xmin
    state.lineplot_xmax = xmax
    state.lineplot_ymin = ymin
    state.lineplot_ymax = ymax

    state.plot_type = 'vectorplot'
    atv_setwindow, state.lineplot_window_id
    erase

endif
  

plottitle = strcompress('Plot of vector [' + $
                        strcompress(string(state.vector_coord1[0]) + ',' + $
                                    string(state.vector_coord1[1]), $
                                    /remove_all) + $
                        '] to [' + $
                        strcompress(string(state.vector_coord2[0]) + ',' + $
                                    string(state.vector_coord2[1]), $
                                    /remove_all) + ']')

plot, vectdist, pixval, $
  xst = 3, yst = 3, psym = 10, $
  title = plottitle, $
  xtitle = 'Vector Distance', $
  ytitle = 'Pixel Value', $
  color = color, xmargin=[15,3], $
  xran = [state.lineplot_xmin, state.lineplot_xmax], $
  yran = [state.lineplot_ymin, state.lineplot_ymax], $
  thick = thick, xthick = thick, ythick = thick, charthick = thick



if (not (keyword_set(ps))) then begin 
  widget_control, state.lineplot_base_id, /clear_events
  atv_resetwindow
endif

end

;--------------------------------------------------------------------

pro atv_surfplot, ps=ps, fullrange=fullrange, newcoord=newcoord

common atv_state
common atv_images

if (keyword_set(ps)) then begin
    thick = 3
    color = 0
endif else begin
    thick = 1
    color = 7
endelse


if (not (keyword_set(ps))) then begin

    newplot = 0
    if (not (xregistered('atv_lineplot', /noshow))) then begin
        atv_lineplot_init
        newplot = 1
    endif
    
    widget_control, state.histbutton_base_id, map=0
    widget_control, state.holdrange_button_id, sensitive=0
    
; set new plot coords if passed from a main window keyboard event
    if (keyword_set(newcoord)) then begin
        plotsize = $
          fix(min([50, state.image_size[0]/2., state.image_size[1]/2.]))
        center = plotsize > state.coord < (state.image_size[0:1] - plotsize) 
        
        shade_image = main_image[center[0]-plotsize:center[0]+plotsize-1, $
                                 center[1]-plotsize:center[1]+plotsize-1]
        
        state.lineplot_xmin = center[0]-plotsize
        state.lineplot_xmax = center[0]+plotsize-1
        state.lineplot_ymin = center[1]-plotsize 
        state.lineplot_ymax = center[1]+plotsize-1
        
        state.plot_coord = state.coord
        
        widget_control, state.lineplot_xmin_id, $
          set_value = state.lineplot_xmin
        widget_control, state.lineplot_xmax_id, $
          set_value = state.lineplot_xmax
        widget_control, state.lineplot_ymin_id, $
          set_value = state.lineplot_ymin
        widget_control, state.lineplot_ymax_id, $
          set_value = state.lineplot_ymax
    endif
    
    if (keyword_set(fullrange)) then begin
        widget_control, state.lineplot_xmin_id, set_value = 0
        widget_control, state.lineplot_xmax_id, $
          set_value = state.image_size[0]-1
        widget_control, state.lineplot_ymin_id, set_value = 0
        widget_control, state.lineplot_ymax_id, $
          set_value = state.image_size[1]-1
    endif

    state.plot_type = 'surfplot'
    atv_setwindow, state.lineplot_window_id
    erase
    
; now get plot coords from the widget box   
    widget_control,state.lineplot_xmin_id, get_value=xmin
    widget_control,state.lineplot_xmax_id, get_value=xmax
    widget_control,state.lineplot_ymin_id, get_value=ymin
    widget_control,state.lineplot_ymax_id, get_value=ymax  

    state.lineplot_xmin = xmin
    state.lineplot_xmax = xmax
    state.lineplot_ymin = ymin
    state.lineplot_ymax = ymax
endif


shade_image =  main_image[state.lineplot_xmin:state.lineplot_xmax, $
                          state.lineplot_ymin:state.lineplot_ymax]
;shades = scaled_image[state.lineplot_xmin:state.lineplot_xmax, $
;                          state.lineplot_ymin:state.lineplot_ymax]

plottitle = $
  strcompress('Surface plot of ' + $
              strcompress('['+string(round(state.lineplot_xmin))+ $
                          ':'+string(round(state.lineplot_xmax))+ $
                          ','+string(round(state.lineplot_ymin))+ $
                          ':'+string(round(state.lineplot_ymax))+ $
                          ']', /remove_all))

xdim = state.lineplot_xmax - state.lineplot_xmin + 1
ydim = state.lineplot_ymax - state.lineplot_ymin + 1

xran = lindgen(xdim) + state.lineplot_xmin
yran = lindgen(ydim) + state.lineplot_ymin

; reload the color table of the main window with default brightness
; and contrast, to make the surface plot come out ok
atv_stretchct, 0.5, 0.5

shade_surf, shade_image, xst=3, yst=3, zst=3, $
  xran, yran, $
  title = plottitle, $
  xtitle = 'X', ytitle = 'Y', ztitle = 'Pixel Value', $
  color = color, charsize=1.5, $
  thick = thick, xthick = thick, ythick = thick, zthick = thick, $
  charthick = thick    ;, shades = shades


if (not (keyword_set(ps))) then begin 
    widget_control, state.lineplot_base_id, /clear_events
    atv_resetwindow
endif

end


;--------------------------------------------------------------------

pro atv_contourplot, ps=ps, fullrange=fullrange, newcoord=newcoord

if (keyword_set(ps)) then begin
    thick = 3
    color = 0
endif else begin
    thick = 1
    color = 7
endelse

common atv_state
common atv_images

if (not (keyword_set(ps))) then begin

    newplot = 0
    if (not (xregistered('atv_lineplot', /noshow))) then begin
        atv_lineplot_init
        newplot = 1
    endif
    
    widget_control, state.histbutton_base_id, map=0
    widget_control, state.holdrange_button_id, sensitive=0
    
    if (keyword_set(newcoord)) then begin
        
        plotsize = $
          fix(min([50, state.image_size[0]/2., state.image_size[1]/2.]))
        center = plotsize > state.coord < (state.image_size[0:1] - plotsize) 
        
        contour_image =  main_image[center[0]-plotsize:center[0]+plotsize-1, $
                                    center[1]-plotsize:center[1]+plotsize-1]
        
        state.lineplot_xmin = center[0]-plotsize
        state.lineplot_xmax = center[0]+plotsize-1
        state.lineplot_ymin = center[1]-plotsize
        state.lineplot_ymax = center[1]+plotsize-1
                
        state.plot_coord = state.coord

        widget_control,state.lineplot_xmin_id, $
          set_value=state.lineplot_xmin
        widget_control,state.lineplot_xmax_id, $
          set_value=state.lineplot_xmax
        widget_control,state.lineplot_ymin_id, $
          set_value=state.lineplot_ymin
        widget_control,state.lineplot_ymax_id, $
          set_value=state.lineplot_ymax
    endif
    
    if (keyword_set(fullrange)) then begin
        widget_control, state.lineplot_xmin_id, set_value = 0
        widget_control, state.lineplot_xmax_id, $
          set_value = state.image_size[0]-1
        widget_control, state.lineplot_ymin_id, set_value = 0
        widget_control, state.lineplot_ymax_id, $
          set_value = state.image_size[1]-1
    endif

    state.plot_type = 'contourplot'
    atv_setwindow, state.lineplot_window_id
    erase
    
; now get plot coords from the widget box   
    widget_control,state.lineplot_xmin_id, get_value=xmin
    widget_control,state.lineplot_xmax_id, get_value=xmax
    widget_control,state.lineplot_ymin_id, get_value=ymin
    widget_control,state.lineplot_ymax_id, get_value=ymax  

    state.lineplot_xmin = xmin
    state.lineplot_xmax = xmax
    state.lineplot_ymin = ymin
    state.lineplot_ymax = ymax   
endif

contour_image =  main_image[state.lineplot_xmin:state.lineplot_xmax, $
                            state.lineplot_ymin:state.lineplot_ymax]


if (state.scaling EQ 1) then begin
    contour_image = alog10(contour_image)
    logflag = 'Log'
endif else begin
    logflag = ''
endelse

plottitle =  $
  strcompress(logflag + $
              ' Contour plot of ' + $
              strcompress('['+string(round(state.lineplot_xmin))+ $
                          ':'+string(round(state.lineplot_xmax))+ $
                          ','+string(round(state.lineplot_ymin))+ $
                          ':'+string(round(state.lineplot_ymax))+ $
                          ']', /remove_all))


xdim = state.lineplot_xmax - state.lineplot_xmin + 1
ydim = state.lineplot_ymax - state.lineplot_ymin + 1

xran = lindgen(xdim) + state.lineplot_xmin
yran = lindgen(ydim) + state.lineplot_ymin

contour, temporary(contour_image), xst=3, yst=3, $
  xran, yran, $
  nlevels = 10, $
  /follow, $
  title = plottitle, $
  xtitle = 'X', ytitle = 'Y', color = color, $
  thick = thick, xthick = thick, ythick = thick, charthick = thick


if (not (keyword_set(ps))) then begin 
  widget_control, state.lineplot_base_id, /clear_events
  atv_resetwindow
endif

end

;----------------------------------------------------------------------

pro atv_histplot, ps=ps, fullrange=fullrange, newcoord=newcoord

common atv_state
common atv_images


if (keyword_set(ps)) then begin
    thick = 3
    color = 0
endif else begin
    thick = 1
    color = 7
endelse

if (not (keyword_set(ps))) then begin

    newplot = 0
    if (not (xregistered('atv_lineplot', /noshow))) then begin
        atv_lineplot_init
        newplot = 1
    endif
    
    widget_control, state.histbutton_base_id, map=1
    widget_control, state.holdrange_button_id, sensitive=0
    
    if (keyword_set(newcoord)) then begin
        
        state.plot_coord = state.coord
        plotsize_x = $
          fix(min([20, state.image_size[0]/2.]))
        plotsize_y = $
          fix(min([20, state.image_size[1]/2.]))
        
; Establish pixel boundaries to histogram
        x1 = (state.plot_coord[0]-plotsize_x) > 0.
        x2 = (state.plot_coord[0]+plotsize_x) < (state.image_size[0]-1)
        y1 = (state.plot_coord[1]-plotsize_y) > 0.
        y2 = (state.plot_coord[1]+plotsize_y) < (state.image_size[1]-1)
        
        widget_control, state.x1_pix_id, set_value=x1
        widget_control, state.x2_pix_id, set_value=x2
        widget_control, state.y1_pix_id, set_value=y1
        widget_control, state.y2_pix_id, set_value=y2
    endif

    state.plot_type = 'histplot'
    atv_setwindow, state.lineplot_window_id
    erase
endif


; get histogram region 
widget_control, state.x1_pix_id, get_value=x1
widget_control, state.x2_pix_id, get_value=x2
widget_control, state.y1_pix_id, get_value=y1
widget_control, state.y2_pix_id, get_value=y2        
hist_image = main_image[x1:x2, y1:y2]

; initialize the binsize if necessary
if (state.binsize EQ 0 OR keyword_set(newcoord)) then begin
    nbins = 50.
    state.binsize = (float(max(hist_image)) - float(min(hist_image)) ) / nbins
    if (abs(state.binsize) GT 10) then $
      state.binsize = fix(state.binsize)
    widget_control, state.histplot_binsize_id, set_value=state.binsize
endif

; Call plothist to create histogram arrays
plothist, hist_image, xhist, yhist, bin=state.binsize, /NaN, /noplot
   
; Only initialize plot window and plot ranges to the min/max ranges
; when histplot window is not already present or plot window is present
; but last plot was not a histplot.  Otherwise, use the values
; currently in the min/max boxes

if (keyword_set(newcoord) OR keyword_set(fullrange)) then begin
    state.lineplot_xmin = min(hist_image)
    state.lineplot_xmax = max(hist_image)
    state.lineplot_ymin = 0.
    state.lineplot_ymax = round(max(yhist) * 1.1)
    widget_control, state.lineplot_xmin_id, set_value = state.lineplot_xmin
    widget_control, state.lineplot_xmax_id, set_value = state.lineplot_xmax
    widget_control, state.lineplot_ymin_id, set_value = state.lineplot_ymin
    widget_control, state.lineplot_ymax_id, set_value = state.lineplot_ymax
endif

widget_control, state.histplot_binsize_id, get_value=binsize
widget_control, state.lineplot_xmin_id, get_value=xmin
widget_control, state.lineplot_xmax_id, get_value=xmax
widget_control, state.lineplot_ymin_id, get_value=ymin
widget_control, state.lineplot_ymax_id, get_value=ymax

state.binsize = binsize
state.lineplot_xmin = xmin
state.lineplot_xmax = xmax
state.lineplot_ymin = ymin
state.lineplot_ymax = ymax

plottitle = $
  strcompress('Histogram plot of ' + $
              strcompress('['+string(round(x1))+ $
                          ':'+string(round(x2))+ $
                          ','+string(round(y1))+ $
                          ':'+string(round(y2))+ $
                          ']', /remove_all))

;Plot histogram with proper ranges
plothist, hist_image, xhist, yhist, bin=state.binsize, /NaN, $
  xtitle='Pixel Value', ytitle='Number', title=plottitle, $
  xran=[state.lineplot_xmin,state.lineplot_xmax], $
  yran=[state.lineplot_ymin,state.lineplot_ymax], $
  xstyle=1, ystyle=1, color=color, $
  thick = thick, xthick = thick, ythick = thick, charthick = thick


if (not (keyword_set(ps))) then begin 
    widget_control, state.lineplot_base_id, /clear_events
    atv_resetwindow
endif

end


;----------------------------------------------------------------------

pro atv_lineplot_event, event

common atv_state
common atv_images

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'lineplot_done': begin
        widget_control, event.top, /destroy
        state.plot_type = ''
    end

    'lineplot_base': begin      ; Resize event
        state.lineplot_size = [event.x, event.y]- state.lineplot_pad
        widget_control, state.lineplot_widget_id, $
          xsize = (state.lineplot_size[0] > state.lineplot_min_size[0]), $
          ysize = (state.lineplot_size[1] > state.lineplot_min_size[1])

        case state.plot_type of
            'rowplot': atv_rowplot
            'colplot': atv_colplot
            'vectorplot': atv_vectorplot
            'histplot': atv_histplot
            'surfplot': atv_surfplot
            'contourplot': atv_contourplot
        endcase
    end

    'lineplot_holdrange': begin
        if (state.holdrange_value eq 1) then state.holdrange_value = 0 $
        else state.holdrange_value = 1
    end
    
    'lineplot_fullrange': begin
       case state.plot_type of
         'rowplot': begin

            widget_control,state.lineplot_xmin_id, $
              set_value=0

            state.lineplot_xmin = 0.0

            widget_control,state.lineplot_xmax_id, $
              set_value=state.image_size[0]

            state.lineplot_xmax = state.image_size[0]

            widget_control,state.lineplot_ymin_id, $
              set_value=min(main_image[*,state.plot_coord[1]])

            state.lineplot_ymin = min(main_image[*,state.plot_coord[1]])

            widget_control,state.lineplot_ymax_id, $
              set_value=max(main_image[*,state.plot_coord[1]])

            state.lineplot_ymax = max(main_image[*,state.plot_coord[1]]) 

            atv_rowplot, /update
         end
         'colplot': begin

            widget_control,state.lineplot_xmin_id, $
              set_value=0

            state.lineplot_xmin = 0.0

            widget_control,state.lineplot_xmax_id, $
              set_value=state.image_size[1]

            state.lineplot_xmax = state.image_size[1]

            widget_control,state.lineplot_ymin_id, $
              set_value=min(main_image[state.plot_coord[0], *])

            state.lineplot_ymin = min(main_image[state.plot_coord[0], *])

            widget_control,state.lineplot_ymax_id, $
              set_value=max(main_image[state.plot_coord[0], *])

            state.lineplot_ymax = max(main_image[state.plot_coord[0], *])

            atv_colplot, /update
         end
         'gaussrowplot': begin

            x2=long((state.plot_coord[0]+10.) < (state.image_size[0]-1.))
            x1=long((state.plot_coord[0]-10.) > 0.)
            y2=long((state.plot_coord[1]+2.) < (state.image_size[1]-1))
            y1=long((state.plot_coord[1]-2.) > 0.)
            x=fltarr(x2-x1+1)
            y=fltarr(x2-x1+1)

            n_x = x2-x1+1
            n_y = y2-y1+1

            for i=0, n_x - 1 do begin
              x[i]=x1+i
              y[i]=total(main_image[x[i],y1:y2])/(n_y)
            endfor

            x_interp=interpol(x,1000)
            y_interp=interpol(y,1000)
            yfit=gaussfit(x_interp,y_interp,a,nterms=4)

            widget_control,state.lineplot_xmin_id, $
              set_value=x[0]

            state.lineplot_xmin = x[0]

            widget_control,state.lineplot_xmax_id, $
              set_value=x[n_x-1]

            state.lineplot_xmax = x[n_x-1]

            widget_control,state.lineplot_ymin_id, $
              set_value=min(y)

            state.lineplot_ymin = min(y)

            widget_control,state.lineplot_ymax_id, $
              set_value=(max(y) > max(yfit))

            state.lineplot_ymax = max(y) > max(yfit)

            atv_gaussrowplot, /update
         end
         'gausscolplot': begin

            x2=long((state.plot_coord[1]+10.) < (state.image_size[1]-1.))
            x1=long((state.plot_coord[1]-10.) > 0.)
            y2=long((state.plot_coord[0]+2.) < (state.image_size[0]-1))
            y1=long((state.plot_coord[0]-2.) > 0.)
            x=fltarr(x2-x1+1)
            y=fltarr(x2-x1+1)

            n_x = x2-x1+1
            n_y = y2-y1+1

            for i=0, n_x - 1 do begin
              x[i]=x1+i
              y[i]=total(main_image[y1:y2,x[i]])/(n_y)
            endfor

            x_interp=interpol(x,1000)
            y_interp=interpol(y,1000)
            yfit=gaussfit(x_interp,y_interp,a,nterms=4)

            widget_control,state.lineplot_xmin_id, $
              set_value=x[0]

            state.lineplot_xmin = x[0]

            widget_control,state.lineplot_xmax_id, $
              set_value=x[n_x-1]

            state.lineplot_xmax = x[n_x-1]

            widget_control,state.lineplot_ymin_id, $
              set_value=min(y)

            state.lineplot_ymin = min(y)

            widget_control,state.lineplot_ymax_id, $
              set_value=(max(y) > max(yfit))

            state.lineplot_ymax = max(y) > max(yfit)

            atv_gausscolplot, /update
         end
         'vectorplot': begin

            d = sqrt((state.vector_coord1[0]-state.vector_coord2[0])^2 + $
                    (state.vector_coord1[1]-state.vector_coord2[1])^2)

            v_d = fix(d + 1)
            dx = (state.vector_coord2[0]-state.vector_coord1[0]) / float(v_d - 1)
            dy = (state.vector_coord2[1]-state.vector_coord1[1]) / float(v_d - 1)

            x = fltarr(v_d)
            y = fltarr(v_d)
            vectdist = indgen(v_d)
            pixval = fltarr(v_d)

            x[0] = state.vector_coord1[0]
            y[0] = state.vector_coord1[1]

            for i = 1, n_elements(x) - 1 do begin
              x[i] = state.vector_coord1[0] + dx * i
              y[i] = state.vector_coord1[1] + dy * i
            endfor

            for j = 0, n_elements(x) - 1 do begin
              col = x[j]
              row = y[j]
              floor_col = floor(col)
              ceil_col = ceil(col)
              floor_row = floor(row)
              ceil_row = ceil(row)
    
              pixval[j] = (total([main_image[floor_col,floor_row], $
                                  main_image[floor_col,ceil_row], $
                                  main_image[ceil_col,floor_row], $
                                  main_image[ceil_col,ceil_row]])) / 4.

            endfor

            widget_control,state.lineplot_xmin_id, set_value=0
            state.lineplot_xmin = 0.0

            widget_control,state.lineplot_xmax_id, set_value=max(vectdist)
            state.lineplot_xmax = max(vectdist)

            widget_control,state.lineplot_ymin_id, set_value=min(pixval)
            state.lineplot_ymin = min(pixval)

            widget_control,state.lineplot_ymax_id, set_value=max(pixval)
            state.lineplot_ymax = max(pixval) 

            atv_vectorplot, /update

         end
         'histplot': begin

            plotsize_x = $
              fix(min([20, state.image_size[0]/2.]))

            plotsize_y = $
              fix(min([20, state.image_size[1]/2.]))

         ; Establish pixel boundaries to histogram
            x1 = (state.plot_coord[0]-plotsize_x) > 0.
            x2 = (state.plot_coord[0]+plotsize_x) < (state.image_size[0]-1)
            y1 = (state.plot_coord[1]-plotsize_y) > 0.
            y2 = (state.plot_coord[1]+plotsize_y) < (state.image_size[1]-1)

         ; Set up histogram pixel array.  User may do rectangular regions.
            hist_image = main_image[x1:x2, y1:y2]

            state.lineplot_xmin = min(hist_image)
            state.lineplot_xmin_orig = state.lineplot_xmin
            state.lineplot_xmax = max(hist_image)
            state.lineplot_xmax_orig = state.lineplot_xmax
            state.lineplot_ymin = 0

            widget_control, state.lineplot_xmin_id, $
              set_value = state.lineplot_xmin

            widget_control, state.lineplot_xmax_id, $
              set_value = state.lineplot_xmax

            widget_control, state.lineplot_ymin_id, $
              set_value = state.lineplot_ymin

            state.binsize = (state.lineplot_xmax - state.lineplot_xmin) * 0.01
            state.binsize = state.binsize > $
              ((state.lineplot_xmax - state.lineplot_xmin) * 1.0e-5)
            state.binsize = fix(state.binsize)

            widget_control, state.x1_pix_id, set_value=x1
            widget_control, state.x2_pix_id, set_value=x2
            widget_control, state.y1_pix_id, set_value=y1
            widget_control, state.y2_pix_id, set_value=y2
            widget_control, state.histplot_binsize_id, set_value=state.binsize

            ;Set lineplot window and erase
            atv_setwindow, state.lineplot_window_id
            erase

            ;Call plothist to create histogram arrays.  Necessary to get 
            ;default ymax
            plothist, hist_image, xhist, yhist, bin=state.binsize, $
              /NaN, /nodata

            state.lineplot_ymax = fix(max(yhist) + 0.05*max(yhist))

            widget_control, state.lineplot_ymax_id, $
              set_value = state.lineplot_ymax            

            atv_histplot, /update

         end
         'surfplot': begin

            plotsize = $
              fix(min([50, state.image_size[0]/2., state.image_size[1]/2.]))
            center = plotsize > state.plot_coord < $
                     (state.image_size[0:1] - plotsize) 

            state.lineplot_xmin = (center[0]-plotsize)
            state.lineplot_xmax = (center[0]+plotsize-1)
            state.lineplot_ymin = (center[1]-plotsize)
            state.lineplot_ymax = (center[1]+plotsize-1)

            widget_control,state.lineplot_xmin_id, $
              set_value=state.lineplot_xmin

            widget_control,state.lineplot_xmax_id, $
              set_value=state.lineplot_xmax

            widget_control,state.lineplot_ymin_id, $
              set_value=state.lineplot_ymin

            widget_control,state.lineplot_ymax_id, $
              set_value=state.lineplot_ymax

            atv_surfplot, /update
        end
         'contourplot': begin

            plotsize = $
              fix(min([50, state.image_size[0]/2., state.image_size[1]/2.]))
            center = plotsize > state.plot_coord < $
                     (state.image_size[0:1] - plotsize) 

            state.lineplot_xmin = (center[0]-plotsize)
            state.lineplot_xmax = (center[0]+plotsize-1)
            state.lineplot_ymin = (center[1]-plotsize)
            state.lineplot_ymax = (center[1]+plotsize-1)

            widget_control,state.lineplot_xmin_id, $
              set_value=state.lineplot_xmin

            widget_control,state.lineplot_xmax_id, $
              set_value=state.lineplot_xmax

            widget_control,state.lineplot_ymin_id, $
              set_value=state.lineplot_ymin

            widget_control,state.lineplot_ymax_id, $
              set_value=state.lineplot_ymax

            atv_contourplot, /update
        end
         'slice3dplot': begin

            widget_control,state.lineplot_xmin_id, $
              set_value=0

            state.lineplot_xmin = 0.0

            widget_control,state.lineplot_xmax_id, $
              set_value=state.image_size[2]

            state.lineplot_xmax = state.image_size[2]

            widget_control,state.lineplot_ymin_id, $
            set_value=min(main_image_stack[state.plot_coord[0], $
                                           state.plot_coord[1], *])

            state.lineplot_ymin = $
              min(main_image_stack[state.plot_coord[0], $
                                   state.plot_coord[1], *])

            widget_control,state.lineplot_ymax_id, $
            set_value=max(main_image_stack[state.plot_coord[0], $
                                           state.plot_coord[1], *])

            state.lineplot_ymax = $
              max(main_image_stack[state.plot_coord[0], $
                                   state.plot_coord[1], *]) 

            atv_slice3dplot, /update
        end
         else:
       endcase
    end
    'lineplot_ps': begin
        fname = strcompress(state.current_dir + 'atv_plot.ps', /remove_all)
        lpforminfo = cmps_form(cancel = canceled, create = create, $
                             parent = state.lineplot_base_id, $
                             /preserve_aspect, $
                             /color, $
                             /nocommon, papersize='Letter', $
                             filename = fname, $
                             button_names = ['Create PS File'])
        
        if (canceled) then return
        if (lpforminfo.filename EQ '') then return
        
        tmp_result = findfile(lpforminfo.filename, count = nfiles)
        
        result = ''
        if (nfiles GT 0) then begin
            mesg = strarr(2)
            mesg[0] = 'Overwrite existing file:'
            tmp_string = strmid(lpforminfo.filename, strpos(lpforminfo.filename, $
                                                          '/') + 1)
            mesg[1] = strcompress(tmp_string + '?', /remove_all)
            result =  dialog_message(mesg, $
                                     /default_no, $
                                     dialog_parent = state.base_id, $
                                     /question)                 
        endif
        
        if (strupcase(result) EQ 'NO') then return
        
        widget_control, /hourglass
        
        screen_device = !d.name
        set_plot, 'ps'
        device, _extra = lpforminfo
        
        case (state.plot_type) of
            'rowplot': atv_rowplot, /ps
            'colplot': atv_colplot, /ps
         'gaussrowplot': begin
            atv_gaussrowplot, /ps
         end
         'gausscolplot': begin
            atv_gausscolplot, /ps
         end
            'vectorplot': atv_vectorplot, /ps
            'histplot': atv_histplot, /ps
         'surfplot': begin
           if (state.lineplot_xmin ge state.lineplot_xmax OR $
               state.lineplot_ymin ge state.lineplot_ymax) then begin
             atv_message, 'XMin and YMin must be less than Xmax and YMax', $
               msgtype='error', /window
             device, /close
             set_plot, screen_device
             return
           endif

            atv_surfplot, /ps
         end

         'contourplot': begin
           if (state.lineplot_xmin ge state.lineplot_xmax OR $
               state.lineplot_ymin ge state.lineplot_ymax) then begin
             atv_message, 'XMin and YMin must be less than Xmax and YMax', $
               msgtype='error', /window
             device, /close
             set_plot, screen_device
             return
           endif

            atv_contourplot, /ps
         end

         'slice3dplot': begin
            atv_slice3dplot, /ps
         end

         else:
        endcase
        
        device, /close
        set_plot, screen_device
        
    end    
    
    'lineplot_newrange': begin
        
        widget_control, state.lineplot_xmin_id, get_value = xmin
        widget_control, state.lineplot_xmax_id, get_value = xmax
        widget_control, state.lineplot_ymin_id, get_value = ymin
        widget_control, state.lineplot_ymax_id, get_value = ymax
        
        ; check plot ranges for validity
        if (state.plot_type EQ 'surfplot' OR $
            state.plot_type EQ 'contourplot') then begin
            
            xmin = fix(round(0 > xmin < (state.image_size[0] - 2)))
            xmax = fix(round(1 > xmax < (state.image_size[0] - 1)))
            ymin = fix(round(0 > ymin < (state.image_size[1] - 2)))
            ymax = fix(round(1 > ymax < (state.image_size[1] - 1)))

            if (event.id EQ state.lineplot_xmin_id) then $
              if (xmin GT xmax) then xmin = xmax-1
            if (event.id EQ state.lineplot_xmax_id) then $
              if (xmax LT xmin) then xmax = xmin+1
            if (event.id EQ state.lineplot_ymin_id) then $
              if (ymin GT ymax) then ymin = ymax-1
            if (event.id EQ state.lineplot_xmax_id) then $
              if (ymax LT ymin) then ymax = ymin+1
            
        endif
        
        state.lineplot_xmin = xmin
        state.lineplot_xmax = xmax
        state.lineplot_ymin = ymin
        state.lineplot_ymax = ymax
        
        widget_control, state.lineplot_xmin_id, set_value = xmin
        widget_control, state.lineplot_xmax_id, set_value = xmax
        widget_control, state.lineplot_ymin_id, set_value = ymin
        widget_control, state.lineplot_ymax_id, set_value = ymax

        case state.plot_type of
            'rowplot': atv_rowplot
            'colplot': atv_colplot
         	'gaussrowplot': atv_gaussrowplot, /update
         	'gausscolplot': atv_gausscolplot, /update
            'vectorplot': atv_vectorplot
            'surfplot': atv_surfplot
            'contourplot': atv_contourplot
         'histplot': begin
            widget_control, state.histplot_binsize_id, get_value = tmp_value
            state.binsize = tmp_value
            atv_histplot, /update
         end
         'surfplot': begin
            atv_surfplot, /update
         end
         'slice3dplot': begin
            atv_slice3dplot, /update
         end
         else:
       endcase
    end

else:
endcase

end

;----------------------------------------------------------------------
pro atv_polarim_event, event

	common atv_state
;	common atv_images
	common atv_pdata

	widget_control, event.id, get_uvalue = uvalue
		nplot = state.polarim_plotindex

	case uvalue of
	    'polarim_display': state.polarim_display = event.value
	    'polarim_highth': state.polarim_highthresh = float(event.value)
	    'polarim_lowth': state.polarim_lowthresh = float(event.value)
	    'magnification': (*(plot_ptr[nplot])).options.magnification = event.value
	    'polarim_offset': (*(plot_ptr[nplot])).options.thetaoffset = event.value
	    'Refresh': atv_refresh
	    'polarim_done': widget_control, event.top, /destroy
	    else:
	endcase
	atv_refresh

end


pro atv_polarim

; aperture photometry front end

common atv_state
common atv_pdata

if state.polarim_present eq 0 then begin
	message,/info,"No polarimetry data is present in ATV!"
	return
endif

if (not (xregistered('atv_polarim', /noshow))) then begin

	nplot = state.polarim_plotindex

    polarim_base = $
      widget_base(/base_align_center, $
                  group_leader = state.base_id, $
                  /column, $
                  title = 'atv Polarimetry', $
                  uvalue = 'polarim_base')
    
    polarim_top_base = widget_base(polarim_base, /row, /base_align_center)

    polarim_data_base1 = widget_base( $
            polarim_top_base, /column, frame=0)

;    polarim_data_base2 = widget_base( $
;            polarim_top_base, /column, frame=0)

;    polarim_draw_base = widget_base( $
 ;           polarim_base, /row, /base_align_center, frame=0)

    polarim_data_base1a = widget_base(polarim_data_base1, /column, frame=1)

    state.polarim_display_id = $
      cw_field(polarim_data_base1a, $
               /long, $
               /return_events, $
               title = 'Overplot vectors?:', $
               uvalue = 'polarim_display', $
               value = state.polarim_display, $
               xsize = 5)

    state.polarim_lowth_id = $
      cw_field(polarim_data_base1a, $
               /float, $
               /return_events, $
               title = 'Minimum Pol Vector: ', $
               uvalue = 'polarim_lowth', $
               value = state.polarim_lowthresh, $
               xsize = 10)

    state.polarim_highth_id = $
      cw_field(polarim_data_base1a, $
               /float, $
               /return_events, $
               title = 'Maximum Pol Vector:', $
               uvalue = 'polarim_highth', $
               value = state.polarim_highthresh, $
               xsize = 10)
    
    
    state.polarim_mag_id = $
      cw_field(polarim_data_base1a, $
               /float, $
               /return_events, $
               title = 'Magnification:', $
               uvalue = 'magnification', $
               value = (*(plot_ptr[nplot])).options.magnification, $
               xsize = 10)

    state.polarim_offset_id = $
      cw_field(polarim_data_base1a, $
               /float, $
               /return_events, $
               title = 'Pos Angle Offset: ', $
               uvalue = 'polarim_offset', $
               value = (*(plot_ptr[nplot])).options.thetaoffset, $
               xsize = 10)

    
    state.photwarning_id = $
      widget_label(polarim_data_base1, $
                   value='-------------------------', $
                   /dynamic_resize, $
                   frame=1)

    photsettings_id = $
      widget_button(polarim_data_base1, $
                    value = 'Refresh', $
                    uvalue = 'Refresh')

    polarim_done = $
      widget_button(polarim_data_base1, $
                    value = 'Done', $
                    uvalue = 'polarim_done')

    widget_control, polarim_base, /realize

    ;widget_control, photzoom_widget_id, get_value=tmp_value
    ;state.photzoom_window_id = tmp_value

    xmanager, 'atv_polarim', polarim_base, /no_block
    
    atv_resetwindow
endif


end


;----------------------------------------------------------------------
;                         help window
;---------------------------------------------------------------------

pro atv_help
common atv_state

h = strarr(235)
i = 0
h[i] =  'ATV HELP'
i = i + 1
h[i] =  ''
i = i + 1
h[i] =  'MENU BAR:'
i = i + 1
h[i] =  'File->ReadFits:         Read in a new fits image from disk'
i = i + 1
h[i] =  'File->WriteFits:          Write out a new fits image to disk (single-plane or entire image)'
i = i + 1
h[i] =  'File->WritePS:          Write a PostScript file of the current display'
i = i + 1
h[i] =  'File->WriteImage:       Write an output png, jpg, or tiff image of the current display'
i = i + 1
h[i] =  'File->GetImage:         Download an archival image based on object name or coordinates'  
i = i + 1
h[i] =  'File->Quit:             Quits atv'
i = i + 1
h[i] =  'ColorMap Menu:          Selects color table'
i = i + 1
h[i] =  'Scaling Menu:           Selects linear, log, or histogram-equalized scaling'
i = i + 1
h[i] =  'Labels->TextLabel:      Brings up a dialog box for text input'
i = i + 1
h[i] =  'Labels->Arrow:            Brings up a dialog box for overplotting arrows'
i = i + 1
h[i] =  'Labels->Contour:        Brings up a dialog box for overplotting contours'
i = i + 1
h[i] =  'Labels->Compass:        Draws a compass (requires WCS info in header)'
i = i + 1
h[i] =  'Labels->Scalebar:       Draws a scale bar (requires WCS info in header)'
i = i + 1
h[i] =  'Labels->Region:           Brings up a dialog box for overplotting regions'
i = i + 1
h[i] =  'Labels->WCS Grid:         Draws a WCS grid on current image'
i = i + 1

h[i] =  'Labels->EraseLast:      Erases the most recent plot label'
i = i + 1
h[i] =  'Labels->EraseAll:       Erases all plot labels'
i = i + 1
h[i] =  'Blink->SetBlink:        Sets the current display to be the blink image'
i = i + 1
h[i] =  '                             for mouse button 1, 2, or 3'
i = i + 1
h[i] =  'Blink->MakeRGB:         Make an RGB truecolor image from the 3 blink channels'
i = i + 1
h[i] =  'Rotate/Zoom->Rotate:    Rotate image clockwise by an arbitrary angle'
i = i + 1
h[i] =  'Rotate/Zoom->90, 180, or 270 deg: rotates clockwise'
i = i + 1
h[i] =  'Rotate/Zoom->Invert:    Inverts image along x, y, or both axes'
i = i + 1
h[i] =  'Zoom->Zoom In:            Zoom in by 2x'
i = i + 1
h[i] =  'Zoom->Zoom Out:           Zoom out by 2x'
i = i + 1
h[i] =  'Rotate/Zoom->1/16x, etc: Sets zoom factor to selected scaling'
i = i + 1
h[i] =  'Zoom->Center:             Recenter image'
i = i + 1
h[i] =  'Zoom->None:               Invert to original image
i = i + 1
h[i] =  'Zoom->Invert X:           Invert the X-axis of the original image'
i = i + 1
h[i] =  'Zoom->Invert Y:           Invert the Y-axis of the original image'
i = i + 1
h[i] =  'Zoom->Invert X&Y:         Invert both the X and Y axes of the original image'
i = i + 1
h[i] =  'Zoom->Rotate:             Rotate image by arbitrary angle'
i = i + 1
h[i] =  'Zoom->0 deg:              Rotate image to original orientation'          
i = i + 1
h[i] =  'Zoom->90 deg:             Rotate original image by 90 degrees'
i = i + 1
h[i] =  'Zoom->180 deg:            Rotate original image by 180 degrees'    
i = i + 1
h[i] =  'Zoom->270 deg:            Rotate original image by 270 degrees'        
i = i + 1

h[i] =  ''
i = i + 1

h[i] =  'ImageInfo->ImageHeader: Display the FITS header, if there is one.'
i = i + 1
h[i] =  'ImageInfo->Photometry:  Brings up photometry window'
i = i + 1
h[i] =  'ImageInfo->Pixel Table:   Brings up a pixel table that tracks as the cursor moves'
i = i + 1
h[i] =  'ImageInfo->Load Regions:  Load in an SAOImage/DS9 region file and overplot on image'
i = i + 1
h[i] =  '                            Region files must be in the following format:'
i = i + 1
h[i] =  ''
i = i + 1
h[i] =  '                              circle( xcenter, ycenter, radius)'
i = i + 1
h[i] =  '                              box( xcenter, ycenter, xwidth, ywidth)'
i = i + 1
h[i] =  '                              ellipse( xcenter, ycenter, xwidth, ywidth, rotation angle)'
i = i + 1
h[i] =  '                              polygon( x1, y1, x2, y2, x3, y3, ...)'
i = i + 1
h[i] =  '                              line( x1, y1, x2, y2)'
i = i + 1
h[i] = ' '
i = i + 1
h[i] =  '                          Coordinates may be specified in pixels or WCS.  Radii and widths'
i = i + 1
h[i] =  '                          are specified in pixels or arcminutes.  For example,'
i = i + 1
h[i] = ' '
i = i + 1
h[i] =  '                              circle( 100.5, 46.3, 10.0)'
i = i + 1
h[i] = ' '
i = i + 1
h[i] =  '                          draws a circle with a radius of 10 pixels, centered at (100.5, 46.3)'
i = i + 1
h[i] =  ' '
i = i + 1
h[i] =  '                              circle(00:47:55.121, -25:22:11.98, 0.567)'
i = i + 1 
h[i] =  ' '
i = i + 1
h[i] =  '                          draws a circle with a radius of 0.567 arcminutes, centered at (00:47:55.121, -25:22:11.98)'
i = i + 1
h[i] =  ' '
i = i + 1
h[i] =  '                          The coordinate system for the region coordinates may be specified by'
i = i + 1
h[i] = ' '
i = i + 1
h[i] =  '                              circle(00:47:55.121, -25:22:11.98, 0.567, J2000)'
i = i + 1
h[i] =  '                              circle(11.97967, -25.36999, 0.567, J2000)'
i = i + 1
h[i] =  '                              circle(00:45:27.846, -25:38:33.51, 0.567, B1950)'
i = i + 1
h[i] =  '                              circle(11.366, -25.6426, 0.567, B1950)'
i = i + 1
h[i] =  '                              circle(98.566, -88.073, 0.567, galactic)'
i = i + 1
h[i] =  '                              circle(0.10622, -27.88563, 0.567, ecliptic)'
i = i + 1
h[i] =  ' '
i = i + 1
h[i] =  '                          If no coordinate system is given and coordinates are in colon-separated WCS format, the'
i = i + 1
h[i] =  '                          native coordinate system is used.'
i = i + 1
h[i] =  ' '
i = i + 1
h[i] =  '                          Region color may be specified for the following colors in the format below:'
i = i + 1
h[i] =  '                          Red, Black, Green, Blue, Cyan, Magenta, Yellow, White'
i = i + 1
h[i] = ' '
i = i + 1
h[i] = '                               circle(100.5, 46.3, 10.0) # color=red
i = i + 1
h[i] = ' '
i = i + 1
h[i] =  '                          Region text may be specified in the following format:'
i = i + 1
h[i] = ' '
i = i + 1
h[i] = '                               circle(100.5, 46.3, 10.0) # text={Text written here}
i = i + 1
h[i] = ' '
i = i + 1
h[i] =  'ImageInfo->RA,dec(J2000): Coordinates displayed are RA,dec (J2000)'
i = i + 1
h[i] =  'ImageInfo->RA,dec(B1950): Coordinates displayed are RA,dec (B1950)'
i = i + 1
h[i] =  'ImageInfo->RA,dec(J2000) deg: Coordinates displayed are RA,dec (J2000) in degrees'
i = i + 1
h[i] =  'ImageInfo->Galactic:      Coordinates displayed are Galactic coordinates'
i = i + 1
h[i] =  'ImageInfo->Ecliptic(J2000): Coordinates displayed are Ecliptic (J2000)'
i = i + 1
h[i] =  'ImageInfo->Native:        Coordinates displayed are those of the image'
i = i + 1
h[i] =  'ImageInfo->Save Regions:  Save currently displayed regions to a SAOImage/DS9 region file'
i = i + 1
h[i] =  'ImageInfo->Archive Image: Dialog to download DSS, 2MASS, or IRAS images into ATV'
i = i + 1
h[i] =  ''
i = i + 1
h[i] =  'CONTROL PANEL ITEMS:'
i = i + 1
h[i] = 'Min:             shows minimum data value displayed; enter new min value here'
i = i + 1
h[i] = 'Max:             shows maximum data value displayed; enter new max value here'
i = i + 1
h[i] = 'Pan Window:      use mouse to drag the image-view box around'
i = i + 1
h[i] = ''
i = i + 1
h[i] = 'MOUSE MODE SELECTOR:'
i = i + 1
h[i] =  'Color:          sets color-stretch mode:'
i = i + 1
h[i] = '                    With mouse button 1 down, drag mouse to change the color stretch.  '
i = i + 1
h[i] = '                    Move vertically to change contrast, and'
i = i + 1
h[i] = '                         horizontally to change brightness.'
i = i + 1 
h[i] = '                    button 2 or 3: center on current position'
i = i + 1
h[i] = 'Zoom:           sets zoom mode:' 
i = i + 1 
h[i] = '                    button1: zoom in & center on current position'
i = i + 1
h[i] = '                    button2: center on current position'
i = i + 1 
h[i] = '                    button3: zoom out & center on current position'
i = i + 1
h[i] = 'Blink:           sets blink mode:'
i = i + 1
h[i] = '                    press mouse button in main window to show blink image'
i = i + 1
h[i] = 'ImExam:          sets ImageExamine mode:'
i = i + 1
h[i] = '                    button 1: photometry'
i = i + 1
h[i] = '                    button 2: center on current position'
i = i + 1
h[i] = '                    button 3: image statistics'
i = i + 1
h[i] = 'Vector:          sets vector mode: click to plot pixel values along a vector'
i = i + 2
h[i] = 'BUTTONS:'
i = i + 1
h[i] = 'Invert:          inverts the current color table'
i = i + 1
h[i] = 'Restretch:       sets min and max to preserve display colors while linearizing the color table'
i = i + 1
h[i] = 'AutoScale:       sets min and max to show data values around image median'
i = i + 1
h[i] = 'FullRange:       sets min and max to show the full data range of the image'
i = i + 1
h[i] = 'ZoomIn:          zooms in by x2'
i = i + 1
h[i] = 'ZoomOut:         zooms out by x2'
i = i + 1
h[i] = 'Zoom1:           sets zoom level to original scale'
i = i + 1
h[i] = 'Center:          centers image on display window'
i = i + 1
h[i] = 'Done:            quits atv'
i = i + 1
h[i] = ''
i = i + 1
h[i] = 'Keyboard commands in display window:'
i = i + 1
h[i] = '    Arrow keys or numeric keypad (with NUM LOCK on) moves cursor'
i = i + 1
h[i] = '    r: row plot'
i = i + 1
h[i] = '    c: column plot'
i = i + 1
h[i] = '    s: surface plot'
i = i + 1
h[i] = '    t: contour plot'
i = i + 1
h[i] = '    p: aperture photometry at current position'
i = i + 1
h[i] = '    i: image statistics at current position'
i = i + 1
h[i] = '    m: cycles through mouse modes'
i = i + 1
h[i] = '    Shift-1,2,3:  sets blink buffer 1, 2, or 3'
i = i + 1
h[i] = '    b,n: When viewing an image cube, move Back or to Next image'
i = i + 1
h[i] = '    +,-: Zoom in or out.'
i = i + 1
h[i] = '    a,f,z,l,g: Change image scale or stretch. '
i = i + 1
h[i] = '    q: quits atv'
i = i + 2
h[i] = 'IDL COMMAND LINE HELP:'
i = i + 1
h[i] =  'To pass an array to atv:'
i = i + 1
h[i] =  '   atv, array_name [, options]'
i = i + 1
h[i] = 'To pass a fits filename to atv:'
i = i + 1
h[i] = '    atv, fitsfile_name [, options] (enclose filename in single quotes) '
i = i + 1
h[i] = 'Command-line options are: '
i = i + 1
h[i]  = '   [,min = min_value] [,max=max_value] [,/linear] [,/log] [,/histeq] [,/asinh]'
i = i + 1
h[i]  = '   [,/block] [,/align] [,/stretch] [,header=header]'
i = i + 2
h[i] = 'To overplot a contour plot on the draw window:'
i = i + 1
h[i] = '    atvcontour, array_name [, options...]'
i = i + 1
h[i] = 'To overplot text on the draw window: '
i = i + 1
h[i] = '    atvxyouts, x, y, text_string [, options]  (enclose string in single quotes)'
i = i + 1
h[i] = 'To overplot points or lines on the current plot:'
i = i + 1
h[i] = '    atvplot, xvector, yvector [, options]'
i = i + 2
h[i] = 'The options for atvcontour, atvxyouts, and atvplot are essentially'
i = i + 1
h[i] =  'the same as those for the idl contour, xyouts, and plot commands,'
i = i + 1
h[i] = 'except that data coordinates are always used.' 
i = i + 1
h[i] = 'The default color for overplots is red.'
i = i + 2
h[i] = 'The lowest 8 entries in the color table are:'
i = i + 1
h[i] = '    0 = black'
i = i + 1
h[i] = '    1 = red'
i = i + 1
h[i] = '    2 = green'
i = i + 1
h[i] = '    3 = blue'
i = i + 1
h[i] = '    4 = cyan'
i = i + 1
h[i] = '    5 = magenta'
i = i + 1
h[i] = '    6 = yellow'
i = i + 1
h[i] = '    7 = white'
i = i + 1
h[i] = '    The top entry in the color table is also reserved for white. '
i = i + 2
h[i] = 'Other commands:'
i = i + 1
h[i] = 'atverase [, N]:       erases all (or last N) plots and text'
i = i + 1
h[i] = 'atvclear: displays a small blank image (can be useful to clear memory)'
i = i + 1
h[i] = 'atv_activate: reanimates a frozen atv if another idl program crashes or hits a stop'
i = i + 1
h[i] = 'atv_shutdown:   quits atv'
i = i + 2
h[i] = 'NOTE: If atv should crash, type atv_shutdown at the idl prompt.'
i = i + 5
h[i] = strcompress('ATV.PRO version '+state.version)
i = i + 1
h[i] = 'For full instructions, or to download the most recent version, go to:'
i = i + 1
h[i] = 'http://www.physics.uci.edu/~barth/atv'


if (not (xregistered('atv_help', /noshow))) then begin

helptitle = strcompress('atv v' + state.version + ' help')

    help_base =  widget_base(group_leader = state.base_id, $
                             /column, $
                             /base_align_right, $
                             title = helptitle, $
                             uvalue = 'help_base')

    help_text = widget_text(help_base, $
                            /scroll, $
                            value = h, $
                            xsize = 85, $
                            ysize = 24)
    
    help_done = widget_button(help_base, $
                              value = 'Done', $
                              uvalue = 'help_done')

    widget_control, help_base, /realize
    xmanager, 'atv_help', help_base, /no_block
    
endif

end

;----------------------------------------------------------------------

pro atv_help_event, event

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'help_done': widget_control, event.top, /destroy
    else:
endcase

end

;----------------------------------------------------------------------
;      Routines for displaying image statistics
;----------------------------------------------------------------------

pro atv_stats_refresh

; Calculate box statistics and update the results

common atv_state
common atv_images

b = round((state.statboxsize - 1) / 2)

xmin = 0 > (state.cursorpos[0] - b) < (state.image_size[0] - 1)
xmax = 0 > (state.cursorpos[0] + b) < (state.image_size[0] - 1)
ymin = 0 > (state.cursorpos[1] - b) < (state.image_size[1] - 1)
ymax = 0 > (state.cursorpos[1] + b) < (state.image_size[1] - 1)

xmin = round(xmin)
xmax = round(xmax)
ymin = round(ymin)
ymax = round(ymax)

cut = float(main_image[xmin:xmax, ymin:ymax])
npix = (xmax - xmin + 1) * (ymax - ymin + 1)

cutmin = min(cut, max=maxx, /nan)
cutmax = maxx
cutmean = mean(cut, /nan, /double)
cutmedian = median(cut)
cutstddev = stddev(cut, /nan, /double)


widget_control, state.statbox_id, set_value=state.statboxsize
widget_control, state.statxcenter_id, set_value = state.cursorpos[0]
widget_control, state.statycenter_id, set_value = state.cursorpos[1]
tmp_string = strcompress('# Pixels in Box:  ' + string(npix))
widget_control, state.stat_npix_id, set_value = tmp_string
tmp_string = strcompress('Min:  ' + string(cutmin))
widget_control, state.statbox_min_id, set_value = tmp_string
tmp_string = strcompress('Max:  ' + string(cutmax))
widget_control, state.statbox_max_id, set_value = tmp_string
tmp_string = strcompress('Mean:  ' + string(cutmean))
widget_control, state.statbox_mean_id, set_value = tmp_string
tmp_string = strcompress('Median:  ' + string(cutmedian))
widget_control, state.statbox_median_id, set_value = tmp_string
tmp_string = strcompress('StdDev:  ' + string(cutstddev))
widget_control, state.statbox_stdev_id, set_value = tmp_string

atv_tvstats

end

;----------------------------------------------------------------------

pro atv_stats_event, event

common atv_state
common atv_images

widget_control, event.id, get_uvalue = uvalue

case uvalue of

    'statbox': begin
        state.statboxsize = long(event.value) > 3
        if ( (state.statboxsize / 2 ) EQ $
             round(state.statboxsize / 2.)) then $
          state.statboxsize = state.statboxsize + 1
        atv_stats_refresh
    end

    'statxcenter': begin
        state.cursorpos[0] = 0 > long(event.value) < (state.image_size[0] - 1)
        atv_stats_refresh
    end

    'statycenter': begin
        state.cursorpos[1] = 0 > long(event.value) < (state.image_size[1] - 1)
        atv_stats_refresh
    end

    'statxyresize': begin
      widget_control, state.stat_xyresize_button_id
      if (state.stat_xyresize eq 1) then state.stat_xyresize = 0 $
      else state.stat_xyresize = 1
    end

    'showstatzoom': begin
        widget_control, state.showstatzoom_id, get_value=val
        case val of
            'Show Region': begin
                widget_control, state.statzoom_widget_id, $
                  xsize=state.statzoom_size, ysize=state.statzoom_size
                widget_control, state.showstatzoom_id, $
                  set_value='Hide Region'
            end
            'Hide Region': begin
                widget_control, state.statzoom_widget_id, $
                  xsize=1, ysize=1
                widget_control, state.showstatzoom_id, $
                  set_value='Show Region'
             end
         endcase
         atv_stats_refresh
    end
  
    'stats_hist': begin

       x1 = 0 > (state.cursorpos[0] - (state.statbox_size/2))
       x2 = n_elements(main_image[*,0])-1 < $ 
            (state.cursorpos[0] + (state.statbox_size/2))
       y1 = 0 > (state.cursorpos[1] - (state.statbox_size/2))
       y2 = n_elements(main_image[0,*])-1 < $
            (state.cursorpos[1] + (state.statbox_size/2))

       if (not (xregistered('atv_lineplot', /noshow))) then begin
         atv_lineplot_init
       endif

       state.plot_coord = state.cursorpos
       widget_control, state.x1_pix_id, set_value=x1
       widget_control, state.x2_pix_id, set_value=x2
       widget_control, state.y1_pix_id, set_value=y1
       widget_control, state.y2_pix_id, set_value=y2
       hist_image = main_image[x1:x2,y1:y2]

       state.lineplot_xmin = min(hist_image)
       state.lineplot_xmax = max(hist_image)
       state.lineplot_ymin = 0.
       state.binsize = (state.lineplot_xmax - state.lineplot_xmin) * 0.01
       state.binsize = state.binsize > $
         ((state.lineplot_xmax - state.lineplot_xmin) * 1.0e-5)
       state.binsize = fix(state.binsize)

       ;Set plot window before calling plothist to get histogram ranges
       atv_setwindow, state.lineplot_window_id
       erase

       plothist, hist_image, xhist, yhist, bin=state.binsize, /NaN, /nodata

       state.lineplot_ymax = max(yhist) + (0.05*max(yhist))

       widget_control, state.lineplot_xmin_id, $
         set_value = state.lineplot_xmin

       widget_control, state.lineplot_xmax_id, $
         set_value = state.lineplot_xmax

       widget_control, state.lineplot_ymin_id, $
         set_value = state.lineplot_ymin

       widget_control, state.lineplot_ymax_id, $
         set_value = state.lineplot_ymax

       widget_control, state.histplot_binsize_id, set_value=state.binsize

       atv_histplot, /update

    end

    'stats_save': begin
       stats_outfile = dialog_pickfile(filter='*.txt', $
                        file='atv_stats.txt', get_path = tmp_dir, $
                        title='Please Select File to Append Stats')

       IF (strcompress(stats_outfile, /remove_all) EQ '') then RETURN

       IF (stats_outfile EQ tmp_dir) then BEGIN
         atv_message, 'Must indicate filename to save.', $
                       msgtype = 'error', /window
         return
       ENDIF

       openw, lun, stats_outfile, /get_lun, /append

       widget_control, state.stat_npix_id, get_value = npix_str
       widget_control, state.statbox_min_id, get_value = minstat_str
       widget_control, state.statbox_max_id, get_value = maxstat_str
       widget_control, state.statbox_mean_id, get_value = meanstat_str
       widget_control, state.statbox_median_id, get_value = medianstat_str
       widget_control, state.statbox_stdev_id, get_value = stdevstat_str

       printf, lun, 'ATV IMAGE BOX STATISTICS--NOTE: IDL Arrays Begin With Index 0'
       printf, lun, '============================================================='
       printf, lun, ''
       printf, lun, 'Image Name: ' + strcompress(state.imagename,/remove_all)
       printf, lun, 'Image Size: ' + strcompress(string(state.image_size[0]) $
                           + ' x ' + string(state.image_size[1]))
 
       printf, lun, 'Image Min: ' + $
                    strcompress(string(state.image_min),/remove_all)
       printf, lun, 'Image Max: ' + $
                    strcompress(string(state.image_max),/remove_all)

       if (state.image_size[2] gt 1) then printf, lun, 'Image Plane: ' + $
                    strcompress(string(state.cur_image_num),/remove_all)

       printf, lun, ''
       printf, lun, 'Selected Box Statistics:'
       printf, lun, '------------------------'
     
       printf, lun, 'X-Center: ' + strcompress(string(state.cursorpos[0]), $
                /remove_all)
       printf, lun, 'Y-Center: ' + strcompress(string(state.cursorpos[1]), $
                /remove_all)
       printf, lun, 'Xmin: ' + $
              strcompress(string(state.cursorpos[0] - state.xstatboxsize/2), $
                /remove_all)
       printf, lun, 'Xmax: ' + $
              strcompress(string(state.cursorpos[0] + state.xstatboxsize/2), $
                /remove_all)
       printf, lun, 'Ymin: ' + $
              strcompress(string(state.cursorpos[1] - state.ystatboxsize/2), $
                /remove_all)
       printf, lun, 'Ymax: ' + $
              strcompress(string(state.cursorpos[1] + state.ystatboxsize/2), $
                /remove_all)

       printf, lun, npix_str
       printf, lun, minstat_str
       printf, lun, maxstat_str
       printf, lun, meanstat_str
       printf, lun, medianstat_str
       printf, lun, stdevstat_str
       printf, lun, ''

       close, lun
       free_lun, lun

    end

    'stats_done': widget_control, event.top, /destroy
    else:
endcase


end

;----------------------------------------------------------------------

pro atv_showstats

; Brings up a widget window for displaying image statistics

common atv_state
common atv_images

common atv_state

state.cursorpos = state.coord

if (not (xregistered('atv_stats', /noshow))) then begin

    stats_base = $
      widget_base(group_leader = state.base_id, $
                  /column, $
                  /base_align_center, $
                  title = 'atv image statistics', $
                  uvalue = 'stats_base')
    state.stats_base_id = stats_base
    
    stats_nbase = widget_base(stats_base, /row, /base_align_center)
    stats_base1 = widget_base(stats_nbase, /column, frame=1)
    stats_base2 = widget_base(stats_nbase, /column)
    stats_base2a = widget_base(stats_base2, /column, frame=1)
    stats_zoombase = widget_base(stats_base, /column)

    tmp_string = strcompress('Image size:  ' + $
                             string(state.image_size[0]) + $
                             ' x ' + $
                             string(state.image_size[1]))

    size_label = widget_label(stats_base1, value = tmp_string)

    tmp_string = strcompress('Image Min:  ' + string(state.image_min))
    min_label= widget_label(stats_base1, value = tmp_string)
    tmp_string = strcompress('Image Max:  ' + string(state.image_max))
    max_label= widget_label(stats_base1, value = tmp_string)

    state.statbox_id = $
      cw_field(stats_base1, $
               /long, $
               /return_events, $
               title = 'Box Size for Stats:', $
               uvalue = 'statbox', $
               value = state.statboxsize, $
               xsize = 5)
    
    state.statxcenter_id = $
      cw_field(stats_base1, $
               /long, $
               /return_events, $
               title = 'Box X Center:', $
               uvalue = 'statxcenter', $
               value = state.cursorpos[0], $ 
               xsize = 5)

    state.statycenter_id = $
      cw_field(stats_base1, $
               /long, $
               /return_events, $
               title = 'Box Y Center:', $
               uvalue = 'statycenter', $
               value = state.cursorpos[1], $ 
               xsize = 5)

    tmp_string = strcompress('# Pixels in Box:  ' + string(100000))
    state.stat_npix_id = widget_label(stats_base2a, value = tmp_string,$
      /align_left)
    tmp_string = strcompress('Min:  ' + '0.0000000000')
    state.statbox_min_id = widget_label(stats_base2a, value = tmp_string,$
      /align_left)
    tmp_string = strcompress('Max:  ' + '0.0000000000')
    state.statbox_max_id = widget_label(stats_base2a, value = tmp_string, $
      /align_left)
    tmp_string = strcompress('Mean:  ' + '0.0000000000')
    state.statbox_mean_id = widget_label(stats_base2a, value = tmp_string, $
      /align_left)
    tmp_string = strcompress('Median:  ' + '0.0000000000')
    state.statbox_median_id = widget_label(stats_base2a, value = tmp_string, $
      /align_left)
    tmp_string = strcompress('StdDev:  ' + '0.0000000000')
    state.statbox_stdev_id = widget_label(stats_base2a, value = tmp_string, $
      /align_left)
    
    state.showstatzoom_id = widget_button(stats_base2, $
          value = 'Hide Region', uvalue = 'showstatzoom')

    stat_hist = widget_button(stats_base2, value = 'Histogram Pixels', $
          uvalue = 'stats_hist')

    stat_save = widget_button(stats_base2, value = 'Save Stats', $
          uvalue = 'stats_save')

    stat_done = $
      widget_button(stats_base2, $
                    value = 'Done', $
                    uvalue = 'stats_done')
    
    state.statzoom_widget_id = widget_draw(stats_zoombase, $
       scr_xsize = state.statzoom_size, scr_ysize = state.statzoom_size)

    widget_control, stats_base, /realize
    
    xmanager, 'atv_stats', stats_base, /no_block
    
    widget_control, state.statzoom_widget_id, get_value = tmp_val
    state.statzoom_window_id = tmp_val

    atv_resetwindow

endif

atv_stats_refresh

end

;---------------------------------------------------------------------


pro atv_tvstats

; Routine to display the zoomed region around a stats point

common atv_state
common atv_images

atv_setwindow, state.statzoom_window_id
erase

x = round(state.cursorpos[0])
y = round(state.cursorpos[1])

boxsize = (state.statboxsize - 1) / 2
xsize = state.statboxsize
ysize = state.statboxsize
image = bytarr(xsize,ysize)

xmin = (0 > (x - boxsize))
xmax = ((x + boxsize) < (state.image_size[0] - 1) )
ymin = (0 > (y - boxsize) )
ymax = ((y + boxsize) < (state.image_size[1] - 1))

startx = abs( (x - boxsize) < 0 )
starty = abs( (y - boxsize) < 0 ) 

image[startx, starty] = scaled_image[xmin:xmax, ymin:ymax]

xs = indgen(xsize) + xmin - startx
ys = indgen(ysize) + ymin - starty

xs_delta = (xs[xsize-1] - xs[0]) / float(xsize - 1.0)
ys_delta = (ys[ysize-1] - ys[0]) / float(ysize - 1.0)
x_ran = [xs[0]-xs_delta/2.0,xs[xsize-1]+xs_delta/2.0]
y_ran = [ys[0]-ys_delta/2.0,ys[ysize-1]+ys_delta/2.0]

dev_width = 0.8 * state.statzoom_size
dev_pos = [0.15 * state.statzoom_size, $
           0.15 * state.statzoom_size, $
           0.95 * state.statzoom_size, $
           0.95 * state.statzoom_size]

x_factor = dev_width / xsize
y_factor = dev_width / ysize
x_offset = (x_factor - 1.0) / x_factor / 2.0
y_offset = (y_factor - 1.0) / y_factor / 2.0
xi = findgen(dev_width) / x_factor - x_offset ;x interp index
yi = findgen(dev_width) / y_factor - y_offset ;y interp index

image = Poly_2D(image, [[0,0],[1.0/x_factor,0]], $
             [[0,1.0/y_factor],[0,0]], $
             0, dev_width, dev_width)

xsize = (size(image))[1]
ysize = (size(image))[2]
out_xs = xi * xs_delta + xs[0]
out_ys = yi * ys_delta + ys[0]

sz = size(image)
xsize = Float(sz[1])       ;image width
ysize = Float(sz[2])       ;image height
dev_width = dev_pos[2] - dev_pos[0] + 1
dev_width = dev_pos[3] - dev_pos[1] + 1

tv, image, /device, dev_pos[0], dev_pos[1], $
  xsize=dev_pos[2]-dev_pos[0], $
  ysize=dev_pos[3]-dev_pos[1]

plot, [0, 1], /noerase, /nodata, xstyle = 1, ystyle = 1, $
  /device, position = dev_pos, color=7, $
  xrange = x_ran, yrange = y_ran, xtitle='X', ytitle='Y'

atv_resetwindow
end

;----------------------------------------------------------------------
;        aperture photometry and radial profile routines
;---------------------------------------------------------------------

pro atv_imcenterf, xcen, ycen

; program to calculate the center of mass of an image around
; the point (x,y), return the answer in (xcen,ycen).
;
; by M. Liu, adapted for inclusion in ATV by AJB
;
; ALGORITHM:
;   1. first finds max pixel value in
;	   a 'bigbox' box around the cursor
;   2. then calculates centroid around the object 
;   3. iterates, recalculating the center of mass 
;      around centroid until the shifts become smaller 
;      than MINSHIFT (0.3 pixels) 

common atv_images
common atv_state

; iteration controls
MINSHIFT = 0.3

; max possible x or y direction shift
MAXSHIFT = 3

; Bug fix 4/16/2000: added call to round to make sure bigbox is an integer
bigbox=round(1.5*state.centerboxsize)

sz = size(main_image)

; box size must be odd
dc = (state.centerboxsize-1)/2
if ( (bigbox / 2 ) EQ round(bigbox / 2.)) then bigbox = bigbox + 1
db = (bigbox-1)/2

; need to start with integers
xx = state.cursorpos[0]
yy = state.cursorpos[1]

; first find max pixel in box around the cursor
x0 = (xx-db) > 0
x1 = (xx+db) < (sz(1)-1)
y0 = (yy-db) > 0
y1 = (yy+db) < (sz(2)-1)
cut = main_image[x0:x1,y0:y1]
cutmax = max(cut)
w=where(cut EQ cutmax)
cutsize = size(cut)
my = (floor(w/cutsize[1]))[0]
mx = (w - my*cutsize[1])[0]

xx = mx + x0
yy = my + y0 
xcen = xx
ycen = yy

; then find centroid 
if  (n_elements(xcen) gt 1) then begin
    xx = round(total(xcen,/nan)/total(finite(xcen))) 
    yy = round(total(ycen,/nan)/total(finite(ycen)))
endif

done = 0
niter = 1
    
;	cut out relevant portion
sz = size(main_image)
x0 = round((xx-dc) > 0)		; need the ()'s
x1 = round((xx+dc) < (sz[1]-1))
y0 = round((yy-dc) > 0)
y1 = round((yy+dc) < (sz[2]-1))
xs = x1 - x0 + 1
ys = y1 - y0 + 1
cut = float(main_image[x0:x1, y0:y1])

; sky subtract before centering
; note that this is a quick and dirty solution, and may cause
; centering problems for poorly behaved data  -- AB, 2/28/07
cut = cut - min(cut)
                                ; find x position of center of mass
cenxx = fltarr(xs, ys, /nozero)
for i = 0L, (xs-1) do $         ; column loop
  cenxx[i, *] = cut[i, *] * i
xcen = total(cenxx,/NaN) / total(cut,/NaN) + x0

                                ; find y position of center of mass
cenyy = fltarr(xs, ys, /nozero)
for i = 0L, (ys-1) do $         ; row loop
  cenyy[*, i] = cut[*, i] * i
ycen = total(cenyy,/NaN) / total(cut,/NaN) + y0

if (finite(xcen,/NAN) or finite(ycen,/NAN)) then $
	state.photwarning = "ERROR: Could not find centroid! Too many NaNs?"
if (abs(xcen-state.cursorpos[0]) gt MAXSHIFT) or $
  (abs(ycen-state.cursorpos[1]) gt MAXSHIFT) then begin
    state.photwarning = 'Warning: Possible mis-centering?'
endif

end

;----------------------------------------------------------------------

function atv_splinefwhm, rad, prof, splrad, splprof

; given a radial profile (counts vs radius) will use
; a spline to extract the FWHM
;
; ALGORITHM
;   finds peak in radial profile, then marches along until finds
;   where radial profile has dropped to half of that,
;   assumes peak value of radial profile is at minimum radius
;
; original version by M. Liu, adapted for ATV by AJB

common atv_state

nrad = n_elements(rad)

; check the peak
w = where(prof eq max(prof),wcnt)
if wcnt lt 1 then begin
	message,/info, "Invalid radii!"
	state.photwarning  = "Radii are invalid! (too large?)"
	return,-1
endif
if float(rad(w[0])) ne min(rad) then begin
state.photwarning = 'Warning: Profile peak is off-center!'
  return,-1
endif

; interpolate radial profile at 50 times as many points
splrad = min(rad) + findgen(nrad*50+1) * (max(rad)-min(rad)) / (nrad*50)
nspl = n_elements(splrad)

; spline the profile
splprof = spline(rad,prof,splrad)

; march along splined profile until cross 0.5*peak value
found = 0
i = 0
repeat begin
  if splprof(i) lt 0.5*max(splprof) then $
	found = 1 $
  else $
	i = i+1
endrep until ((found) or (i eq nspl))

if (i lt 2) or (i eq nspl) then begin
state.photwarning = 'Warning: Unable to measure FWHM!'
  return,-1
endif

; now interpolate across the 2 points straddling the 0.5*peak
fwhm = splrad(i)+splrad(i-1)

return,fwhm
end

;-----------------------------------------------------------------------

pro atv_radplotf, x, y, fwhm

; Program to calculate radial profile of an image
; given aperture location, range of sizes, and inner and 
; outer radius for sky subtraction annulus.  Calculates sky by
; median.
; 
; original version by M. Liu, adapted for inclusion in ATV by AJB

common atv_state
common atv_images

; set defaults
inrad = 0.5*sqrt(2)
outrad = round(state.outersky * 1.2)
drad=1.
insky = outrad+drad
outsky = insky+drad+20.

; initialize arrays
inrad = float(inrad)
outrad = float(outrad)
drad = float(drad)
nrad = ceil((outrad-inrad)/drad) + 1
out = fltarr(nrad,12)

; extract relevant image subset (may be rectangular), translate coord origin,
;   bounded by edges of image
;   (there must be a cute IDL way to do this neater)
sz = size(main_image)
x0 = floor(x-outsky) 
x1 = ceil(x+outsky)   ; one pixel too many?
y0 = floor(y-outsky) 
y1 = ceil(y+outsky)
x0 = x0 > 0.0
x1 = x1 < (sz[1]-1)
y0 = y0 > 0.0
y1 = y1 < (sz[2]-1)
nx = x1 - x0 + 1
ny = y1 - y0 + 1

; trim the image, translate coords
img = main_image[x0:x1,y0:y1]
xcen = x - x0
ycen = y - y0

; for debugging, can make some masks showing different regions
skyimg = fltarr(nx,ny)			; don't use /nozero!!
photimg = fltarr(nx,ny)			; don't use /nozero!!

; makes an array of (distance)^2 from center of aperture
;   where distance is the radial or the semi-major axis distance.
;   based on DIST_CIRCLE and DIST_ELLIPSE in Goddard IDL package,
;   but deals with rectangular image sections
distsq = fltarr(nx,ny,/nozero)

xx = findgen(nx)
yy = findgen(ny)
x2 = (xx - xcen)^(2.0)
y2 = (yy - ycen)^(2.0)
for i = 0L,(ny-1) do $          ; row loop
  distsq[*,i] = x2 + y2(i)

; get sky level by masking and then medianing remaining pixels
; note use of "gt" to avoid picking same pixels as flux aperture
ns = 0
msky = 0.0
errsky = 0.0

in2 = insky^(2.0)
out2 = outsky^(2.0)
if (in2 LT max(distsq)) then begin
    w = where((distsq gt in2) and (distsq le out2),ns)
    skyann = img[w] 
endif else begin
    w = where(distsq EQ distsq)
    skyann = img[w]
    state.photwarning = 'Not enough pixels in sky!'
endelse

msky = median(skyann)
errsky = stddev(skyann,/NaN)
skyimg[w] = -5.0
photimg = skyimg

errsky2 = errsky * errsky

out[*,8] = msky
out[*,9] = ns
out[*,10]= errsky

; now loop through photometry radii, finding the total flux, differential
;	flux, and differential average pixel value along with 1 sigma scatter
; 	relies on the fact the output array is full of zeroes
for i = 0,nrad-1 do begin
    
    dr = drad
    if i eq 0 then begin
        rin =  0.0
        rout = inrad
        rin2 = -0.01
    endif else begin
        rin = inrad + drad *(i-1)	
        rout = (rin + drad) < outrad
        rin2 = rin*rin
    endelse
    rout2 = rout*rout
    
; 	get flux and pixel stats in annulus, wary of counting pixels twice
;	checking if necessary if there are pixels in the sector
    w = where(distsq gt rin2 and distsq le rout2,np)
    
    pfrac = 1.0                 ; fraction of pixels in each annulus used
    
    if np gt 0 then begin
        ann = img[w]
		; deal with NaN pixels
		pfrac = total(finite(ann))/n_elements(ann) 
        dflux = total(ann,/NaN) * 1./pfrac
        dnpix = np
        dnet = dflux - (dnpix * msky) * 1./pfrac
        davg = dnet / (dnpix * 1./pfrac)
        if (np gt 1) and (total(finite(ann)) gt 1) then dsig = stddev(ann,/NaN) else dsig = 0.00
        
;		std dev in each annulus including sky sub error
        derr = sqrt(dsig*dsig + errsky2)
        
        photimg[w] = rout2
        
        out[i,0] = (rout+rin)/2.0
        out[i,1] = out[i-1>0,1] + dflux
        out[i,2] = out[i-1>0,2] + dnet
        out[i,3] = out[i-1>0,3] + dnpix
        out[i,4] = dflux
        out[i,5] = dnpix
        out[i,6] = davg
        out[i,7] = dsig
        out[i,11] = derr
    endif else if (i ne 0) then begin
        out[i,0]= rout
        out[i,1:3] = out[i-1,1:3]
        out[i, 4:7] = 0.0
        out[i,11] = 0.0
    endif else begin
        out[i, 0] = rout
    endelse
    
endfor

; fill radpts array after done with differential photometry
w = where(distsq ge 0.0 and distsq le outrad*outrad)
radpts = dblarr(2,n_elements(w))
radpts[0,*] = sqrt(distsq[w])
radpts[1,*] = img[w]

; compute FWHM via spline interpolation of radial profile
fwhm = atv_splinefwhm(out[*,0],out[*,6])

; plot the results

if n_elements(radpts(1, *)) gt 100 then pp = 3 else pp = 1

yminpoint = msky
ymaxpoint = max(radpts[1,*])
blankspace=0.08
ymin = yminpoint - blankspace*(ymaxpoint - yminpoint)
ymax = ymaxpoint + blankspace*(ymaxpoint - yminpoint)

if keyword_set(state.photplotmode) then ymin = min(radpts[1, *],/nan)

plot, radpts[0, *], radpts[1, *], /nodata, xtitle = 'Radius (pixels)', $
  ytitle = 'Counts', color=7, charsize=1.2, yrange = [ymin,ymax], yst=1, $
  ylog = keyword_set(state.photplotmode)
oplot, radpts[0, *], radpts[1, *], psym = pp, color=6
oploterror, out[*, 0], out[*, 6]+out[*, 8], $
  out[*, 11]/sqrt(out[*, 5]), psym=-4, color=7, errcolor=7

end

;-----------------------------------------------------------------------

pro atv_apphot_refresh

; Do aperture photometry using idlastro daophot routines.

common atv_state
common atv_images

state.photwarning = 'Warnings: None.'

; Center on the object position nearest to the cursor
if (state.centerboxsize GT 0) then begin
    atv_imcenterf, x, y
endif else begin   ; no centering
    x = state.cursorpos[0]
    y = state.cursorpos[1]
endelse

; Make sure that object position is on the image
x = 0 > x < (state.image_size[0] - 1)
y = 0 > y < (state.image_size[1] - 1)

if ((x - state.outersky) LT 0) OR $
  ((x + state.outersky) GT (state.image_size[0] - 1)) OR $
  ((y - state.outersky) LT 0) OR $
  ((y + state.outersky) GT (state.image_size[1] - 1)) then $
  state.photwarning = 'Warning: Sky apertures fall outside image!'

; Condition to test whether phot aperture is off the image
if (x LT state.aprad) OR $
  ((state.image_size[0] - x) LT state.aprad) OR $
  (y LT state.aprad) OR $
  ((state.image_size[1] - y) LT state.aprad) then begin
    flux = -1.
    state.photwarning = 'Warning: Aperture Outside Image Border!'
endif
    
phpadu = 1.0                    ; don't convert counts to electrons
apr = [state.aprad]
skyrad = [state.innersky, state.outersky]
; Assume that all pixel values are good data
badpix = [state.image_min-1, state.image_max+1]  

if (state.skytype EQ 1) then begin    ; calculate median sky value

    xmin = (x - state.outersky) > 0
    xmax = (xmin + (2 * state.outersky + 1)) < (state.image_size[0] - 1)
    ymin = (y - state.outersky) > 0
    ymax = (ymin + (2 * state.outersky + 1)) < (state.image_size[1] - 1)
    
    small_image = main_image[xmin:xmax, ymin:ymax]
    nx = (size(small_image))[1]
    ny = (size(small_image))[2]
    i = lindgen(nx)#(lonarr(ny)+1)
    j = (lonarr(nx)+1)#lindgen(ny)
    xc = x - xmin
    yc = y - ymin
    
    w = where( (((i - xc)^2 + (j - yc)^2) GE state.innersky^2) AND $
               (((i - xc)^2 + (j - yc)^2) LE state.outersky^2),  nw)
    
    if ((x - state.outersky) LT 0) OR $
      ((x + state.outersky) GT (state.image_size[0] - 1)) OR $
      ((y - state.outersky) LT 0) OR $
      ((y + state.outersky) GT (state.image_size[1] - 1)) then $
      state.photwarning = 'Warning: Sky apertures fall outside image!'
    
    if (nw GT 0) then  begin
        skyval = median(small_image(w)) 
    endif else begin
        skyval = -1
        state.photwarning = 'Warning: No pixels in sky!'
    endelse
endif

; Do the photometry now
case state.skytype of
    0: aper, main_image, [x], [y], flux, errap, sky, skyerr, phpadu, apr, $
      skyrad, badpix, flux=abs(state.magunits-1), /silent
    1: aper, main_image, [x], [y], flux, errap, sky, skyerr, phpadu, apr, $
      skyrad, badpix, flux=abs(state.magunits-1), /silent, $
      setskyval = skyval
    2: aper, main_image, [x], [y], flux, errap, sky, skyerr, phpadu, apr, $
      skyrad, badpix, flux=abs(state.magunits-1), /silent, $
      setskyval = 0
endcase

flux = flux[0]
sky = sky[0]

if (flux EQ 99.999) then begin
    state.photwarning = 'Warning: Error in computing flux!'
    flux = -1.0
endif

if (state.magunits EQ 1) then begin    ; apply zeropoint
    flux = flux + state.photzpt - 25.0
endif

; Run atv_radplotf and plot the results

atv_setwindow, state.radplot_window_id
atv_radplotf, x, y, fwhm

if state.photplotmode eq 0 then begin
	; linear
	range = !y.crange 
	labels = range[0]+(range[1]-range[0])*[0.92,0.85,0.75]
endif else begin 
	; logarithmic
	range = 10^!y.crange
	labels = 10^(!y.crange[0]+(!y.crange[1]-!y.crange[0])*[0.92,0.85,0.75])
endelse
if (not (keyword_set(ps))) then begin
	color1 = 2
	color2 = 4
	color3 = 5
	color4 = 1
endif else begin
	color1 = 0
	color2 = 0
	color3 = 0
	color4 = 0
endelse

; overplot the phot apertures on radial plot
plots, [state.aprad, state.aprad], !y.crange, line = 1, color=2, thick=2, psym=0

ymin = !y.crange(0)
ymax = !y.crange(1)
ypos = ymin + 0.85*(ymax-ymin)
xyouts, /data, state.aprad, ypos, ' aprad', $
  color=2, charsize=1.5
if (state.skytype NE 2) then begin
    plots, [state.innersky,state.innersky], range, $
      line = 1, color=color2, thick=2, psym=0
    xyouts, /data, state.innersky, labels[1], ' insky', $
      color=color2, charsize=1.5
    plots, [state.outersky,state.outersky], range, $
      line = 1, color=color3, thick=2, psym=0
    xyouts, /data, state.outersky * 0.82, labels[2], ' outsky', $
      color=color3, charsize=1.5
endif
plots, !x.crange, [sky, sky], color=color4, thick=2, psym=0, line = 2
xyouts, /data, state.innersky + (0.1*(state.outersky-state.innersky)), $
  sky+0.06*(!y.crange[1] - sky), 'sky level', color=color4, charsize=1.5

atv_resetwindow

; output the results

case state.magunits of
    0: begin
         fluxstr = 'Object counts: '
         errstr = 'Counts error: '
     end
    1: begin
         fluxstr = 'Magnitude: '
         errstr = 'Mag error: '
     end
 else:
endcase
  
state.centerpos = [x, y]

tmp_string = string(state.cursorpos[0], state.cursorpos[1], $
                    format = '("Cursor position:  x=",i4,"  y=",i4)' )
tmp_string1 = string(state.centerpos[0], state.centerpos[1], $
                    format = '("Object centroid:  x=",f6.1,"  y=",f6.1)' )
tmp_string2 = strcompress(fluxstr+string(flux, format = '(g12.6)' ))
tmp_string3 = string(sky, format = '("Sky level: ",g12.6)' )
tmp_string4 = string(fwhm, format='("FWHM (pix): ",g7.3)' )
tmp_string5 = strcompress(errstr+string(errap, format = '(g12.6)' ))

widget_control, state.centerbox_id, set_value = state.centerboxsize
widget_control, state.cursorpos_id, set_value = tmp_string
widget_control, state.centerpos_id, set_value = tmp_string1
widget_control, state.radius_id, set_value = state.aprad 
widget_control, state.outersky_id, set_value = state.outersky
widget_control, state.innersky_id, set_value = state.innersky
widget_control, state.skyresult_id, set_value = tmp_string3
widget_control, state.photresult_id, set_value = tmp_string2
widget_control, state.photerror_id, set_value = tmp_string5
widget_control, state.fwhm_id, set_value = tmp_string4
widget_control, state.photwarning_id, set_value=state.photwarning

; Uncomment next lines if you want atv to output the WCS coords of 
; the centroid for the photometry object:
;if (state.wcstype EQ 'angle') then begin
;    xy2ad, state.centerpos[0], state.centerpos[1], *(state.astr_ptr), $
;      clon, clat
;    wcsstring = atv_wcsstring(clon, clat, (*state.astr_ptr).ctype,  $
;                state.equinox, state.display_coord_sys, state.display_equinox)
;    print, 'Centroid WCS coords: ', wcsstring
;endif

atv_tvphot

atv_resetwindow
end

;----------------------------------------------------------------------

pro atv_tvphot

; Routine to display the zoomed region around a photometry point,
; with circles showing the photometric apterture and sky radii.

common atv_state
common atv_images

atv_setwindow, state.photzoom_window_id
erase

x = round(state.centerpos[0])
y = round(state.centerpos[1])

boxsize = round(state.outersky * 1.2)
xsize = (2 * boxsize) + 1
ysize = (2 * boxsize) + 1
image = bytarr(xsize,ysize)

xmin = (0 > (x - boxsize))
xmax = ((x + boxsize) < (state.image_size[0] - 1) )
ymin = (0 > (y - boxsize) )
ymax = ((y + boxsize) < (state.image_size[1] - 1))

startx = abs( (x - boxsize) < 0 )
starty = abs( (y - boxsize) < 0 ) 

image[startx, starty] = scaled_image[xmin:xmax, ymin:ymax]

xs = indgen(xsize) + xmin - startx
ys = indgen(ysize) + ymin - starty

xs_delta = (xs[xsize-1] - xs[0]) / float(xsize - 1.0)
ys_delta = (ys[ysize-1] - ys[0]) / float(ysize - 1.0)
x_ran = [xs[0]-xs_delta/2.0,xs[xsize-1]+xs_delta/2.0]
y_ran = [ys[0]-ys_delta/2.0,ys[ysize-1]+ys_delta/2.0]

dev_width = 0.8 * state.photzoom_size
dev_pos = [0.15 * state.photzoom_size, $
           0.15 * state.photzoom_size, $
           0.95 * state.photzoom_size, $
           0.95 * state.photzoom_size]

x_factor = dev_width / xsize
y_factor = dev_width / ysize
x_offset = (x_factor - 1.0) / x_factor / 2.0
y_offset = (y_factor - 1.0) / y_factor / 2.0
xi = findgen(dev_width) / x_factor - x_offset ;x interp index
yi = findgen(dev_width) / y_factor - y_offset ;y interp index

image = Poly_2D(image, [[0,0],[1.0/x_factor,0]], $
             [[0,1.0/y_factor],[0,0]], $
             0, dev_width, dev_width)

xsize = (size(image))[1]
ysize = (size(image))[2]
out_xs = xi * xs_delta + xs[0]
out_ys = yi * ys_delta + ys[0]

sz = size(image)
xsize = Float(sz[1])       ;image width
ysize = Float(sz[2])       ;image height
dev_width = dev_pos[2] - dev_pos[0] + 1
dev_width = dev_pos[3] - dev_pos[1] + 1


tv, image, /device, dev_pos[0], dev_pos[1], $
  xsize=dev_pos[2]-dev_pos[0], $
  ysize=dev_pos[3]-dev_pos[1]

plot, [0, 1], /noerase, /nodata, xstyle = 1, ystyle = 1, $
  /device, position = dev_pos, color=7, $
  xrange = x_ran, yrange = y_ran

tvcircle, /data, state.aprad, state.centerpos[0], state.centerpos[1], $
  color=2, thick=2, psym=0
if (state.skytype NE 2) then begin
    tvcircle, /data, state.innersky, state.centerpos[0], state.centerpos[1], $
      color=4, thick=2, psym=0
    tvcircle, /data, state.outersky, state.centerpos[0], state.centerpos[1], $
      color=5, thick=2, psym=0
endif

atv_resetwindow
end

;----------------------------------------------------------------------

pro atv_apphot_event, event

common atv_state
common atv_images

widget_control, event.id, get_uvalue = uvalue

case uvalue of

    'centerbox': begin
        if (event.value EQ 0) then begin
            state.centerboxsize = 0
        endif else begin
            state.centerboxsize = long(event.value) > 3
            if ( (state.centerboxsize / 2 ) EQ $
                 round(state.centerboxsize / 2.)) then $
              state.centerboxsize = state.centerboxsize + 1
        endelse
        atv_apphot_refresh
    end
        
    'radius': begin
        state.aprad = 1.0 > event.value < state.innersky
        atv_apphot_refresh
    end

    'innersky': begin
        state.innersky = state.aprad > event.value < (state.outersky - 1)
        state.innersky = 2 > state.innersky
        if (state.outersky EQ state.innersky + 1) then $
          state.outersky = state.outersky + 1
        atv_apphot_refresh
    end

    'outersky': begin
        state.outersky = event.value > (state.innersky + 2)
        atv_apphot_refresh
    end

    'showradplot': begin
        widget_control, state.showradplot_id, get_value=val
        case val of
            'Show Radial Profile': begin
                ysize = 350 < (state.screen_ysize - 350)
                widget_control, state.radplot_widget_id, $
                  xsize=500, ysize=ysize
                widget_control, state.showradplot_id, $
                  set_value='Hide Radial Profile'
            end
            'Hide Radial Profile': begin
                widget_control, state.radplot_widget_id, $
                  xsize=1, ysize=1
                widget_control, state.showradplot_id, $
                  set_value='Show Radial Profile'
             end
         endcase
         atv_apphot_refresh
    end

    'magunits': begin
        state.magunits = event.value
        atv_apphot_refresh
    end

    'photsettings': atv_apphot_settings

    'radplot_stats_save': begin
        radplot_stats_outfile = dialog_pickfile(filter='*.txt', $
                        file='atv_phot.txt', get_path = tmp_dir, $
                        title='Please Select File to Append Photometry Stats')

        IF (strcompress(radplot_stats_outfile, /remove_all) EQ '') then RETURN

        IF (radplot_stats_outfile EQ tmp_dir) then BEGIN
          atv_message, 'Must indicate filename to save.', $
                       msgtype = 'error', /window
          return
        ENDIF

        openw, lun, radplot_stats_outfile, /get_lun, /append

        widget_control, state.cursorpos_id, get_value = cursorpos_str
        widget_control, state.centerbox_id, get_value = centerbox_str
        widget_control, state.centerpos_id, get_value = centerpos_str
        widget_control, state.radius_id, get_value = radius_str
        widget_control, state.innersky_id, get_value = innersky_str
        widget_control, state.outersky_id, get_value = outersky_str
        widget_control, state.fwhm_id ,get_value = fwhm_str
        widget_control, state.skyresult_id, get_value = skyresult_str
        widget_control, state.photresult_id, get_value = objectcounts_str
        widget_control, state.photerror_id, get_value = error_str

        printf, lun, 'ATV PHOTOMETRY RESULTS--NOTE: IDL Arrays Begin With Index 0'
        printf, lun, '============================================================='
        if (state.image_size[2] gt 1) then printf, lun, 'Image Plane: ' + $
                     strcompress(string(state.cur_image_num),/remove_all)
        printf, lun, strcompress(cursorpos_str)
        printf, lun, 'Centering box size (pix): ' + $
                      strcompress(string(centerbox_str),/remove_all)
        printf, lun, strcompress(centerpos_str)
        printf, lun, 'Aperture radius: ' + $
                     strcompress(string(radius_str), /remove_all)
        printf, lun, 'Inner sky radius: ' + $
                     strcompress(string(innersky_str), /remove_all)
        printf, lun, 'Outer sky radius: ' + $
                     strcompress(string(outersky_str), /remove_all)
        printf, lun, strcompress(fwhm_str)
        printf, lun, strcompress(skyresult_str)
        printf, lun, objectcounts_str
        printf, lun, error_str
        printf, lun, ''

        close, lun
        free_lun, lun
    end

    'apphot_ps': begin

        fname = strcompress(state.current_dir + 'atv_phot.ps', /remove_all)
        forminfo = cmps_form(cancel = canceled, create = create, $
                     /preserve_aspect, $
                     /color, $
                     /nocommon, papersize='Letter', $
                     filename = fname, $
                     button_names = ['Create PS File'])

        if (canceled) then return
        if (forminfo.filename EQ '') then return

        tmp_result = findfile(forminfo.filename, count = nfiles)

        result = ''
        if (nfiles GT 0) then begin
          mesg = strarr(2)
          mesg[0] = 'Overwrite existing file:'
          tmp_string = strmid(forminfo.filename, strpos(forminfo.filename, $
                              '/') + 1)
          mesg[1] = strcompress(tmp_string + '?', /remove_all)
          result =  dialog_message(mesg, $
                             /default_no, $
                             dialog_parent = state.base_id, $
                             /question)                 
        endif

        if (strupcase(result) EQ 'NO') then return
    
        widget_control, /hourglass

        screen_device = !d.name

        set_plot, 'ps'
        device, _extra = forminfo

        atv_apphot_refresh, /ps

        device, /close
        set_plot, screen_device

    end

    'apphot_done': widget_control, event.top, /destroy
    else:
endcase

end

;----------------------------------------------------------------------

pro atv_apphot_settings

; Routine to get user input on various photometry settings

common atv_state

skyline = strcompress('0, button, IDLPhot Sky|Median Sky|No Sky Subtraction,'+$
                      'exclusive,' + $
                      'label_left=Select Sky Algorithm: , set_value = ' + $
                      string(state.skytype))

magline = strcompress('0, button, Counts|Magnitudes, exclusive,' + $
                      'label_left = Select Output Units: , set_value =' + $
                      string(state.magunits))

plotline = strcompress('0, button, Linear|Log, exclusive,' + $
                      'label_left = Select Plot Scale: , set_value =' + $
                      string(state.photplotmode))

zptline = strcompress('0, float,'+string(state.photzpt) + $
                      ',label_left = Magnitude Zeropoint:,' + $
                      'width = 12')

formdesc = [skyline, $
            magline, $
			plotline,$
            zptline, $
            '0, label, [Magnitude = zeropoint - 2.5 log (counts)]', $
            '0, button, Apply Settings, quit', $
            '0, button, Cancel, quit']

textform = cw_form(formdesc, /column, $
                   title = 'atv photometry settings')

if (textform.tag6 EQ 1) then return ; cancelled

state.skytype = textform.tag0
state.magunits = textform.tag1
state.photplotmode = textform.tag2
state.photzpt = textform.tag3

atv_apphot_refresh

end

;----------------------------------------------------------------------

pro atv_apphot

; aperture photometry front end

common atv_state

state.cursorpos = state.coord

if (not (xregistered('atv_apphot', /noshow))) then begin

    apphot_base = $
      widget_base(/base_align_center, $
                  group_leader = state.base_id, $
                  /column, $
                  title = 'atv aperture photometry', $
                  uvalue = 'apphot_base')
    
    apphot_top_base = widget_base(apphot_base, /row, /base_align_center)

    apphot_data_base1 = widget_base( $
            apphot_top_base, /column, frame=0)

    apphot_data_base2 = widget_base( $
            apphot_top_base, /column, frame=0)

    apphot_draw_base = widget_base( $
            apphot_base, /row, /base_align_center, frame=0)

    apphot_data_base1a = widget_base(apphot_data_base1, /column, frame=1)
    tmp_string = $
      string(1000, 1000, $
             format = '("Cursor position:  x=",i4,"  y=",i4)' )

    state.cursorpos_id = $
      widget_label(apphot_data_base1a, $
                   value = tmp_string, $
                   uvalue = 'cursorpos')

    state.centerbox_id = $
      cw_field(apphot_data_base1a, $
               /long, $
               /return_events, $
               title = 'Centering box size (pix):', $
               uvalue = 'centerbox', $
               value = state.centerboxsize, $
               xsize = 5)
    
    tmp_string1 = $
      string(99999.0, 99999.0, $
             format = '("Object centroid:  x=",f7.1,"  y=",f7.1)' )
    
    state.centerpos_id = $
      widget_label(apphot_data_base1a, $
                   value = tmp_string1, $
                   uvalue = 'centerpos')

    state.radius_id = $
      cw_field(apphot_data_base1a, $
               /floating, $
               /return_events, $
               title = 'Aperture radius:', $
               uvalue = 'radius', $
               value = state.aprad, $
               xsize = 8)
    
    state.innersky_id = $
      cw_field(apphot_data_base1a, $
               /floating, $
               /return_events, $
               title = 'Inner sky radius:', $
               uvalue = 'innersky', $
               value = state.innersky, $
               xsize = 8)
    
    state.outersky_id = $
      cw_field(apphot_data_base1a, $
               /floating, $
               /return_events, $
               title = 'Outer sky radius:', $
               uvalue = 'outersky', $
               value = state.outersky, $
               xsize = 8)
    
    photzoom_widget_id = widget_draw( $
         apphot_data_base2, $
         scr_xsize=state.photzoom_size, scr_ysize=state.photzoom_size)

    tmp_string4 = string(0.0, format='("FWHM (pix): ",g7.3)' )
    state.fwhm_id = widget_label(apphot_data_base2, $
                                 value=tmp_string4, $
                                 uvalue='fwhm')

    tmp_string3 = string(10000000.00, $
                         format = '("Sky level: ",g12.6)' )
    
    state.skyresult_id = $
      widget_label(apphot_data_base2, $
                   value = tmp_string3, $
                   uvalue = 'skyresult')
    
    tmp_string2 = string(1000000000.00, $
                         format = '("Object counts: ",g12.6)' )
    
    state.photresult_id = $
      widget_label(apphot_data_base2, $
                   value = tmp_string2, $
                   uvalue = 'photresult', $
                   /frame)
    tmp_string5 = string(1000000000.00, $
                         format = '("Counts error: ",g12.6)' )

    state.photerror_id = $
      widget_label(apphot_data_base2, $
                   value = tmp_string5)


    state.photwarning_id = $
      widget_label(apphot_data_base1, $
                   value='-------------------------', $
                   /dynamic_resize, $
                   frame=1)

    photsettings_id = $
      widget_button(apphot_data_base1, $
                    value = 'Photometry Settings...', $
                    uvalue = 'photsettings')

    radplot_log_save = $
      widget_button(apphot_data_base1, $
                    value = 'Save Photometry Stats', $
                    uvalue = 'radplot_stats_save')

    state.showradplot_id = $
      widget_button(apphot_data_base1, $
                    value = 'Hide Radial Profile', $
                    uvalue = 'showradplot')
    
    state.radplot_widget_id = widget_draw( $
         apphot_draw_base, scr_xsize=500, $
           scr_ysize=(350 < (state.screen_ysize - 350)))

    apphot_ps = $
      widget_button(apphot_data_base2, $
                    value = 'Create Profile PS', $
                    uvalue = 'apphot_ps')

    apphot_done = $
      widget_button(apphot_data_base2, $
                    value = 'Done', $
                    uvalue = 'apphot_done')

    widget_control, apphot_base, /realize

    widget_control, photzoom_widget_id, get_value=tmp_value
    state.photzoom_window_id = tmp_value
    widget_control, state.radplot_widget_id, get_value=tmp_value
    state.radplot_window_id = tmp_value

    xmanager, 'atv_apphot', apphot_base, /no_block
    
    atv_resetwindow
endif

atv_apphot_refresh

end
;--------------------------------------------------------------------
;    atv main program.  needs to be last in order to compile.
;---------------------------------------------------------------------

; Main program routine for ATV.  If there is no current ATV session,
; then run atv_startup to create the widgets.  If ATV already exists,
; then display the new image to the current ATV window.

pro atv, image, $
         min = minimum, $
         max = maximum, $
         autoscale = autoscale,  $
         linear = linear, $
         log = log, $
         histeq = histeq, $
         block = block, $
         align = align, $
         stretch = stretch, $
         header = header,$
		 stokesq = stokesq, stokesu = stokesu, _extra=extra, names=names

common atv_state
common atv_images

if (long(strmid(!version.release,0,1)) LT 6) then begin
    print, 'ATV requires IDL version 6.0 or greater.'
    retall
endif

if (not(keyword_set(block))) then block = 0 else block = 1

newimage = 0

if ( (n_params() EQ 0) AND (xregistered('atv', /noshow))) then begin
    print, 'USAGE: atv, array_name OR fitsfile'
    print, '         [,min = min_value] [,max=max_value] '
    print, '         [,/linear] [,/log] [,/histeq] [,/block]'
    print, '         [,/align] [,/stretch] [,header=header]'
    return
endif

if (!d.name NE 'X' AND !d.name NE 'WIN' AND !d.name NE 'MAC') then begin
    print, 'Graphics device must be set to X, WIN, or MAC for ATV to work.'
	return
endif

; remove old common block variables from prior invocations
delvarx,image_names

; If image is a filename, read in the file
if ( (n_params() NE 0) AND (size(image, /tname) EQ 'STRING')) then begin
    ifexists = findfile(image, count=count)
    if (count EQ 0) then begin
        print, 'ERROR: File not found!'
		return
    endif else begin
		needtoreadfile = 1 ; don't read now; wait until after atv_init
    endelse
endif else begin
	; Check for existence of array
	if ( (n_params() NE 0) AND (size(image, /tname) NE 'STRING') AND $
	   (size(image, /tname) EQ 'UNDEFINED')) then begin
	   	 print, 'ERROR: Data array does not exist!'
		return
	endif
endelse

; Before starting up atv, get the user's external window id.  We can't
; use the atv_getwindow routine yet because we haven't run atv
; startup.  A subtle issue: atv_resetwindow won't work the first time
; through because xmanager doesn't get called until the end of this
; routine.  So we have to deal with the external window explicitly in
; this routine.
if (not (xregistered('atv', /noshow))) then begin
   userwindow = !d.window
   atv_startup
   align = 0B     ; align, stretch keywords make no sense if we are
   stretch = 0B   ; just starting up. 

; Startup message, if desired   
;   print
;   msgstring = strcompress('ATV ' + state.version + ' starting. ')
;   print, msgstring  
;   print

endif


if (n_elements(align) EQ 0) then align = state.default_align
if (n_elements(stretch) EQ 0) then stretch = state.default_stretch

;  Now we can actually read the array in. 
;  This has to go *after* atv_startup since that initializes the
;  common block variables used by readmultifits. -MDP 2004-05-03
;  
; If image is a filename, read in the file
if keyword_set(needtoreadfile) then begin
        atv_readmultifits, fitsfilename=image, newimage=newimage
endif

; Check for existence of array
if ( (n_params() NE 0) AND (size(image, /tname) NE 'STRING') AND $
   (size(image, /tname) EQ 'UNDEFINED')) then begin
    print, 'ERROR: Data array does not exist!'
	atv_shutdown
endif

if (size(image,/tname) eq "COMPLEX") then begin
	print,"ATV can't handle complex images. Try again with abs()?"
	atv_shutdown
	stop
endif

; need to do this *before* opening images, so the title_extras will get set
; appropriately below.
if keyword_set(names) then image_names=names

; If user has passed atv a data array, read it into main_image.
if ( (n_params() NE 0) AND (size(image, /tname) NE 'STRING') AND $
   (size(image, /tname) NE 'UNDEFINED')) then begin
; Make sure it's a 2-d or 3-d array
    if ( (size(image))[0] NE 2 AND (size(image))[0] NE 3) then begin
        print, 'ERROR: Input data must be a 2-d or 3-d array!'    
	    if (size(image))[0] eq 4 then begin
			  print, 'Reforming a 4D array to 3D...'
			  sz = size(image)
			  image = reform(image, sz[1], sz[2], sz[3]*sz[4])
		  endif else begin
		  	atv_shutdown
			return
		  endelse
    endif else begin
        if (size(image))[0] EQ 2 THEN begin
            main_image = image
            newimage = 1
            state.image_size = [(size(main_image_stack))[1:2], 1]
            state.imagename = ''
            state.title_extras = ''
			image_names = '' ; no image names. 
            atv_setheader, header
            
            widget_control, state.curimnum_base_id,map=0,xsize=1,ysize=1


        endif else begin ; case of 3-d stack of images [x,y,n]
            main_image_stack = image
            main_image = main_image_stack[*, *, 0]
            state.image_size = (size(main_image_stack))[1:3]
            state.cur_image_num = 0
            newimage = 1
            state.imagename = ''
			if n_elements(image_names) lt 1 then begin
				image_names = "Plane "+strtrim(string(indgen(state.image_size[2])))
			endif
			state.title_extras = image_names[state.cur_image_num]
            atv_setheader, header

            widget_control,state.curimnum_base_id,map=1, $
              xsize=state.draw_window_size[0],ysize=45

            widget_control, state.curimnum_text_id, sensitive = 1, $
                     set_value = 0
            widget_control, state.curimnum_slidebar_id, sensitive = 1, $
                     set_value = 0, set_slider_min = 0, $
                     set_slider_max = state.image_size[2]-1
        endelse
          ;Reset image rotation angle to 0 and inversion to none
           state.rot_angle = 0.
           state.invert_image = 'none'
        if (state.firstimage EQ 1) then begin
            align = 0
            stretch = 0
        endif
    endelse
 endif
 

;   Define default startup image 
if (n_elements(main_image) LE 1) then begin
    ;main_image = cos(((findgen(500)- 250)*2.9) # ((findgen(500)-250)*1.5))
    gridsize = 512
    centerpix = 256
; square pattern:
    x = ((findgen(gridsize) # replicate(1, gridsize)) - centerpix + 0.001)*0.5
    y = ((replicate(1, gridsize) # findgen(gridsize)) - centerpix + 0.001)*0.5
    main_image = (sin(x)/x) * (sin(y)/y)
    state.min_value = 0.0
    state.max_value = 1.0
    stretch = 1
    autoscale = 0
    imagename = ''
    newimage = 2             ; flag for startup image
    atv_setheader, ''
    state.title_extras = 'firstimage'
endif


if (newimage GE 1) then begin  
; skip this part if new image is invalid or if user selected 'cancel'
; in dialog box
    atv_getstats, align=align
    
    display_image = 0

    if n_elements(minimum) GT 0 then begin
        state.min_value = minimum
    endif
    
    if n_elements(maximum) GT 0 then begin 
        state.max_value = maximum
    endif
    
    if state.min_value GE state.max_value then begin
        state.min_value = state.max_value - 1.
    endif
    
    if (keyword_set(linear)) then state.scaling = 0
    if (keyword_set(log))    then state.scaling = 1
    if (keyword_set(histeq)) then state.scaling = 2
    if (keyword_set(asinh))  then state.scaling = 3
    
; Perform autoscale if current stretch invalid or stretch keyword
; not set, or if this is the first image
if n_elements(stretch) eq 0 then stretch=1 ; MDP addition March 2007 to stop crashes
    IF (state.min_value EQ state.max_value) OR (stretch EQ 0) THEN BEGIN 

       if (keyword_set(autoscale) OR $
           ((state.default_autoscale EQ 1) AND (n_elements(minimum) EQ 0) $
            AND (n_elements(maximum) EQ 0)) ) $
         then atv_autoscale
    ENDIF 

	; MDP 2007-06-15 don't autoscale if max or min set explicitly
    if (state.firstimage EQ 1 AND newimage EQ 1) and (~keyword_set(minimum) and ~keyword_set(maximum)) then atv_autoscale
    if (newimage EQ 1) then state.firstimage = 0  ; now have a real image
        
    atv_set_minmax
    
    IF ((NOT keyword_set(align)) AND state.default_align EQ 0) THEN BEGIN 
       state.zoom_level = 0
       state.zoom_factor = 1.0
    ENDIF 

    atv_displayall
    atv_settitle
    
    atv_resetwindow
	if state.autozoom and not(keyword_set(align)) then atv_autozoom
endif
if keyword_set(stokesq) and keyword_set(stokesu) then $
	atvpol,stokesq,stokesu,/noxmcheck,_extra=extra

state.block = block

; Register the widget with xmanager if it's not already registered
if (not(xregistered('atv', /noshow))) then begin
    nb = abs(block - 1)
    xmanager, 'atv', state.base_id, no_block = nb, cleanup = 'atv_shutdown'
    wset, userwindow
    ; if blocking mode is set, then when the procedure reaches this
    ; line atv has already been terminated.  If non-blocking, then
    ; the procedure continues below.  If blocking, then the state
    ; structure doesn't exist any more so don't set active window.
    if (block EQ 0) then state.active_window_id = userwindow
endif




end





