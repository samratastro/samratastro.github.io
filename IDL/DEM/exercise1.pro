;files=fltarr(6)
f1= file_search('/home/samrat/Work/samrat_dem/94/*.fits')
f2= file_search('/home/samrat/Work/samrat_dem/131/*.fits')
f3= file_search('/home/samrat/Work/samrat_dem/171/*.fits')
f4= file_search('/home/samrat/Work/samrat_dem/193/*.fits')
f5= file_search('/home/samrat/Work/samrat_dem/211/*.fits')
f6= file_search('/home/samrat/Work/samrat_dem/335/*.fits')

imgnum=2

for i=0,imgnum-1 do begin
tic,/profiler
files=[f1(i),f2(i),f3(i),f4(i),f5(i),f6(i)]
; Author: Mark Cheung
; Purpose: Sample script for running sparse DEM inversion (see Cheung
; et al. 2015) for details.
; Revision history: 2015-10-20 First version
;                   2015-10-23 Added note about lgT axis
print, files
; Restore some prepared AIA level 1.5 data (i.e. already been aia_preped)
;restgen, file=aiadatafile, struct=s
mreadfits,files,oindex,odata
 for j=0,5 do begin
    odata[*,*,j] = odata[*,*,j]/oindex[j].exptime ; exposure time normalization ---> DN/s
 endfor
; Initialize solver. 
; This step builds up the response functions (which are
; time-dependent) and the basis functions.
; If running inversions over a set of data spanning over a few hours
; or even days, it is not necessary to re-initialize

;IF (0) THEN BEGIN 
; As discussed in the Appendix of Cheung et al. (2015), the inversion
; can push some EM into temperature bins if the lgT axis goes below
; logT ~ 5.7 (see Fig. 18). It is suggested the user check the
; dependence of the solution by varying lgTmin.
lgTmin = 5.7   ; minimum for lgT axis for inversion 
dlgT   = 0.1 ; width of lgT bin
nlgT   = 21 ; number of lgT bins
aia_sparse_em_init, timedepend =oindex[0].date_obs, /evenorm, use_lgtaxis=findgen(nlgT)*dlgT+lgTmin,bases_sigmas=[0.0,0.1,0.2,0.6]

lgtaxis = aia_sparse_em_lgtaxis()
; ENDIF

; We will use the data in s.img to invert. s.img is a stack of level
; 1.5, exposure time normalized AIA pixel values. 
exptimestr = '['+strjoin(string(oindex.exptime,format='(F8.5)'),',')+']'

; This defines the tolerance function for the inversion. 
; y denotes the values in DN/s/pixel (i.e. level 1.5, exposure time
; normalized)

; aia_bp_estimate_error gives the estimated uncertainty for a given
; DN / pixel for a given AIA channel. Since y contains DN/pixel/s,
; we have to multiply by exposure time first to pass
; aia_bp_estimate_error. Then we will divide the output of
; aia_bp_estimate_error by the exposure time again.
; If one wants to include uncertainties in atomic data, suggest to add
; the /temp keyword to aia_bp_estimate_error()
tolfunc = 'aia_bp_estimate_error(y*'+exptimestr+', [94,131,171,193,211,335], num_images='+strtrim(string(1.0^2,format='(I4)'),2)+')/'+exptimestr

; Do DEM solve. 
; Note!!! The solver assumes the third dimension is arranged according to [94,131,171,193,211,335]
aia_sparse_em_solve, odata, tolfunc=tolfunc, tolfac=1.4, oem=emcube, status=status, coeff=coeff
; emcube contains the emission measure contained in each lgt bin
; status contains a mask indicating whether a solution was
; found. status = 0 means a solution was found within the given constriants.
; If you want to actual coefficients of the basis functions
; (i.e. x_i's of Kx = y), then add coeff=coeff to the call to aia_sparse_em_solve

; Now show results
obs_time = strmid(oindex[0].T_rec,11,2)+strmid(oindex[0].T_rec,14,2)+strmid(oindex[0].T_rec,17,2)
obs_time_p=strmid(oindex[0].T_rec,11,2)+':'+strmid(oindex[0].T_rec,14,2)+':'+strmid(oindex[0].T_rec,17,2)
dispimage= aia_sparse_em_display(em=emcube)
window,xs=1380, ys=768,TITLE='EM Map for AIA_94_T'+obs_time,/pixmap
tv,dispimage,/true
;write_png,'/home/samrat/Plasmoid_event/dem_output/image_files/dem_map_'+obs_time+'.png',TVRD(/true)
save,/variables,filename='/home/samrat/Work/samrat_dem/sav_files/dem_map_'+obs_time+'.sav'
print, 'saving file number  ',i+1
emcube=0
coeff=0
DISPIMAGE=0
FILES=0
ODATA=0
OINDEX=0
TOLFUNC=0
toc
endfor 

end





