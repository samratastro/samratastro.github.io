;goto, iia

x_cut=[320,500]
y_cut=[300,450]

pixscale = 0.6

org=[3000+x_cut[0]-2048,1500+y_cut[0]-2048]*pixscale 
scale=[pixscale,pixscale]
filenames=file_search('/home/samrat/Plasmoid_event/dem_output/sav_files_crop/*.sav')
;f1=file_search('/home/samrat/Plasmoid_event/dem_output/*.fits')

savfilenums=1979

for p=0,savfilenums-1 do begin
tic,/profiler
;mreadfits, f1(i), hdr, img
restore, filenames[p]
lgtem = fltarr((1+x_cut[1]-x_cut[0]),(1+y_cut[1]-y_cut[0]))
ntbins = size(lgtaxis,/dim)

EMsum = fltarr((1+x_cut[1]-x_cut[0]),(1+y_cut[1]-y_cut[0]))

EMsum = total(emcube,3)

for a=0,(1+x_cut[1]-x_cut[0])-1 do begin
for b=0,(1+y_cut[1]-y_cut[0])-1 do begin
lgtem(a,b)= total(lgtaxis(*)*emcube(a,b,*)/EMsum(a,b))
endfor
endfor
;img1=img[150:850,100:800]
;lgtems=lgtem[150:850,100:800]
window, xs=1100, ys=1000,/pixmap
;!P.MULTI = [0, 2, 1]
;aia_lct, wave=131, /load
;plot_image, img1^.20,background=cgcolor('white'), color=cgcolor('black'),xticklen=-0.01,yticklen=-0.01,charthick=2,charsize=2,position=[0.07,0.1,0.46,0.9],origin=orgn, scale=scl, xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)',title='AIA_131_T_'+string(obs_time_p)
loadct,5
plot_image, lgtem, background=fsc_color('white'), color = 0,xticklen=-0.01,yticklen=-0.01,charthick=2,charsize=2,position=[0.1,0.1,0.8,0.9],origin=org, scale=scale, xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)',title='Temperature Map (t='+obs_time+')'
colorbar, /vertical, range=[lgtaxis[0],lgtaxis[20]],color=0, position=[0.9,0.1,0.95,0.9], title='logT',charthick=2,charsize=2,format='(F4.1)';divisions=(lgtaxis[20]-lgtaxis[0])/21
;save,emsum,lgtem,obs_time,obs_time_p,filename='/home/samrat/Plasmoid_event/dem_output/Tmap/sav/Tmap'+obs_time+'.sav'
write_png,'/home/samrat/Plasmoid_event/dem_output/Tmap/Tmap'+obs_time+'.png',TVRD(/true)
lgtem=0
EMsum=0
hdr=0
img=0
img1=0
EMsum=0
lgtems=0
emcube=0
coeff=0
print,'Saving image Number = ',p+1
toc
endfor

end

;stop

;lgtems=lgtem[200:900,100:800]
;im = IMAGE(lgtems, RGB_TABLE=5, LAYOUT=[2,1,1],$
;   ;POSITION=[0.10,0.05,0.95,0.9], $
;   FONT_COLOR='Blue', FONT_SIZE=16, $
;   TITLE='Temprature_Maps')
;c = COLORBAR(TARGET=im, ORIENTATION=1, RANGE = [lgtaxis[0],lgtaxis[20]], $
; ;  POSITION=[0.10,0.05,0.13,0.9], 
; FONT_SIZE=8,$
;  TITLE='lgt')
;  aia_lct, wave=131, /load
;  ;loadct,4,/silent
; im = IMAGE(img^0.25,/CURRENT, LAYOUT=[2,1,2],$
;  ; POSITION=[0.10,0.05,0.95,0.9], $
;   FONT_COLOR='Blue', FONT_SIZE=16, $
;   TITLE='AIA_131') 

