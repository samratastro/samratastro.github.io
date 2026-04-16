;.compile aia_sparse_em_init
;restore,'/home/samrat/Plasmoid_event/dem_output/sav_files_high/dem_map_161011.sav',/v

xrange=[320,500]
yrange=[300,450]

pxscale = 0.6
origin=[3000+xrange[0]-2048,1500+yrange[0]-2048]*pxscale 
scl=[pxscale,pxscale]




fl=file_search('/home/samrat/Plasmoid_event/dem_output/sav_files_crop/*.sav')


loadct,5
window,0,xs=3000,ys=1700,/pixmap

filenums=n_elements(fl)

for j1=0,filenums-1 do begin

restore, fl[j1]

em1=total(emcube[*,*,0:3],3)
em2=total(emcube[*,*,3:7],3)
em3=total(emcube[*,*,7:11],3)
em4=total(emcube[*,*,11:14],3)
em5=total(emcube[*,*,14:17],3)
em6=total(emcube[*,*,17:20],3)




!P.Multi = [0, 3, 2]

plot_image, em1, background=fsc_color('white'), color = 0, position=[0.04,0.6,0.32,0.97],xticklen=-0.01,yticklen=-0.01,charthick=3,charsize=3,origin=origin, scale=scl, xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)',title='Log T=[5.7-6.0] (t='+obs_time+')' 
;colorbar, /horizontal,color=0, position=[0.06,0.52,0.31,0.54], title='EM [10^26 cm^-5]',charthick=2,charsize=2,format='(F5.1)';divisions=(lgtaxis[20]-lgtaxis[0])/21

plot_image, em2, background=fsc_color('white'), color = 0, position=[0.37,0.6,0.65,0.97],xticklen=-0.01,yticklen=-0.01,charthick=3,charsize=3,origin=origin, scale=scl, xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)',title='Log T=[6.0-6.4] (t='+obs_time+')' 
colorbar, /horizontal,color=0, position=[0.3,0.52,0.71,0.54], title='EM [10^26 cm^-5]',charthick=2,charsize=3,format='(F5.1)';divisions=(lgtaxis[20]-lgtaxis[0])/21

plot_image, em3, background=fsc_color('white'), color = 0, position=[0.7,0.6,0.98,0.97],xticklen=-0.01,yticklen=-0.01,charthick=3,charsize=3,origin=origin, scale=scl, xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)',title='Log T=[6.4-6.8] (t='+obs_time+')'
;colorbar, /horizontal,color=0, position=[0.72,0.52,0.97,0.54], title='EM [10^26 cm^-5]',charthick=2,charsize=2,format='(F5.1)';divisions=(lgtaxis[20]-lgtaxis[0])/21


plot_image, em4, background=fsc_color('white'), color = 0, position=[0.04,0.11,0.32,0.48],xticklen=-0.01,yticklen=-0.01,charthick=3,charsize=3,origin=origin, scale=scl, xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)',title='Log T=[6.8-7.1] (t='+obs_time+')'
;colorbar, /horizontal,color=0, position=[0.06,0.03,0.31,0.05], title='EM [10^26 cm^-5]',charthick=2,charsize=2,format='(F5.1)';divisions=(lgtaxis[20]-lgtaxis[0])/21

plot_image, em5, background=fsc_color('white'), color = 0, position=[0.37,0.11,0.65,0.48],xticklen=-0.01,yticklen=-0.01,charthick=3,charsize=3,origin=origin, scale=scl, xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)',title='Log T=[7.1-7.4] (t='+obs_time+')'
;colorbar, /horizontal,color=0, position=[0.4,0.03,0.64,0.05], title='EM [10^26 cm^-5]',charthick=2,charsize=2,format='(F5.1)';divisions=(lgtaxis[20]-lgtaxis[0])/21

plot_image, em6, background=fsc_color('white'), color = 0, position=[0.7,0.11,0.98,0.48],xticklen=-0.01,yticklen=-0.01,charthick=3,charsize=3,origin=origin, scale=scl, xtitle='Solar X (arcsec)', ytitle='Solar Y (arcsec)',title='Log T=[7.4-7.7] (t='+obs_time+')'
;colorbar, /horizontal,color=0, position=[0.73,0.03,0.97,0.05], title='EM [10^26 cm^-5]',charthick=2,charsize=2,format='(F5.1)';divisions=(lgtaxis[20]-lgtaxis[0])/21


write_png,'/home/samrat/Plasmoid_event/dem_output/em_images/em_map'+obs_time+'.png',TVRD(/true)

print, "Saving image number ", j1+1
endfor

goto, jump



lgemrange = [26, 29.5] 
  ;loadct, 22, /silent
  k= 1.0000000e+26
data = bytscl(alog10(total(emcube[*,*,4:7],3)*k),min=lgemrange[0],max=lgemrange[1])
;data = alog10(total(EMCUBE[*,*,1:3],3)*k)
;        caption = 'EM in lgt=['+$
;                strtrim(string(5.75,format='(F5.2)'),2)+','+$
;                strtrim(string(6.05,format='(F5.2)'),2)+']'

im = IMAGE(data, RGB_TABLE=22,$
   POSITION=[0.10,0.05,0.95,0.9], $
   FONT_COLOR='Blue', FONT_SIZE=16, $
   TITLE='DEM')
   
;arr =  fltarr(6,1001) 
;mini =  lgemrange[0]
;maxm =  lgemrange[1]
;step=(maxm-mini)/1001

;for i = 0,6-1 do begin 
;arr (i,*) = mini + findgen(1001) * step
;endfor
;im1=IMAGE(arr, RGB_TABLE=22,LAYOUT=[2,1,2],$
;   POSITION=[0.22,0.05,0.29,0.9], $
;   FONT_COLOR='Blue', FONT_SIZE=12, $
;   TITLE='EM')

 
c = COLORBAR(TARGET=im, ORIENTATION=1, RANGE = lgemrange, $
   POSITION=[0.10,0.05,0.13,0.9], FONT_SIZE=8,$
   TITLE='EM')
c.TEXTPOS = 0
c.TICKDIR = 1
c.BORDER_ON = 1
c.COLOR = 'Blue'
c.FONT_STYLE = 'Italic'




jump:

END



