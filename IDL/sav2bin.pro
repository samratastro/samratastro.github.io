pro sav2bin

indir ='/Users/samratsen/extrapolation/output/ar13753/extrapolation/15_07_2024_h1210/'
restore,indir+'Bnfff_ar13753_15_07_2024_h1210_vectorb_784_512_512.sav',/v

   openw,1,indir+'bx.bin'
   writeu,1, bvx
   close,1

   openw,1,indir+'by.bin'
   writeu,1, bvy
   close,1

   openw,1,indir+'bz.bin'
  writeu,1,bvz
 close,1

;restore,'jx.sav',/v
;   openw,1,'jx.bin'
;  writeu,1,jx
; close,1

;restore,'jy.sav',/v
;   openw,1,'jy.bin'
;  writeu,1,jy
; close,1

;restore,'jz.sav',/v
;   openw,1,'jz.bin'
;  writeu,1,jz
; close,1

;restore,indir+'slqtcube.sav',/v
;   openw,1,indir+'slq.bin'
;  writeu,1,slqt
; close,1
;
; restore,indir+'twcube.sav',/v
;   openw,1,indir+'tw.bin'
;  writeu,1,tw
; close,1

stop
end
