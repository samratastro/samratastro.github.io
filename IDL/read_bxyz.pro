pro read_bxyz

fn = './bfield_ar12280_256_256_256.sav'
nx = 256
ny = 256
nz = 256
restore, fn
Bx = BXC
By = BYC
Bz = BZC


openw, var_lun, 'mag_3d.dat',/get_lun
for k=0,nz-1 do begin
  for j=0,ny-1 do begin
    for i=0,nx-1 do begin
       printf,var_lun,format='(3e25.16)', Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    endfor
  endfor
endfor
close, var_lun
free_lun,var_lun


end
