#!/bin/bash

cd /Applications/vapor.app/Contents/MacOS

indir="/Users/samratsen/extrapolation/output/ar13753/extrapolation/15_07_2024_h1210/"
echo $indir
./vdccreate -dimension 784x512x512 -vars3d bx:by:bz $indir/ar13753.vdc

./raw2vdc -ts 0 -type float64 -varname bx $indir/ar13753.vdc $indir/bx.bin
./raw2vdc -ts 0 -type float64 -varname by $indir/ar13753.vdc $indir/by.bin
./raw2vdc -ts 0 -type float64 -varname bz $indir/ar13753.vdc $indir/bz.bin
###./raw2vdc -ts 0 -type float64 -varname slqt $indir/ar13753.vdc $indir/slq.bin
###./raw2vdc -ts 0 -type float64 -varname tw $indir/ar13753.vdc $indir/tw.bin

cd $indir
