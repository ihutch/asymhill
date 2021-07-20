make AsymHill
make electronholeorbits

# Illustrative explanation plot
./AsymHill -p.09 -e.06 -f -g1.5,1.,.3 -g-1.5,1.,.7 -x -w-3
mv plot0001.ps fiplots.ps
epstopdf fiplots.ps
mv plot0002.ps explot.ps
epstopdf explot.ps
mv plot0004.ps convplot.ps
epstopdf convplot.ps
mv plot0008.ps forceplot.ps
epstopdf forceplot.ps
mv plot0006.ps ftrapped.ps
epstopdf ftrapped.ps

./AsymHill -p.09 -e.06 -g1.5,1.,.3 -g-1.5,1.,.7 -s7 -x -w-3
mv plot0001.ps dynamic.ps
epstopdf dynamic.ps

./AsymHill -p.4 -e.00 -g.75,1.,.3 -s0 -m -w-3
mv plot0001.ps multiden.ps
epstopdf multiden.ps

./AsymHill -p.4  -g1.1972,1.,.7 -g-1.1972,.4472,.3 -w-3 -a
mv plot0001.ps phiscaled.ps
epstopdf phiscaled.ps
mv plot0002.ps scaled.ps
epstopdf scaled.ps

sed -i -e 's/nphi=200/nphi=600/' AsymHill.f90
make AsymHill
./AsymHill -p.09 -e.06 -f -g1.5,1.,.3 -g-1.5,1.,.7 -x -w-3
mv plot0007.ps phiofmodx.ps
epstopdf phiofmodx.ps
sed -i -e 's/nphi=600/nphi=200/' AsymHill.f90
make AsymHill

./electronholeorbits
mv plot0001.ps electronholeorbits.ps
epstopdf  electronholeorbits.ps

