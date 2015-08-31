set encoding iso_8859_1
set xrange [1:250]
set yrange [0.83:1.05]
set autoscale xfix
set ter png enh interlace size 2400,840 font 'Nimbus,40'
set xlabel 'Position'
set ylabel 'Reliability           ' 
set ytics nomirror 0.9,0.1,1
set out '/big/software/TOPCONS2/topcons2_webserver/test/out2//seq_0///Topcons/topcons.large.png'
set tmargin 1.3
set lmargin 11.5
set rmargin 6.5
set label 'TOPCONS' font 'Nimbus,42' at screen 0.022,0.775
set object 1 rect from 38.5,1.02 to 49.5,1.021625 fc rgb "red" fs noborder
set object 2 rect from 110.5,1.02 to 115.5,1.021625 fc rgb "red" fs noborder
set object 3 rect from 163.5,1.02 to 179.5,1.021625 fc rgb "red" fs noborder
set object 4 rect from 233.5,1.02 to 250.5,1.021625 fc rgb "red" fs noborder
set object 5 rect from 0.5,1.027125 to 17.5,1.02875 fc rgb "blue" fs noborder
set object 6 rect from 70.5,1.027125 to 89.5,1.02875 fc rgb "blue" fs noborder
set object 7 rect from 136.5,1.027125 to 142.5,1.02875 fc rgb "blue" fs noborder
set object 8 rect from 200.5,1.027125 to 212.5,1.02875 fc rgb "blue" fs noborder
set object 9 rect from 17.5,1.02 to 38.5,1.02875 fc rgb "white"
set object 10 rect from 49.5,1.02 to 70.5,1.02875 fc rgb "grey" fs noborder
set object 11 rect from 89.5,1.02 to 110.5,1.02875 fc rgb "white"
set object 12 rect from 115.5,1.02 to 136.5,1.02875 fc rgb "grey" fs noborder
set object 13 rect from 142.5,1.02 to 163.5,1.02875 fc rgb "white"
set object 14 rect from 179.5,1.02 to 200.5,1.02875 fc rgb "grey" fs noborder
set object 15 rect from 212.5,1.02 to 233.5,1.02875 fc rgb "white"
set object 16 rect from 17.5,1.02 to 38.5,1.02875 fc rgb "white"
set object 17 rect from 49.5,1.02 to 70.5,1.02875 fc rgb "grey" fs noborder
set object 18 rect from 89.5,1.02 to 110.5,1.02875 fc rgb "white"
set object 19 rect from 115.5,1.02 to 136.5,1.02875 fc rgb "grey" fs noborder
set object 20 rect from 142.5,1.02 to 163.5,1.02875 fc rgb "white"
set object 21 rect from 179.5,1.02 to 200.5,1.02875 fc rgb "grey" fs noborder
set object 22 rect from 212.5,1.02 to 233.5,1.02875 fc rgb "white"
plot '/big/software/TOPCONS2/topcons2_webserver/test/out2//seq_0///Topcons/reliability.final' w l t '' lc rgb "black" lw 4
exit
