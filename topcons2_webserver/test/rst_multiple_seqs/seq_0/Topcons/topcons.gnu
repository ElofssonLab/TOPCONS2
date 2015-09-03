set encoding iso_8859_1
set xrange [1:359]
set yrange [0.83:1.15]
set autoscale xfix
set ter png enh interlace size 2400,840 font 'Nimbus,40'
set xlabel 'Position'
set ylabel 'Reliability           ' 
set ytics nomirror 0.7,0.1,1
set out '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_0///Topcons/topcons.large.png'
set tmargin 1.3
set lmargin 11.5
set rmargin 6.5
set label 'TOPCONS' font 'Nimbus,42' at screen 0.022,0.775
set object 1 rect from 70.5,1.06 to 84.5,1.064875 fc rgb "red" fs noborder
set object 2 rect from 147.5,1.06 to 159.5,1.064875 fc rgb "red" fs noborder
set object 3 rect from 234.5,1.06 to 248.5,1.064875 fc rgb "red" fs noborder
set object 4 rect from 317.5,1.06 to 359.5,1.064875 fc rgb "red" fs noborder
set object 5 rect from 0.5,1.081375 to 49.5,1.08625 fc rgb "blue" fs noborder
set object 6 rect from 105.5,1.081375 to 126.5,1.08625 fc rgb "blue" fs noborder
set object 7 rect from 180.5,1.081375 to 213.5,1.08625 fc rgb "blue" fs noborder
set object 8 rect from 269.5,1.081375 to 296.5,1.08625 fc rgb "blue" fs noborder
set object 9 rect from 49.5,1.06 to 70.5,1.08625 fc rgb "white"
set object 10 rect from 84.5,1.06 to 105.5,1.08625 fc rgb "grey" fs noborder
set object 11 rect from 126.5,1.06 to 147.5,1.08625 fc rgb "white"
set object 12 rect from 159.5,1.06 to 180.5,1.08625 fc rgb "grey" fs noborder
set object 13 rect from 213.5,1.06 to 234.5,1.08625 fc rgb "white"
set object 14 rect from 248.5,1.06 to 269.5,1.08625 fc rgb "grey" fs noborder
set object 15 rect from 296.5,1.06 to 317.5,1.08625 fc rgb "white"
set object 16 rect from 49.5,1.06 to 70.5,1.08625 fc rgb "white"
set object 17 rect from 84.5,1.06 to 105.5,1.08625 fc rgb "grey" fs noborder
set object 18 rect from 126.5,1.06 to 147.5,1.08625 fc rgb "white"
set object 19 rect from 159.5,1.06 to 180.5,1.08625 fc rgb "grey" fs noborder
set object 20 rect from 213.5,1.06 to 234.5,1.08625 fc rgb "white"
set object 21 rect from 248.5,1.06 to 269.5,1.08625 fc rgb "grey" fs noborder
set object 22 rect from 296.5,1.06 to 317.5,1.08625 fc rgb "white"
plot '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_0///Topcons/reliability.final' w l t '' lc rgb "black" lw 4
exit
