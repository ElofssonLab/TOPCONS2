set encoding iso_8859_1
set xrange [1:334]
set yrange [0.83:1.25]
set autoscale xfix
set ter png enh interlace size 2400,840 font 'Nimbus,40'
set xlabel 'Position'
set ylabel 'Reliability           ' 
set ytics nomirror 0.5,0.1,1
set out '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_2///Topcons/topcons.large.png'
set tmargin 1.3
set lmargin 11.5
set rmargin 6.5
set label 'TOPCONS' font 'Nimbus,42' at screen 0.022,0.775
set object 1 rect from 67.5,1.1 to 81.5,1.108125 fc rgb "red" fs noborder
set object 2 rect from 136.5,1.1 to 154.5,1.108125 fc rgb "red" fs noborder
set object 3 rect from 219.5,1.1 to 249.5,1.108125 fc rgb "red" fs noborder
set object 4 rect from 303.5,1.1 to 334.5,1.108125 fc rgb "red" fs noborder
set object 5 rect from 0.5,1.135625 to 46.5,1.14375 fc rgb "blue" fs noborder
set object 6 rect from 102.5,1.135625 to 115.5,1.14375 fc rgb "blue" fs noborder
set object 7 rect from 175.5,1.135625 to 198.5,1.14375 fc rgb "blue" fs noborder
set object 8 rect from 270.5,1.135625 to 282.5,1.14375 fc rgb "blue" fs noborder
set object 9 rect from 46.5,1.1 to 67.5,1.14375 fc rgb "white"
set object 10 rect from 81.5,1.1 to 102.5,1.14375 fc rgb "grey" fs noborder
set object 11 rect from 115.5,1.1 to 136.5,1.14375 fc rgb "white"
set object 12 rect from 154.5,1.1 to 175.5,1.14375 fc rgb "grey" fs noborder
set object 13 rect from 198.5,1.1 to 219.5,1.14375 fc rgb "white"
set object 14 rect from 249.5,1.1 to 270.5,1.14375 fc rgb "grey" fs noborder
set object 15 rect from 282.5,1.1 to 303.5,1.14375 fc rgb "white"
set object 16 rect from 46.5,1.1 to 67.5,1.14375 fc rgb "white"
set object 17 rect from 81.5,1.1 to 102.5,1.14375 fc rgb "grey" fs noborder
set object 18 rect from 115.5,1.1 to 136.5,1.14375 fc rgb "white"
set object 19 rect from 154.5,1.1 to 175.5,1.14375 fc rgb "grey" fs noborder
set object 20 rect from 198.5,1.1 to 219.5,1.14375 fc rgb "white"
set object 21 rect from 249.5,1.1 to 270.5,1.14375 fc rgb "grey" fs noborder
set object 22 rect from 282.5,1.1 to 303.5,1.14375 fc rgb "white"
plot '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_2///Topcons/reliability.final' w l t '' lc rgb "black" lw 4
exit
