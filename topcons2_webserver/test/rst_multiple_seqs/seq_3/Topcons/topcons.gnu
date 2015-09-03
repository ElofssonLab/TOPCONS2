set encoding iso_8859_1
set xrange [1:330]
set yrange [0.83:1.2]
set autoscale xfix
set ter png enh interlace size 2400,840 font 'Nimbus,40'
set xlabel 'Position'
set ylabel 'Reliability           ' 
set ytics nomirror 0.6,0.1,1
set out '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_3///Topcons/topcons.large.png'
set tmargin 1.3
set lmargin 11.5
set rmargin 6.5
set label 'TOPCONS' font 'Nimbus,42' at screen 0.022,0.775
set object 1 rect from 63.5,1.08 to 75.5,1.0865 fc rgb "red" fs noborder
set object 2 rect from 132.5,1.08 to 152.5,1.0865 fc rgb "red" fs noborder
set object 3 rect from 215.5,1.08 to 247.5,1.0865 fc rgb "red" fs noborder
set object 4 rect from 299.5,1.08 to 330.5,1.0865 fc rgb "red" fs noborder
set object 5 rect from 0.5,1.1085 to 42.5,1.115 fc rgb "blue" fs noborder
set object 6 rect from 96.5,1.1085 to 111.5,1.115 fc rgb "blue" fs noborder
set object 7 rect from 173.5,1.1085 to 194.5,1.115 fc rgb "blue" fs noborder
set object 8 rect from 268.5,1.1085 to 278.5,1.115 fc rgb "blue" fs noborder
set object 9 rect from 42.5,1.08 to 63.5,1.115 fc rgb "white"
set object 10 rect from 75.5,1.08 to 96.5,1.115 fc rgb "grey" fs noborder
set object 11 rect from 111.5,1.08 to 132.5,1.115 fc rgb "white"
set object 12 rect from 152.5,1.08 to 173.5,1.115 fc rgb "grey" fs noborder
set object 13 rect from 194.5,1.08 to 215.5,1.115 fc rgb "white"
set object 14 rect from 247.5,1.08 to 268.5,1.115 fc rgb "grey" fs noborder
set object 15 rect from 278.5,1.08 to 299.5,1.115 fc rgb "white"
set object 16 rect from 42.5,1.08 to 63.5,1.115 fc rgb "white"
set object 17 rect from 75.5,1.08 to 96.5,1.115 fc rgb "grey" fs noborder
set object 18 rect from 111.5,1.08 to 132.5,1.115 fc rgb "white"
set object 19 rect from 152.5,1.08 to 173.5,1.115 fc rgb "grey" fs noborder
set object 20 rect from 194.5,1.08 to 215.5,1.115 fc rgb "white"
set object 21 rect from 247.5,1.08 to 268.5,1.115 fc rgb "grey" fs noborder
set object 22 rect from 278.5,1.08 to 299.5,1.115 fc rgb "white"
plot '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_3///Topcons/reliability.final' w l t '' lc rgb "black" lw 4
exit
