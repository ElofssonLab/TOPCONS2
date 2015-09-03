set encoding iso_8859_1
set xrange [1:338]
set yrange [0.83:1.1]
set autoscale xfix
set ter png enh interlace size 2400,840 font 'Nimbus,40'
set xlabel 'Position'
set ylabel 'Reliability           ' 
set ytics nomirror 0.8,0.1,1
set out '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_1///Topcons/topcons.large.png'
set tmargin 1.3
set lmargin 11.5
set rmargin 6.5
set label 'TOPCONS' font 'Nimbus,42' at screen 0.022,0.775
set object 1 rect from 43.5,1.04 to 55.5,1.04325 fc rgb "red" fs noborder
set object 2 rect from 114.5,1.04 to 136.5,1.04325 fc rgb "red" fs noborder
set object 3 rect from 202.5,1.04 to 244.5,1.04325 fc rgb "red" fs noborder
set object 4 rect from 297.5,1.04 to 338.5,1.04325 fc rgb "red" fs noborder
set object 5 rect from 0.5,1.05425 to 22.5,1.0575 fc rgb "blue" fs noborder
set object 6 rect from 76.5,1.05425 to 93.5,1.0575 fc rgb "blue" fs noborder
set object 7 rect from 157.5,1.05425 to 181.5,1.0575 fc rgb "blue" fs noborder
set object 8 rect from 265.5,1.05425 to 276.5,1.0575 fc rgb "blue" fs noborder
set object 9 rect from 22.5,1.04 to 43.5,1.0575 fc rgb "white"
set object 10 rect from 55.5,1.04 to 76.5,1.0575 fc rgb "grey" fs noborder
set object 11 rect from 93.5,1.04 to 114.5,1.0575 fc rgb "white"
set object 12 rect from 136.5,1.04 to 157.5,1.0575 fc rgb "grey" fs noborder
set object 13 rect from 181.5,1.04 to 202.5,1.0575 fc rgb "white"
set object 14 rect from 244.5,1.04 to 265.5,1.0575 fc rgb "grey" fs noborder
set object 15 rect from 276.5,1.04 to 297.5,1.0575 fc rgb "white"
set object 16 rect from 22.5,1.04 to 43.5,1.0575 fc rgb "white"
set object 17 rect from 55.5,1.04 to 76.5,1.0575 fc rgb "grey" fs noborder
set object 18 rect from 93.5,1.04 to 114.5,1.0575 fc rgb "white"
set object 19 rect from 136.5,1.04 to 157.5,1.0575 fc rgb "grey" fs noborder
set object 20 rect from 181.5,1.04 to 202.5,1.0575 fc rgb "white"
set object 21 rect from 244.5,1.04 to 265.5,1.0575 fc rgb "grey" fs noborder
set object 22 rect from 276.5,1.04 to 297.5,1.0575 fc rgb "white"
plot '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_1///Topcons/reliability.final' w l t '' lc rgb "black" lw 4
exit
