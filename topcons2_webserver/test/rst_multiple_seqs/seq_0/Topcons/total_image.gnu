set style line 11 lt 1 lw 1 lc rgb "blue"
set encoding iso_8859_1
set yrange [0:50]
set xrange [1:359]
set y2range [-2:36.6]
set autoscale xfix
set ter png enh interlace size 2400,1680 font 'Nimbus,40'
set y2label '{/Symbol D}G (kcal/mol)                                             ' tc lt 3
set ytics scale 1,0.0 nomirror ("2lnlA" 26.9 0, "SPOCTOPUS" 32.9 0, "SCAMPI" 35.9 0, "PolyPhobius" 38.9 0, "Philius" 41.9 0, "OCTOPUS" 44.9 0, "TOPCONS" 47.9 0)
set y2tics nomirror -3,2,15
set out '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_0///Topcons/total_image.large.png'
set lmargin 13.5
set rmargin 6.5
set tmargin 1.3
set object 1 rect from screen 0.19,0.986 to screen 0.21,0.992 fc rgb "red" fs noborder
set label 'Inside' font 'Nimbus,30' at screen 0.215,0.982
set object 2 rect from screen 0.28,0.986 to screen 0.30,0.992 fc rgb "blue" fs noborder
set label 'Outside' font 'Nimbus,30' at screen 0.305,0.982
set object 3 rect from screen 0.38,0.978 to screen 0.40,1 fc rgb "grey" fs noborder
set label 'TM-helix (IN->OUT)' font 'Nimbus,30' at screen 0.405,0.982
set object 4 rect from screen 0.57,0.978 to screen 0.59,1 fc rgb "white"
set label 'TM-helix (OUT->IN)' font 'Nimbus,30' at screen 0.595,0.982
set object 5 rect from screen 0.76,0.978 to screen 0.78,1 fc rgb "black"
set label 'Signal peptide' font 'Nimbus,30' at screen 0.785,0.982
set ytic scale 0
set object 6 rect from 71.5,26.4 to 85.5,26.6 fc rgb "red" fs noborder
set object 7 rect from 139.5,26.4 to 163.5,26.6 fc rgb "red" fs noborder
set object 8 rect from 229.5,26.4 to 253.5,26.6 fc rgb "red" fs noborder
set object 9 rect from 316.5,26.4 to 349.5,26.6 fc rgb "red" fs noborder
set object 10 rect from 350.5,26.4 to 359.5,26.6 fc rgb "red" fs noborder
set object 11 rect from 9.5,27.2 to 24.5,27.4 fc rgb "blue" fs noborder
set object 12 rect from 25.5,27.2 to 49.5,27.4 fc rgb "blue" fs noborder
set object 13 rect from 105.5,27.2 to 117.5,27.4 fc rgb "blue" fs noborder
set object 14 rect from 183.5,27.2 to 207.5,27.4 fc rgb "blue" fs noborder
set object 15 rect from 275.5,27.2 to 294.5,27.4 fc rgb "blue" fs noborder
set object 16 rect from 49.5,26.4 to 71.5,27.4 fc rgb "white"
set object 17 rect from 85.5,26.4 to 105.5,27.4 fc rgb "grey" fs noborder
set object 18 rect from 117.5,26.4 to 139.5,27.4 fc rgb "white"
set object 19 rect from 163.5,26.4 to 183.5,27.4 fc rgb "grey" fs noborder
set object 20 rect from 207.5,26.4 to 229.5,27.4 fc rgb "white"
set object 21 rect from 253.5,26.4 to 275.5,27.4 fc rgb "grey" fs noborder
set object 22 rect from 294.5,26.4 to 316.5,27.4 fc rgb "white"
set object 23 rect from 49.5,26.4 to 71.5,27.4 fc rgb "white"
set object 24 rect from 85.5,26.4 to 105.5,27.4 fc rgb "grey" fs noborder
set object 25 rect from 117.5,26.4 to 139.5,27.4 fc rgb "white"
set object 26 rect from 163.5,26.4 to 183.5,27.4 fc rgb "grey" fs noborder
set object 27 rect from 207.5,26.4 to 229.5,27.4 fc rgb "white"
set object 28 rect from 253.5,26.4 to 275.5,27.4 fc rgb "grey" fs noborder
set object 29 rect from 294.5,26.4 to 316.5,27.4 fc rgb "white"
set object 30 rect from 68.5,32.4 to 84.5,32.6 fc rgb "red" fs noborder
set object 31 rect from 145.5,32.4 to 159.5,32.6 fc rgb "red" fs noborder
set object 32 rect from 233.5,32.4 to 248.5,32.6 fc rgb "red" fs noborder
set object 33 rect from 318.5,32.4 to 359.5,32.6 fc rgb "red" fs noborder
set object 34 rect from 0.5,33.2 to 47.5,33.4 fc rgb "blue" fs noborder
set object 35 rect from 115.5,33.2 to 124.5,33.4 fc rgb "blue" fs noborder
set object 36 rect from 180.5,33.2 to 212.5,33.4 fc rgb "blue" fs noborder
set object 37 rect from 269.5,33.2 to 297.5,33.4 fc rgb "blue" fs noborder
set object 38 rect from 47.5,32.4 to 68.5,33.4 fc rgb "white"
set object 39 rect from 84.5,32.4 to 115.5,33.4 fc rgb "grey" fs noborder
set object 40 rect from 124.5,32.4 to 145.5,33.4 fc rgb "white"
set object 41 rect from 159.5,32.4 to 180.5,33.4 fc rgb "grey" fs noborder
set object 42 rect from 212.5,32.4 to 233.5,33.4 fc rgb "white"
set object 43 rect from 248.5,32.4 to 269.5,33.4 fc rgb "grey" fs noborder
set object 44 rect from 297.5,32.4 to 318.5,33.4 fc rgb "white"
set object 45 rect from 47.5,32.4 to 68.5,33.4 fc rgb "white"
set object 46 rect from 84.5,32.4 to 115.5,33.4 fc rgb "grey" fs noborder
set object 47 rect from 124.5,32.4 to 145.5,33.4 fc rgb "white"
set object 48 rect from 159.5,32.4 to 180.5,33.4 fc rgb "grey" fs noborder
set object 49 rect from 212.5,32.4 to 233.5,33.4 fc rgb "white"
set object 50 rect from 248.5,32.4 to 269.5,33.4 fc rgb "grey" fs noborder
set object 51 rect from 297.5,32.4 to 318.5,33.4 fc rgb "white"
set object 52 rect from 72.5,35.4 to 84.5,35.6 fc rgb "red" fs noborder
set object 53 rect from 147.5,35.4 to 162.5,35.6 fc rgb "red" fs noborder
set object 54 rect from 234.5,35.4 to 250.5,35.6 fc rgb "red" fs noborder
set object 55 rect from 316.5,35.4 to 359.5,35.6 fc rgb "red" fs noborder
set object 56 rect from 0.5,36.2 to 51.5,36.4 fc rgb "blue" fs noborder
set object 57 rect from 105.5,36.2 to 126.5,36.4 fc rgb "blue" fs noborder
set object 58 rect from 183.5,36.2 to 213.5,36.4 fc rgb "blue" fs noborder
set object 59 rect from 271.5,36.2 to 295.5,36.4 fc rgb "blue" fs noborder
set object 60 rect from 51.5,35.4 to 72.5,36.4 fc rgb "white"
set object 61 rect from 84.5,35.4 to 105.5,36.4 fc rgb "grey" fs noborder
set object 62 rect from 126.5,35.4 to 147.5,36.4 fc rgb "white"
set object 63 rect from 162.5,35.4 to 183.5,36.4 fc rgb "grey" fs noborder
set object 64 rect from 213.5,35.4 to 234.5,36.4 fc rgb "white"
set object 65 rect from 250.5,35.4 to 271.5,36.4 fc rgb "grey" fs noborder
set object 66 rect from 295.5,35.4 to 316.5,36.4 fc rgb "white"
set object 67 rect from 51.5,35.4 to 72.5,36.4 fc rgb "white"
set object 68 rect from 84.5,35.4 to 105.5,36.4 fc rgb "grey" fs noborder
set object 69 rect from 126.5,35.4 to 147.5,36.4 fc rgb "white"
set object 70 rect from 162.5,35.4 to 183.5,36.4 fc rgb "grey" fs noborder
set object 71 rect from 213.5,35.4 to 234.5,36.4 fc rgb "white"
set object 72 rect from 250.5,35.4 to 271.5,36.4 fc rgb "grey" fs noborder
set object 73 rect from 295.5,35.4 to 316.5,36.4 fc rgb "white"
set object 74 rect from 73.5,38.4 to 84.5,38.6 fc rgb "red" fs noborder
set object 75 rect from 148.5,38.4 to 159.5,38.6 fc rgb "red" fs noborder
set object 76 rect from 237.5,38.4 to 249.5,38.6 fc rgb "red" fs noborder
set object 77 rect from 316.5,38.4 to 359.5,38.6 fc rgb "red" fs noborder
set object 78 rect from 0.5,39.2 to 50.5,39.4 fc rgb "blue" fs noborder
set object 79 rect from 105.5,39.2 to 126.5,39.4 fc rgb "blue" fs noborder
set object 80 rect from 182.5,39.2 to 214.5,39.4 fc rgb "blue" fs noborder
set object 81 rect from 270.5,39.2 to 296.5,39.4 fc rgb "blue" fs noborder
set object 82 rect from 50.5,38.4 to 73.5,39.4 fc rgb "white"
set object 83 rect from 84.5,38.4 to 105.5,39.4 fc rgb "grey" fs noborder
set object 84 rect from 126.5,38.4 to 148.5,39.4 fc rgb "white"
set object 85 rect from 159.5,38.4 to 182.5,39.4 fc rgb "grey" fs noborder
set object 86 rect from 214.5,38.4 to 237.5,39.4 fc rgb "white"
set object 87 rect from 249.5,38.4 to 270.5,39.4 fc rgb "grey" fs noborder
set object 88 rect from 296.5,38.4 to 316.5,39.4 fc rgb "white"
set object 89 rect from 50.5,38.4 to 73.5,39.4 fc rgb "white"
set object 90 rect from 84.5,38.4 to 105.5,39.4 fc rgb "grey" fs noborder
set object 91 rect from 126.5,38.4 to 148.5,39.4 fc rgb "white"
set object 92 rect from 159.5,38.4 to 182.5,39.4 fc rgb "grey" fs noborder
set object 93 rect from 214.5,38.4 to 237.5,39.4 fc rgb "white"
set object 94 rect from 249.5,38.4 to 270.5,39.4 fc rgb "grey" fs noborder
set object 95 rect from 296.5,38.4 to 316.5,39.4 fc rgb "white"
set object 96 rect from 74.5,41.4 to 83.5,41.6 fc rgb "red" fs noborder
set object 97 rect from 147.5,41.4 to 159.5,41.6 fc rgb "red" fs noborder
set object 98 rect from 238.5,41.4 to 247.5,41.6 fc rgb "red" fs noborder
set object 99 rect from 317.5,41.4 to 359.5,41.6 fc rgb "red" fs noborder
set object 100 rect from 0.5,42.2 to 49.5,42.4 fc rgb "blue" fs noborder
set object 101 rect from 105.5,42.2 to 126.5,42.4 fc rgb "blue" fs noborder
set object 102 rect from 182.5,42.2 to 215.5,42.4 fc rgb "blue" fs noborder
set object 103 rect from 270.5,42.2 to 300.5,42.4 fc rgb "blue" fs noborder
set object 104 rect from 49.5,41.4 to 74.5,42.4 fc rgb "white"
set object 105 rect from 83.5,41.4 to 105.5,42.4 fc rgb "grey" fs noborder
set object 106 rect from 126.5,41.4 to 147.5,42.4 fc rgb "white"
set object 107 rect from 159.5,41.4 to 182.5,42.4 fc rgb "grey" fs noborder
set object 108 rect from 215.5,41.4 to 238.5,42.4 fc rgb "white"
set object 109 rect from 247.5,41.4 to 270.5,42.4 fc rgb "grey" fs noborder
set object 110 rect from 300.5,41.4 to 317.5,42.4 fc rgb "white"
set object 111 rect from 49.5,41.4 to 74.5,42.4 fc rgb "white"
set object 112 rect from 83.5,41.4 to 105.5,42.4 fc rgb "grey" fs noborder
set object 113 rect from 126.5,41.4 to 147.5,42.4 fc rgb "white"
set object 114 rect from 159.5,41.4 to 182.5,42.4 fc rgb "grey" fs noborder
set object 115 rect from 215.5,41.4 to 238.5,42.4 fc rgb "white"
set object 116 rect from 247.5,41.4 to 270.5,42.4 fc rgb "grey" fs noborder
set object 117 rect from 300.5,41.4 to 317.5,42.4 fc rgb "white"
set object 118 rect from 68.5,44.4 to 84.5,44.6 fc rgb "red" fs noborder
set object 119 rect from 145.5,44.4 to 159.5,44.6 fc rgb "red" fs noborder
set object 120 rect from 233.5,44.4 to 248.5,44.6 fc rgb "red" fs noborder
set object 121 rect from 318.5,44.4 to 359.5,44.6 fc rgb "red" fs noborder
set object 122 rect from 0.5,45.2 to 47.5,45.4 fc rgb "blue" fs noborder
set object 123 rect from 115.5,45.2 to 124.5,45.4 fc rgb "blue" fs noborder
set object 124 rect from 180.5,45.2 to 212.5,45.4 fc rgb "blue" fs noborder
set object 125 rect from 269.5,45.2 to 297.5,45.4 fc rgb "blue" fs noborder
set object 126 rect from 47.5,44.4 to 68.5,45.4 fc rgb "white"
set object 127 rect from 84.5,44.4 to 115.5,45.4 fc rgb "grey" fs noborder
set object 128 rect from 124.5,44.4 to 145.5,45.4 fc rgb "white"
set object 129 rect from 159.5,44.4 to 180.5,45.4 fc rgb "grey" fs noborder
set object 130 rect from 212.5,44.4 to 233.5,45.4 fc rgb "white"
set object 131 rect from 248.5,44.4 to 269.5,45.4 fc rgb "grey" fs noborder
set object 132 rect from 297.5,44.4 to 318.5,45.4 fc rgb "white"
set object 133 rect from 47.5,44.4 to 68.5,45.4 fc rgb "white"
set object 134 rect from 84.5,44.4 to 115.5,45.4 fc rgb "grey" fs noborder
set object 135 rect from 124.5,44.4 to 145.5,45.4 fc rgb "white"
set object 136 rect from 159.5,44.4 to 180.5,45.4 fc rgb "grey" fs noborder
set object 137 rect from 212.5,44.4 to 233.5,45.4 fc rgb "white"
set object 138 rect from 248.5,44.4 to 269.5,45.4 fc rgb "grey" fs noborder
set object 139 rect from 297.5,44.4 to 318.5,45.4 fc rgb "white"
set object 140 rect from 70.5,47.4 to 84.5,47.6 fc rgb "red" fs noborder
set object 141 rect from 147.5,47.4 to 159.5,47.6 fc rgb "red" fs noborder
set object 142 rect from 234.5,47.4 to 248.5,47.6 fc rgb "red" fs noborder
set object 143 rect from 317.5,47.4 to 359.5,47.6 fc rgb "red" fs noborder
set object 144 rect from 0.5,48.2 to 49.5,48.4 fc rgb "blue" fs noborder
set object 145 rect from 105.5,48.2 to 126.5,48.4 fc rgb "blue" fs noborder
set object 146 rect from 180.5,48.2 to 213.5,48.4 fc rgb "blue" fs noborder
set object 147 rect from 269.5,48.2 to 296.5,48.4 fc rgb "blue" fs noborder
set object 148 rect from 49.5,47.4 to 70.5,48.4 fc rgb "white"
set object 149 rect from 84.5,47.4 to 105.5,48.4 fc rgb "grey" fs noborder
set object 150 rect from 126.5,47.4 to 147.5,48.4 fc rgb "white"
set object 151 rect from 159.5,47.4 to 180.5,48.4 fc rgb "grey" fs noborder
set object 152 rect from 213.5,47.4 to 234.5,48.4 fc rgb "white"
set object 153 rect from 248.5,47.4 to 269.5,48.4 fc rgb "grey" fs noborder
set object 154 rect from 296.5,47.4 to 317.5,48.4 fc rgb "white"
set object 155 rect from 49.5,47.4 to 70.5,48.4 fc rgb "white"
set object 156 rect from 84.5,47.4 to 105.5,48.4 fc rgb "grey" fs noborder
set object 157 rect from 126.5,47.4 to 147.5,48.4 fc rgb "white"
set object 158 rect from 159.5,47.4 to 180.5,48.4 fc rgb "grey" fs noborder
set object 159 rect from 213.5,47.4 to 234.5,48.4 fc rgb "white"
set object 160 rect from 248.5,47.4 to 269.5,48.4 fc rgb "grey" fs noborder
set object 161 rect from 296.5,47.4 to 317.5,48.4 fc rgb "white"
plot '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_0///DG1.txt' axes x1y2 w l t '' lt 3 lw 4
exit
