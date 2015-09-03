set style line 11 lt 1 lw 1 lc rgb "blue"
set encoding iso_8859_1
set yrange [0:50]
set xrange [1:334]
set y2range [-2:28.8]
set autoscale xfix
set ter png enh interlace size 2400,1680 font 'Nimbus,40'
set y2label '{/Symbol D}G (kcal/mol)                                             ' tc lt 3
set ytics scale 1,0.0 nomirror ("2y02A" 26.9 0, "SPOCTOPUS" 32.9 0, "SCAMPI" 35.9 0, "PolyPhobius" 38.9 0, "Philius" 41.9 0, "OCTOPUS" 44.9 0, "TOPCONS" 47.9 0)
set y2tics nomirror -2,2,12
set out '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_2///Topcons/total_image.large.png'
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
set object 6 rect from 66.5,26.4 to 84.5,26.6 fc rgb "red" fs noborder
set object 7 rect from 132.5,26.4 to 160.5,26.6 fc rgb "red" fs noborder
set object 8 rect from 216.5,26.4 to 254.5,26.6 fc rgb "red" fs noborder
set object 9 rect from 303.5,26.4 to 321.5,26.6 fc rgb "red" fs noborder
set object 10 rect from 11.5,27.2 to 35.5,27.4 fc rgb "blue" fs noborder
set object 11 rect from 36.5,27.2 to 48.5,27.4 fc rgb "blue" fs noborder
set object 12 rect from 107.5,27.2 to 112.5,27.4 fc rgb "blue" fs noborder
set object 13 rect from 179.5,27.2 to 189.5,27.4 fc rgb "blue" fs noborder
set object 14 rect from 274.5,27.2 to 280.5,27.4 fc rgb "blue" fs noborder
set object 15 rect from 48.5,26.4 to 66.5,27.4 fc rgb "white"
set object 16 rect from 84.5,26.4 to 95.5,27.4 fc rgb "grey" fs noborder
set object 17 rect from 96.5,26.4 to 107.5,27.4 fc rgb "grey" fs noborder
set object 18 rect from 112.5,26.4 to 132.5,27.4 fc rgb "white"
set object 19 rect from 160.5,26.4 to 179.5,27.4 fc rgb "grey" fs noborder
set object 20 rect from 189.5,26.4 to 195.5,27.4 fc rgb "white"
set object 21 rect from 196.5,26.4 to 216.5,27.4 fc rgb "white"
set object 22 rect from 254.5,26.4 to 274.5,27.4 fc rgb "grey" fs noborder
set object 23 rect from 280.5,26.4 to 288.5,27.4 fc rgb "white"
set object 24 rect from 289.5,26.4 to 302.5,27.4 fc rgb "white"
set object 25 rect from 48.5,26.4 to 66.5,27.4 fc rgb "white"
set object 26 rect from 84.5,26.4 to 107.5,27.4 fc rgb "grey" fs noborder
set object 27 rect from 112.5,26.4 to 132.5,27.4 fc rgb "white"
set object 28 rect from 160.5,26.4 to 179.5,27.4 fc rgb "grey" fs noborder
set object 29 rect from 189.5,26.4 to 216.5,27.4 fc rgb "white"
set object 30 rect from 254.5,26.4 to 274.5,27.4 fc rgb "grey" fs noborder
set object 31 rect from 280.5,26.4 to 302.5,27.4 fc rgb "white"
set object 32 rect from 67.5,32.4 to 81.5,32.6 fc rgb "red" fs noborder
set object 33 rect from 137.5,32.4 to 154.5,32.6 fc rgb "red" fs noborder
set object 34 rect from 219.5,32.4 to 249.5,32.6 fc rgb "red" fs noborder
set object 35 rect from 0.5,33.2 to 46.5,33.4 fc rgb "blue" fs noborder
set object 36 rect from 102.5,33.2 to 116.5,33.4 fc rgb "blue" fs noborder
set object 37 rect from 175.5,33.2 to 198.5,33.4 fc rgb "blue" fs noborder
set object 38 rect from 270.5,33.2 to 334.5,33.4 fc rgb "blue" fs noborder
set object 39 rect from 46.5,32.4 to 67.5,33.4 fc rgb "white"
set object 40 rect from 81.5,32.4 to 102.5,33.4 fc rgb "grey" fs noborder
set object 41 rect from 116.5,32.4 to 137.5,33.4 fc rgb "white"
set object 42 rect from 154.5,32.4 to 175.5,33.4 fc rgb "grey" fs noborder
set object 43 rect from 198.5,32.4 to 219.5,33.4 fc rgb "white"
set object 44 rect from 249.5,32.4 to 270.5,33.4 fc rgb "grey" fs noborder
set object 45 rect from 46.5,32.4 to 67.5,33.4 fc rgb "white"
set object 46 rect from 81.5,32.4 to 102.5,33.4 fc rgb "grey" fs noborder
set object 47 rect from 116.5,32.4 to 137.5,33.4 fc rgb "white"
set object 48 rect from 154.5,32.4 to 175.5,33.4 fc rgb "grey" fs noborder
set object 49 rect from 198.5,32.4 to 219.5,33.4 fc rgb "white"
set object 50 rect from 249.5,32.4 to 270.5,33.4 fc rgb "grey" fs noborder
set object 51 rect from 70.5,35.4 to 78.5,35.6 fc rgb "red" fs noborder
set object 52 rect from 134.5,35.4 to 158.5,35.6 fc rgb "red" fs noborder
set object 53 rect from 218.5,35.4 to 251.5,35.6 fc rgb "red" fs noborder
set object 54 rect from 303.5,35.4 to 334.5,35.6 fc rgb "red" fs noborder
set object 55 rect from 0.5,36.2 to 49.5,36.4 fc rgb "blue" fs noborder
set object 56 rect from 99.5,36.2 to 113.5,36.4 fc rgb "blue" fs noborder
set object 57 rect from 179.5,36.2 to 197.5,36.4 fc rgb "blue" fs noborder
set object 58 rect from 272.5,36.2 to 282.5,36.4 fc rgb "blue" fs noborder
set object 59 rect from 49.5,35.4 to 70.5,36.4 fc rgb "white"
set object 60 rect from 78.5,35.4 to 99.5,36.4 fc rgb "grey" fs noborder
set object 61 rect from 113.5,35.4 to 134.5,36.4 fc rgb "white"
set object 62 rect from 158.5,35.4 to 179.5,36.4 fc rgb "grey" fs noborder
set object 63 rect from 197.5,35.4 to 218.5,36.4 fc rgb "white"
set object 64 rect from 251.5,35.4 to 272.5,36.4 fc rgb "grey" fs noborder
set object 65 rect from 282.5,35.4 to 303.5,36.4 fc rgb "white"
set object 66 rect from 49.5,35.4 to 70.5,36.4 fc rgb "white"
set object 67 rect from 78.5,35.4 to 99.5,36.4 fc rgb "grey" fs noborder
set object 68 rect from 113.5,35.4 to 134.5,36.4 fc rgb "white"
set object 69 rect from 158.5,35.4 to 179.5,36.4 fc rgb "grey" fs noborder
set object 70 rect from 197.5,35.4 to 218.5,36.4 fc rgb "white"
set object 71 rect from 251.5,35.4 to 272.5,36.4 fc rgb "grey" fs noborder
set object 72 rect from 282.5,35.4 to 303.5,36.4 fc rgb "white"
set object 73 rect from 70.5,38.4 to 79.5,38.6 fc rgb "red" fs noborder
set object 74 rect from 136.5,38.4 to 156.5,38.6 fc rgb "red" fs noborder
set object 75 rect from 219.5,38.4 to 249.5,38.6 fc rgb "red" fs noborder
set object 76 rect from 303.5,38.4 to 334.5,38.6 fc rgb "red" fs noborder
set object 77 rect from 0.5,39.2 to 46.5,39.4 fc rgb "blue" fs noborder
set object 78 rect from 104.5,39.2 to 113.5,39.4 fc rgb "blue" fs noborder
set object 79 rect from 178.5,39.2 to 198.5,39.4 fc rgb "blue" fs noborder
set object 80 rect from 270.5,39.2 to 281.5,39.4 fc rgb "blue" fs noborder
set object 81 rect from 46.5,38.4 to 70.5,39.4 fc rgb "white"
set object 82 rect from 79.5,38.4 to 104.5,39.4 fc rgb "grey" fs noborder
set object 83 rect from 113.5,38.4 to 136.5,39.4 fc rgb "white"
set object 84 rect from 156.5,38.4 to 178.5,39.4 fc rgb "grey" fs noborder
set object 85 rect from 198.5,38.4 to 219.5,39.4 fc rgb "white"
set object 86 rect from 249.5,38.4 to 270.5,39.4 fc rgb "grey" fs noborder
set object 87 rect from 281.5,38.4 to 303.5,39.4 fc rgb "white"
set object 88 rect from 46.5,38.4 to 70.5,39.4 fc rgb "white"
set object 89 rect from 79.5,38.4 to 104.5,39.4 fc rgb "grey" fs noborder
set object 90 rect from 113.5,38.4 to 136.5,39.4 fc rgb "white"
set object 91 rect from 156.5,38.4 to 178.5,39.4 fc rgb "grey" fs noborder
set object 92 rect from 198.5,38.4 to 219.5,39.4 fc rgb "white"
set object 93 rect from 249.5,38.4 to 270.5,39.4 fc rgb "grey" fs noborder
set object 94 rect from 281.5,38.4 to 303.5,39.4 fc rgb "white"
set object 95 rect from 0.5,41.4 to 79.5,41.6 fc rgb "red" fs noborder
set object 96 rect from 136.5,41.4 to 153.5,41.6 fc rgb "red" fs noborder
set object 97 rect from 220.5,41.4 to 249.5,41.6 fc rgb "red" fs noborder
set object 98 rect from 301.5,41.4 to 334.5,41.6 fc rgb "red" fs noborder
set object 99 rect from 104.5,42.2 to 113.5,42.4 fc rgb "blue" fs noborder
set object 100 rect from 177.5,42.2 to 198.5,42.4 fc rgb "blue" fs noborder
set object 101 rect from 273.5,42.2 to 280.5,42.4 fc rgb "blue" fs noborder
set object 102 rect from 79.5,41.4 to 104.5,42.4 fc rgb "grey" fs noborder
set object 103 rect from 113.5,41.4 to 136.5,42.4 fc rgb "white"
set object 104 rect from 153.5,41.4 to 177.5,42.4 fc rgb "grey" fs noborder
set object 105 rect from 198.5,41.4 to 220.5,42.4 fc rgb "white"
set object 106 rect from 249.5,41.4 to 273.5,42.4 fc rgb "grey" fs noborder
set object 107 rect from 280.5,41.4 to 301.5,42.4 fc rgb "white"
set object 108 rect from 79.5,41.4 to 104.5,42.4 fc rgb "grey" fs noborder
set object 109 rect from 113.5,41.4 to 136.5,42.4 fc rgb "white"
set object 110 rect from 153.5,41.4 to 177.5,42.4 fc rgb "grey" fs noborder
set object 111 rect from 198.5,41.4 to 220.5,42.4 fc rgb "white"
set object 112 rect from 249.5,41.4 to 273.5,42.4 fc rgb "grey" fs noborder
set object 113 rect from 280.5,41.4 to 301.5,42.4 fc rgb "white"
set object 114 rect from 67.5,44.4 to 81.5,44.6 fc rgb "red" fs noborder
set object 115 rect from 137.5,44.4 to 154.5,44.6 fc rgb "red" fs noborder
set object 116 rect from 219.5,44.4 to 249.5,44.6 fc rgb "red" fs noborder
set object 117 rect from 0.5,45.2 to 46.5,45.4 fc rgb "blue" fs noborder
set object 118 rect from 102.5,45.2 to 116.5,45.4 fc rgb "blue" fs noborder
set object 119 rect from 175.5,45.2 to 198.5,45.4 fc rgb "blue" fs noborder
set object 120 rect from 270.5,45.2 to 334.5,45.4 fc rgb "blue" fs noborder
set object 121 rect from 46.5,44.4 to 67.5,45.4 fc rgb "white"
set object 122 rect from 81.5,44.4 to 102.5,45.4 fc rgb "grey" fs noborder
set object 123 rect from 116.5,44.4 to 137.5,45.4 fc rgb "white"
set object 124 rect from 154.5,44.4 to 175.5,45.4 fc rgb "grey" fs noborder
set object 125 rect from 198.5,44.4 to 219.5,45.4 fc rgb "white"
set object 126 rect from 249.5,44.4 to 270.5,45.4 fc rgb "grey" fs noborder
set object 127 rect from 46.5,44.4 to 67.5,45.4 fc rgb "white"
set object 128 rect from 81.5,44.4 to 102.5,45.4 fc rgb "grey" fs noborder
set object 129 rect from 116.5,44.4 to 137.5,45.4 fc rgb "white"
set object 130 rect from 154.5,44.4 to 175.5,45.4 fc rgb "grey" fs noborder
set object 131 rect from 198.5,44.4 to 219.5,45.4 fc rgb "white"
set object 132 rect from 249.5,44.4 to 270.5,45.4 fc rgb "grey" fs noborder
set object 133 rect from 67.5,47.4 to 81.5,47.6 fc rgb "red" fs noborder
set object 134 rect from 136.5,47.4 to 154.5,47.6 fc rgb "red" fs noborder
set object 135 rect from 219.5,47.4 to 249.5,47.6 fc rgb "red" fs noborder
set object 136 rect from 303.5,47.4 to 334.5,47.6 fc rgb "red" fs noborder
set object 137 rect from 0.5,48.2 to 46.5,48.4 fc rgb "blue" fs noborder
set object 138 rect from 102.5,48.2 to 115.5,48.4 fc rgb "blue" fs noborder
set object 139 rect from 175.5,48.2 to 198.5,48.4 fc rgb "blue" fs noborder
set object 140 rect from 270.5,48.2 to 282.5,48.4 fc rgb "blue" fs noborder
set object 141 rect from 46.5,47.4 to 67.5,48.4 fc rgb "white"
set object 142 rect from 81.5,47.4 to 102.5,48.4 fc rgb "grey" fs noborder
set object 143 rect from 115.5,47.4 to 136.5,48.4 fc rgb "white"
set object 144 rect from 154.5,47.4 to 175.5,48.4 fc rgb "grey" fs noborder
set object 145 rect from 198.5,47.4 to 219.5,48.4 fc rgb "white"
set object 146 rect from 249.5,47.4 to 270.5,48.4 fc rgb "grey" fs noborder
set object 147 rect from 282.5,47.4 to 303.5,48.4 fc rgb "white"
set object 148 rect from 46.5,47.4 to 67.5,48.4 fc rgb "white"
set object 149 rect from 81.5,47.4 to 102.5,48.4 fc rgb "grey" fs noborder
set object 150 rect from 115.5,47.4 to 136.5,48.4 fc rgb "white"
set object 151 rect from 154.5,47.4 to 175.5,48.4 fc rgb "grey" fs noborder
set object 152 rect from 198.5,47.4 to 219.5,48.4 fc rgb "white"
set object 153 rect from 249.5,47.4 to 270.5,48.4 fc rgb "grey" fs noborder
set object 154 rect from 282.5,47.4 to 303.5,48.4 fc rgb "white"
plot '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_2///DG1.txt' axes x1y2 w l t '' lt 3 lw 4
exit
