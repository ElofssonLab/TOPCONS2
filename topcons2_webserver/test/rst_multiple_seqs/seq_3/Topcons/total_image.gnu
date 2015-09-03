set style line 11 lt 1 lw 1 lc rgb "blue"
set encoding iso_8859_1
set yrange [0:50]
set xrange [1:330]
set y2range [-2:28.8]
set autoscale xfix
set ter png enh interlace size 2400,1680 font 'Nimbus,40'
set y2label '{/Symbol D}G (kcal/mol)                                             ' tc lt 3
set ytics scale 1,0.0 nomirror ("4dajA" 26.9 0, "SPOCTOPUS" 32.9 0, "SCAMPI" 35.9 0, "PolyPhobius" 38.9 0, "Philius" 41.9 0, "OCTOPUS" 44.9 0, "TOPCONS" 47.9 0)
set y2tics nomirror -2,2,12
set out '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_3///Topcons/total_image.large.png'
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
set object 6 rect from 63.5,26.4 to 79.5,26.6 fc rgb "red" fs noborder
set object 7 rect from 129.5,26.4 to 153.5,26.6 fc rgb "red" fs noborder
set object 8 rect from 212.5,26.4 to 250.5,26.6 fc rgb "red" fs noborder
set object 9 rect from 298.5,26.4 to 319.5,26.6 fc rgb "red" fs noborder
set object 10 rect from 3.5,27.2 to 36.5,27.4 fc rgb "blue" fs noborder
set object 11 rect from 42.5,27.2 to 46.5,27.4 fc rgb "blue" fs noborder
set object 12 rect from 101.5,27.2 to 108.5,27.4 fc rgb "blue" fs noborder
set object 13 rect from 174.5,27.2 to 192.5,27.4 fc rgb "blue" fs noborder
set object 14 rect from 271.5,27.2 to 277.5,27.4 fc rgb "blue" fs noborder
set object 15 rect from 46.5,26.4 to 63.5,27.4 fc rgb "white"
set object 16 rect from 79.5,26.4 to 101.5,27.4 fc rgb "grey" fs noborder
set object 17 rect from 108.5,26.4 to 129.5,27.4 fc rgb "white"
set object 18 rect from 153.5,26.4 to 174.5,27.4 fc rgb "grey" fs noborder
set object 19 rect from 192.5,26.4 to 212.5,27.4 fc rgb "white"
set object 20 rect from 250.5,26.4 to 271.5,27.4 fc rgb "grey" fs noborder
set object 21 rect from 277.5,26.4 to 298.5,27.4 fc rgb "white"
set object 22 rect from 46.5,26.4 to 63.5,27.4 fc rgb "white"
set object 23 rect from 79.5,26.4 to 101.5,27.4 fc rgb "grey" fs noborder
set object 24 rect from 108.5,26.4 to 129.5,27.4 fc rgb "white"
set object 25 rect from 153.5,26.4 to 174.5,27.4 fc rgb "grey" fs noborder
set object 26 rect from 192.5,26.4 to 212.5,27.4 fc rgb "white"
set object 27 rect from 250.5,26.4 to 271.5,27.4 fc rgb "grey" fs noborder
set object 28 rect from 277.5,26.4 to 298.5,27.4 fc rgb "white"
set object 29 rect from 63.5,32.4 to 76.5,32.6 fc rgb "red" fs noborder
set object 30 rect from 133.5,32.4 to 150.5,32.6 fc rgb "red" fs noborder
set object 31 rect from 215.5,32.4 to 247.5,32.6 fc rgb "red" fs noborder
set object 32 rect from 0.5,33.2 to 42.5,33.4 fc rgb "blue" fs noborder
set object 33 rect from 97.5,33.2 to 112.5,33.4 fc rgb "blue" fs noborder
set object 34 rect from 171.5,33.2 to 194.5,33.4 fc rgb "blue" fs noborder
set object 35 rect from 268.5,33.2 to 330.5,33.4 fc rgb "blue" fs noborder
set object 36 rect from 42.5,32.4 to 63.5,33.4 fc rgb "white"
set object 37 rect from 76.5,32.4 to 97.5,33.4 fc rgb "grey" fs noborder
set object 38 rect from 112.5,32.4 to 133.5,33.4 fc rgb "white"
set object 39 rect from 150.5,32.4 to 171.5,33.4 fc rgb "grey" fs noborder
set object 40 rect from 194.5,32.4 to 215.5,33.4 fc rgb "white"
set object 41 rect from 247.5,32.4 to 268.5,33.4 fc rgb "grey" fs noborder
set object 42 rect from 42.5,32.4 to 63.5,33.4 fc rgb "white"
set object 43 rect from 76.5,32.4 to 97.5,33.4 fc rgb "grey" fs noborder
set object 44 rect from 112.5,32.4 to 133.5,33.4 fc rgb "white"
set object 45 rect from 150.5,32.4 to 171.5,33.4 fc rgb "grey" fs noborder
set object 46 rect from 194.5,32.4 to 215.5,33.4 fc rgb "white"
set object 47 rect from 247.5,32.4 to 268.5,33.4 fc rgb "grey" fs noborder
set object 48 rect from 66.5,35.4 to 74.5,35.6 fc rgb "red" fs noborder
set object 49 rect from 130.5,35.4 to 154.5,35.6 fc rgb "red" fs noborder
set object 50 rect from 214.5,35.4 to 247.5,35.6 fc rgb "red" fs noborder
set object 51 rect from 299.5,35.4 to 330.5,35.6 fc rgb "red" fs noborder
set object 52 rect from 0.5,36.2 to 45.5,36.4 fc rgb "blue" fs noborder
set object 53 rect from 95.5,36.2 to 109.5,36.4 fc rgb "blue" fs noborder
set object 54 rect from 175.5,36.2 to 193.5,36.4 fc rgb "blue" fs noborder
set object 55 rect from 268.5,36.2 to 278.5,36.4 fc rgb "blue" fs noborder
set object 56 rect from 45.5,35.4 to 66.5,36.4 fc rgb "white"
set object 57 rect from 74.5,35.4 to 95.5,36.4 fc rgb "grey" fs noborder
set object 58 rect from 109.5,35.4 to 130.5,36.4 fc rgb "white"
set object 59 rect from 154.5,35.4 to 175.5,36.4 fc rgb "grey" fs noborder
set object 60 rect from 193.5,35.4 to 214.5,36.4 fc rgb "white"
set object 61 rect from 247.5,35.4 to 268.5,36.4 fc rgb "grey" fs noborder
set object 62 rect from 278.5,35.4 to 299.5,36.4 fc rgb "white"
set object 63 rect from 45.5,35.4 to 66.5,36.4 fc rgb "white"
set object 64 rect from 74.5,35.4 to 95.5,36.4 fc rgb "grey" fs noborder
set object 65 rect from 109.5,35.4 to 130.5,36.4 fc rgb "white"
set object 66 rect from 154.5,35.4 to 175.5,36.4 fc rgb "grey" fs noborder
set object 67 rect from 193.5,35.4 to 214.5,36.4 fc rgb "white"
set object 68 rect from 247.5,35.4 to 268.5,36.4 fc rgb "grey" fs noborder
set object 69 rect from 278.5,35.4 to 299.5,36.4 fc rgb "white"
set object 70 rect from 66.5,38.4 to 75.5,38.6 fc rgb "red" fs noborder
set object 71 rect from 132.5,38.4 to 152.5,38.6 fc rgb "red" fs noborder
set object 72 rect from 215.5,38.4 to 245.5,38.6 fc rgb "red" fs noborder
set object 73 rect from 299.5,38.4 to 330.5,38.6 fc rgb "red" fs noborder
set object 74 rect from 0.5,39.2 to 42.5,39.4 fc rgb "blue" fs noborder
set object 75 rect from 99.5,39.2 to 109.5,39.4 fc rgb "blue" fs noborder
set object 76 rect from 174.5,39.2 to 193.5,39.4 fc rgb "blue" fs noborder
set object 77 rect from 268.5,39.2 to 277.5,39.4 fc rgb "blue" fs noborder
set object 78 rect from 42.5,38.4 to 66.5,39.4 fc rgb "white"
set object 79 rect from 75.5,38.4 to 99.5,39.4 fc rgb "grey" fs noborder
set object 80 rect from 109.5,38.4 to 132.5,39.4 fc rgb "white"
set object 81 rect from 152.5,38.4 to 174.5,39.4 fc rgb "grey" fs noborder
set object 82 rect from 193.5,38.4 to 215.5,39.4 fc rgb "white"
set object 83 rect from 245.5,38.4 to 268.5,39.4 fc rgb "grey" fs noborder
set object 84 rect from 277.5,38.4 to 299.5,39.4 fc rgb "white"
set object 85 rect from 42.5,38.4 to 66.5,39.4 fc rgb "white"
set object 86 rect from 75.5,38.4 to 99.5,39.4 fc rgb "grey" fs noborder
set object 87 rect from 109.5,38.4 to 132.5,39.4 fc rgb "white"
set object 88 rect from 152.5,38.4 to 174.5,39.4 fc rgb "grey" fs noborder
set object 89 rect from 193.5,38.4 to 215.5,39.4 fc rgb "white"
set object 90 rect from 245.5,38.4 to 268.5,39.4 fc rgb "grey" fs noborder
set object 91 rect from 277.5,38.4 to 299.5,39.4 fc rgb "white"
set object 92 rect from 0.5,41.4 to 8.5,41.6 fc rgb "red" fs noborder
set object 93 rect from 66.5,41.4 to 75.5,41.6 fc rgb "red" fs noborder
set object 94 rect from 132.5,41.4 to 152.5,41.6 fc rgb "red" fs noborder
set object 95 rect from 214.5,41.4 to 244.5,41.6 fc rgb "red" fs noborder
set object 96 rect from 299.5,41.4 to 330.5,41.6 fc rgb "red" fs noborder
set object 97 rect from 24.5,42.2 to 44.5,42.4 fc rgb "blue" fs noborder
set object 98 rect from 99.5,42.2 to 109.5,42.4 fc rgb "blue" fs noborder
set object 99 rect from 175.5,42.2 to 194.5,42.4 fc rgb "blue" fs noborder
set object 100 rect from 269.5,42.2 to 276.5,42.4 fc rgb "blue" fs noborder
set object 101 rect from 8.5,41.4 to 24.5,42.4 fc rgb "grey" fs noborder
set object 102 rect from 44.5,41.4 to 66.5,42.4 fc rgb "white"
set object 103 rect from 75.5,41.4 to 99.5,42.4 fc rgb "grey" fs noborder
set object 104 rect from 109.5,41.4 to 132.5,42.4 fc rgb "white"
set object 105 rect from 152.5,41.4 to 175.5,42.4 fc rgb "grey" fs noborder
set object 106 rect from 194.5,41.4 to 214.5,42.4 fc rgb "white"
set object 107 rect from 244.5,41.4 to 269.5,42.4 fc rgb "grey" fs noborder
set object 108 rect from 276.5,41.4 to 299.5,42.4 fc rgb "white"
set object 109 rect from 8.5,41.4 to 24.5,42.4 fc rgb "grey" fs noborder
set object 110 rect from 44.5,41.4 to 66.5,42.4 fc rgb "white"
set object 111 rect from 75.5,41.4 to 99.5,42.4 fc rgb "grey" fs noborder
set object 112 rect from 109.5,41.4 to 132.5,42.4 fc rgb "white"
set object 113 rect from 152.5,41.4 to 175.5,42.4 fc rgb "grey" fs noborder
set object 114 rect from 194.5,41.4 to 214.5,42.4 fc rgb "white"
set object 115 rect from 244.5,41.4 to 269.5,42.4 fc rgb "grey" fs noborder
set object 116 rect from 276.5,41.4 to 299.5,42.4 fc rgb "white"
set object 117 rect from 63.5,44.4 to 74.5,44.6 fc rgb "red" fs noborder
set object 118 rect from 133.5,44.4 to 150.5,44.6 fc rgb "red" fs noborder
set object 119 rect from 215.5,44.4 to 247.5,44.6 fc rgb "red" fs noborder
set object 120 rect from 0.5,45.2 to 42.5,45.4 fc rgb "blue" fs noborder
set object 121 rect from 95.5,45.2 to 112.5,45.4 fc rgb "blue" fs noborder
set object 122 rect from 171.5,45.2 to 194.5,45.4 fc rgb "blue" fs noborder
set object 123 rect from 268.5,45.2 to 330.5,45.4 fc rgb "blue" fs noborder
set object 124 rect from 42.5,44.4 to 63.5,45.4 fc rgb "white"
set object 125 rect from 74.5,44.4 to 95.5,45.4 fc rgb "grey" fs noborder
set object 126 rect from 112.5,44.4 to 133.5,45.4 fc rgb "white"
set object 127 rect from 150.5,44.4 to 171.5,45.4 fc rgb "grey" fs noborder
set object 128 rect from 194.5,44.4 to 215.5,45.4 fc rgb "white"
set object 129 rect from 247.5,44.4 to 268.5,45.4 fc rgb "grey" fs noborder
set object 130 rect from 42.5,44.4 to 63.5,45.4 fc rgb "white"
set object 131 rect from 74.5,44.4 to 95.5,45.4 fc rgb "grey" fs noborder
set object 132 rect from 112.5,44.4 to 133.5,45.4 fc rgb "white"
set object 133 rect from 150.5,44.4 to 171.5,45.4 fc rgb "grey" fs noborder
set object 134 rect from 194.5,44.4 to 215.5,45.4 fc rgb "white"
set object 135 rect from 247.5,44.4 to 268.5,45.4 fc rgb "grey" fs noborder
set object 136 rect from 63.5,47.4 to 75.5,47.6 fc rgb "red" fs noborder
set object 137 rect from 132.5,47.4 to 152.5,47.6 fc rgb "red" fs noborder
set object 138 rect from 215.5,47.4 to 247.5,47.6 fc rgb "red" fs noborder
set object 139 rect from 299.5,47.4 to 330.5,47.6 fc rgb "red" fs noborder
set object 140 rect from 0.5,48.2 to 42.5,48.4 fc rgb "blue" fs noborder
set object 141 rect from 96.5,48.2 to 111.5,48.4 fc rgb "blue" fs noborder
set object 142 rect from 173.5,48.2 to 194.5,48.4 fc rgb "blue" fs noborder
set object 143 rect from 268.5,48.2 to 278.5,48.4 fc rgb "blue" fs noborder
set object 144 rect from 42.5,47.4 to 63.5,48.4 fc rgb "white"
set object 145 rect from 75.5,47.4 to 96.5,48.4 fc rgb "grey" fs noborder
set object 146 rect from 111.5,47.4 to 132.5,48.4 fc rgb "white"
set object 147 rect from 152.5,47.4 to 173.5,48.4 fc rgb "grey" fs noborder
set object 148 rect from 194.5,47.4 to 215.5,48.4 fc rgb "white"
set object 149 rect from 247.5,47.4 to 268.5,48.4 fc rgb "grey" fs noborder
set object 150 rect from 278.5,47.4 to 299.5,48.4 fc rgb "white"
set object 151 rect from 42.5,47.4 to 63.5,48.4 fc rgb "white"
set object 152 rect from 75.5,47.4 to 96.5,48.4 fc rgb "grey" fs noborder
set object 153 rect from 111.5,47.4 to 132.5,48.4 fc rgb "white"
set object 154 rect from 152.5,47.4 to 173.5,48.4 fc rgb "grey" fs noborder
set object 155 rect from 194.5,47.4 to 215.5,48.4 fc rgb "white"
set object 156 rect from 247.5,47.4 to 268.5,48.4 fc rgb "grey" fs noborder
set object 157 rect from 278.5,47.4 to 299.5,48.4 fc rgb "white"
plot '/big/tmp/tmp1/TOPCONS2/topcons2_webserver/test/rst2//seq_3///DG1.txt' axes x1y2 w l t '' lt 3 lw 4
exit
