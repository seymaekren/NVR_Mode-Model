# Rigid cluster decomposition coloring script for PyMol
#
# Created by Dan Farrell, Brandon Hespenheide.
# Department of Physics and Astronomy
# Biophysics Theory Group
# Arizona State University
############################################################

from pymol import cmd
from pymol.cgo import *

bg_color white
# script to make a movie of the FRODA generated structures.

from glob import glob
filelist = glob ("1EF1-NMA.7.model6_froda_*.pdb") 
filelist.sort()
cmd.load("1EF1-NMA.7.model6_RCD.pdb", "1EF1-NMA.7.model6_froda")
for file in filelist: cmd.load( file, "1EF1-NMA.7.model6_froda") 
show lines, 1EF1-NMA.7.model6_froda
color black
color 0x0000b2, ( b > 0.99 and b < 1.01)
color 0x57eb0f, ( b > 1.99 and b < 2.01)
# Draw hbonds and hydrophoic tethers as pymol distance objects
set dash_gap, 0.1
distance hbonds = id 3 , id 13
distance hbonds = id 57 , id 1078
distance hbonds = id 93 , id 1114
distance hbonds = id 162 , id 120
distance hbonds = id 340 , id 301
distance hbonds = id 354 , id 315
distance hbonds = id 373 , id 326
distance hbonds = id 388 , id 338
distance hbonds = id 402 , id 352
distance hbonds = id 418 , id 371
distance hbonds = id 440 , id 386
distance hbonds = id 450 , id 400
distance hbonds = id 472 , id 416
distance hbonds = id 508 , id 438
distance hbonds = id 611 , id 1173
distance hbonds = id 640 , id 1130
distance hbonds = id 719 , id 961
distance hbonds = id 886 , id 286
distance hbonds = id 928 , id 884
distance hbonds = id 928 , id 903
distance hbonds = id 982 , id 698
distance hbonds = id 1061 , id 36
distance hbonds = id 1084 , id 1032
distance hbonds = id 1097 , id 75
distance hbonds = id 1116 , id 655
distance hbonds = id 1132 , id 113
distance hbonds = id 1151 , id 621
distance hbonds = id 1194 , id 595
distance hbonds = id 1207 , id 609
color red, hbonds
hide labels, hbonds
disable hbonds
distance hydrophobics = id 10 , id 1073
distance hydrophobics = id 13 , id 217
distance hydrophobics = id 43 , id 214
distance hydrophobics = id 43 , id 210
distance hydrophobics = id 47 , id 217
distance hydrophobics = id 47 , id 894
distance hydrophobics = id 50 , id 1109
distance hydrophobics = id 50 , id 86
distance hydrophobics = id 60 , id 1054
distance hydrophobics = id 66 , id 175
distance hydrophobics = id 68 , id 96
distance hydrophobics = id 80 , id 172
distance hydrophobics = id 82 , id 411
distance hydrophobics = id 122 , id 407
distance hydrophobics = id 122 , id 1135
distance hydrophobics = id 217 , id 277
distance hydrophobics = id 217 , id 894
distance hydrophobics = id 231 , id 263
distance hydrophobics = id 281 , id 331
distance hydrophobics = id 294 , id 331
distance hydrophobics = id 308 , id 803
distance hydrophobics = id 311 , id 939
distance hydrophobics = id 311 , id 807
distance hydrophobics = id 614 , id 1200
distance hydrophobics = id 629 , id 803
distance hydrophobics = id 643 , id 1121
distance hydrophobics = id 660 , id 1103
distance hydrophobics = id 705 , id 952
distance hydrophobics = id 712 , id 759
distance hydrophobics = id 734 , id 1002
distance hydrophobics = id 734 , id 1005
distance hydrophobics = id 803 , id 1105
distance hydrophobics = id 807 , id 898
distance hydrophobics = id 851 , id 939
distance hydrophobics = id 879 , id 919
distance hydrophobics = id 894 , id 1073
distance hydrophobics = id 972 , id 1069
distance hydrophobics = id 1125 , id 1157
color green, hydrophobics
hide labels, hydrophobics
disable hydrophobics
# Rigid Cluster 1 has 77 atoms.
create RC1, ( b > 0.99 and b < 1.01)
show sticks, RC1
set line_width = 3, RC1
color 0x0000b2, RC1

# Rigid Cluster 2 has 21 atoms.
create RC2, ( b > 1.99 and b < 2.01)
show sticks, RC2
set line_width = 3, RC2
color 0x57eb0f, RC2

# Rigid Cluster BIN2
create BIN2, ( b > 2.99 and b < 8.01)
show sticks, BIN2
set line_width = 3, BIN2
color gray, BIN2
disable BIN2

# Rigid Cluster BIN1
create BIN1, ( b > 8.99 and b < 371.01)
show sticks, BIN1
set line_width = 3, BIN1
color gray, BIN1
disable BIN1

