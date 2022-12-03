mol default style {VDW}
#mol default material {Diffuse}
display projection orthographic
# position the stage and axes
axes location off
# stage location off
#mol new  E:/VMD/triJanus/*.pdb
set list [ eval exec ls [glob *.pdb] ]
#load pdb E:/test/tri56.pdb
for { set i 0 } { $i < [llength $list] } { incr i } {
mol new [lindex $list $i] 
mol off $i
set sel [atomselect $i "name N"]
$sel set radius 0.30
set sel [atomselect $i "name O"]
$sel set radius 0.2975
## change the color of background and Atoms
color Display Background white
light 2 on
display depthcue off
#color Display BackgroundTop white
#color Display BackgroundBot blue2
color Name N blue2
color Name O orange
mol modselect 0 $i "name N" 
mol addrep $i 
mol modselect 1 $i "name O"
#display height 800
display resize 1280 1024
#display shadows on
#pbc join res -border 5
pbc box -center origin -style tubes -width 0.7 -color gray
#D:\\Program Files\\POV-Ray for Windows v3.6\\bin\\
#render POV3 1.pov +H500 +W400 -I%s -O%s.tga +D +X +A +FT
#povray +H500 +W400 -I%s -O%s.tga +D +X +A +FT
}
