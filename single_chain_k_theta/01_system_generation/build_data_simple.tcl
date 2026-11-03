# build_data_simple.tcl
# Requires VMD with the TopoTools plugin.
#
# Usage:
#   vmd -dispdev text -e build_data_simple.tcl -args <n> <pdbfile> <outprefix>
#
# Same SAW backbone as the blocklength study, but here every bead is
# identical: one atom type, charge = 0. Only the angle stiffness k
# is varied later, per LAMMPS run -- topology here is the same for
# every k, so this is written ONCE and reused for all of them.

package require topotools

set n       [lindex $argv 0]
set pdbfile [lindex $argv 1]
set outpref [lindex $argv 2]

mol new $pdbfile waitfor all

# explicitly set the periodic box -- must match generate_saw.py's --box
set boxlen 100.0
molinfo top set a $boxlen
molinfo top set b $boxlen
molinfo top set c $boxlen
molinfo top set alpha 90.0
molinfo top set beta  90.0
molinfo top set gamma 90.0

set sel [atomselect top "all"]
set natoms [$sel num]
if {$natoms != $n} {
    puts "ERROR: pdb has $natoms atoms, expected $n. Aborting."
    exit 1
}

# every bead: same type, no charge
$sel set type    1
$sel set charge  0.0
$sel set resname "BEA"
$sel set name    "CA"
$sel set mass    1.0

# bonds: linear chain, bead i -- bead i+1 (explicit type 1)
set bondlist {}
for {set i 0} {$i < [expr {$n - 1}]} {incr i} {
    lappend bondlist [list $i [expr {$i + 1}] 1]
}
topo setbondlist $bondlist

# angles: auto-derive from bond topology (needed for angle_style harmonic)
topo guessangles

topo writelammpsdata "${outpref}.data" full
animate write psf "${outpref}.psf"

puts "Wrote ${outpref}.data and ${outpref}.psf"
exit
