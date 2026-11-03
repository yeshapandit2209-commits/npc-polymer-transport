# build_data.tcl
# Requires VMD with the TopoTools plugin.
#
# Usage:
#   vmd -dispdev text -e build_data.tcl -args <n> <blocklength> <pdbfile> <outprefix>
#
# Reads the SAW backbone (chain.pdb, same for every blocklength),
# assigns charge type in alternating blocks of size <blocklength>
# (type 1 = positive, +1e; type 2 = negative, -1e), adds linear bonds,
# and writes a LAMMPS full-style data file + PSF for this blocklength.

package require topotools

set n        [lindex $argv 0]
set blocklen [lindex $argv 1]
set pdbfile  [lindex $argv 2]
set outpref  [lindex $argv 3]

if {[expr {$n % (2*$blocklen)}] != 0} {
    puts "WARNING: n ($n) is not an even multiple of 2*blocklength ($blocklen)."
    puts "         Net charge may not come out exactly zero -- check the printed sum below."
}

mol new $pdbfile waitfor all

# explicitly set the periodic box -- must match generate_saw.py's --box
# (without this, writelammpsdata falls back to a bogus/near-zero box,
# which is what caused the "bond extent > half box" crash)
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

# assign type/charge per bead in blocks of $blocklen, alternating + / -
set types {}
set charges {}
for {set i 0} {$i < $n} {incr i} {
    set block_idx [expr {int($i / $blocklen)}]
    if {$block_idx % 2 == 0} {
        lappend types 1
        lappend charges 1.0
    } else {
        lappend types 2
        lappend charges -1.0
    }
}

$sel set type    $types
$sel set charge  $charges
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

# sanity check: net charge should be 0
set qsum 0.0
foreach q $charges { set qsum [expr {$qsum + $q}] }
puts "blocklength=$blocklen : net charge = $qsum"

topo writelammpsdata "${outpref}_bl${blocklen}.data" full
animate write psf "${outpref}_bl${blocklen}.psf"

puts "Wrote ${outpref}_bl${blocklen}.data and ${outpref}_bl${blocklen}.psf"
exit
