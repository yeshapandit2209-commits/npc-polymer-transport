package require topotools
package require pbctools

# --------------------------------------------------
# PBC
# --------------------------------------------------

pbc set {150 150 150}

# --------------------------------------------------
# Parameters
# --------------------------------------------------

set numChains 88
set chainLength 40

# No rigidLength needed here: build_cyl_system.py already bakes the
# RRRFFF rigid/flex pattern into each atom's name (C_rigid/C_flex) when
# it writes the coordinate file, so the "Atom types" section below just
# reads that back off atom name -- it doesn't need to know the pattern.

# --------------------------------------------------
# Clear topology
# --------------------------------------------------

topo clearbonds

set sel [atomselect top all]
$sel set mass 1.0
$sel set charge 0.0
$sel delete

# --------------------------------------------------
# Polymer atoms
# --------------------------------------------------

set polymer [atomselect top "name C_anchor C_rigid C_flex A"]

set polymerAtoms [$polymer get index]
set numPolymerAtoms [$polymer num]

puts "Found $numPolymerAtoms polymer atoms"

$polymer delete

if {$numPolymerAtoms != [expr {$numChains*$chainLength}]} {

    puts "ERROR:"
    puts "Expected [expr {$numChains*$chainLength}] polymer atoms"
    puts "Found $numPolymerAtoms"

    return
}

# --------------------------------------------------
# Chain labels
# --------------------------------------------------

set letters {
A B C D E F G H I J K L M N O P Q R S T
U V W X Y Z
AA AB AC AD AE AF AG AH AI AJ AK AL AM AN AO AP AQ AR AS AT
AU AV AW AX AY AZ
BA BB BC BD BE BF BG BH BI BJ BK BL BM BN BO BP BQ BR BS BT
BU BV BW BX BY BZ
CA CB CC CD CE CF CG CH CI CJ CK CL CM CN CO CP CQ CR CS CT
CU CV CW CX CY CZ
DA DB DC DD DE DF DG DH DI DJ DK DL DM DN DO DP DQ DR DS DT
}

# --------------------------------------------------
# Assign chain names
# --------------------------------------------------

for {set chainID 0} {$chainID < $numChains} {incr chainID} {

    set offset [expr {$chainID*$chainLength}]
    set chainName [lindex $letters $chainID]

    for {set i 0} {$i < $chainLength} {incr i} {

        set atomIndex [lindex $polymerAtoms [expr {$offset+$i}]]

        set sel [atomselect top "index $atomIndex"]

        $sel set resname $chainName

        $sel delete
    }
}

puts "Chain names assigned"

# --------------------------------------------------
# Atom types
# --------------------------------------------------

set sel [atomselect top "name C_anchor"]
$sel set type 1
$sel delete

set sel [atomselect top "name C_rigid"]
$sel set type 2
$sel delete

set sel [atomselect top "name C_flex"]
$sel set type 3
$sel delete

set sel [atomselect top "name X"]
$sel set type 4
$sel delete

set sel [atomselect top "name A"]
$sel set type 5
$sel delete

puts "Atom types assigned"

# --------------------------------------------------
# Polymer selection for bonds
# --------------------------------------------------

set polymer [atomselect top "name C_anchor C_rigid C_flex A"]

set chain_atoms [$polymer get index]

$polymer delete

# --------------------------------------------------
# Backbone bonds
# --------------------------------------------------

puts "Creating backbone bonds..."

for {set chainID 0} {$chainID < $numChains} {incr chainID} {

    set offset [expr {$chainID*$chainLength}]

    for {set i 0} {$i < [expr {$chainLength-1}]} {incr i} {

        set atom1 [lindex $chain_atoms [expr {$offset+$i}]]
        set atom2 [lindex $chain_atoms [expr {$offset+$i+1}]]

        topo addbond $atom1 $atom2
    }
}

puts "Backbone bonds created"

# --------------------------------------------------
# Generate angles
# --------------------------------------------------

puts "Generating angles..."

topo guessangles

puts "Angles generated"

# --------------------------------------------------
# Diagnostics
# --------------------------------------------------

puts "Anchors : [[atomselect top {name C_anchor}] num]"
puts "Rigid   : [[atomselect top {name C_rigid}] num]"
puts "Flex    : [[atomselect top {name C_flex}] num]"
puts "A       : [[atomselect top {name A}] num]"
puts "Wall    : [[atomselect top {name X}] num]"

# --------------------------------------------------
# Output
# --------------------------------------------------

puts "Writing PDB..."

set all_atoms [atomselect top all]

$all_atoms writepdb 100_s80.pdb

puts "Writing PSF..."

animate write psf 100_s80.psf

puts "Writing LAMMPS data..."

topo writelammpsdata 100_s80.data full

puts "Done."

$all_atoms delete