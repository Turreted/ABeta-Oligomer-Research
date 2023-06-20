# Load the protein structure
set structure_file "example-data/relaxed_model_3_multimer_v3_pred_3.pdb"
set color_file "example-data/result_model_3_multimer_v3_pred_3_vmd.csv"

# load in new file and set its molid
mol new $structure_file
set sel [atomselect top "all"]

# Get the number of residues in the protein structure
set selCA [atomselect top "name CA"]
set num_residues [$selCA num]

# Define the array of values for each residue (0-100)
set color_values {}

# read color values from a csv file
set fileHandle [open $color_file r]
while {[gets $fileHandle value] != -1} {
    lappend color_values $value
}

close $fileHandle

puts "Residues: $num_residues"
puts "Color values: $color_values"

set scaled_colors {}

for {set i 0} {$i < $num_residues} {incr i} {
    set residue_value [lindex $color_values $i]

    # scale the color to 0-100 then scale that to R-B
    set scaled_value $residue_value

    # for each atom in the residue, set it to the scaled value.
    # NOTE: lists are 0-indexed and the residues are 1-indexed
    set sel_atoms [atomselect top "resid [expr $i+1]"]
    set num_res [$sel_atoms num]

    for {set j 0} {$j < $num_res} {incr j} {
        lappend scaled_colors $scaled_value
    }

    puts "Residue $i: $scaled_value"
}

# set beta channel of protein to confidence values
$sel set beta $scaled_colors

# Center and zoom the view on the protein
display resetview
