#!/usr/local/bin/expect

set fd [ open /dev/dsp1.3 r ]
binary scan [ read $fd 768 ] i192 offsets
close $fd


set offsets
set offlist [split $offsets]
set offsets1 {}
set offsets2 {}
set offsets3 {}
set offsets4 {}
set offsets_full {}
set channels {}

for {set counter 1} {$counter < 13} {incr counter} {
	lappend channels $counter
	lappend channels $counter
}

for {set i 0} {$i < 24} {incr i} {
	lappend offsets1 [ lindex $offlist $i ]
	lappend offsets2 [ lindex $offlist [ expr $i+48 ] ]
	lappend offsets3 [ lindex $offlist [ expr $i+96 ] ]
	lappend offsets4 [ lindex $offlist [ expr $i+144 ] ]
}

set formatStr {%15s%20s%20s%17s%17s}
puts [format $formatStr "Channel" "Voltage Offsets 1" "Voltage Offsets 2" "Power Offsets 1" "Power Offsets 2"]

puts [format $formatStr "---------------" "-------------------" "-------------------" "---------------" "---------------"]
foreach channelsValue $channels offsets1Value $offsets1 offsets2Value $offsets2 offsets3Value $offsets3 offsets4Value $offsets4 {
    puts [format $formatStr $channelsValue $offsets1Value $offsets2Value $offsets3Value $offsets4Value]
}

