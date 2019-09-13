#!/usr/local/bin/expect

# Number of channels the BOLODSP FPGA module was built for

proc get_nchan {} {
	set fp [open /etc/acq400/0/aggregator r]
	set aggstr [gets $fp]
	close $fp
	regexp {sites=([1-6,]+)} $aggstr nc sites
	set sn [split $sites ,]
	return [expr [llength $sn] * 8 ]
}

set nchan [get_nchan]

# Each channel has I, Q, PI, PQ offsets
set nvals [ expr { $nchan * 4 } ]
# Each offset is 4 bytes
set nbytes [ expr { $nvals * 4 } ]
set fd [ open /dev/dsp1.3 r ]
binary scan [ read $fd $nbytes ] i${nvals} offsets
close $fd


set offsets
set offlist [split $offsets]
set offsets1 {}
set offsets2 {}
set offsets3 {}
set offsets4 {}
set offsets_full {}
set channels {}

for {set i 0} {$i < $nchan} {incr i} {
    lappend channels [expr $i + 1]
    lappend offsets1 [ lindex $offlist [ expr {2 * $i} ] ]
    lappend offsets2 [ lindex $offlist [ expr {2 * $i + 1} ] ]
    lappend offsets3 [ lindex $offlist [ expr {2 * ($i + $nchan)} ] ]
    lappend offsets4 [ lindex $offlist [ expr {2 * ($i + $nchan) + 1} ] ]
}

set formatStr {%15s%20s%20s%17s%17s}
puts [format $formatStr "Channel" "Voltage Offsets I" "Voltage Offsets Q" "Power Offsets I" "Power Offsets Q"]

puts [format $formatStr "---------------" "-------------------" "-------------------" "---------------" "---------------"]
foreach channelsValue $channels offsets1Value $offsets1 offsets2Value $offsets2 offsets3Value $offsets3 offsets4Value $offsets4 {
    puts [format $formatStr $channelsValue $offsets1Value $offsets2Value $offsets3Value $offsets4Value]
}

