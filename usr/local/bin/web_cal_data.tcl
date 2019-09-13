#!/usr/local/bin/expect

set file [ open "/tmp/calibfit.log" "r" ]

set data [ read $file ]
set formatted_data [ split $data ]

set ch {}
set i0 {}
set q0 {}
set sens {}
set tau {}
set status {}

set formatted_data [lsearch -all -inline -not -exact $formatted_data "calibration_data"]
#set formatted_data [lreplace $formatted_data $idx $idx]


for {set i 0} {$i <= [  llength $formatted_data ]} {incr i 5} {
	catch {
        	lappend ch [ lindex $formatted_data $i ]
	        lappend i0   [format {%.3f} [ lindex $formatted_data [ expr $i+1 ] ] ]
        	lappend q0   [format {%.3f} [ lindex $formatted_data [ expr $i+2 ] ] ]
		set _sens    [format {%.3f} [ lindex $formatted_data [ expr $i+3 ] ] ]
	        lappend sens $_sens
	        lappend tau  [format {%.3f} [ lindex $formatted_data [ expr $i+4 ] ] ]
		if { $_sens >= 1 && $_sens <= 10 } {
			lappend status "PASS"
		} else {
			lappend status "FAIL"
		}
	}
}

set formatStr {%15s%20s%20s%17s%17s%17s}
puts [format $formatStr "Channel" "i0" "q0" "Sensitivity" "Cooling time" "Status" ]

puts [format $formatStr "---------------" "-------------------" "-------------------" "---------------" "---------------" "-----------------"]

foreach ch $ch i0 $i0 q0 $q0 sens $sens tau $tau status $status {
    puts [format $formatStr $ch $i0 $q0 $sens $tau $status]
}


