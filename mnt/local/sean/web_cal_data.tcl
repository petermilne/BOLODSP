#!/usr/local/bin/expect

set file [ open "/tmp/calibfit.log" "r" ]

set data [ read $file ]
set formatted_data [ split $data ]

set ch {}
set i0 {}
set q0 {}
set sens {}

set formatted_data [lsearch -all -inline -not -exact $formatted_data "load_offset_channel.tcl"]
#set formatted_data [lreplace $formatted_data $idx $idx]


for {set i 0} {$i <= [  llength $formatted_data ]} {incr i 4} {
        lappend ch [ lindex $formatted_data $i ]
        lappend i0 [ lindex $formatted_data [ expr $i+1 ] ]
        lappend q0 [ lindex $formatted_data [ expr $i+2 ] ]
        lappend sens [ lindex $formatted_data [ expr $i+3 ] ]

}

set formatStr {%15s%20s%20s%17s}
puts [format $formatStr "Channel" "i0" "q0" "Sensitivity" ]

puts [format $formatStr "---------------" "-------------------" "-------------------" "---------------"]
foreach ch $ch i0 $i0 q0 $q0 sens $sens { 
    puts [format $formatStr $ch $i0 $q0 $sens ]
}


