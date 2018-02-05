#!/usr/local/bin/expect
puts "Checking that 4 arguments are passed to this file"
if { $argc < 4 } {
    puts "Usage: $argv0 <channel> <i0> <q0> <sens>"
    puts "Loads the offsets for voltage and power for a given channel"
    exit
}

puts "Confirmed: 4 args have been passed"

# Sadly, there is currently no seeking support for the offset
# coefficients file. So the only way to load one channel's
# offset at a time is to read the whole offsets RAM, change
# the channel in question and then write the whole RAM back
# to the FPGA with the previous contents of the RAM plus the
# modified channel

puts "opening file dsp1.3 in read mode"
set fd [ open /dev/dsp1.3 r ]
puts "opened file in read mode"

# 48 channels, I, Q, PI, PQ per channel, 4 bytes per value
binary scan [ read $fd 768 ] i192 offsets


close $fd
puts "fd closed"

# Channels numbered from 1, but list indexed from 0

set channel [ expr { [ lindex $argv 0 ] - 1 } ]
puts "channel value = $channel"

set i0 [ lindex $argv 1 ]
puts "i0 value = $i0"

set q0 [ lindex $argv 2 ]
puts "q0 value = $q0"

set sens [ lindex $argv 3 ]
puts "sens value = $sens"

set scale [ expr { 1.25/2**25*20/18 } ]
puts "scale value = $scale"

set i0_int [ expr { int($i0/$scale) } ]
puts "i0_int value = $i0_int"

set q0_int [ expr { int($q0/$scale) } ]
puts "q0_int value = $q0_int"

set pscale [ expr { 1.25/2**19*20/18 } ]
puts "pscale value = $pscale"

set pi0_int [ expr { int(($i0/$sens)/$pscale) } ]
puts "pi0_int value = $pi0_int"

set pq0_int [ expr { int(($q0/$sens)/$pscale) } ]
puts "pq0_int value = $pq0_int" 

# Replace the I0, Q0, PI0 and PQ0 elements of the offsets
# list with updated values



set i0_index [ expr { $channel*2 } ]
puts "i0_index = $i0_index"

set q0_index [ expr $i0_index+1 ]
puts "q0_index = $i0_index"

# The power offsets come after the voltage offsets, so
# after 2*nchannels values

set nchannels 48
set pi0_index [ expr { 2*($nchannels+$channel) } ]
set pq0_index [ expr { $pi0_index+1 } ]

set offsets [ lreplace $offsets $i0_index $q0_index $i0_int $q0_int ]
puts "replaced index $i0_index with $i0_int and replacing index $q0_index with $q0_int , both from array: $offsets"
puts ""
set array_size [llength $offsets]
puts "array size = $array_size"

set offsets [ lreplace $offsets $pi0_index $pq0_index $pi0_int $pq0_int ]
puts "replaced index $pi0_index with $pi0_int and replacing index $pq0_index with $pq0_int , both from array: $offsets"


set fd [ open /dev/dsp1.3 w ]

puts -nonewline $fd [ binary format i* $offsets ]


flush $fd

close $fd
