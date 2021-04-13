#!/usr/local/bin/expect
if { $argc < 4 } {
    puts "Usage: $argv0 <channel> <i0> <q0> <sens>"
    puts "Loads the offsets for voltage and power for a given channel"
    exit
}
# Sadly, there is currently no seeking support for the offset
# coefficients file. So the only way to load one channel's
# offset at a time is to read the whole offsets RAM, change
# the channel in question and then write the whole RAM back
# to the FPGA with the previous contents of the RAM plus the
# modified channel
set fd [ open /dev/dsp1.3 r ]
# 48 channels, I, Q, PI, PQ per channel, 4 bytes per value
binary scan [ read $fd 768 ] i192 offsets
close $fd
# Channels numbered from 1, but list indexed from 0
set channel [ expr { [ lindex $argv 0 ] - 1 } ]
set i0 [ lindex $argv 1 ]
set q0 [ lindex $argv 2 ]
set sens [ lindex $argv 3 ]
set scale [ expr { 1.25/2**25*20/18 } ]
set i0_int [ expr { int($i0/$scale) } ]
set q0_int [ expr { int($q0/$scale) } ]
set pscale [ expr { 1.25/2**19*20/18 } ]
set pi0_int [ expr { $sens == 0 ? 0 : int(($i0/$sens)/$pscale) } ]
set pq0_int [ expr { $sens == 0 ? 0 : int(($q0/$sens)/$pscale) } ]
# Replace the I0, Q0, PI0 and PQ0 elements of the offsets
# list with updated values
set i0_index [ expr { $channel*2 } ]
set q0_index [ expr $i0_index+1 ]
# The power offsets come after the voltage offsets, so
# after 2*nchannels values
set nchannels 48
set pi0_index [ expr { 2*($nchannels+$channel) } ]
set pq0_index [ expr { $pi0_index+1 } ]
set offsets [ lreplace $offsets $i0_index $q0_index $i0_int $q0_int ]
set offsets [ lreplace $offsets $pi0_index $pq0_index $pi0_int $pq0_int ]
set fd [ open /dev/dsp1.3 w ]
puts -nonewline $fd [ binary format i* $offsets ]
flush $fd
close $fd
