#!/usr/local/bin/expect
puts ""
puts "Diagnostics - Displaying RAM Content"
puts ""

set fd [ open /dev/dsp1.3 r ]

binary scan [ read $fd 768 ] i192 offsets

close $fd

set offsets
puts "array: $offsets"
puts ""
