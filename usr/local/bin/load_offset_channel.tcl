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

# 48 channels, I, Q, PI, PQ per channel, 4 bytes per value
# Why 48? MAX 6 x BOLO8 = 48 channels in one box. One size fits all..

set MAXCHAN 48

proc read_bin {fn bdat} {
	upvar $bdat _bdat
	set fd [ open $fn r ]
	binary scan [ read $fd 768 ] i192 _bdat
	close $fd
}

proc write_bin {fn bdat} {
	upvar $bdat _bdat
	set fd [ open $fn w ]
	puts -nonewline $fd [ binary format i* $_bdat ]
	close $fd
}

proc read_knob {knob} {
	set fd [ open $knob r ]
	set kv [ string trim [ read $fd ] ]
	close $fd
	return $kv
}

proc get_vgain {site physchan} {
	# Calculate the maximum voltage by reading ADC gain registers.
	# This is significantly faster than using EPICS PV access, which is
	# necessary if this script is to run for 48 channels in a reasonable
	# amount of time.
	set a0 [read_knob /dev/bolo8/${site}/ADC_${physchan}_A0 ]
	set a1 [read_knob /dev/bolo8/${site}/ADC_${physchan}_A1 ]
	set vg [ expr { 10.0 / 2**"0b${a1}${a0}" } ]
	return $vg
}

proc get_excite_scale {} {
	set dac_excite_amp_knob [ read_knob /etc/acq400/14/DAC_EXCITE_AMP ]
	set excite_lut [ dict create 3 18.0 2 15.0 1 9.0 0 1.0 ];  # From wavegen.vhd
	set excite_amp [ dict get $excite_lut $dac_excite_amp_knob ]
	set scale [ expr { 20.0 / $excite_amp } ]
	return $scale
}

set ch [ lindex $argv 0 ]

read_bin /dev/dsp1.3 offsets
# for post mortem examine before
write_bin /tmp/b8_offsets${ch}.1 offsets

# Channels numbered from 1, but list indexed from 0
set ch0 [ expr $ch - 1 ]
set i0 [ lindex $argv 1 ]
set q0 [ lindex $argv 2 ]
set sens [ lindex $argv 3 ]
# Calculate the volts-to-counts scaling, depends on the gain
set site [ expr {int($ch0 / 8) + 1} ]
set physchan [ expr {int($ch0 % 8) + 1} ]

set vgain [ get_vgain $site $physchan ]
set excscale [ get_excite_scale ]

# Note, I and Q scales use 2**25 rather than 2**24 like A because
# the offset correction works on I=A/2 cos(phi) and Q=-A/2 sin(phi).
set scale [ expr { $vgain*$excscale/2**25 } ]
set i0_int [ expr { int($i0/$scale) } ]
set q0_int [ expr { int($q0/$scale) } ]

# Ditto with power: use 2**19 instead of 2**18.
set pscale [ expr { $vgain*$excscale/2**19 } ]
set pi0_int [ expr { $sens == 0 ? 0 : int(($i0/$sens)/$pscale) } ]
set pq0_int [ expr { $sens == 0 ? 0 : int(($q0/$sens)/$pscale) } ]
# Replace the I0, Q0, PI0 and PQ0 elements of the offsets
# list with updated values
set i0_index [ expr { $ch0*2 } ]
set q0_index [ expr $i0_index+1 ]
# The power offsets come after the voltage offsets, so
# after 2*MAXCHAN values

set pi0_index [ expr { 2*($MAXCHAN+$ch0) } ]
set pq0_index [ expr { $pi0_index+1 } ]
set offsets [ lreplace $offsets $i0_index $q0_index $i0_int $q0_int ]
set offsets [ lreplace $offsets $pi0_index $pq0_index $pi0_int $pq0_int ]

write_bin /dev/dsp1.3 offsets
# for post mortem examine before
write_bin /tmp/b8_offsets${ch}.2 offsets


