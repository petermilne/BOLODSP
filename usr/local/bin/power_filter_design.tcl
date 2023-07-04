#!/usr/bin/expect

# This script takes the cooling time and sensitivity from a calibration,
# and designs a deconvolution filter for a given channel.
# The filter is G = F + tau*dF/dt, where F is a low-pass windowed sinc
# FIR filter
# It reads the requested filter bandwidth from one of the DSP site knobs

# Number of coefficients in the filter
set NTAPS 500
# Bandwidth, i.e. cutoff frequency, read from file
set fil [ open /etc/acq400/14/FILTER_BANDWIDTH ]
set FC [ read $fil ]
close $fil
# Sampling rate gives Nyquist frequency
set FS 250000.0
# Requested stop band attenuation
set ATTEN 60.0
# Binary point for fixed point conversion
set BINARY_POINT 16

# Source the functions used to generate a normalised windowed-sinc low pass
# filter, with a Kaiser-Bessel window
source [file dirname [info script]]/firwin_kb.tcl

# Get the sensitivity and cooling time as command line arguments, as well as
# the channel number
if { $argc < 3 } {
    puts "Usage: $argv0 <channel> <sens> <tau>"
    exit
}
set channel [lindex $argv 0]
set sens [lindex $argv 1]
set tau [lindex $argv 2]

proc differentiate_filter {taps_norm fs} {
    # Calculates the time derivative of the filter taps_norm, for a given
    # sample rate fs, by doing dfdt[i] = diff(taps_norm)/dt. The difference
    # is calculated using second-order central finite difference on most of the
    # list, and second order forward and backward finite differences on the
    # first and last elements respectively
    # To speed up the calculation, use 1/dt = fs
    set dfdt {}
    # First element is second order forward finite difference
    lappend dfdt [ expr { (-1.5 * [lindex $taps_norm 0] \
                              + 2 * [lindex $taps_norm 1] \
                              - 0.5 * [lindex $taps_norm 2]) * $fs } ]
    # Central differences for all elements except first and last
    foreach tpre [ lrange $taps_norm 0 end-2 ] tpost [ lrange $taps_norm 2 end ] {
        lappend dfdt [ expr { 0.5 * ($tpost - $tpre) * $fs } ]
    }
    # Final element is second order backward finite difference
    lappend dfdt [ expr { (0.5 * [lindex $taps_norm end-2] \
                              - 2 * [lindex $taps_norm end-1] \
                              + 1.5 * [lindex $taps_norm end]) * $fs } ]
    return $dfdt
}

proc create_deconvolution_filter {taps_norm dfdt sens tau} {
    # Creates the actual deconvolution filter: G = 1/sens(F + tau*dFdt),
    # where f is given by the taps_norm filter coefficients
    set filt [ list ]
    set invsens [ expr { 1.0/$sens } ]
    foreach f $taps_norm df $dfdt {
	lappend filt [ expr { $invsens * ($f + $tau * $df) } ]
    }
    return $filt
}

proc write_bin {fn bdat} {
	upvar $bdat _bdat
	set fd [ open $fn w ]
	puts -nonewline $fd [ binary format i* $_bdat ]
	close $fd
}

proc output_coef_data_2015 {filt channel binary_point} {
    # Converts filt to fixed point, scaling by 2**binary_point,
    # and formats the filter to be written to the FPGA
    set filt_fixed [ list ]
    foreach f $filt {
	lappend filt_fixed [ expr { int($f * 2**$binary_point) } ]
    }
    # The coefficients need to be reordered for the power FIR compiler
    set filt_fixed_reordered [ list {*}[ lrange $filt_fixed 400 499 ] \
				   {*}[ lrange $filt_fixed 300 399 ] \
				   {*}[ lrange $filt_fixed 200 299 ] \
				   {*}[ lrange $filt_fixed 100 199 ] \
				   {*}[ lrange $filt_fixed 0 99 ]]
    # The filter number (effectively the channel number) needs to be pre-prended
    # to the coefficient set before reload. Channels are indexed from 0 on the FPGA,
    # but from 1 in software
    set filt_out [linsert $filt_fixed_reordered 0 [expr { $channel-1 } ]]
    # Write to the power coefficients device file
    #write_bin /dev/dsp1.2 filt_out
    write_bin "/tmp/power_filter$channel" filt_out
}

proc output_coef_data {filt channel binary_point} {
    # Converts filt to fixed point, scaling by 2**binary_point,
    # and formats the filter to be written to the FPGA
    set filt_fixed [ list ]
    foreach f $filt {
	lappend filt_fixed [ expr { int($f * 2**$binary_point) } ]
    }
    # The coefficients need to be reordered for the power FIR compiler
    set filt_fixed_reordered [ list {*}[ lrange $filt_fixed 375 499 ] \
				    {*}[ lrange $filt_fixed 250 374 ] \
				    {*}[ lrange $filt_fixed 125 249 ] \				   
				    {*}[ lrange $filt_fixed   0 124 ]]
    # The filter number (effectively the channel number) needs to be pre-prended
    # to the coefficient set before reload. Channels are indexed from 0 on the FPGA,
    # but from 1 in software
    set filt_out [linsert $filt_fixed_reordered 0 [expr { $channel-1 } ]]
    # Write to the power coefficients device file
    #write_bin /dev/dsp1.2 filt_out
    write_bin "/tmp/power_filter$channel" filt_out
}

proc produce_filter_data {NTAPS FC FS ATTEN BINARY_POINT sens tau channel} {
    if { $sens == 0 } {
        # A failed calibration precludes designing a useful filter.
        # Explicitly zero the output to make clear it shouldn't be used.
        set filt [ lrepeat $NTAPS 0 ]
    } else {
        # Produce a windowed-sinc low pass filter
        set taps_norm [ windowed_sinc $NTAPS $FC $FS $ATTEN]
        # Differentiate the filter
        set dfdt [ differentiate_filter $taps_norm $FS ]
        # Now create the deconvolution filter
        set filt [ create_deconvolution_filter $taps_norm $dfdt $sens $tau ]
    }
    # Finally, output the data
    output_coef_data $filt $channel $BINARY_POINT
}

produce_filter_data $NTAPS $FC $FS $ATTEN $BINARY_POINT $sens $tau $channel
