#!/usr/local/bin/tclsh8.6

# This script designs a Kaiser-Bessel windowed sinc FIR filter kernel
# The output is a series of 4-byte signed integers to be redirected to
# the DSP module's filter reload address space.
# It reads the required filter bandwidth from a DSP site knob

# Number of coefficients in the filter
set NTAPS 400
# Bandwidth, i.e. cutoff frequency, read from file
set fil [ open /etc/acq400/14/FILTER_BANDWIDTH ]
set FC [ read $fil ]
close $fil
# Sampling rate gives Nyquist frequency
set FS 250000.0
# Filter stop-band attenuation
set ATTEN 60.0
# Binary point for fixed point conversion
set BINARY_POINT 22

# Source the script containing the windowed-sinc-generating functions. It
# must be in the same directory as this script
source [file dirname [info script]]/firwin_kb.tcl

proc fixed_point_windowed_sinc {ntaps fc fs atten binary_point} {
    # Designs a low-pass Kaiser-Bessel windowed sinc FIR filter, which is
    # normalised and scaled by 2**binary_point, and returned as integers
    set taps_norm [ windowed_sinc $ntaps $fc $fs $atten ]
    set taps_int [ list ]
    foreach t $taps_norm {
        lappend taps_int [ expr { int($t * 2**$binary_point) } ]
    }
    return $taps_int
}

proc output_filter_data {taps_int} {
    # Reorder the data in to the Vivado-specified reload order
    set taps_reordered [ list {*}[lrange $taps_int 100 199] \
                             {*}[lrange $taps_int 0 99]]
    # Write the data as binary 4-byte signed integers to the quadrature
    # coefficients device file
    set fd [ open /dev/dsp1.1 w ]
    puts -nonewline $fd [ binary format i* $taps_reordered ]
    close $fd
}

proc design_and_load_filter {NTAPS FC FS ATTEN BINARY_POINT} {
    set taps_int [fixed_point_windowed_sinc $NTAPS $FC $FS $ATTEN $BINARY_POINT]
    output_filter_data $taps_int
}

design_and_load_filter $NTAPS $FC $FS $ATTEN $BINARY_POINT
