# This TCL module contains functions to generate a Kaiser-Bessel-windowed sinc
# FIR filter. The filter can then be converted to fixed point to be loaded
# onto the FPGA, or can be used to design a smooth-and-differentiate
# deconvolution filter for real-time power calculations

proc factorial {n} {
    # Computes n!
    # 0! = 1! = 1
    expr { $n < 2 ? 1 : $n*[factorial [incr n -1]] }
}

proc I_0 {x} {
    # Computes the modified Bessel function of the first kind, order 0, for a
    # given value x. The mathematical calculation uses an ininite sum, which we
    # truncate at the point where additional increments are less than the
    # eventual fixed point precision of the filter.
    set k 0
    set tot 0.0
    set tol [expr { 2.0**-17 }]; # Fixed point coeffs are 17 bits wide
    while true {
        set nextval [ expr { (0.25*$x*$x)**$k / [ factorial $k ]**2 } ]
        if [ expr { abs($nextval) < $tol }] {
            break
        } else {
            set tot [ expr { $tot + $nextval } ]
            incr k
        }
    }
    return $tot
}

proc kaiser_beta {atten} {
    # Calculates the Kaiser beta parameter for a given stop-band attenuation
    # in dB
    expr {
          $atten > 50 ? 0.1102*($atten - 8.7) :
          21 <= $atten <= 50 ? 0.5842*($atten - 21.0)**0.4 + 0.07886*($atten - 21.0) :
          0.0
      }
}

proc kaiser_window {ntaps beta} {
    # Designs a kaiser window with a given number of taps ntaps, and beta
    # parameter which trades off roll-off width with stop band attentuation
    set window_half [list]
    set m [ expr { $ntaps - 1 } ]; # Filter order
    set alpha [expr { $m / 2.0 }]
    set i0_beta [ I_0 $beta ]
    set i0_beta_inv [ expr { 1.0 / $i0_beta } ]
    # The Kaiser window is symmetric, so only calculate half the values
    set nhalf [ expr { $ntaps / 2 } ]
    for { set i 0 } { $i < $nhalf } { incr i } {
        lappend window_half [ expr { [ I_0 [expr {$beta*sqrt(1.0 - (($i - $alpha)/$alpha)**2)} ]] * $i0_beta_inv } ]
    }
    set window [ list {*}$window_half {*}[lreverse $window_half ]]
    return $window
}

proc sinc {ntaps fcs} {
    # Designs a sinc filter (no window) of length ntaps and normalised cutoff
    # frequency fcs = (cutoff frequency/sample frequency)
    set taps [list]
    set pi [ expr { atan(1) * 4 } ]
    set m [ expr { $ntaps - 1 } ]; # Filter order
    for { set i 0 } { $i < $ntaps } { incr i } {
	set imhalf [ expr { ($i-$m/2.0) } ]
        if { $i != [ expr { $ntaps/2 } ] } {
            set twopiioverm [ expr { 2.0*$pi*$i/$m } ]
            lappend taps [ expr { sin(2.0*$pi*$fcs*$imhalf)/$imhalf } ]
        } else {
            lappend taps [ expr { 2.0*$pi*$fcs } ]
        }
    }
    return $taps
}

proc windowed_sinc {ntaps fc fs atten} {
    # Designs a low-pass windowed sinc filter, using a Kaiser-Bessel window
    # The filter has ntaps coefficients, cutoff frequency fc, sample frequency
    # fs, and stop-band attenuation atten. The filter has unity gain
    # Design the sinc part first
    set fcs [ expr double($fc)/$fs ]
    set sinctaps [ sinc $ntaps $fcs ]
    # And the window
    set beta [ kaiser_beta $atten ]
    set window [ kaiser_window $ntaps $beta ]
    # Now create the windowed sinc taps
    set taps [list]
    foreach si $sinctaps wi $window {
        lappend taps [ expr { $si * $wi } ]
    }
    # Normalise the filter to unity gain (i.e. sum of all coefficients = 1)
    set norm_factor [ expr { 1.0 / [ ::tcl::mathop::+ {*}$taps ] } ]
    set taps_norm [ list ]
    foreach t $taps {
	lappend taps_norm [ expr { $t * $norm_factor } ]
    }
    return $taps_norm
}
