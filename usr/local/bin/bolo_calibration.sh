#!/bin/sh

# This script runs a calibration: soft triggers the device, waits 0.1s,
# then sets the bias voltage, waits for heating and then sets the bias voltage
# back to reverse. It assumes that a transient capture has already been set up
# and armed to recover the data

# Read calibration parameters from knobs stored here:
KNOBDIR=/etc/acq400/14/
cal_en=$(cat ${KNOBDIR}/CAL_EN)
tcool=$(cat ${KNOBDIR}/TCOOL)
theat=$(cat ${KNOBDIR}/THEAT)
vbias=$(cat ${KNOBDIR}/VBIAS)
cal_delay=$(cat ${KNOBDIR}/CAL_DELAY)
diode_drop_v=$(cat ${KNOBDIR}/DIODE_DROP_V)

# Exit if we don't want to perform a calibration
[ $cal_en -eq 0 ] && exit

# Reset the DSP system to initialise counters etc.
echo 1 > /dev/dsp1/DSP_RESET
sleep 0.001
echo 0 > /dev/dsp1/DSP_RESET
sleep 1
# Tell the FPGA this is a calibration run
echo 1 > /dev/dsp1/CALIBRATION
# Reverse bias the offset DAC to begin with
osdac_reverse=0x03a000
echo $osdac_reverse > /dev/dsp1/OSDAC_REG_DATA
# Empiracally measured drop in voltage across diodes
diode_drop=$(expect -c "puts [ expr {${diode_drop_v}/15.0*2**15}]")
echo $diode_drop > /dev/dsp1/DIODE_DROP
# Calculate the required offset DAC setting based on the bias voltage: 
# Voltage is differential, on a 16-bit signed DAC with 7.5V reference
osdac_value=$(expect -c "puts [expr {int(${vbias}/15.0*2**15+$diode_drop)}]")
osdac_hex=$(printf %04x $osdac_value)
# Start data stream
soft_trigger 1
# Wait for 0.1s to allow an offset measurement to be made
sleep 0.1
# Set the bias voltage
# N.B. The top byte of the offset-DAC register should be 0x03, as this will
# be translated into the broadcast address of the offset-DAC on the FPGA
echo 0x03${osdac_hex} > /dev/dsp1/OSDAC_REG_DATA
# Wait for heating to equilibrium
sleep $theat
# Turn off the heating and reverse bias the offset DAC
echo $osdac_reverse > /dev/dsp1/OSDAC_REG_DATA
# Wait for the cooling to finish
sleep $tcool
# Unset the trigger, ready for next run
soft_trigger 0
# Unset the calibration run setting too. It'll be re-set to 1 the next
# time this script runs anyway
# Wait until acquisition has finished first, else VDC/current will be replaced
# with power in the data stream
while true
do
    acquire_idle=$(set.site 0 transient_state | awk '{print $1}')
    [ "$acquire_idle" -eq 0 ] && break
    sleep 1
done
echo 0 > /dev/dsp1/CALIBRATION
