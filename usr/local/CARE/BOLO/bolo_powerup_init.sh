#!/bin/sh
# Set up the DSP
set.fdrive 18e3

CUSTOM=/mnt/local/sysconfig/bolo.sh
if [ -e $CUSTOM ]; then
	source $CUSTOM
	case $CALIBFIT_EXCITEV in
	18.0)
		DAC_EXCITE_AMP=3;;
	15.0)
		DAC_EXITE_AMP=2;;
	12.0)
		DAC_EXITE_AMP=1;;
	9.0)
		DAC_EXCITE_AMP=0;;
	*)
		echo "WARNING: CALIBFIT_EXCITEV not specified, setting 9.0"
		DAC_EXCITE_AMP=0;;
	esac

	set.site 14 DAC_EXCITE_AMP=$DAC_EXCITE_AMP
fi
# Turn all the gains up to 1V25, since bolometer signals are always small
hostname=$(hostname)
echo "Setting gains to 1V25"
# Read which sites are active from the distributor - we set the distributor
# to control all BOLO sites in bolodsp.init
sites=$(set.site 0 distributor | awk '{print $2}' | grep -o '[1-6]')
for site in $sites; do
    ao420_init $site
    for channel in $(seq 1 8); do
        set.site $site B8:GAIN:$channel 1V2  
    done
done

# Set internal rising edge trigger
set.site 1 trg=1,1,1
# 1MHz sample clock set with role - ROLE=fpmaster for front panel clk
set.site 0 sync_role ${ROLE:-master} 1M ${FIN:-} TRG:DX=d1 CLKDIV=1


# Reset offset DACs. They seem to be confused by the switch to sideport control,
# or possibly initialising EPICS
echo 1 | tee /dev/bolo8/*/OS_DAC_RESETn > /dev/null
echo 0 | tee /dev/bolo8/*/OS_DAC_RESETn > /dev/null

set.sys /sys/module/acq420fmc/parameters/FIFERR 0

# Start a stream to set the DAC output going, to Ohmically heat the sensor.
# It will continue even after the stream is stopped
/usr/local/bin/streamtonowhered start
sleep 1
soft_trigger
sleep 1
/usr/local/bin/streamtonowhered stop
#/usr/local/bin/web_diagnostics_ram

