#!/bin/sh
# Wait until the filters are ready for another reload sequence, then
# reset the DSP system by writing to the reset register and then
# clearing it
wait_for_filters
# Now we can safely reset
set.site 14 DSP_RESET=1
sleep 0.1
set.site 14 DSP_RESET=0
