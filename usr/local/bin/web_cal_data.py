#!/usr/bin/python

"""
Formats the calibration data for display in the BOLO_CAL web page.
"""

import os

SENS_MIN = float(os.getenv("SENS_MIN", 1))
SENS_MAX = float(os.getenv("SENS_MAX", 15))

format_str = "{0:<15s}{1:<20s}{2:<20s}{3:<17s}{4:<17s}{5:<17s}"

print(format_str.format("Channel", "I0", "Q0", "Sensitivity", "Cooling time", "Status"))
print("-------------- ------------------- ------------------- ---------------- "
      "---------------- ----------------")

with open("/tmp/calibfit.log", "r") as f:
    for line in f:
        ch, sens, tau, ioff, qoff = line.split()
        status = "PASS" if SENS_MIN <= float(sens) <= SENS_MAX else "FAIL"
        print(format_str.format(ch, ioff, qoff, sens, tau, status))
