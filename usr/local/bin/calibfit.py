#!/usr/bin/python

"""
Calculate calibration constants from acquired calibration data.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import reduce
from itertools import repeat
import logging
from collections import namedtuple
import math
from operator import mul
import os
from pathlib import Path
import numpy as np


# Logging configuration.
logger = logging.getLogger("calibfit")
logger.setLevel(os.environ.get("LOGLEVEL", "INFO"))
logger.addHandler(logging.StreamHandler())
include_traceback = (logger.level == logging.DEBUG)


# Raw data scale factors.
# Amplitude scale depends on EXCITEV and CORDIC scale factor Z_n.
# n is nbits - 2 according to Xilinx, so we need Z_30.
Z_30 = reduce(mul, (math.sqrt(1 + 2**(-2 * i)) for i in range(1, 31)))
EXCITEV = float(os.getenv("CALIBFIT_EXCITEV", "18"))
GAIN_PV = {"1V2": 1.25, "2V5": 2.5, "5V0": 5.0, "10V": 10.0}
AMP_SCALE = 1 / 2**24 * 20 / EXCITEV / Z_30  # Multiplied by vgain
PHASE_SCALE = 2**-29
VDC_SCALE = 1 / 2**15  # Multiplied by vgain
IDC_SCALE = 128 / 100 / 2**12 * 25 / 3 / 1000
SAMPLE_RATE = 1e4


# Number of samples to skip due to FPGA filter warmup.
SKIP = 20


# Some more descriptive containers for data returned from functions.
ChannelData = namedtuple("ChannelData", ["amplitude", "phase", "vdc", "idc", "time"])
LinearFit = namedtuple("LinearFit", ["c0", "tau"])
FitParams = namedtuple("FitParams", ["sens", "tau", "Ioff", "Qoff"])


def read_channel(channel: int, nsamples: int, gainpv: str, fileroot: Path) -> ChannelData:
    """
    Read the amplitude, phase and DC bias data for a given channel.

    :param channel: the channel number, from 1 to NCHAN
    :param nsamples: the number of samples to read
    :param gainpv: the B8:GAIN PV setting for the channel to be read
    :param fileroot: the root directory where the data is stored
    :return: a tuple (amplitude, phase, vdc, idc, time) of arrays
    """
    site = (channel - 1) // 8 + 1
    sitechannel = (channel - 1) % 8 + 1
    ampchannel = (sitechannel - 1) * 3 + 1
    phasechannel = (sitechannel - 1) * 3 + 2
    dcchannel = (sitechannel - 1) * 3 + 3
    datadir = fileroot / f"{site}"
    ampfile = datadir / f"{ampchannel:02d}"
    phasefile = datadir / f"{phasechannel:02d}"
    dcfile = datadir / f"{dcchannel:02d}"
    # np.fromfile only gained support for Path objects in 1.17: use str for 1.16.
    amplitude = np.fromfile(str(ampfile), np.int32, nsamples)
    phase = np.fromfile(str(phasefile), np.int32, nsamples)
    dcdata = np.fromfile(str(dcfile), np.int32, nsamples)
    vdc = (dcdata // 2**16)
    idc = (dcdata.view(np.uint32) % 2**16)
    time = np.arange(nsamples)
    # Trim off samples from incomplete filter output.
    amplitude = amplitude[SKIP:]
    phase = phase[SKIP:]
    vdc = vdc[SKIP:]
    idc = idc[SKIP:]
    time = time[SKIP:]
    # Convert to physical units. Use float32 to save memory and speed up arithmetic.
    vgain = GAIN_PV[gainpv]
    amplitude = (amplitude * vgain * AMP_SCALE).astype(np.float32)
    phase = (phase * PHASE_SCALE).astype(np.float32)
    vdc = (vdc * vgain * VDC_SCALE).astype(np.float32)
    idc = (idc * IDC_SCALE) .astype(np.float32)
    time = (time / SAMPLE_RATE).astype(np.float32)
    return ChannelData(amplitude, phase, vdc, idc, time)


def get_preheating_indices(vdc: np.ndarray,
                           heating_threshold: float,
                           cooling_threshold: float) -> np.ndarray:
    """
    Return the array indices corresponding to the period before heating.

    Heating is when vdc > heating threshold. The pre-heating phase is
    useful for measuring the I and Q voltage offsets, which will be the
    same as the cooling phase.

    :param vdc: the DC voltage data for a single channel
    :param heating_threshold: DC voltage above which to consider heating
    :param cooling_threshold: DC voltage below which to consider no heating
    :return: the indices in vdc for the pre-heating period
    """
    heating_indices = get_heating_indices(vdc, heating_threshold)
    nonheating_indices = np.flatnonzero(vdc < cooling_threshold)
    heating_start = heating_indices.min()
    preheating_indices = nonheating_indices[nonheating_indices < heating_start]
    return preheating_indices


def get_heating_indices(vdc: np.ndarray, heating_threshold: float) -> np.ndarray:
    """
    Return the array indices corresponding to heating.

    Heating is when vdc > heating_threshold.

    :param vdc: the DC voltage data for a single channel
    :param heating_threshold: DC voltage above which to consider heating
    :return: the indices in vdc for the heating period
    """
    indices = np.flatnonzero(vdc > heating_threshold)
    if indices.size == 0:
        raise RuntimeError("No heating measured: reduce heating threshold.")
    return indices


def get_cooling_indices(vdc: np.ndarray,
                        cooling_threshold: float,
                        heating_indices: np.ndarray) -> np.ndarray:
    """
    Return the array indices corresponding to cooling.

    Cooling is when vdc < cooling threshold, and after heating.

    :param vdc: the DC voltage data for a single channel
    :param cooling_threshld: DC voltage below which to consider cooling
    :param heating_indices: the indices in vdc corresponding to heating
    :return: the indices in vdc for the cooling period
    """
    non_heating_indices = np.flatnonzero(vdc < cooling_threshold)
    last_heating = heating_indices.max()
    cooling_indices = non_heating_indices[non_heating_indices > last_heating]
    if cooling_indices.size == 0:
        raise RuntimeError("No cooling measured. "
                           "Check foil is not broken, or reduce cooling threshold.")
    return cooling_indices


def calculate_offset(voltage: np.ndarray) -> float:
    """
    Calculate the offset (voltage at zero power).

    :param voltage: voltage data (I or Q) from the twait part of the data
    :return: the voltage offset
    """
    return voltage.mean().item()


def calculate_heating_power(vheat: np.ndarray, iheat: np.ndarray) -> float:
    """
    Calculate the applied Ohmic heating power.

    :param vheat: the DC voltage during heating.
    :param iheat: the DC current during heating.
    :return: the heating power, averaged over the heating period.
    """
    # Two resistors are heated in the circuit.
    pheat = 2 * vheat * iheat
    return pheat.mean().item()


def fit_cooling(vcool: np.ndarray, tvec: np.ndarray) -> LinearFit:
    """
    Fit the cooling curve to the function V = C0 exp(-t/tau).

    :param vcool: the voltage trace (I or Q) during cooling,
                  with offset subtracted
    :param tvec: the time vector
    :return: a tuple containing the C0 and tau calibration constants

    For computational efficiency, the fit is actually a linear fit to
    the natural log of the voltage: this means only positive values of
    the voltage are used for the fit.
    """
    # Only fit data within the first couple of cooling times: after this the
    # signals are so small that offset errors and noise influence the result.
    vstart = vcool[SKIP:SKIP + 50].mean()
    vthreshold = math.exp(-2) * vstart
    indices = np.flatnonzero(vcool > vthreshold)
    logvcool = np.log(vcool[indices])
    tvecpos = tvec[indices]
    t = tvecpos - tvecpos[0]
    fit = np.polynomial.polynomial.polyfit(t, logvcool, 1)
    c0 = math.exp(fit[0])
    tau = -1 / fit[1]
    logger.debug("Fit paramters %f, %f", c0, 1 / tau)
    return LinearFit(c0, tau)


def calculate_sensitivity(c0_I: float, c0_Q: float, pheat: float) -> float:
    """
    Calculate the sensitivity of a single channel.

    :param c0_I: the C0 value from fitting the I cooling curve
    :param c0_Q: the C0 value from fitting the Q cooling curve
    :param pheat: the heating power applied
    :return: the sensitivity in V/W
    """
    return 2 * math.hypot(c0_I, c0_Q) / pheat


def calculate_cooling_time(tau_I: float, tau_Q: float) -> float:
    """
    Calculate the cooling time of a single channel.

    The cooling time is the average of the fitted decay constants from
    the I and Q cooling curves.

    :param tau_I: the decay constant from the I cooling curve
    :param tau_Q: the decay constant from the Q cooling curve
    :return: the cooling time
    """
    return (tau_I + tau_Q) / 2


def calibrate_single_channel(channel: int,
                             options: argparse.Namespace) -> FitParams:
    """
    Calculate calibration constants for a single channel.

    :param channel: the channel number, from 1 to NCHAN
    :param options: the parsed command-line options
    :return: a tuple (sensitivity, cooling time, ioffset, qoffset)
    """
    try:
        amp, phase, vdc, idc, time = read_channel(channel, options.nsamples,
                                                  options.gainpv, options.root_dir)
        preheating_indices = get_preheating_indices(vdc, options.heating_threshold,
                                                    options.cooling_threshold)
        heating_indices = get_heating_indices(vdc, options.heating_threshold)
        cooling_indices = get_cooling_indices(vdc, options.cooling_threshold,
                                              heating_indices)
        tcool = time[cooling_indices]
        I = 0.5 * amp * np.cos(phase)
        Q = -0.5 * amp * np.sin(phase)
        Ioff = calculate_offset(I[preheating_indices])
        Qoff = calculate_offset(Q[preheating_indices])
        Icool = (I - Ioff)[cooling_indices]
        Qcool = (Q - Qoff)[cooling_indices]
        # For the logarithmic fit to work, both Icool and Qcool must be
        # linearly decreasing. Depending on the phase, the trace may need to be
        # inverted.
        if Icool[:100].mean() < Icool[-100:].mean():
            Icool = -Icool
        if Qcool[:100].mean() < Qcool[-100:].mean():
            Qcool = -Qcool
        c0_I, tau_I = fit_cooling(Icool, tcool)
        c0_Q, tau_Q = fit_cooling(Qcool, tcool)
        iheat = idc[heating_indices] - idc[cooling_indices].mean()
        vheat = vdc[heating_indices]
        pheat = calculate_heating_power(vheat, iheat)
        sens = calculate_sensitivity(c0_I, c0_Q, pheat)
        tau = calculate_cooling_time(tau_I, tau_Q)
        if not np.all(np.isfinite([sens, tau])):
            raise RuntimeError("sens and/or tau are inf or nan.")
    except Exception as e:
        logger.exception("Fit failed for channel %s: %s", channel, e,
                         exc_info=include_traceback)
        sens, tau, Ioff, Qoff = 0, 0, 0, 0
    return FitParams(sens, tau, Ioff, Qoff)


def main():
    """
    Run the calibration program.
    """
    parser = argparse.ArgumentParser(
        description="Calculate calibration constants for foils from calibration data"
    )
    parser.add_argument("-n", "--nsamples", type=int, default=100000,
                        help="Number of samples to read per channel")
    parser.add_argument("-C", "--cooling_threshold", type=float, default=0.01,
                        help="Voltage below which the foil is assumed to be cooling")
    parser.add_argument("-H", "--heating_threshold", type=float, default=0.95,
                        help="Voltage above which the foil is assumed to be heating")
    parser.add_argument("-d", "--root_dir", type=Path, default="/dev/acq400/data",
                        help="Root directory where the data is stored")
    parser.add_argument("-t", "--terse", action="store_true",
                        help="Omit headers from output")
    parser.add_argument("-g", "--gainpv", choices=list(GAIN_PV),
                        help="B8:GAIN value used for all channels to calibrate.")
    parser.add_argument("channels", type=int, nargs="+",
                        help="Channels to calibrate")
    options = parser.parse_args()
    with ProcessPoolExecutor() as ex:
        calibrations = list(ex.map(calibrate_single_channel,
                                   options.channels, repeat(options)))
    if not options.terse:
        print("Channel   Sens      Tau       Ioff          Qoff")
    for channel, (sens, tau, ioff, qoff) in zip(options.channels, calibrations):
        print(f"{channel: <10d}{sens: <14.5g}{tau: <14.5g}{ioff: <14.5g}{qoff:.5g}")


if __name__ == "__main__":
    main()
