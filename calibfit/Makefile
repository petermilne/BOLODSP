# Makefile for compiling the calibration fitting routine for a Busybox Linux
# lmfit is statically compiled into the final executable, so no extra library is needed
# Ensure that arm-xilinx-linux-gnueabi-{gcc g++} are in PATH before calling

CROSS_COMPILE=arm-linux-gnueabihf-
CURR_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
LMDIR=$(CURR_DIR)/lmfit
CFLAGS=-Wall -std=c99 -O3 -mcpu=cortex-a9 -mfpu=neon -ffast-math
CXXFLAGS=-Wall -std=c++11 -O3 -mcpu=cortex-a9 -mfpu=neon -ffast-math -I$(LMDIR)/lib

all: calibfit

lmcurve.o: $(LMDIR)/lib/lmcurve.c
	$(CROSS_COMPILE)gcc -c $(CFLAGS) $(LMDIR)/lib/lmcurve.c

lmmin.o: $(LMDIR)/lib/lmmin.c
	$(CROSS_COMPILE)gcc -c $(CFLAGS) $(LMDIR)/lib/lmmin.c

calibfit: calibfit.cxx lmcurve.o lmmin.o
	$(CROSS_COMPILE)g++ $(CXXFLAGS) -o calibfit calibfit.cxx lmmin.o lmcurve.o

clean:
	rm -f calibfit
	rm -f *.o

.PHONY: clean
