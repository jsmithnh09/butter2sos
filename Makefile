#=================================
# Butterworth design makefile
#=================================

SRC_DIR = $(CURDIR)/include
BUILD_DIR = $(CURDIR)/bin
CORE_FILES = $(SRC_DIR)/butter2sos_design.c $(SRC_DIR)/butter2sos_design.h
CC = gcc
CFLAGS = -std=c99 -I$(SRC_DIR)
DLLFLAGS = -shared
SBAND_TARGET = butter2sos
MBAND_TARGET = butterband
LIB_TARGET = butterlib

# single band design LPF/HPF/APF
sband:
	@echo #== building single-band ==#
	$(CC) $(CFLAGS) $(CORE_FILES) $(SRC_DIR)/$(SBAND_TARGET).c -o $(BUILD_DIR)/$(SBAND_TARGET).exe

# multi-band (bandpass/bandstop) design.
mband:
	@echo #== building multi-band ==#
	$(CC) $(CFLAGS) $(CORE_FILES) $(SRC_DIR)/$(MBAND_TARGET).c -o $(BUILD_DIR)/$(MBAND_TARGET).exe

# DLL generation for Python/Julia/MEX interfacing.
lib:
	@echo #== building library ===#
	$(CC) $(CFLAGS) $(CORE_FILES) $(DLLFLAGS) -o $(BUILD_DIR)/$(LIB_TARGET).dll

libcopy:
	@echo #== copying DLL into relevant sub-project directories ==#
	copy source $(BUILD_DIR)/$(LIB_TARGET).dll destination $(CURDIR)/python/pybutter/src
	copy source $(BUILD_DIR)/$(LIB_TARGET).dll destination $(CURDIR)/julia

.PHONY: clean

all: sband mband lib libcopy

# cleaning is failing on Windows...likely need the Windows "rm" alternative.
clean:
	@echo #== cleaning up ==#
	rm -f $(BUILD_DIR)/*.exe
	rm -f $(BUILD_DIR)/*.dll
