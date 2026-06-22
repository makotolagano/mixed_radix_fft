VERILATOR ?= verilator
TOP       ?= tb_mr_fft_preadder

SRC_DIR   := rtl/src
TB_DIR    := rtl/tb
BUILD_DIR := build
VCD_FILE  := $(BUILD_DIR)/$(TOP).vcd

RTL_SRCS  := $(wildcard $(SRC_DIR)/*.sv)
TB_SRC    := $(TB_DIR)/$(TOP).sv
SV_SRCS   := $(RTL_SRCS) $(TB_SRC)

.PHONY: all build sim waves clean tree

all: sim

build:
	$(VERILATOR) --binary --timing --trace -trace-structs -Wall --Wno-fatal \
		--top-module $(TOP) \
		--Mdir $(BUILD_DIR) \
		$(SV_SRCS)

sim:
	./$(BUILD_DIR)/V$(TOP)

waves:
	gtkwave $(VCD_FILE)

clean:
	rm -rf $(BUILD_DIR) obj_dir

tree:
	@echo "RTL files:" && ls -1 $(SRC_DIR)/*.sv
	@echo "TB files:" && ls -1 $(TB_DIR)/*.sv
