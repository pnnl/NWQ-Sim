include compile_config.mk

XACC_DIR = $(HOME)/.xacc
XACC_LIBS = -I$(XACC_DIR)/include/xacc -I$(XACC_DIR)/include/cppmicroservices4 -I$(XACC_DIR)/include/quantum/gate -L$(XACC_DIR)/lib -lxacc -lxacc-quantum-gate -lCppMicroServices -lxacc-pauli

TARGET = ../bell.cpp

SRC = $(TARGET) ../nwq_accelerator.cpp 

.PHONY: clean

all: xacc_nwq

xacc_nwq: $(SRC)
	$(CC) $(CC_FLAGS) $(CC_LIBS) -I../../include/ $(XACC_LIBS) $(SRC) -o $@

clean:
	rm -f xacc_nwq
