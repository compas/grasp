# Determine the location of the root of the repository
export GRASP := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))

# If present, load the Make.user file which may contain user-defined overrides
# to environment variables.
MAKE_USER_FILE := $(GRASP)/Make.user
ifeq (exists, $(shell [ -e $(MAKE_USER_FILE) ] && echo exists ))
include $(MAKE_USER_FILE)
endif

# Variables affecting the GRASP build
export FC_FLAGS ?= -O2 -fno-automatic
export FC_LD ?=
export FC_MPI ?= mpifort
export LAPACK_LIBS ?= -llapack -lblas
export FC_MPIFLAGS ?= $(FC_FLAGS)
export FC_MPILD ?= $(FC_LD)

LIBRARIES = libmod lib9290 libdvd90 libmcp90 librang90 mpi90
APPLICATIONS = HF jjgen90 rangular90 rbiotransform90 rci90 rcsfgenerate90 \
	rcsfzerofirst90 rmcdhf90 rnucleus90 rtransition90_mpi sms90 jj2lsj90  \
	rangular90_mpi rbiotransform90_mpi rci90_mpi rcsfinteract90 rhfs90    \
	rmcdhf90_mpi  rtransition90  rwfnestimate90

LIBRARY_TARGETS = $(foreach library,$(LIBRARIES),src/lib/$(library))
APPLICATION_TARGETS = $(foreach application,$(APPLICATIONS),src/appl/$(application))

.PHONY: all lib appl tool $(LIBRARY_TARGETS) $(APPLICATION_TARGETS)
all: lib appl tool
appl: $(APPLICATION_TARGETS)
lib: $(LIBRARY_TARGETS)
$(LIBRARY_TARGETS): src/lib/%:
	@echo "Building: $@"
	$(MAKE) -C $@
$(APPLICATION_TARGETS): src/appl/%: lib
	@echo "Building: $@"
	$(MAKE) -C $@
tool: lib
	@echo "Building: src/tool"
	$(MAKE) -C src/tool

LIBRARY_CLEAN_TARGETS = $(foreach library,$(LIBRARIES),clean/lib/$(library))
APPLICATION_CLEAN_TARGETS = $(foreach application,$(APPLICATIONS),clean/appl/$(application))
.PHONY: clean clean/lib clean/appl clean/tool $(LIBRARY_CLEAN_TARGETS) $(APPLICATION_CLEAN_TARGETS)
clean: clean/lib clean/appl clean/tool $(LIBRARY_CLEAN_TARGETS) $(APPLICATION_CLEAN_TARGETS)
clean/lib: $(LIBRARY_CLEAN_TARGETS)
	rm -vf $(GRASP)/lib/*.a
$(LIBRARY_CLEAN_TARGETS): clean/lib/%:
	$(MAKE) -C src/lib/$* clean
clean/appl: $(APPLICATION_CLEAN_TARGETS)
	rm -vf $(GRASP)/bin/*
$(APPLICATION_CLEAN_TARGETS): clean/appl/%:
	$(MAKE) -C src/appl/$* clean
clean/tool:
	$(MAKE) -C src/tool clean
