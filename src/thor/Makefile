# ---------------------------------------------------------------------------
# Makefile for the THOR Fortran routines (hst1pass, thor_go_gaia, thor_go_nogaia)
#
# The Fortran compiler is auto-detected so this works on any machine that has
# a reasonably modern Fortran compiler installed (gfortran, ifort/ifx, flang,
# etc.). You can always override the choice explicitly:
#
#     make FC=ifort
#     make FC=/opt/intel/bin/ifx
#
# The two .F sources use cpp directives (#define / #include); the uppercase
# .F extension makes every common compiler run the preprocessor automatically.
# ---------------------------------------------------------------------------

# Auto-detect a Fortran compiler only if the user has not set FC. `FC` is a
# built-in Make variable whose default value is "f77", so `origin FC` is
# "default" until someone overrides it on the command line or environment.
ifeq ($(origin FC),default)
  FC := $(shell \
    for c in gfortran ifx ifort flang-new flang g77 f77; do \
      command -v "$$c" >/dev/null 2>&1 && { echo "$$c"; exit 0; }; \
    done; \
    IFS=:; for dir in $$PATH; do \
      for exe in "$$dir"/gfortran-mp-[0-9]* "$$dir"/gfortran-[0-9]*; do \
        [ -x "$$exe" ] && { echo "$$exe"; exit 0; }; \
      done; \
    done)
endif

# Optional flags. Left empty by default so we stay compiler-agnostic; users
# may pass their own, e.g.  make FFLAGS="-O2".
FFLAGS ?=

EXES := hst1pass thor_go_gaia thor_go_nogaia

.PHONY: all clean check
.DEFAULT_GOAL := all

all: check $(EXES)

# Fail early with a helpful message if no compiler was found.
check:
	@if [ -z "$(strip $(FC))" ]; then \
	  echo "ERROR: no Fortran compiler found on PATH."; \
	  echo "       Install one (e.g. gfortran) or run: make FC=<your-compiler>"; \
	  exit 1; \
	fi
	@echo "Using Fortran compiler: $(FC)"

hst1pass: hst1pass.F
	$(FC) $(FFLAGS) $< -o $@

# Both thor_go variants build the same way; only the source differs
# (thor_go_gaia.F uses a local Gaia_*.dB file, thor_go_nogaia.F does not).
thor_go_%: thor_go_%.F
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -f $(EXES)
