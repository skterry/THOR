# THOR
Terry Hubble Observations of Roman (fields)

THOR is a Data Reduction Pipeline (DRP) for parallel HST Wide-field Camera 3 (WFC3) and Advanced Camera for Surveys (ACS) images taken as part 
of program GO-17776 (Terry et al.), a Precursor Survey of the Roman Galactic Bulge Time Domain Survey (GBTDS) Fields.

### HAMRR
Hubble Advanced Mining Routine for Roman (HAMRR; pronounced "Hammer") is a subroutine which performs a cone search on the HST catalog and returns a list 
of stars with measured photometry, astrometry, color-magnitude diagrams luminosity functions, image cutouts, and more.

## Installation

After cloning, run the following from the top-level directory:

```bash
pip install -e .
```

This single command performs the full setup:

1. **Compiles the Fortran sources** in `src/thor/` (`hst1pass.F`, `thor_go.F`) into
   the bare executables `hst1pass` and `thor_go`. The Fortran compiler is
   auto-detected (gfortran, ifort/ifx, flang, ...); you can force a specific one
   with `make FC=<compiler>` in `src/thor/`, or by setting `FC` before installing.
2. **Downloads the two large (~1 GB each) data archives** from Google Drive into
   `data/`, with a progress bar for each.

Both steps are idempotent — already-built executables and already-downloaded
files are skipped, so re-running `pip install -e .` is cheap. Set
`THOR_SKIP_SETUP=1` to skip the compile/download steps entirely.

> **Requirements:** a Fortran compiler and `make` on your `PATH`. The Python
> build dependencies (`gdown`, `tqdm`) are installed automatically by pip.
