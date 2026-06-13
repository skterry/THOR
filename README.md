### THOR
Terry Hubble Observations of Roman (THOR) is a Data Reduction Pipeline (DRP) for parallel HST Wide-field Camera 3 (WFC3) and Advanced Camera for Surveys (ACS) images taken as part 
of program GO-17776 (Terry et al.), a Precursor Survey of the Roman Galactic Bulge Time Domain Survey (GBTDS) Fields.

### HAMRR
Hubble Advanced Mining Routine for Roman (HAMRR; pronounced "Hammer") is a subroutine which performs a cone search on the HST catalog and returns a list 
of stars with measured photometry, astrometry, color-magnitude diagrams luminosity functions, image cutouts, and more.

## Installation

Installing THOR/HAMRR is straightforward. From a terminal, type:

```bash
git clone https://github.com/skterry/THOR.git
cd THOR
pip install -e .
```

This install command performs the full setup:

1. **Compiles the Fortran sources** in `src/thor/` (`hst1pass.F`, `thor_go.F`).
   Most Fortran compilers are supported (gfortran, ifort/ifx, flang, ...); however 
   you can force a specific one with `make FC=<compiler>` in `src/thor/`, or by setting `FC` before installing.
2. **Downloads two large data files (~1 GB each)** into `data/`. These are:
   
   a. A sample of images (from field HD_138) to test the THOR data reduction. `data/field_HD138/`
   
   b. The full catalog of sources detected in the HST bulge survey. `data/thor_hst_wfc3_acs_bulge_v0.2_cat.fits.zip`
      note: this is a shared-risk catalog and has not been peer-reviewed or published (as of June 2026).

### Running THOR
Reducing the HST bulge survey data with THOR requires two steps:

1. Reduce the ACS and WFC3 images:

   ```bash
   reduce_acs.src
   reduce_wfc3.src
   ```

2. Collate the reduced data:

   ```bash
   collate_thor.src
   ```

These scripts can be run from inside the directory `example/thor_HD138/` to learn the reduction process. If you have downloaded
external images from the survey and wish to reduce them with THOR, you will need to modify the `*.src` scripts.

### Running HAMRR
Using HAMRR to query the HST bulge survey catalog(s) is relatively straightforward:
```bash
python path/to/hamrr.py
```
HAMRR can be run in two modes: 
1. Interactive-mode via command line interface (CLI).

2. Parameter-file mode: If the file named `params.in` exists in the directory where `hamrr.py` is run, the code will read all 
pre-defined arguments written in the param file.

Required arguments for `hamrr.py`:
- `catalog_path`: Path to THOR-generated catalog (e.g. `../../data/hlsp_thor_hst_wfc3_acs_bulge_early-release_v1.0_cat.fits.zip`)
- `ra`: Right Ascension sky coordinate (`Decimal degrees, HMS/DMS, or Colon-separated`)
- `dec`: Declination sky coordinate (`Decimal degrees, HMS/DMS, or Colon-separated`)
- `radius_arcsec`: Search radius (`arcseconds`)

Optional arguments for `hamrr.py`:
- `save_results`: Save the result table to a file? (`yes/no`) default: no
- `save_format`: File format for saved table (`fits/ecsv/votable`) default: none
- `plot_cmd`: Save the V-I, I color-magnitude diagram (CMD)? (`yes/no`) default: no
- `plot_lf`: Save a luminosity function (LF)? (`yes/no`) default: no
- `lf_band`: Passband for LF (`F606W/F814W`) default: F814W
- `lf_normalise`: Normalize the LF by sky area? (`yes/no`) default: none
- `image_cut`: Save cutout images around target? (`yes/no`) default: no
- `imagecut_size`: Pixel size for the cutout (pix). default: pixel size cutout will match cone search size from `radius_arcsec` argument above
- `output_dir`: Output directory name. default: `output/`


If you use data from the early-release catalog in `data/hlsp_thor_hst_wfc3_acs_bulge_early-release_v1.0_cat.fits.zip`
please cite the following work:

    @article{terry2026hst,
            title={An HST Wide-field Survey of the Galactic Bulge: Overview, Strategy, and First Results},
            author={Terry, Sean K and Anderson, Jay and Beichman, Charles A and Bennett, David P and Bhattacharya, Aparna and Beaulieu, Jean-Philippe and Gaudi, B Scott and Green, Joel and Huston, Macy J and Lu, Jessica R and others},
            journal={The Astrophysical Journal Letters},
            volume={1003},
            number={1},
            pages={L1},
            year={2026},
            publisher={The American Astronomical Society}
            }

If you use data from the full (shared-risk) catalog in `data/thor_hst_wfc3_acs_bulge_v0.2_cat.fits.zip`
please cite the following work:

    @article{terry2026thor,
            title={THOR and HAMRR},
            author={Terry, Sean K and Anderson, Jay},
            journal={arXiv preprint arXiv:}
            year={2026},
            }




