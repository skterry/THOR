"""
image_cut.py — MAST Image Search and Cutout Module
Part of the Hubble Advanced Mining Routine for Roman (HAMRR) package.
Authors: Sean K. Terry, Claude (Anthropic)

Provides get_mast_cutout(), which:
  1. Queries the MAST archive for HST observations from a given program that
     cover a specified sky position.
  2. Downloads the best available drizzled science image for each filter
     (F606W and F814W).
  3. Creates a square pixel cutout centred on the target coordinates for
     each filter.
  4. Saves each cutout as a FITS file in the specified output directory,
     named by filter (e.g. cutout_ra..._dec..._f606w.fits).
  5. Deletes each full-frame download after its cutout is saved.
"""

import os
import urllib.request
import re

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astroquery.mast import Observations

import pdb


# =============================================================================
# Public function
# =============================================================================

def get_mast_cutout(ra: float, dec: float, program_id: int = 17776,
                    size: int = 200, output_dir: str = "output") -> None:
    """
    Query MAST for HST observations from the specified program ID that cover
    the specified sky coordinates. For each available filter (F606W and F814W),
    downloads the best drizzled science image, creates a size×size pixel cutout
    centred on the given coordinates, saves it, then deletes the full-frame
    download to conserve disk space.

    Cutout filenames encode the filter, e.g.:
        cutout_ra266.405_dec-29.007_f606w.fits
        cutout_ra266.405_dec-29.007_f814w.fits

    Args:
        ra         : Right Ascension of the cutout centre (decimal degrees).
        dec        : Declination of the cutout centre (decimal degrees).
        program_id : HST program ID to query (default: 17776).
        size       : Cutout size in pixels — square (default: 200).
        output_dir : Directory in which to save the cutout FITS files (default: "output").

    Returns:
        None — saves cutout FITS files and prints progress to stdout.
    """
    print(f"\nQuerying MAST for HST program {program_id} observations "
          f"covering RA={ra:.3f}°, DEC={dec:.3f}°...")

    coord = SkyCoord(ra, dec, unit='deg', frame='icrs')
    _print_coord_summary(ra, dec)

    # ------------------------------------------------------------------
    # 1. Query MAST
    # ------------------------------------------------------------------
    # Radius is kept tight (~3 arcmin) to avoid pulling in neighbouring HST
    # fields from the same program — ACS/WFC and WFC3/UVIS pointings are
    # typically separated by more than this, so only the field that genuinely
    # contains the target coordinates should be returned.
    obs_table = Observations.query_criteria(
        proposal_id=program_id,
        dataproduct_type="image",
        obs_collection="HST",
        coordinates=coord,
        radius=0.0014 * u.deg, #high-precision to avoid confusion.
    )

    if len(obs_table) == 0:
        print(f"  No observations found for program {program_id}.")
        return

    # Keep only HST observations
    hst_obs = obs_table[obs_table['obs_collection'] == 'HST']
    if len(hst_obs) == 0:
        print(f"  No HST observations found for program {program_id}.")
        return

    # Prefer WFC3/UVIS or ACS/WFC when available
    preferred_instruments = ['WFC3/UVIS', 'ACS/WFC']
    if 'instrument_name' in hst_obs.colnames:
        relevant_obs = hst_obs[np.isin(hst_obs['instrument_name'], preferred_instruments)]
    else:
        relevant_obs = hst_obs

    if len(relevant_obs) == 0:
        relevant_obs = hst_obs   # fall back to all instruments

    #print(f"  Found {len(relevant_obs)} relevant HST observation(s).")

    # ------------------------------------------------------------------
    # 2. Filter observations to the closest field by angular separation
    # ------------------------------------------------------------------
    # Each filter (F606W, F814W) is a separate observation with its own
    # obsid in MAST. Rather than picking a single closest observation, we
    # keep all observations within the field, fetch products for all of
    # them, and pool everything together so _select_products_by_filter
    # can find the right _drc file for each filter independently.
    if 's_ra' in relevant_obs.colnames and 's_dec' in relevant_obs.colnames:
        obs_coords = SkyCoord(
            ra=np.array(relevant_obs['s_ra'].data.astype(float)),
            dec=np.array(relevant_obs['s_dec'].data.astype(float)),
            unit='deg', frame='icrs',
        )
        separations = coord.separation(obs_coords)
        # Keep all observations within 3 arcmin of the closest pointing
        # centre — this captures both filter observations for the same
        # field while still excluding genuinely different fields.
        min_sep = separations.min()
        field_mask = separations <= (min_sep + 3 * u.arcmin)
        field_obs = relevant_obs[field_mask]
    else:
        field_obs = relevant_obs   # fallback if pointing columns are absent

    obs_ids_str = ', '.join(str(o['obs_id']) for o in field_obs)
    #print(f"  Using {len(field_obs)} observation(s) for this field: {obs_ids_str}")

    # ------------------------------------------------------------------
    # 3. Fetch product list once for all field observations, then restrict
    #    strictly to ifhl*_drc.fits files
    # ------------------------------------------------------------------
    # get_product_list() accepts a whole table of observations and makes a
    # single API call, which is much faster than looping per obsid.
    # We then immediately discard everything that isn't an ifhl*_drc.fits
    # file — these are the only products we ever need.
    all_products = Observations.get_product_list(field_obs)

    pattern = re.compile(r'(ifhl|jfhl)[a-z0-9]{2}_drc\.fits$', re.IGNORECASE)

    science_products = all_products[
    [bool(pattern.search(f)) for f in all_products['productFilename']]
    ]

    if len(science_products) == 0:
        print("  No ifhl?? or jfhl?? *_drc.fits products found for this field.")
        return

    print("  Found ifhl*_drc.fits products:")
    #for i, prod in enumerate(science_products):
        #print(f"    [{i}] {prod['productFilename']}")

    # Select one product per filter (F606W and F814W)
    filter_products = _select_products_by_filter(science_products)

    if not filter_products:
        print("  No suitable F606W or F814W ifhl*_drc.fits products found.")
        return

    # ------------------------------------------------------------------
    # 4–5. For each filter: download, validate, cutout, save, delete
    # ------------------------------------------------------------------
    download_dir = './mast_downloads'
    os.makedirs(download_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    for filter_tag, product in filter_products.items():
        print(f"\n  Processing filter: {filter_tag.upper()}")
        print(f"  Downloading product: {product['productFilename']}")

        url = f"https://mast.stsci.edu/api/v0.1/Download/file?uri={product['dataURI']}"
        local_path = os.path.join(download_dir, product['productFilename'])
        print(f"  Downloading from: {url}")
        urllib.request.urlretrieve(url, local_path)

        # Validate: MAST returns an HTML error page for proprietary data
        # instead of a FITS file, which would crash fits.open() confusingly.
        if not _is_valid_fits(local_path):
            print(f"  Error: Downloaded file for {filter_tag.upper()} is not a valid FITS file.")
            print("  This may mean the data are still proprietary or the download failed.")
            print("  Check your MAST credentials (~/.mast token) or try a different program ID.")
            os.remove(local_path)
            continue   # try next filter

        try:
            with fits.open(local_path) as hdul:
                sci_ext = _find_sci_extension(hdul)
                data        = hdul[sci_ext].data
                header      = hdul[sci_ext].header
                pri_header  = hdul[0].header
                wcs         = WCS(header)

                cutout = Cutout2D(data, coord, size, wcs=wcs, mode='partial')

                cutout_header = cutout.wcs.to_header()
                for key in ['FILTER', 'EXPTIME', 'INSTRUME', 'TELESCOP', 'DATE-OBS']:
                    if key in header:
                        cutout_header[key] = (header[key], header.comments[key])
                    elif key in pri_header:
                        cutout_header[key] = (pri_header[key], pri_header.comments[key])
                cutout_header['FILTER'] = (filter_tag.upper(),
                                           'Filter')

                cutout_hdu = fits.PrimaryHDU(data=cutout.data, header=cutout_header)

                out_name = os.path.join(
                    output_dir,
                    f"cutout_ra{ra:.3f}_dec{dec:.3f}_{filter_tag}.fits",
                )
                cutout_hdu.writeto(out_name, overwrite=True)
                print(f"  Cutout saved → {out_name}")

            # Delete full-frame download now the cutout is safely on disk.
            os.remove(local_path) #delete or don't delete the full-frame DRC images.

        except Exception as e:
            # Best-effort cleanup on failure
            if os.path.exists(local_path):
                try:
                    os.remove(local_path)
                    print(f"  Deleted full-frame download after error: {local_path}")
                except OSError:
                    pass
            if "do not overlap" in str(e):
                print(f"  Warning: Target coordinates not in {filter_tag.upper()} image bounds.")
                print(f"  (Coordinates are outside the field of view of "
                      f"{product['productFilename']})")
            else:
                print(f"  Error creating {filter_tag.upper()} cutout: {e}")
            print("  You may need to:")
            print("  - Verify the target coordinates are correct")
            print("  - Adjust the search radius or cone search parameters")
            print("  - The observation found may not fully cover your exact coordinates")
            print("  - Verify the target coordinates are within the field of view.")

    # Remove the download directory now that all full-frame files have been
    # deleted — it will be empty at this point.
    try:
        os.rmdir(download_dir)
    except OSError:
        pass   # non-empty or already removed — leave it alone


# =============================================================================
# Private helpers
# =============================================================================

def _print_coord_summary(ra: float, dec: float) -> None:
    """Print RA/DEC in both decimal and sexagesimal for verification."""
    # RA → HMS
    total_h = ra / 15.0
    h  = int(total_h)
    rm = (total_h - h) * 60
    m  = int(rm)
    s  = (rm - m) * 60
    ra_hms = f"{h:02d}h {m:02d}m {s:05.2f}s"

    # DEC → DMS
    sign    = "-" if dec < 0 else "+"
    abs_dec = abs(dec)
    d  = int(abs_dec)
    rm = (abs_dec - d) * 60
    am = int(rm)
    as_ = (rm - am) * 60
    dec_dms = f"{sign}{d:02d}° {am:02d}' {as_:05.2f}\""

    print(f"  Coordinates: RA  = {ra:.6f}°  ({ra_hms})")
    print(f"               DEC = {dec:+.6f}°  ({dec_dms})")


def _select_products_by_filter(science_products) -> dict:
    """
    Select one _drc science product per filter (F606W and F814W) from a
    pre-filtered MAST product list that already contains only _drc.fits files.

    For each filter the first matching product is taken (the product list
    coming in is already restricted to _drc.fits files by the caller).

    Returns:
        dict mapping filter tag ('f606w', 'f814w') to the chosen product row.
        Only filters for which a product was found are included.
    """
    filters = ['f606w', 'f814w']
    result  = {}

    for filt in filters:
        candidates = [p for p in science_products
                      if filt in p['productFilename'].lower()]

        if not candidates:
            print(f"  No _drc.fits product found for filter {filt.upper()} — skipping.")
            continue

        result[filt] = candidates[0]

    return result


def _find_sci_extension(hdul) -> int:
    """
    Return the index of the science data extension in an open FITS HDUList.
    Looks for an extension named 'SCI'; falls back to extension 1, or 0 for
    simple single-extension files.
    """
    for i, hdu in enumerate(hdul):
        if hdu.data is not None and hdu.header.get('EXTNAME') == 'SCI':
            return i
    return 1 if len(hdul) > 1 else 0


def _is_valid_fits(path: str) -> bool:
    """Return True if the file begins with the FITS magic keyword 'SIMPLE'."""
    try:
        with open(path, 'rb') as f:
            return f.read(6).startswith(b'SIMPLE')
    except OSError:
        return False
